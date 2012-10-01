#include <staden_config.h>

#include <tcl.h>
#include <assert.h>
#include <ctype.h>

#include <tg_gio.h>
#include "misc.h"
#include "export_contigs.h"
#include "gap_cli_arg.h"
#include "newgap_structs.h"
#include "list_proc.h"
#include "xalloc.h"
#include "dna_utils.h"
#include "consensus.h"
#include "sam_index.h"
#include "dstring.h"
#include "tagdb.h"
#include "bam.h"

/* Sequence formats */
#define FORMAT_FASTA 0
#define FORMAT_FASTQ 1
#define FORMAT_SAM   2
#define FORMAT_BAF   3
#define FORMAT_ACE   4
#define FORMAT_CAF   5
#define FORMAT_BAM   6

/* Annotation formats */
#define FORMAT_GFF   10


static int export_contigs(GapIO *io, int cc, contig_list_t *cv, int format,
			  char *fn, int fixmates, int depad);
static int export_tags(GapIO *io, int cc, contig_list_t *cv, int format,
		       int consensus, int unpadded, char *fn);

/* ------------------------------------------------------------------------
 * Tcl interfaces
 */
typedef struct {
    GapIO *io;
    char  *inlist;
    char  *format;
    char  *outfile;
    int    fixmates;
    int    depad;
} ec_arg;

int tcl_export_contigs(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[])
{
    int rargc, format_code, res;
    contig_list_t *rargv;
    ec_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,     offsetof(ec_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,     offsetof(ec_arg, inlist)},
	{"-format",	ARG_STR, 1, "baf",    offsetof(ec_arg, format)},
	{"-outfile",    ARG_STR, 1, "out.baf",offsetof(ec_arg, outfile)},
	{"-fixmates",   ARG_INT, 1, "0",      offsetof(ec_arg, fixmates)},
	{"-depad",      ARG_INT, 1, "1",      offsetof(ec_arg, depad)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (0 == strcmp(args.format, "fasta"))
	format_code = FORMAT_FASTA;
    else if( 0 == strcmp(args.format, "fastq"))
	format_code = FORMAT_FASTQ;
    else if( 0 == strcmp(args.format, "sam"))
	format_code = FORMAT_SAM;
    else if( 0 == strcmp(args.format, "bam"))
	format_code = FORMAT_BAM;
    else if( 0 == strcmp(args.format, "baf"))
	format_code = FORMAT_BAF;
    else if( 0 == strcmp(args.format, "ace"))
	format_code = FORMAT_ACE;
    else if( 0 == strcmp(args.format, "caf"))
	format_code = FORMAT_CAF;
    else
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    res = export_contigs(args.io, rargc, rargv, format_code, args.outfile,
			 args.fixmates, args.depad);

    free(rargv);

    return res == 0 ? TCL_OK : -1;
}

typedef struct {
    GapIO *io;
    char *inlist;
    char *format;
    char *outfile;
    int   consensus; /* map sequence tags to the consensus */
    int   unpadded;
} et_arg;

int tcl_export_tags(ClientData clientData, Tcl_Interp *interp,
		    int objc, Tcl_Obj *CONST objv[])
{
    int rargc, format_code, res;
    contig_list_t *rargv;
    et_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,     offsetof(et_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,     offsetof(et_arg, inlist)},
	{"-format",	ARG_STR, 1, "gff",    offsetof(et_arg, format)},
	{"-outfile",    ARG_STR, 1, "out.gff",offsetof(et_arg, outfile)},
	{"-consensus",  ARG_INT, 1, "1",      offsetof(et_arg, consensus)},
	{"-unpadded",   ARG_INT, 1, "1",      offsetof(et_arg, unpadded)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (0 == strcmp(args.format, "gff"))
	format_code = FORMAT_GFF;
    else
	return TCL_ERROR;

    active_list_contigs_extended(args.io, args.inlist, &rargc, &rargv);

    res = export_tags(args.io, rargc, rargv, format_code,
		      args.consensus, args.unpadded, args.outfile);

    free(rargv);

    return res == 0 ? TCL_OK : -1;
}

/* ------------------------------------------------------------------------
 * export_contigs implementation
 */


/*
 * A FIFO queue. We stack up sequences in here until we have sufficient
 * data to punt the entire thing out including all the annotations.
 */
typedef struct fifo {
    struct fifo *next;
    struct fifo *prev;
    rangec_t r; /* Make this a generic payload? */
    int clipped_pos;
} fifo_t;

typedef struct {
    fifo_t *head;
    fifo_t *tail;
} fifo_queue_t;

static fifo_queue_t *fifo_queue_create(void) {
    fifo_queue_t *f = malloc(sizeof(*f));

    if (!f)
	return NULL;

    f->head = NULL;
    f->tail = NULL;

    return f;
}

static void fifo_queue_destroy(fifo_queue_t *f) {
    fifo_t *i, *n;

    if (!f)
	return;

    for (i = f->head; i; i = n) {
	n = i->next;
	free(i);
    }
    free(f);
}

static fifo_t *fifo_queue_push(fifo_queue_t *f, rangec_t *r) {
    fifo_t *i = malloc(sizeof(*i));

    i->r    = *r;
    i->next = NULL;
    i->prev = f->tail;

    if (!f->head)
	f->head = i;
    if (f->tail)
	f->tail->next = i;
    f->tail = i;

    return i;
}

static fifo_t *fifo_queue_head(fifo_queue_t *f) {
    return f->head;
}

static fifo_t *fifo_queue_pop(fifo_queue_t *f) {
    fifo_t *i, *n;

    i = f->head;
    n = i->next;

    if (f->head)
	f->head = f->head->next;

    if (f->head == NULL)
	f->tail = NULL;

    if (n)
	n->prev = NULL;

    return i;
}

/* Generalised swap the place of two fifo queue members (untested) */
static void fifo_queue_swap(fifo_queue_t *f, fifo_t *a, fifo_t *b) {
    fifo_t *ap = a->prev, *an = a->next;
    fifo_t *bp = b->prev, *bn = b->next;

    if (ap != b) {
	b->prev = ap;
	if (ap)
	    ap->next = b;
    } else {
	a->next = b;
	b->prev = a;
    }
    
    if (an != b) {
	b->next = an;
	if (an)
	    an->prev = b;
    } else {
	a->prev = b;
	b->next = a;
    }

    if (bp != a) {
	a->prev = bp;
	if (bp)
	    bp->next = a;
    } else {
	a->prev = b;
	b->next = a;
    }

    if (bn != a) {
	a->next = bn;
	if (bn)
	    bn->prev = a;
    } else {
	a->next = b;
	b->prev = a;
    }

    if (f->head == a)
	f->head = b;
    else if (f->head == b)
	f->head = a;

    if (f->tail == a)
	f->tail = b;
    else if (f->tail == b)
	f->tail = a;
}

/*
 * Computes and returns a false name for when the sequence has no read name.
 * This is just the record number of the sequence or it's pair (whatever is
 * lowest), optionally with a suffix.
 */
static char *false_name(GapIO *io, seq_t *s, int suffix, int *name_len) {
    static char false_name[1024];
    tg_rec p = sequence_get_pair(io, s);

    if (suffix)
	sprintf(false_name, "seq_%"PRIrec"%s",
		s->rec < p ? p : s->rec,
		s->rec < p ? ".f" : ".r");
    else
	sprintf(false_name, "seq_%"PRIrec,
		s->rec < p ? p : s->rec);

    if (name_len)
	*name_len = strlen(false_name);

    return false_name;
}

static int export_header_sam(GapIO *io, bam_file_t *bf,
			     int cc, contig_list_t *cv,
			     int fixmates) {
    int i;
    dstring_t *ds;

    /* Construct a header string */
    ds = dstring_create(NULL);
    if (!ds)
	return -1;

    /* Generated in sorted order, adhering to version 1.4 */
    dstring_appendf(ds, "@HD\tVN:1.4\tSO:coordinate\n");

    /* Inefficient as we have to loop twice - here and when outputting reads */
    if (fixmates) {
	/* When listing mate pairs we need to list all seqs due to spanning
	 * reads. We can't predict which will be mentioned either.
	 */
	for (i = 0; i < io->db->Ncontigs; i++) {
	    contig_t *c;
	    int len;
	    tg_rec cnum;

	    cnum = arr(tg_rec, io->contig_order, i);

	    if (NULL == (c = cache_search(io, GT_Contig, cnum)))
		return -1;
	    
	    len = c->end;
	    if (c->start <= 0)
		len += 1-c->start;
	    dstring_appendf(ds, "@SQ\tSN:%s\tLN:%d\n",  c->name, len);
	}
    } else {
	for (i = 0; i < cc; i++) {
	    contig_t *c;
	    int len;

	    if (NULL == (c = cache_search(io, GT_Contig, cv[i].contig)))
		return -1;

	    len = c->end;
	    if (c->start <= 0)
		len += 1-c->start;
	    dstring_appendf(ds, "@SQ\tSN:%s\tLN:%d\n",  c->name, len);
	}
    }

    /* Libraries - well read-groups in SAM. */
    for (i = 0; i < io->db->Nlibraries; i++) {
	tg_rec lrec = arr(tg_rec, io->library, i);
	library_t *lib = cache_search(io, GT_Library, lrec);
	if (lib->name) {
	dstring_appendf(ds, "@RG\tID:%s\tSM:%s\tLB:%s\n",
		    lib->name, "unknown", lib->name);
	} else {
	    char buf[100];
	    sprintf(buf, "rg#%"PRIrec, lib->rec);
	    dstring_appendf(ds, "@RG\tID:SM:%s\t%s\tLB:%s\n",
			    buf, "unknown", buf);
	}
    }

    /* Copy it into the bam header and parse it to populate refs[] array */
    if (bf->header)
	free(bf->header);
    bf->header = strdup(dstring_str(ds));
    bf->header_len = strlen(bf->header);
    
    dstring_destroy(ds);

    if (-1 == bam_parse_header(bf))
	return -1;

    bam_write_header(bf);

    return 0;
}

static uint32_t *sam_padded_cigar(char *seq, int left, int right, int olen,
				  int *ncigar_p) {
    static uint32_t *cigar = NULL;
    static int cg_alloc = 0;
    int i, j, ncigar = 0;
    int oplen = 0, cglen = 0;
    enum cigar_op op = BAM_CSOFT_CLIP;

    for (i = j = 0; i < olen; i++) {
	//printf("%c %c\n", s->seq[i], cons[pos+i-c->start]);
	if (ncigar >= cg_alloc) {
	    cg_alloc += 32;
	    cigar = realloc(cigar, cg_alloc*4);
	}

	if (seq[i] == '*') {
	    if (op != BAM_CSOFT_CLIP || i == left) {
		if (op != BAM_CDEL && oplen > 0) {
		    cigar[ncigar++] = (oplen<<4) + op;
		    cglen += oplen;
		    oplen = 0;
		}
		op = BAM_CDEL;
		oplen++;
	    }
	    /* else, gap in soft-clips.
	     * Right now this causes tview to fail an assertion though.
	     */
	} else if (i < left) {
	    if (op != BAM_CSOFT_CLIP && oplen) {
		cigar[ncigar++] = (oplen<<4) + op;
		if (op != BAM_CDEL) cglen += oplen;
		oplen = 0;
	    }
	    op = BAM_CSOFT_CLIP;
	    oplen++;
	} else if (i >= left && i < right) {
	    if (op != BAM_CMATCH && oplen) {
		cigar[ncigar++] = (oplen<<4) + op;
		if (op != BAM_CDEL) cglen += oplen;
		oplen = 0;
	    }
	    op = BAM_CMATCH;
	    oplen++;
	} else {
	    if (op != BAM_CSOFT_CLIP && oplen) {
		/* Switching from deletion to BAM_CSOFT_CLIP breaks tview */
		if (op != BAM_CDEL) {
		    cigar[ncigar++] = (oplen<<4) + op;
		    if (op != BAM_CDEL) cglen += oplen;
		}
		oplen = 0;
	    }
	    op = BAM_CSOFT_CLIP;
	    oplen++;
	}

	j++;
    }
    if (oplen > 0) {
	if (ncigar >= cg_alloc) {
	    cg_alloc += 32;
	    cigar = realloc(cigar, cg_alloc*4);
	}

	cigar[ncigar++] = (oplen<<4) + op;
	if (op != BAM_CDEL) cglen += oplen;
    }

    //assert(cglen == len);

    /* Sam cannot handle gaps at the end of sequences */
    assert(ncigar);
    if ((cigar[ncigar-1] & 15) == BAM_CDEL) {
	ncigar--;
    }

    *ncigar_p = ncigar;
    return cigar;
}

static uint32_t *sam_depadded_cigar(char *seq, int left, int right, int olen,
				    char *cons, int *ncigar_p) {
    static uint32_t *cigar = NULL;
    static int cg_alloc = 0;
    int i, j;
    int oplen = 0, cglen = 0;
    enum cigar_op op = BAM_CSOFT_CLIP, last_op = -1;
    int ncigar = 0, last_ncigar = 0, last_oplen = 0;

    for (i = j = 0; i < olen; i++) {
	//printf("%2d %c %c\n", i, seq[i], 
	//       i >= left && i < right ? cons[j] : '.');
	if (ncigar >= cg_alloc) {
	    cg_alloc += 32;
	    cigar = realloc(cigar, cg_alloc*4);
	}

	if (seq[i] == '*' && i >= left && i < right && cons[j] == '*') {
	    /* gap in both so skip unless we're neighbouring I/D */
	    if (op != BAM_CPAD && oplen) {
		/* May need to merge next op with previous op */
		last_op    = op;
		last_oplen = oplen;
		last_ncigar= ncigar;
		cigar[ncigar++] = (oplen<<4) + op;
		if (op != BAM_CDEL) cglen += oplen;
		oplen = 0;
	    }
	    op = BAM_CPAD;
	    oplen++;
	} else if (seq[i] == '*' && i >= left & i < right) {
	    /* gap is seq only */
	    //if (op != BAM_CSOFT_CLIP) {
		if (op != BAM_CDEL && oplen > 0) {
		    if (op == BAM_CPAD && last_op == BAM_CDEL) {
			oplen = last_oplen;
			ncigar = last_ncigar;
		    } else {
			cigar[ncigar++] = (oplen<<4) + op;
			cglen += oplen;
			oplen = 0;
		    }
		}
		op = BAM_CDEL;
		oplen++;
	     //}
	} else if (seq[i] == '*') {
	    /* else, gap in soft-clips - remove as SAM cannot represent this */
	} else if (i < left) {
	    if (op != BAM_CSOFT_CLIP && oplen) {
		cigar[ncigar++] = (oplen<<4) + op;
		if (op != BAM_CDEL) cglen += oplen;
		oplen = 0;
	    }
	    op = BAM_CSOFT_CLIP;
	    oplen++;
	} else if (i >= left && i < right) {
	    if (cons[j] == '*') {
		/* Insertion */
		if (op != BAM_CINS && oplen) {
		    cigar[ncigar++] = (oplen<<4) + op;
		    if (op != BAM_CDEL) cglen += oplen;
		    oplen = 0;
		}
		/* Handle seqs starting in xPyI */
		if (op == BAM_CSOFT_CLIP) {
		    int plen = 0;
		    while (cons[j-plen-1] == '*')
			plen++;
		    if (plen)
			cigar[ncigar++] = (plen<<4) + BAM_CPAD;
		}
		op = BAM_CINS;
		oplen++;
	    } else {
		if (op != BAM_CMATCH && oplen) {
		    if (op == BAM_CPAD && last_op == BAM_CMATCH) {
			oplen = last_oplen;
			ncigar= last_ncigar;
		    } else {
			cigar[ncigar++] = (oplen<<4) + op;
			if (op != BAM_CDEL) cglen += oplen;
			oplen = 0;
		    }
		}
		op = BAM_CMATCH;
		oplen++;
	    }
	} else {
	    if (op != BAM_CSOFT_CLIP && oplen) {
		/* Switching from deletion to BAM_CSOFT_CLIP breaks tview */
		if (op != BAM_CDEL) {
		    cigar[ncigar++] = (oplen<<4) + op;
		    if (op != BAM_CDEL) cglen += oplen;
		}
		oplen = 0;
	    }
	    op = BAM_CSOFT_CLIP;
	    oplen++;
	}

	if (i >= left)
	    j++;
    }
    if (oplen > 0) {
	if (ncigar >= cg_alloc) {
	    cg_alloc += 32;
	    cigar = realloc(cigar, cg_alloc*4);
	}

	cigar[ncigar++] = (oplen<<4) + op;
	if (op != BAM_CDEL) cglen += oplen;
    }

    //assert(cglen == len);

    /* Sam cannot handle gaps at the end of sequences */
    assert(ncigar);
    if ((cigar[ncigar-1] & 15) == BAM_CDEL) {
	ncigar--;
    }

    *ncigar_p = ncigar;
    return cigar;
}

/*
 * Exports a single consensus tag as a fake sam sequence.
 */
static void sam_export_cons_tag(GapIO *io, bam_file_t *bf, fifo_t *fi,
				contig_t *c, int offset,
				int depad, char *cons,
				int *npad, int *pad_to) {
    static unsigned char *bam = NULL;
    static size_t bam_alloc = 0;
    int bam_idx, bam_len;
    static dstring_t *ds = NULL;
    uint32_t cigar;
    anno_ele_t *a = (anno_ele_t *)cache_search(io, GT_AnnoEle, fi->r.rec);
    int start, end, i;
    char type[5];

#define BAM_EXTEND(sz)					      \
    do {						      \
	while (bam_alloc < bam_idx + (sz)) {		      \
	    bam_alloc *= 2;				      \
	    bam = realloc(bam, bam_alloc);		      \
	    if (bam == NULL) {				      \
		abort();				      \
		fprintf(stderr, "Error allocating memory\n"); \
	    }						      \
        }						      \
    } while (0)

#define BAM_AUX_Z(type, str, len);           \
    do {				     \
	size_t l = (len);		     \
	BAM_EXTEND(l+4);		     \
	bam[bam_idx++] = type[0];	     \
	bam[bam_idx++] = type[1];	     \
	bam[bam_idx++] = 'Z';		     \
	memcpy(&bam[bam_idx], (str), l);     \
	bam_idx += l;			     \
	bam[bam_idx++] = 0;                  \
    } while (0)

    /* Initial allocation of bam string and dstring */
    bam_len = 2 + 9*36 + 4 + 1 + 1;
    if (bam_alloc < bam_len) {
	bam_alloc = bam_len;
	bam = realloc(bam, bam_alloc);
    }

    if (ds)
	dstring_empty(ds);
    else
	ds = dstring_create(NULL);

    /* Depad coords */
    start = fi->r.start;
    end   = fi->r.end;

    if (depad) {
	int np2;

	/* pos is currently padded */
	for (i = *pad_to; i < start; i++) {
	    if (cons[i - (c->start)] == '*')
		(*npad)++;
	}
	*pad_to = i;

	start -= *npad;

	np2 = *npad;
	for (; i < end; i++) {
	    if (cons[i - c->start] == '*')
		np2++;
	}
	end -= np2;
    }


    /* sometimes 784? bottom strand ones? */
    /* Basic SAM record */
    cigar = ((end-start+1)<<4) | BAM_CMATCH;
    bam_idx = bam_construct_seq(bam, bam_alloc,
				"*", 1,
				768,
				bam_name2ref(bf, c->name),
				start + offset,
				start, end,
				255,
				1, &cigar,
				-1, 0, 0,
				0, "", "");
    //dstring_appendf(ds, "*\t768\t%s\t%d\t255\t%dM\t*\t0\t0\t*\t*\t",
    //                c->name, start + offset, end - start + 1);
    
    /* Annotation itself. */
    dstring_appendf(ds, "%c;", a->direction);
    dstring_append_hex_encoded(ds, type2str(fi->r.mqual, type), ";|");
    if (a->comment && *a->comment) {
	dstring_append(ds, ";Note=");
	dstring_append_hex_encoded(ds, a->comment, ";|");
    }

    BAM_AUX_Z("CT", dstring_str(ds), dstring_length(ds));

    //dstring_append_char(ds, '\n');
    //fputs(dstring_str(ds), fp);

    BAM_EXTEND(1); bam[bam_idx] = 0; // terminate aux list.
    ((bam_seq_t *)bam)->blk_size = &bam[bam_idx] -
	(unsigned char *)&((bam_seq_t *)bam)->ref;
    
    bam_put_seq(bf, bam);
}

/*
 * Exports a single sam sequence, marrying up annotations to the appropriate
 * sequences and adding these in auxillary tags.
 *
 * We use the new format PT:Z:<tag_list> for these.
 *
 * Depad indicates whether we wish to output in depadded consensus coordinates.
 * If true, we use cons to identify where the pads are. For efficiency we
 * keep track of the number of pads seen so far up to base 'pad_to', and
 * update these fields as we go (with an assumption data is in sorted order).
 */
static void sam_export_seq(GapIO *io, bam_file_t *bf,
			   fifo_t *fi, fifo_queue_t *tq,
			   int fixmates, tg_rec crec, contig_t *c, int offset,
			   int depad, char *cons, int *npad, int *pad_to) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, fi->r.rec), *sorig = s;
    int len = s->len < 0 ? -s->len : s->len, olen = len;
    int lenQ = 0;

    char *tname, *cp;
    int i, j, flag, tname_len, pos;
    int left, right;
    int iend, isize;
    char *mate_ref;
    library_t *lib = NULL;
    tg_rec last_lrec = -1;
    int mqual;
    char rg_buf[1024], *aux_ptr;
    int first_tag = 1;
    int *depad_map = NULL;
    static unsigned char *bam = NULL;
    static size_t bam_alloc = 0;
    int bam_idx = 0, bam_len, start, end;
    int ref_id, mate_ref_id;
    fifo_t *last, *ti;

    /* Sequence and quality buffers */
    static int qalloc = 0;
    static char *Q = NULL, *S = NULL;

    /* Cigar encoding */
    uint32_t *cigar = NULL, *cigar_tmp;
    int ncigar, ncigar_tmp;

    /* A single SAM line */
    static dstring_t *ds = NULL;

    /*--- Reserve enough memory */
    if (len > qalloc) {
	qalloc = len+1; /* Minimum of 1 byte */
	Q = realloc(Q, qalloc);
	S = realloc(S, qalloc);
    }

    
    /*--- Complement sequence if needed */
    if ((s->len < 0) ^ fi->r.comp) {
	s = dup_seq(s);
	complement_seq_t(s);
    }

    /*--- Depad sequence/qual */
    for (i = j = 0; i < len; i++) {
	int v = s->conf[i];
	if (v < 0) v = 0;
	if (v > 93) v = 93;
	Q[j] = v;
	if (s->conf[i])
	    lenQ = 1;

	if (s->seq[i] == '-')
	    S[j++] = 'n';
	else if (s->seq[i] != '*')
	    S[j++] = s->seq[i];
	/* else gap */
    }
    len = j;
    if (lenQ) {
	lenQ = len;
    } else { 
	lenQ = len;
	memset(Q, 255, len);
    }


    /*--- Best guess at template name */
    tname = s->name;
    tname_len = s->template_name_len;
    if (tname_len > s->name_len)
	tname_len = s->name_len;

    if (tname_len == s->name_len || tname_len == 0) {
	if (cp = strchr(s->name, '/'))
	    tname_len = (int)(cp-s->name);
	else
	    tname_len = s->name_len;
    }

    /*--- Compute SAM flag field */
    flag = 0;
	
    /*
     * Is this correct? Does "mapped in a proper pair" mean they
     * are consistent, for whatever meaning "consistent" has in this
     * context?
     */
    if (fi->r.flags & GRANGE_FLAG_TYPE_PAIRED) {
	flag |= 0x01; /* paired */
	if (!fi->r.pair_rec) {
	    flag |= 0x08; /* mate unmapped */
	} else {
	    /*
	     * Mate may still be unmapped, need to check.
	     * See "fixmates" code below.
	     */
	}
    }

    /* strand */
    if ((sorig->len < 0) ^ fi->r.comp)  flag |= 0x10; /* reverse strand */

    /*
     * 1st/2nd in pair. This is not quite the same as the fwd/rev flags
     * we store I think as we may have fwd/fwd pairs under some
     * protocols. Is this true?
     */
    if (flag & 0x01) {
	/* if ((fi->r.flags & GRANGE_FLAG_END_MASK) !=
	   (fi->r.flags & GRANGE_FLAG_PEND_MASK)) */ {
	       /* Paired data with differing flags, so it's likely valid */
	       if ((fi->r.flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_REV)
		   flag |= 0x80;
	       else
		   flag |= 0x40;

	       /*
		* FIXME: need to check read-group library to know what the
		* expected orientation is for a proper pair.
		* This also doesn't take into account distance, or even
		* which comes first. Ie  --> <--  vs  <-- -->
		*/
	       if ((fi->r.flags & GRANGE_FLAG_END_MASK) !=
		   (fi->r.flags & GRANGE_FLAG_PEND_MASK)) {
		   /* Opposite orientations */
		   //flag |= 0x02; /* proper pair */
	       }

	       /*
		* Mate complemented flag isn't bullet proof. If the bin
		* it is in is complemented then we need to negate the
		* flag, but we can't tell that without finding the mate's
		* bin and checking which is slow. Take best guess?
		*/
	       if (((fi->r.flags & GRANGE_FLAG_COMP2) != 0) ^ (fi->r.comp))
		   flag |= 0x20; /* mate complemented */
	   }
    }

    /* Find depadded coord */
    pos = fi->r.start;
    left  = s->left-1;
    right = s->right;
    pos += left;
	
    if (depad) {
	/* pos is currently padded */
	for (i = *pad_to; i < pos; i++) {
	    if (cons[i - (c->start)] == '*')
		(*npad)++;
	}
	*pad_to = i;

	pos -= *npad;
    }

    /*--- Generate cigar string. */
    if (depad)
	cigar = sam_depadded_cigar(s->seq, left, right, olen,
				   &cons[pos - c->start + *npad],
				   &ncigar);
    else
	cigar = sam_padded_cigar(s->seq, left, right, olen, &ncigar);
    
    /*--- Compute start and end value for bin */
    start = end =  pos + offset;
    for (i = 0; i < ncigar; i++) {
	switch(cigar[i]&15) {
	case BAM_CMATCH:
	case BAM_CDEL:
	case BAM_CREF_SKIP:
	case BAM_CBASE_MATCH:
	case BAM_CBASE_MISMATCH:
	    end += cigar[i]>>4;
	}
    }


    /*--- Compute insert sizes */
    if (fixmates && fi->r.pair_rec) {
	tg_rec other_c;
	int other_st, other_en, other_dir, comp;
	seq_t *pair_s;
	range_t pair_r;

	sequence_get_position2(io, fi->r.pair_rec, &other_c, &other_st,
			       &other_en, &other_dir, NULL, &pair_r, &pair_s);

	comp = (pair_s->len >= 0) ^ other_dir;
	cache_decr(io, pair_s);
	    
	iend = other_st < other_en ? other_st : other_en;

	if (crec == other_c) {
	    int ll, rr;
	    mate_ref = "=";
	    ll = (sorig->len >= 0) ? pos : pos - sorig->len - 1;
	    rr = comp ? other_st-1 : other_en+1;
	    isize = rr-ll;
	} else {
	    contig_t *oc = cache_search(io, GT_Contig, other_c);
	    mate_ref = oc->name;
	    isize = 0;
	}

	if (!comp)
	    flag |=  0x20; /* strand of mate */

	/* FIXME: Can also check here if proper pair, based on
	 * library type.
	 */

	if ((pair_r.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISUMSEQ) {
	    flag |=  0x08; /* mate unmapped */
	    flag &= ~0x02; /* proper pair */
	    isize = 0;
	}
    } else {
	iend = 0;
	isize = 0;
	mate_ref = "*";
    }

    if ((fi->r.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISUMSEQ) {
	flag |=  0x04; /* query unmapped */
	flag &= ~0x02; /* proper pair */
	mqual = 0;
	ncigar_tmp = 0;
	cigar_tmp = NULL;
	isize = 0;
    } else {
	mqual = fi->r.mqual;
	ncigar_tmp = ncigar;
	cigar_tmp = cigar;
    }

    if (ds)
	dstring_empty(ds);
    else
	ds = dstring_create(NULL);


    /* Create BAM struct */
    bam_len = tname_len + len + (len+1)/2 + 9*36 + ncigar_tmp*4;

    /* Initial allocation of bam string */
    if (bam_alloc < bam_len) {
	bam_alloc = bam_len;
	bam = realloc(bam, bam_alloc);
    }

    ref_id = bam_name2ref(bf, c->name);
    mate_ref_id = (mate_ref[0] == '=' && mate_ref[1] == '\0')
	? ref_id : bam_name2ref(bf, mate_ref);
    bam_idx = bam_construct_seq(bam, bam_alloc,
				tname, tname_len,
				flag,
				ref_id,
				pos+offset,
				start, end,
				mqual,
				ncigar_tmp, cigar_tmp,
				mate_ref_id,
				iend,
				isize,
				len, S, Q);

    /*--- Aux strings */
    if (s->parent_type == GT_Library) {
	if (last_lrec != s->parent_rec) {
	    last_lrec  = s->parent_rec;
	    lib = cache_search(io, GT_Library, s->parent_rec);
	}
	if (lib->name)
	    sprintf(rg_buf, "%s", lib->name);
	else
	    sprintf(rg_buf, "rg#%"PRIrec, lib->rec);

	BAM_AUX_Z("RG", rg_buf, strlen(rg_buf));
    }

    if (tname_len != s->name_len) {
	sprintf(rg_buf, "\tFS:Z:%.*s",
		s->name_len - tname_len,
		s->name + tname_len);

	BAM_AUX_Z("FS", s->name + tname_len, s->name_len - tname_len);
    }

    if (s->aux_len) {
	BAM_EXTEND(s->aux_len);
	memcpy(&bam[bam_idx], s->sam_aux, s->aux_len);
	bam_idx += s->aux_len;
    }

    /*--- Attach tags for this sequence too */
    if (fifo_queue_head(tq)) {
	char *d;
	int i = 0, j = 0, l = ABS(s->len), n = 0;
	int op, op_len = 0;
	depad_map = malloc(l * sizeof(*depad_map));
	if (!depad_map)
	    return;

	d = s->seq;
	
	for (i = 0; i < l; i++) {
	    if (op_len == 0) {
		op = cigar[n] & 15;
		op_len = cigar[n] >> 4;
		n++;
	    }

	    depad_map[i] = j;
	    if (d[i] != '*' || op == BAM_CPAD || (!depad && op == BAM_CDEL)) {
		//printf("%2d %c/%c %2d  %d%c\n", i+1, d[i], S[j], j+1, op_len, op);
		op_len--;
		j++;
	    } else {
		//printf("%2d %c/- %2d  -\n", i+1, d[i], j+1);
	    }
	}
    }
    for (last = NULL, ti = fifo_queue_head(tq); ti;) {
	anno_ele_t *a;
	char type[5];
	int st, en;

	if (ti->r.start > fi->r.end)
	    break;

	if (ti->r.pair_rec != fi->r.rec &&
	    (ti->r.flags & GRANGE_FLAG_TAG_SEQ)) {
	    last = ti;
	    ti = ti->next;
	    continue;
	}

	/* Tag is for this seq. */
	a = cache_search(io, GT_AnnoEle, ti->r.rec);
	if (!(ti->r.flags & GRANGE_FLAG_TAG_SEQ)) {
	    last = ti;
	    ti = ti->next;
	    continue;
	}

	
	if (first_tag) {
	    //dstring_append(ds, "\tPT:Z:");
	    first_tag = 0;
	} else {
	    dstring_append(ds, "|");
	}
	    
	if (ti->r.start < fi->r.start)
	    ti->r.start = fi->r.start;
	if (ti->r.start > fi->r.end)
	    ti->r.start = fi->r.end;
	if (ti->r.end < fi->r.start)
	    ti->r.end = fi->r.start;
	if (ti->r.end > fi->r.end)
	    ti->r.end = fi->r.end;

	if ((s->len >= 0) ^ fi->r.comp) {
	    st = ti->r.start - (fi->r.start-1);
	    en = ti->r.end - (fi->r.start-1);
	} else {
	    st = ABS(s->len)+1 - (ti->r.end - (fi->r.start-1));
	    en = ABS(s->len)+1 - (ti->r.start - (fi->r.start-1));
	}
	//printf("%d..%d => %d..%d in seq len %d\n", 
	//       st, en, depad_map[st-1]+1, depad_map[en-1]+1, s->len);
	if (depad_map)
	    dstring_appendf(ds, "%d;%d;%c;",
			    depad_map[st-1]+1, depad_map[en-1]+1,
			    a->direction);
	else
	    dstring_appendf(ds, "%d;%d;%c;", st, en, a->direction);

	dstring_append_hex_encoded(ds, type2str(ti->r.mqual, type), ";|");
	if (a->comment && *a->comment) {
	    dstring_append(ds, ";Note=");
	    dstring_append_hex_encoded(ds, a->comment, ";|");
	}

	/* Remove ti from the list (NB fifo is wrong abstract type) */
	if (last) {
	    last->next = ti->next;
	    if (tq->tail == ti)
		tq->tail = last;
	    free(ti);
	    ti = last->next;
	} else {
	    fifo_queue_pop(tq);
	    free(ti);
	    ti = fifo_queue_head(tq);
	}
    }
    if (!first_tag) {
	BAM_AUX_Z("PT", dstring_str(ds), dstring_length(ds));
    }
    if (depad_map)
	free(depad_map);


    /*--- Finally output it */
    //dstring_append_char(ds, '\n');
    //fputs(dstring_str(ds), fp);

    /*--- Finally output it */
    BAM_EXTEND(1); bam[bam_idx] = 0; // terminate aux list.
    ((bam_seq_t *)bam)->blk_size = &bam[bam_idx] -
	(unsigned char *)&((bam_seq_t *)bam)->ref;

    bam_put_seq(bf, bam);

    if (s != sorig)
	free(s);

}

static int export_contig_sam(GapIO *io, bam_file_t *bf,
			     tg_rec crec, int start, int end,
			     int fixmates, int depad) {
    contig_iterator *ci;
    rangec_t *r;
    contig_t *c;
    char *Q = NULL, *S = NULL;
    char *cigar = NULL;
    int offset, ustart, uend;
    char *cons = NULL;
    /* Seq fragment and tag queues */
    fifo_queue_t *fq = fifo_queue_create(), *tq = fifo_queue_create();
    fifo_t *fi;
    int last_start = 0;
    int npads = 0, pad_to;
    int expanded_start, expanded_end;

    c = (contig_t *)cache_search(io, GT_Contig, crec);
    if (!c)
	return -1;
    cache_incr(io, c);

    iterator_expand_range(io, crec, start, end,
			  &expanded_start, &expanded_end);

    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
				 expanded_start, expanded_end,
				 GRANGE_FLAG_ISANY);

    pad_to = c->start;

    if (depad) {
	int first = c->start;
	int last = c->end;
	int len = last-first+1;

	cons = (char *)xmalloc(len+1);
	cons[0] = 'N'; /* To allow cons[?-1] */
	calculate_consensus_simple(io, crec, first, last, cons+1, NULL);
    }

    /* Sam can only have coordinates from 1 onwards, so shift if needed */
    consensus_valid_range(io, crec, &ustart, &uend);
    offset = ustart < 1
	? 1-ustart
	: 0;

    while (r = contig_iter_next(io, ci)) {
	/* Add new items to fifo */
	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    seq_t *s;
	    fifo_t *fi, *fl;

	    if (r->end < start || r->start > end)
		continue;

	    fi = fifo_queue_push(fq, r);

	    /* Compute clipped start */
	    s = cache_search(io, GT_Seq, fi->r.rec);
	    if ((s->len >= 0) ^ fi->r.comp) {
		fi->clipped_pos = fi->r.start + s->left-1;
	    } else {
		fi->clipped_pos = fi->r.start + (ABS(s->len) - s->right);
	    }

	    /* Sort by clipped left-end */
	    while ((fl = fi->prev) && fi->clipped_pos < fl->clipped_pos) {
		fifo_queue_swap(fq, fi, fl);
	    }

	    last_start = r->start;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    if (r->flags & GRANGE_FLAG_TAG_SEQ) {
		fifo_queue_push(tq, r);
	    } else {
		fifo_t *fi, *fl;

		fi = fifo_queue_push(fq, r); /* consensus tag */
		fi->clipped_pos = r->start;

		/* Sort by clipped left-end */
		while ((fl = fi->prev) && fi->clipped_pos < fl->clipped_pos) {
		    fifo_queue_swap(fq, fi, fl);
		}

		last_start = r->start;
	    }
	}

	/* And pop off when they're history */
	while (fi = fifo_queue_head(fq)) {
	    if (fi->r.end < last_start) {
		fifo_queue_pop(fq);
		if ((fi->r.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
		    sam_export_cons_tag(io, bf, fi, c, offset,
					depad, cons+1, &npads, &pad_to);
		} else {
		    sam_export_seq(io, bf, fi, tq, fixmates, crec, c, offset,
				   depad, cons+1, &npads, &pad_to);
		}
		free(fi);
	    } else {
		break;
	    }
	}
    }

    /* Flush the rest of queues */
    while (fi = fifo_queue_head(fq)) {
	fifo_queue_pop(fq);
	if ((fi->r.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    sam_export_cons_tag(io, bf, fi, c, offset,
				depad, cons+1, &npads, &pad_to);
	} else {
	    sam_export_seq(io, bf, fi, tq, fixmates, crec, c, offset,
			   depad, cons+1, &npads, &pad_to);
	}
	free(fi);
    }

    fifo_queue_destroy(fq);
    fifo_queue_destroy(tq);

    contig_iter_del(ci);
    cache_decr(io, c);

    if (Q) free(Q);
    if (S) free(S);
    if (cigar) free(cigar);
    if (cons) free(cons);

    return 0;
}

/* #define FASTQ_COMPLEMENTED */
static int export_contig_fastq(GapIO *io, FILE *fp,
			       tg_rec crec, int start, int end) {
    contig_iterator *ci = contig_iter_new(io, crec, 0, CITER_FIRST,
					  start, end);
    rangec_t *r;
    int qalloc = 0;
    char *q = NULL;
    char *S = NULL;
    char *name;
    int name_len;
    
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	int len = s->len < 0 ? -s->len : s->len;
	int i, j;

#ifdef FASTQ_COMPLEMENTED
	seq_t *origs = s;
	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}
#endif

	if (len > qalloc) {
	    qalloc = len;
	    q = realloc(q, qalloc);
	    S = realloc(S, qalloc);
	}
	for (i = j = 0; i < len; i++) {
	    int v;

	    if (s->seq[i] == '*') /* Skip pads */
		continue;

	    v = '!' + s->conf[i];
	    if (v < '!') v = '!';
	    if (v > 255) v = 255;
	    q[j] = v;
	    S[j++] = s->seq[i];
	}
	len = j;

	if (s->name_len) {
	    name_len = s->name_len;
	    name = s->name;
	} else {
	    name = false_name(io, s, 1, &name_len);
	}

#ifdef FASTQ_COMPLEMENTED
	fprintf(fp, "@%.*s %d\n%.*s\n+\n%.*s\n",
		name_len, name, s != origs, len, S, len, q);

	if (s != origs)
	    free(s);
#else
	fprintf(fp, "@%.*s\n%.*s\n+\n%.*s\n",
		name_len, name, len, S, len, q);
#endif
    }
    contig_iter_del(ci);

    if (q)
	free(q);
    if (S)
	free(S);

    return 0;
}

static int export_contig_fasta(GapIO *io, FILE *fp,
			       tg_rec crec, int start, int end) {
    contig_iterator *ci = contig_iter_new(io, crec, 0, CITER_FIRST,
					  start, end);
    rangec_t *r;
    char *name;
    int name_len;
    
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	int len = s->len < 0 ? -s->len : s->len;

	if (s->name_len) {
	    name_len = s->name_len;
	    name = s->name;
	} else {
	    name = false_name(io, s, 1, &name_len);
	}

	fprintf(fp, ">%.*s\n%.*s\n", name_len, name, len, s->seq);
    }
    contig_iter_del(ci);

    return 0;
}

static int baf_export_seq(GapIO *io, FILE *fp, fifo_t *fi, fifo_queue_t *tq) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, fi->r.rec);
    int len = s->len < 0 ? -s->len : s->len;
    int i;
    static int qalloc = 0;
    static char *q = NULL, *S = NULL;
    char *name, *cp;
    int name_len;
    fifo_t *last, *ti;

    /* Output the sequence */
    if (len > qalloc) {
	qalloc = len;
	q = realloc(q, qalloc);
	S = realloc(S, qalloc);
    }
    for (i = 0; i < len; i++) {
	int v = '!' + s->conf[i];
	if (v < '!') v = '!';
	if (v > 255) v = 255;
	q[i] = v;

	if (s->seq[i] == '-')
	    S[i] = 'n';
	else if (s->seq[i] == '*')
	    S[i] = '-';
	else
	    S[i] = s->seq[i];
    }

    if (s->name_len) {
	name_len = s->name_len;
	name = s->name;
    } else {
	name = false_name(io, s, 1, &name_len);
    }
    fprintf(fp, "RD=%.*s\n", name_len, name);
    fprintf(fp, "DR=%d\n", (s->len >= 0) ^ fi->r.comp ? 1 : -1);
    fprintf(fp,"AP=%d\n", 
	    (s->len >= 0) ^ fi->r.comp
	    ? fi->r.start + s->left-1
	    : fi->r.start + (ABS(s->len) - s->right));
    fprintf(fp, "QL=%d\n", s->left);
    fprintf(fp, "QR=%d\n", s->right);
    fprintf(fp, "PR=%d\n",
	    (s->flags & SEQ_END_MASK) == SEQ_END_FWD ? 0 : 1);
    if (s->template_name_len == s->name_len) {
	/* Best guess at template name */
	if ((cp = strchr(name, '.')) ||
	    (cp = strchr(name, '/')))
	    fprintf(fp, "TN=%.*s\n", (int)(cp-name), name);
    } else {
	fprintf(fp, "TN=%.*s\n", s->template_name_len, name);
    }
    if (s->trace_name_len)
	fprintf(fp,"TR=%.*s\n", s->trace_name_len, s->trace_name);
    if (s->alignment_len)
	fprintf(fp, "AL=%.*s\n", s->alignment_len, s->alignment);
    if (fi->r.mqual != 255)
	fprintf(fp, "MQ=%d\n", fi->r.mqual);
    fprintf(fp, "SQ=%.*s\n", len, S);
    fprintf(fp, "FQ=%.*s\n\n", len, q);


    /* Find tags on this sequence too */
    for (last = NULL, ti = fifo_queue_head(tq); ti;) {
	anno_ele_t *a;
	char type[5];

	if (ti->r.start > fi->r.end)
	    break;

	if (ti->r.pair_rec != fi->r.rec) {
	    last = ti;
	    ti = ti->next;
	    continue;
	}

	/* Tag is for this seq. */
	a = cache_search(io, GT_AnnoEle, ti->r.rec);
	fprintf(fp, "AN=%s\n", type2str(ti->r.mqual, type));
	if ((s->len >= 0) ^ fi->r.comp) {
	    fprintf(fp, "LO=%d\n", ti->r.start - (fi->r.start-1));
	} else {
	    fprintf(fp, "LO=%d\n", ABS(s->len)+1 - (ti->r.end - (fi->r.start-1)));
	}
	fprintf(fp, "LL=%d\n", ti->r.end - ti->r.start+1);
	if (a->comment && *a->comment) {
	    char *escaped = escape_C_string(a->comment);
	    fprintf(fp, "TX=%s\n", escaped);
	    free(escaped);
	}
	fprintf(fp, "\n");

	/* Remove ti from the list (NB fifo is wrong abstract type) */
	if (last) {
	    last->next = ti->next;
	    if (tq->tail == ti)
		tq->tail = last;
	    free(ti);
	    ti = last->next;
	} else {
	    fifo_queue_pop(tq);
	    free(ti);
	    ti = fifo_queue_head(tq);
	}
    }

    return 0;
}

static int export_contig_baf(GapIO *io, FILE *fp,
			     tg_rec crec, int start, int end) {
    contig_iterator *ci;
    rangec_t *r;
    char *q = NULL, *S = NULL;
    contig_t *c;
    fifo_queue_t *fq = fifo_queue_create(), *tq = fifo_queue_create();
    fifo_t *fi;
    int last_start = 0;
    int expanded_start, expanded_end;

    iterator_expand_range(io, crec, start, end,
			  &expanded_start, &expanded_end);
    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
				 expanded_start, expanded_end,
				 GRANGE_FLAG_ISANY);

    /* Contig record */
    c = (contig_t *)cache_search(io, GT_Contig, crec);
    fprintf(fp, "CO=%s\nLN=%d\n\n", c->name, c->end-c->start+1);

    /* Seq/tag records */
    while (r = contig_iter_next(io, ci)) {
	/* Add new items to fifo */
	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    if (r->end < start || r->start > end)
		continue;

	    fifo_queue_push(fq, r);
	    last_start = r->start;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    if (!(r->flags & GRANGE_FLAG_TAG_SEQ)) {
		/* Contig tag */
		anno_ele_t *a = cache_search(io, GT_AnnoEle, r->rec);
		char type[5];

		fprintf(fp, "AN=%s\nAT=C\n", type2str(r->mqual, type));
		fprintf(fp, "LO=@%d\n", r->start);
		if (r->start != r->end)
		    fprintf(fp, "LL=%d\n", r->end - r->start + 1);

		if (a->comment && *a->comment) {
		    char *escaped = escape_C_string(a->comment);
		    fprintf(fp, "TX=%s\n", escaped);
		    free(escaped);
		}
		fprintf(fp, "\n");

	    } else {
		/* Technically this should be a list and not a fifo */
		fifo_queue_push(tq, r);
	    }
	}

	/* And pop off when they're history */
	while (fi = fifo_queue_head(fq)) {
	    if (fi->r.end < last_start) {
		fifo_queue_pop(fq);
		baf_export_seq(io, fp, fi, tq);
		free(fi);
	    } else {
		break;
	    }
	}
    }

    /* Flush the rest of queues */
    while (fi = fifo_queue_head(fq)) {
	fifo_queue_pop(fq);
	baf_export_seq(io, fp, fi, tq);
	free(fi);
    }

    fifo_queue_destroy(fq);
    fifo_queue_destroy(tq);
    contig_iter_del(ci);

    if (q) free(q);
    if (S) free(S);

    return 0;
}

static int caf_export_seq(GapIO *io, FILE *fp, fifo_t *fi, fifo_queue_t *tq,
			  dstring_t *ds) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, fi->r.rec);
    int len = s->len < 0 ? -s->len : s->len;
    int i, wrap;
    char *name, *cp, name_buf[1024];
    fifo_t *last, *ti;

    /* Output the Sequence record */
    if (s->name_len) {
	memcpy(name_buf, s->name, s->name_len);
	name = name_buf;
	name[s->name_len] = 0;
	/* Best guess - not sure what else we can do atm. */
	if ((cp = strchr(name, '.')) == NULL) {
	    sprintf(name+strlen(name), ".%"PRIrec, s->rec);
	}
    } else {
	name = false_name(io, s, 1, NULL);
    }
    if ((s->len >= 0) ^ fi->r.comp) {
	dstring_appendf(ds, "Assembled_from %s %d %d %d %d\n",
			name,
			fi->r.start + s->left-1,
			fi->r.start + s->right-1,
			s->left, s->right);
    } else {
	dstring_appendf(ds, "Assembled_from %s %d %d %d %d\n",
			name,
			fi->r.start + (ABS(s->len) - s->left),
			fi->r.start + (ABS(s->len) - s->right),
			s->left, s->right);
    }

    fprintf(fp, "Sequence : %s\nIs_read\nPadded\n", name);
    fprintf(fp, "Clipping QUAL %d %d\n", s->left, s->right);
    fprintf(fp, "Strand %s\n", 
	    (s->flags & SEQ_END_MASK) == SEQ_END_FWD ? "Forward" : "Reverse");

    fprintf(fp, "SCF_File %.*s\n", 
	    s->trace_name_len ? s->trace_name_len : (int)strlen(name),
	    s->trace_name_len ? s->trace_name     : name);

    //fprintf(fp, "Insert_size %d %d\n", fixme, fixme);

    if (s->template_name_len == s->name_len) {
	/* Best guess at template name */
	if ((cp = strchr(name, '.')) ||
	    (cp = strchr(name, '/')))
	    fprintf(fp, "Template %.*s\n", (int)(cp-name), name);
    } else {
	    fprintf(fp, "Template %.*s\n", s->template_name_len, name);
    }

    /* Insert size and ligation information if available */
    if (s->parent_type == GT_Library) {
	double mean, sd;
	int count;
	int min = 0, max = 0;
	library_t *lib;

	lib = cache_search(io, GT_Library, s->parent_rec);
	get_library_stats(io, lib->rec, &mean, &sd, NULL, &count);

	/* Parse libraries generated solely from CAF Insert_size records */
	if (strncmp(lib->name, "ins_size=", 9) == 0) {
	    if (sscanf(lib->name, "ins_size=%d-%d", &min, &max) != 2) {
		fprintf(stderr, "Unexpected library name format '%s'\n", 
			lib->name);
		min = max = 0;
	    }
	} else {
	    fprintf(fp, "Ligation_no %s\n", lib->name);
	}

	if ((min == 0 && max == 0) || count < 100) {
	    min = MAX(0, mean - 3*sd);
	    max = MAX(0, mean + 3*sd);
	}
	fprintf(fp, "Insert_size %d %d\n", min, max);
    }

    //fprintf(fp, "Seq_vec SVEC %d %d\n", fixme, fixme);


    /* Find tags on this sequence too */
    for (last = NULL, ti = fifo_queue_head(tq); ti;) {
	anno_ele_t *a;
	char type[5];
	int tpos; 

	if (ti->r.start > fi->r.end)
	    break;

	if (ti->r.pair_rec != fi->r.rec) {
	    last = ti;
	    ti = ti->next;
	    continue;
	}

	/* Tag is for this seq. */
	a = cache_search(io, GT_AnnoEle, ti->r.rec);
	tpos = (s->len >= 0) ^ fi->r.comp
	    ? ti->r.start - (fi->r.start-1)
	    : ABS(s->len)+1 - (ti->r.end - (fi->r.start-1));

	fprintf(fp, "Tag %s %d %d",
		type2str(ti->r.mqual, type),
		tpos,
		tpos + ti->r.end - ti->r.start);
	if (a->comment && *a->comment) {
	    char *escaped = escape_C_string(a->comment);
	    fprintf(fp, " \"%s\"\n", escaped);
	    free(escaped);
	} else {
	    fprintf(fp, "\n");
	}

	/* Remove ti from the list (NB fifo is wrong abstract type) */
	if (last) {
	    last->next = ti->next;
	    if (tq->tail == ti)
		tq->tail = last;
	    free(ti);
	    ti = last->next;
	} else {
	    fifo_queue_pop(tq);
	    free(ti);
	    ti = fifo_queue_head(tq);
	}
    }

    /* And finally sequence / quality records */
    fprintf(fp, "\nDNA : %s\n", name);
    for (i = 0; i < len; i += 60) {
	int j,l = MIN(60, len-i);
	for (j=0; j<l; j++) {
	    switch (s->seq[i+j]) {
	    case '*': fputc('-', fp); break;
	    case '-': fputc('N', fp); break;
	    default:  fputc(s->seq[i+j], fp); break;
	    }
	}
	fputc('\n', fp);
    }
    fprintf(fp, "\nBaseQuality : %s\n", name);
    for (wrap = i = 0; i < len; i++) {
	wrap += fprintf(fp, "%d ", s->conf[i]);
	if (wrap > 60) {
	    wrap = 0;
	    fprintf(fp, "\n");
	}
    }
    if (wrap)
	fprintf(fp, "\n");
    fprintf(fp, "\n");

    return 0;
}

static int export_contig_caf(GapIO *io, FILE *fp,
			     tg_rec crec, int start, int end) {
    contig_iterator *ci;
    rangec_t *r;
    char *q = NULL, *S = NULL;
    contig_t *c;
    fifo_queue_t *fq = fifo_queue_create(), *tq = fifo_queue_create(),
	*cq = fifo_queue_create();
    fifo_t *fi;
    int last_start = 0;
    dstring_t *ds = dstring_create(NULL);
    int first_base, last_base, len;
    char *cons;
    float *qual;
    int i, wrap;
    int expanded_start, expanded_end;

    iterator_expand_range(io, crec, start, end,
			  &expanded_start, &expanded_end);

    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
				 expanded_start, expanded_end,
				 GRANGE_FLAG_ISANY);

    /* Seq/tag records */
    while (r = contig_iter_next(io, ci)) {
	/* Add new items to fifo */
	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    if (r->end < start || r->start > end)
		continue;

	    fifo_queue_push(fq, r);
	    last_start = r->start;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    if (!(r->flags & GRANGE_FLAG_TAG_SEQ)) {
		/* Contig tag - delay until later */
		fifo_queue_push(cq, r);
	    } else {
		/* Technically this should be a list and not a fifo */
		fifo_queue_push(tq, r);
	    }
	}

	/* And pop off when they're history */
	while (fi = fifo_queue_head(fq)) {
	    if (fi->r.end < last_start) {
		fifo_queue_pop(fq);
		caf_export_seq(io, fp, fi, tq, ds);
		free(fi);
	    } else {
		break;
	    }
	}
    }

    /* Flush the rest of queues */
    while (fi = fifo_queue_head(fq)) {
	fifo_queue_pop(fq);
	caf_export_seq(io, fp, fi, tq, ds);
	free(fi);
    }

    fifo_queue_destroy(fq);
    fifo_queue_destroy(tq);
    contig_iter_del(ci);

    
    /* Now also write out the contig record */
    c = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, c);
    fprintf(fp, "Sequence : %s\nIs_contig\nPadded\n", c->name);
    fputs(dstring_str(ds), fp);

    /* And the consensus tags */
    while (fi = fifo_queue_head(cq)) {
	fifo_queue_pop(cq);
	anno_ele_t *a;
	char type[5];

	a = cache_search(io, GT_AnnoEle, fi->r.rec);
	fprintf(fp, "Tag %s %d %d",
		type2str(fi->r.mqual, type),
		fi->r.start, fi->r.end);

	if (a->comment && *a->comment) {
	    char *escaped = escape_C_string(a->comment);
	    fprintf(fp, " \"%s\"\n", escaped);
	    free(escaped);
	} else {
	    fprintf(fp, "\n");
	}
	free(fi);
    }
    fifo_queue_destroy(cq);
    fprintf(fp, "\n");


    /* Plus DNA and BaseQuality parts */
    consensus_valid_range(io, crec, &first_base, &last_base);
    len = last_base - first_base + 1;
    cons = (char *)xmalloc(len);
    qual = (float *)xmalloc(len * sizeof(*qual));
    calculate_consensus_simple(io, crec, first_base, last_base, cons, qual);

    fprintf(fp, "DNA : %s\n", c->name);
    for (i = 0; i < len; i += 60) {
	int j,l = MIN(60, len-i);
	for (j=0; j<l; j++) {
	    fputc(cons[i+j] != '*' ? cons[i+j] : '-', fp);
	}
	fputc('\n', fp);
    }

    fprintf(fp, "\nBaseQuality : %s\n", c->name);
    for (wrap = i = 0; i < len; i++) {
	int qv = qual[i]+.5;
	wrap += fprintf(fp, "%d ", qv);
	if (wrap > 60) {
	    wrap = 0;
	    fprintf(fp, "\n");
	}
    }
    if (wrap)
	fprintf(fp, "\n");
    fprintf(fp, "\n");

    free(cons);
    free(qual);

    cache_decr(io, c);
    if (q) free(q);
    if (S) free(S);
    dstring_destroy(ds);

    return 0;
}

static int export_header_ace(GapIO *io, FILE *fp,
			     int cc, contig_list_t *cv) {
    int nreads = 0, i;

    /* Inefficient as we have to loop twice - here and when outputting reads */
    for (i = 0; i < cc; i++) {
	contig_iterator *ci = contig_iter_new(io, cv[i].contig, 0, CITER_FIRST,
					      cv[i].start, cv[i].end);
	while (contig_iter_next(io, ci))
	    nreads++;

	contig_iter_del(ci);
    }

    fprintf(fp, "AS %d %d\n", cc, nreads);
    return 0;
}

/* NB: Ignores start and end */
static int export_contig_ace(GapIO *io, FILE *fp,
			     tg_rec crec, int start, int end) {
    contig_iterator *ci;
    rangec_t *r;
    int qalloc = 0;
    char *q = NULL, *S = NULL;
    contig_t *c;
    int nreads = 0;
    int i, j, slen, nBS, last, first_base, last_base;
    char *cons;
    float *qual;

    /* Contig header - inefficient! */
    c = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    consensus_valid_range(io, crec, &first_base, &last_base);

    ci = contig_iter_new(io, crec, 0, CITER_FIRST, c->start, c->end);
    last = 1;
    nBS = 0;
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	nreads++;
	if (r->start + s->right-1 > last) {
	    nBS++;
	    last = r->start + s->right-1;
	}
    }
    slen = last_base - first_base + 1;

    fprintf(fp, "\nCO %s %d %d %d U\n",
	    c->name, slen, nreads, nBS);


    /* Contig sequence */
    cons = (char *)xmalloc(slen);
    qual = (float *)xmalloc(slen * sizeof(*qual));
    calculate_consensus_simple(io, crec, first_base, last_base, cons, qual);
    for (i = 0; i < slen; i += 50) {
	fprintf(fp, "%.*s\n", slen-i > 50 ? 50 : slen-i, &cons[i]);
    }

    /* Contig Quality */
    fprintf(fp, "\nBQ\n");
    for (i = j = 0; i < slen; i++) {
	int qv;
	if (cons[i] == '*')
	    continue;

	qv = qual[i]+0.5;
	if (qv < 0) qv = 0;
	if (qv > 97) qv = 97;
	fprintf(fp, " %d", qv);

	if (++j == 150) {
	    j = 0;
	    fprintf(fp, "\n");
	}
    }
    fprintf(fp, "\n");

    /* Contig AF records - inefficient again */
    fprintf(fp, "\n");
    contig_iter_del(ci);
    ci = contig_iter_new(io, crec, 0, CITER_FIRST, c->start, c->end);
    first_base--;
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	fprintf(fp, "AF %s.%c %c %d\n",
		s->name,
		(r->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD
		    ?'f' :'r',
		"UC"[(s->len < 0) ^ r->comp], r->start - first_base);
    }

    /* Contig BS records - ignore for now? */
    contig_iter_del(ci);
    ci = contig_iter_new(io, crec, 0, CITER_FIRST, c->start, c->end);
    last = 0;
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	if (r->start + s->right-1 - first_base > last) {
	    fprintf(fp, "BS %d %d %s.%c\n",
		    MAX(last+1, r->start + s->left-1 - first_base),
		    r->start + s->right-1 - first_base, s->name,
		    (r->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD
		        ?'f' :'r');
	    last = r->start + s->right-1 - first_base;
	}
    }

    /* Seq records themselves */
    contig_iter_del(ci);
    ci = contig_iter_new(io, crec, 0, CITER_FIRST, c->start, c->end);
    while (r = contig_iter_next(io, ci)) {
	seq_t *origs, *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	int len = s->len < 0 ? -s->len : s->len;

	origs = s;
	if (len > qalloc) {
	    qalloc = len;
	    q = realloc(q, qalloc);
	    S = realloc(S, qalloc);
	}

	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	for (i = 0; i < len; i++) {
	    int v = '!' + s->conf[i];
	    if (v < '!') v = '!';
	    if (v > 255) v = 255;
	    q[i] = v;

	    if (s->seq[i] == '-')
		S[i] = 'n';
	    else if (s->seq[i] == '*')
		S[i] = '-';
	    else
		S[i] = s->seq[i];
	}

	fprintf(fp, "\nRD %.*s.%c %d 0 0\n",
		s->name_len, s->name,
		(r->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD
		    ?'f' :'r',
		len);
	for (i = 0; i < len; i += 50) {
	    fprintf(fp, "%.*s\n", len-i > 50 ? 50 : len-i, &s->seq[i]);
	}

	fprintf(fp, "\nQA %d %d %d %d\n",
		s->left, s->right, s->left, s->right);

	if (s->template_name_len != s->name_len) {
	    fprintf(fp, "DS CHROMAT_FILE: %s PHD_FILE: %s.phd.1 TEMPLATE: %.*s TIME: Thu Jan 01 00:00:01 GMT 1970\n",
		    s->name, s->name, s->template_name_len, s->name);
	} else {
	    fprintf(fp, "DS CHROMAT_FILE: %s PHD_FILE: %s.phd.1 TIME: Thu Jan 01 00:00:01 GMT 1970\n", s->name, s->name);
	}

	if (s != origs)
	    free(s);
    }
    contig_iter_del(ci);

    cache_decr(io, c);

    if (q) free(q);
    if (S) free(S);

    return 0;
}

static int export_contigs(GapIO *io, int cc, contig_list_t *cv, int format,
			  char *fn, int fixmates, int depad) {
    int i;
    FILE *fp = NULL;
    bam_file_t *bf = NULL;
    
    switch (format) {
    case FORMAT_SAM:
	if (NULL == (bf = bam_open(fn, "w"))) {
	    perror(fn);
	    return -1;
	}
	break;

    case FORMAT_BAM:
	if (NULL == (bf = bam_open(fn, "wb"))) {
	    perror(fn);
	    return -1;
	}
	break;

    default:
	if (0 == strcmp(fn, "-")) {
	    fp = stdout;
	} else {
	    if (NULL == (fp = fopen(fn, "w"))) {
		perror(fn);
		return -1;
	    }
	}
    }

    /* Header */
    switch (format) {
    case FORMAT_ACE:
	export_header_ace(io, fp, cc, cv);
	break;

    case FORMAT_SAM:
    case FORMAT_BAM:
	export_header_sam(io, bf, cc, cv, fixmates);
	break;
    }

    /* Per contig */
    for (i = 0; i < cc; i++) {
	switch (format) {
	case FORMAT_SAM:
	case FORMAT_BAM:
	    export_contig_sam(io, bf, cv[i].contig, cv[i].start, cv[i].end,
			      fixmates, depad);
	    break;

	case FORMAT_FASTQ:
	    export_contig_fastq(io, fp, cv[i].contig, cv[i].start, cv[i].end);
	    break;

	case FORMAT_FASTA:
	    export_contig_fasta(io, fp, cv[i].contig, cv[i].start, cv[i].end);
	    break;

	case FORMAT_BAF:
	    export_contig_baf(io, fp, cv[i].contig, cv[i].start, cv[i].end);
	    break;

	case FORMAT_ACE:
	    export_contig_ace(io, fp, cv[i].contig, cv[i].start, cv[i].end);
	    break;

	case FORMAT_CAF:
	    export_contig_caf(io, fp, cv[i].contig, cv[i].start, cv[i].end);
	    break;

	default:
	    verror(ERR_WARN, "export_contigs", "Unknown format code %d",
		   format);
	    fclose(fp);
	    return 1;
	}
    }

    switch (format) {
    case FORMAT_SAM:
    case FORMAT_BAM:
	bam_close(bf);
	break;

    default:
	if (fp != stdout)
	    fclose(fp);
    }

    return 0;
}

/* ------------------------------------------------------------------------
 * export_tags implementation
 */

/*
 * Breaks an ACD format tag down into key=value constituents.
 *
 * Fills out the key and value strings (with lengths) and returns the
 * next str value to pass into this function, or NULL when no more
 * key=value pairs have been found.
 *
 * On input key_len and val_len are the size of the supplied key and val
 * buffers. On return they hold the filled lengths.
 */
static char *parse_acd_tag(char *str,
			   char *key, size_t *key_len,
			   char *val, size_t *val_len) {
    char *cp;
    size_t len;

    /* Key */
    for (len = 0, cp = str; len < *key_len && *cp && *cp != '='; len++)
	*key++ = *cp++;
    if (!*cp)
	return NULL;
    *key_len = len;
    str = ++cp;
	
    /* Value */
    for (len = 0, cp = str; len < *val_len && *cp && *cp != '\n'; len++) {
	if (*cp == '\\') {
	    if (*++cp == 'n') {
		*val++ = '\n';
		cp++;
	    } else {
		*val++ = *cp++;
	    }
	} else {
	    *val++ = *cp++;
	}
    }

    *val_len = len;
    str = ++cp;

    return str;
}

static int export_header_tags_gff(GapIO *io, FILE *fp,
				  int cc, contig_list_t *cv) {
    fprintf(fp, "##gff-version 3\n");
    return 0;
}

/*
 * Example gff output:
 *
 * (obj prog type start end  score strand frame attributes)
 * seq1 gap5 COMM 100   200  .     .      .     arbitrary free text?
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int export_tags_gff(GapIO *io, FILE *fp,
			   tg_rec crec, int start, int end,
			   int consensus, int unpadded) {
    contig_iterator *ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
						  start, end,
						  GRANGE_FLAG_ISANNO);
    rangec_t *r;
    contig_t *c;
    char *con = NULL;
    int *map = NULL;
    int db;
    char type[5];

    c = cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    /* Generate a padded to unpadded mapping table */
    if (unpadded) {
	int i, np;
	int cstart;

	if (NULL == (con = malloc(c->end - c->start + 2))) {
	    cache_decr(io, c);
	    return -1;
	}
	if (NULL == (map = malloc((c->end - c->start + 2) *
				  sizeof(int)))) {
	    cache_decr(io, c);
	    free(con);
	    return -1;
	}
	calculate_consensus_simple(io, crec, c->start, c->end, con, NULL);

	if (-1 == consensus_valid_range(io, crec, &cstart, NULL))
	    cstart = c->start;

	for (np = 0, i = c->start; i <= c->end; i++) {
	    map[i-c->start] = i-cstart+1 - np;
	    if (con[i-c->start] == '*')
		np++;
	}
    }

    /* Export */
    while (r = contig_iter_next(io, ci)) {
	anno_ele_t *a = cache_search(io, GT_AnnoEle, r->rec);
	char *name, *ename, *etype;
	int st, en;
	char score[100];
	char phase;
	char gff_type[1024];
	static dstring_t *ds = NULL;

	if (ds)
	    dstring_empty(ds);
	else
	    ds = dstring_create(NULL);

	cache_incr(io, a);

	if (!consensus && (r->flags & GRANGE_FLAG_TAG_SEQ)) {
	    int seq_start, seq_end;
	    seq_t *s;

	    if (r->flags & GRANGE_FLAG_TAG_SEQ) {
		/* Seq tag */
		sequence_get_position(io, r->pair_rec, NULL,
				      &seq_start, &seq_end, NULL);

		st = r->start - seq_start;
		en = r->end - seq_start;

		s = cache_search(io, GT_Seq, r->pair_rec);
		name = s->name;

		if (unpadded) {
		    int i, len = s->len < 0 ? -s->len : s->len;
		    int npads;
		    int dupped = 0;

		    /*--- Complement sequence if needed */
		    if ((s->len < 0) ^ r->comp) {
			s = dup_seq(s);
			complement_seq_t(s);
			dupped = 1;
		    }

		    /*--- Compensate for padding in sequence */
		    for (i = npads = 0; i < len && i < st; i++) {
			if (s->seq[i] == '*')
			    npads++;
		    }
		    st -= npads;
		    for (; i < len && i < en; i++) {
			if (s->seq[i] == '*')
			    npads++;
		    }
		    en -= npads;

		    if (dupped)
			free(s);
		}

	    } else {
		/* Consensus tag */
		name = c->name;

		st = r->start;
		en = r->end;

		if (unpadded) {
		    if (st >= c->start && st <= c->end)
			st = map[st - c->start];
		    if (en >= c->start && en <= c->end)
			en = map[en - c->start];
		}
	    }
	} else {
	    /* Tag as if on consensus, regardless whether it already is */
	    name = c->name;

	    st = r->start;
	    en = r->end;

	    /* Depad it too */
	    if (unpadded) {
		if (st >= c->start && st <= c->end)
		    st = map[st - c->start];
		if (en >= c->start && en <= c->end)
		    en = map[en - c->start];
	    }
	}

	assert(st <= en);
	
	score[0] = '.'; score[1] = 0;
	phase = '.';
	strcpy(gff_type, "remark");

	/*
	 * ACDtags can be broken down properly into key=value pairs.
	 * Everything else is a free-for-all and is simply shoved into
	 * the Note= field.
	 */
	ename = escape_hex_string(name, "!#&'(),/\\=<>{}[]`~");
	etype = escape_hex_string(type2str(r->mqual, type), ",=;");
	
	dstring_appendf(ds, "type=%s", etype);

	db = idToIndex(type);
	if (tag_db[db].default_text &&
	    strncmp(tag_db[db].default_text, "#!acdtag:", 9) == 0) {
	    char *cp = a->comment ? a->comment : "";
	    char key[1024], val[8192];
	    size_t key_len, val_len;

	    printf("ACD TAG with comm %s\n", a->comment ? a->comment : "(null)");
	    while (cp = parse_acd_tag(cp,
				      key, (key_len=1024, &key_len),
				      val, (val_len=8192, &val_len))) {
		printf("Key='%.*s' val='%.*s'\n",
		       key_len, key, val_len, val);
		if (key_len == 5 && strncmp(key, "score", 5) == 0) {
		    if ((val_len > 0) &&
			(*val == '+' || *val == '-' || isdigit(*val))) {
			sprintf(score, "%.10g", atof(val));
		    }
		} else if (key_len == 5 && strncmp(key, "phase", 5) == 0) {
		    if (val_len > 0 && *val >= '0' && *val <= '2')
			phase = *val;
		} else if (key_len == 8 && strncmp(key, "gff_type", 8) == 0) {
		    size_t l = MIN(1023, val_len);
		    strncpy(gff_type, val, l);
		    gff_type[l] = 0;
		} else if (key_len == 11 &&
			   strncmp(key, "gff_attribs", 11) == 0) {
		    char *k0, *k1, *v0, *v1;
		    /* GFF attribute key=value pairs packed into a single
		     * ACDtag text item.
		     */
		    for (v1 = val;;) {
			char tmp_k, tmp_v, *esc_k, *esc_v;

			for (k0 = k1 =v1; val_len && *k1 != '=' && *k1 != '\n'; k1++, val_len--)
			    ;
			if (k1 == k0)
			    break;

			if (*k1 != '=') {
			    v0 = k0;
			    v1 = k1;
			    k0 = "note";
			} else {
			    val_len--;
			    for (v0 = v1 =k1+1; val_len && *v1 != '\n'; v1++, val_len--)
				;
			}

			tmp_k = *k1; tmp_v = *v1;
			*k1 = 0; *v1 = 0;
			esc_k = escape_hex_string(k0, ",=;");
			esc_v = escape_hex_string(v0, ",=;");
			*k1 = tmp_k; *v1 = tmp_v;
			dstring_appendf(ds, ";%s=%s", esc_k, esc_v);
			 
			if (!val_len || *v1++ != '\n')
			    break;
			val_len--;
		    }
		} else {
		    /* Everything else is just additional key=value attribs */
		    char tmp = val[val_len], *escaped;
		    val[val_len] = 0;

		    dstring_append_char(ds, ';');
		    escaped = escape_hex_string(val, ",=;");
		    dstring_appendf(ds, "%.*s=%s", key_len, key, escaped);
		    free(escaped);

		    val[val_len] = tmp;
		}
	    }
	} else {
	    if (a->comment && *a->comment) {
		char *escaped = escape_hex_string(a->comment, ",=;");
		dstring_appendf(ds, ";Gap5=1;Note=%s", escaped);
		free(escaped);
	    }
	}

	fprintf(fp, "%s\tgap5\t%s\t%d\t%d\t%s\t%c\t%c\t%s\n",
		ename, gff_type, st, en, score, a->direction,
		phase, dstring_str(ds));

	free(ename);
	free(etype);

	cache_decr(io, a);
    }
    contig_iter_del(ci);

    if (con)
	free(con);
    if (map)
	free(map);

    cache_decr(io, c);

    return 0;
}

static int export_tags(GapIO *io, int cc, contig_list_t *cv, int format,
		       int consensus, int unpadded, char *fn) {
    int i;
    FILE *fp;
    
    if (NULL == (fp = fopen(fn, "w"))) {
	perror(fn);
	return -1;
    }

    /* Header */
    switch (format) {
    case FORMAT_GFF:
	export_header_tags_gff(io, fp, cc, cv);
	break;
    }

    /* Per contig */
    for (i = 0; i < cc; i++) {
	switch (format) {
	case FORMAT_GFF:
	    export_tags_gff(io, fp, cv[i].contig, cv[i].start, cv[i].end,
			    consensus, unpadded);
	    break;

	default:
	    verror(ERR_WARN, "export_tags", "Unknown format code %d",
		   format);
	    fclose(fp);
	    return 1;
	}
    }

    fclose(fp);

    return 0;
}
