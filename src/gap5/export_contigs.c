#include <tcl.h>
#include <assert.h>

#include <tg_gio.h>
#include "export_contigs.h"
#include "gap_cli_arg.h"
#include "newgap_structs.h"
#include "list_proc.h"
#include "xalloc.h"
#include "dna_utils.h"
#include "consensus.h"
#include "sam_index.h"
#include "dstring.h"

/* Sequence formats */
#define FORMAT_FASTA 0
#define FORMAT_FASTQ 1
#define FORMAT_SAM   2
#define FORMAT_BAF   3
#define FORMAT_ACE   4
#define FORMAT_CAF   5

/* Annotation formats */
#define FORMAT_GFF   10


static int export_contigs(GapIO *io, int cc, contig_list_t *cv, int format,
			  char *fn, int fixmates);
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
			 args.fixmates);

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

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    res = export_tags(args.io, rargc, rargv, format_code,
		      args.consensus, args.unpadded, args.outfile);

    free(rargv);

    return res == 0 ? TCL_OK : -1;
}

/* ------------------------------------------------------------------------
 * export_contigs implementation
 */

/*
 * Allocates and returns an escaped version of str. This relaces quotes,
 * newlines, and other non-printable characters with backslashed versions of
 * them in a C string style formatting.
 *
 * Returns malloced string on success
 *         NULL on failure.
 */
static char *escape_C_string(char *str) {
    size_t l = strlen(str);
    size_t new_l = l*1.1+10;
    char *new = malloc(new_l);
    size_t oi, ni;
    static char type[256];
    static int type_init = 0;

    /* A once-only lookup table to speed up the loop below */
    if (!type_init) {
	int i;
	
	for (i = 0; i < 256; i++) {
	    if (isprint(i) && i != '"' && i != '\\') {
		/* directly printable */
		type[i] = 0;
	    } else {
		switch(i) {
		    /* backslash single-char */
		case '"':
		case '\\':
		    type[i] = i;
		    break;
		case '\n':
		    type[i] = 'n';
		    break;
		case '\r':
		    type[i] = 'r';
		    break;
		case '\t':
		    type[i] = 't';
		    break;
		case '\a':
		    type[i] = 'a';
		    break;

		default:
		    type[i] = 1; /* octal escape */
		}
	    }
	}
	type_init = 1;
    }


    if (!new)
	return NULL;

    for (oi = ni = 0; oi < l; oi++) {
	char c = str[oi];

	/* Make enough room */
	if (ni + 5 >= new_l) {
	    new_l = new_l * 1.2 + 10;
	    if (NULL == (new = realloc(new, new_l)))
		return NULL;
	}

	switch(type[(unsigned char)c]) {
	case 0:
	    new[ni++] = c;
	    break;
	    
	case 1:
	    sprintf(&new[ni], "%03o", c);
	    ni+=3;
	    break;

	default:
	    new[ni++] = '\\';
	    new[ni++] = type[(unsigned char)c];
	}
    }
    new[ni++] = 0;

    return new;
}

/*
 * Computes and returns a false name for when the sequence has no read name.
 * This is just the record number of the sequence or it's pair (whatever is
 * lowest), optionally with a suffix.
 */
static char *false_name(GapIO *io, seq_t *s, int suffix, int *name_len) {
    static char false_name[1024];
    int p = sequence_get_pair(io, s);

    if (suffix)
	sprintf(false_name, "seq_%d%s",
		s->rec < p ? p : s->rec,
		s->rec < p ? ".f" : ".r");
    else
	sprintf(false_name, "seq_%d",
		s->rec < p ? p : s->rec);

    if (name_len)
	*name_len = strlen(false_name);

    return false_name;
}

static int export_header_sam(GapIO *io, FILE *fp,
			     int cc, contig_list_t *cv) {
    int nreads = 0, i;

    /* Inefficient as we have to loop twice - here and when outputting reads */
    for (i = 0; i < cc; i++) {
	contig_t *c;
	int len;

	if (NULL == (c = cache_search(io, GT_Contig, cv[i].contig)))
	    return -1;

	len = c->end;
	if (c->start <= 0)
	    len += 1-c->start;
	fprintf(fp, "@SQ\tSN:%s\tLN:%d\n",  c->name, len);
    }

    /* Libraries - well read-groups in SAM. */
    for (i = 0; i < io->db->Nlibraries; i++) {
	int lrec = arr(GCardinal, io->library, i);
	library_t *lib = cache_search(io, GT_Library, lrec);
	if (lib->name) {
	    fprintf(fp, "@RG\tID:%s\tLB:%s\n", lib->name, lib->name);
	} else {
	    char buf[100];
	    sprintf(buf, "rg#%d", lib->rec);
	    fprintf(fp, "@RG\tID:%s\tLB:%s\n", buf, buf);
	}
    }

    return 0;
}

static int export_contig_sam(GapIO *io, FILE *fp,
			     int crec, int start, int end,
			     int fixmates) {
    contig_iterator *ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
						  start, end,
						  GRANGE_FLAG_ISANY);
    rangec_t *r;
    contig_t *c;
    int qalloc = 0;
    char *Q = NULL, *S = NULL;
    int cg_alloc = 0;
    char *cigar = NULL;
    int offset, ustart, uend;
#if 0
    char *cons;
    int last_start;
#endif

    c = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

#if 0
    {
	int first = c->start;
	int last = c->end;
	int len = last-first+1;
	last_start = c->start;

	cons = (char *)xmalloc(len);
	calculate_consensus_simple(io, crec, first, last, cons, NULL);
    }
#endif

    /* Sam can only have coordinates from 1 onwards, so shift if needed */
    consensus_valid_range(io, crec, &ustart, &uend);
    offset = ustart < 1
	? 1-ustart
	: 0;

    while (r = contig_iter_next(io, ci)) {
	seq_t *s, *sorig;
	char *tname, *cp;
	int len, olen;
	int i, j, flag, tname_len, pos;
	int left, right, op, oplen, cglen;
	int iend, isize;
	char *mate_ref;
	library_t *lib = NULL;
	int last_lrec = -1;
	int mqual;
	char *cigar_tmp;

	if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ &&
	    (r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISUMSEQ)
	    continue;

	sorig = s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	olen = len = s->len < 0 ? -s->len : s->len;

	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	if (len > qalloc) {
	    qalloc = len;
	    Q = realloc(Q, qalloc);
	    S = realloc(S, qalloc);
	}

	/* Best guess at template name */
	tname = s->name;
	tname_len = s->name_len;

	if (tname_len == 0)
	    tname = false_name(io, s, 0, &tname_len);

	if (cp = strchr(s->name, '/'))
	    tname_len = cp-s->name;

	/* Depad sequence/qual */
	for (i = j = 0; i < len; i++) {
	    int v = '!' + s->conf[i];
	    if (v < '!') v = '!';
	    if (v > 255) v = 255;
	    Q[j] = v;

	    if (s->seq[i] == '-')
		S[j++] = 'n';
	    else if (s->seq[i] != '*')
		S[j++] = s->seq[i];
	    /* else gap */
	}
	len = j;

	flag = 0;
	
	/*
	 * Is this correct? Does "mapped in a proper pair" mean they
	 * are consistent, for whatever meaning "consistent" has in this
	 * context?
	 */
	if (r->pair_rec) flag |= 0x03;

	/* strand */
	if (sorig->len < 0)  flag |= 0x10;

	/*
	 * 1st/2nd in pair. This is not quite the same as the fwd/rev flags
	 * we store I think as we may have fwd/fwd pairs under some
	 * protocols. Is this true?
	 */
	if (r->pair_rec) {
	    if ((r->flags & GRANGE_FLAG_END_MASK) !=
		(r->flags & GRANGE_FLAG_PEND_MASK)) {
		/* Paired data with differing flags, so it's likely valid */
		if ((r->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_REV)
		    flag |= 0x40;
		else
		    flag |= 0x80;
	    }
	}

	pos = r->start;

	//puts(s->name);

	/* Generate cigar string. */
	cp = cigar;
	left  = s->left-1;
	right = s->right;
	pos += left;
	
	op = 'S'; oplen = 0;
	cglen = 0;
#if 1
	for (i = j = 0; i < olen; i++) {
	    //printf("%c %c\n", s->seq[i], cons[pos+i-c->start]);
	    if (cp - cigar + 50 > cg_alloc) {
		ptrdiff_t d = cp-cigar;
		cg_alloc += 100;
		cigar = realloc(cigar, cg_alloc);
		cp = cigar + d;
	    }

	    if (s->seq[i] == '*') {
		if (op != 'S') {
		    if (op != 'D' && oplen > 0) {
			cp += sprintf(cp, "%d%c", oplen, op);
			cglen += oplen;
			oplen = 0;
		    }
		    op = 'D';
		    oplen++;
		}
		/* else, gap in soft-clips.
		 * Right now this causes tview to fail an assertion though.
		 */
	    } else if (i < left) {
		if (op != 'S' && oplen) {
		    cp += sprintf(cp, "%d%c", oplen, op);
		    if (op != 'D') cglen += oplen;
		    oplen = 0;
		}
		op = 'S';
		oplen++;
	    } else if (i >= left && i < right) {
		if (op != 'M' && oplen) {
		    cp += sprintf(cp, "%d%c", oplen, op);
		    if (op != 'D') cglen += oplen;
		    oplen = 0;
		}
		op = 'M';
		oplen++;
	    } else {
		if (op != 'S' && oplen) {
		    /* Switching from deletion to 'S' breaks tview */
		    if (op != 'D') {
			cp += sprintf(cp, "%d%c", oplen, op);
			if (op != 'D') cglen += oplen;
		    }
		    oplen = 0;
		}
		op = 'S';
		oplen++;
	    }

	    j++;
	}
	if (oplen > 0) {
	    if (cp - cigar + 50 > cg_alloc) {
		ptrdiff_t d = cp-cigar;
		cg_alloc += 100;
		cigar = realloc(cigar, cg_alloc);
		cp = cigar + d;
	    }

	    cp += sprintf(cp, "%d%c", oplen, op);
	    if (op != 'D') cglen += oplen;
	}
#endif

#if 0
	for (i = last_start; i<pos; i++) {
	    if (cons[i-c->start] == '*')
		offset--;
	}
	last_start = pos;
	if (left) {
	    cp += sprintf(cp, "%dS", left);
	    cglen = left;
	}
	for (i = left; i < right; i++) {
	    char cb = cons[pos+i-c->start];
	    char sb = s->seq[i];

	    if (cp - cigar + 50 > cg_alloc) {
		ptrdiff_t d = cp-cigar;
		cg_alloc += 100;
		cigar = realloc(cigar, cg_alloc);
		cp = cigar + d;
	    }

	    if (cb == '*' && sb == '*') {
		continue;
	    }

	    if (cb == '*') {
		if (op != 'I') {
		    if (oplen > 0)
			cp += sprintf(cp, "%d%c", oplen, op);
		    if (op != 'D')
			cglen += oplen;
		    oplen = 0;
		    op = cglen ? 'I' : 'S';
		}
		oplen++;

	    } else if (sb == '*') {
		if (op != 'D') {
		    if (oplen > 0)
			cp += sprintf(cp, "%d%c", oplen, op);
		    cglen += oplen;
		    oplen = 0;
		    op = 'D';
		}
		oplen++;

	    } else {
		if (op != 'M') {
		    if (oplen > 0)
			cp += sprintf(cp, "%d%c", oplen, op);
		    if (op != 'D')
			cglen += oplen;
		    oplen = 0;
		    op = 'M';
		}
		oplen++;
	    }
	}
	if (oplen > 0) {
	    if (op != 'D')
		cp += sprintf(cp, "%d%c", oplen, op);
	    cglen += oplen;
	}
	if (olen-right) {
	    cp += sprintf(cp, "%dS", olen-right);
	    cglen += olen-right;
	}
	//puts(cigar);
#endif

	assert(cglen == len);

	/* Sam cannot handle gaps at the end of sequences */
	if (*(cp-1) == 'D') {
	    cp-=2;
	    while(*cp >= '0' && *cp <= '9')
		cp--;
	    cp++;
	    *cp = 0;
	}

	/* Compute insert sizes */
	if (fixmates && r->pair_rec) {
	    int other_c, other_st, other_en, other_dir, comp;
	    seq_t *pair_s;
	    range_t pair_r;

	    sequence_get_position2(io, r->pair_rec, &other_c, &other_st,
				   &other_en, &other_dir, &pair_r, &pair_s);

	    comp = (pair_s->len >= 0) ^ other_dir;
	    cache_decr(io, pair_s);
	    
	    iend = other_st < other_en ? other_st : other_en;

	    if (crec == other_c) {
		int l, r;
		mate_ref = "=";
		l = (sorig->len >= 0) ? pos : pos - sorig->len - 1;
		r = comp ? other_st-1 : other_en+1;
		isize = r-l;
	    } else {
		contig_t *oc = cache_search(io, GT_Contig, other_c);
		mate_ref = oc->name;
		isize = 0;
	    }

	    if (!comp)
		flag |= 0x20; /* strand of mate */
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

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISUMSEQ) {
	    flag |=  0x04; /* query unmapped */
	    flag &= ~0x02; /* proper pair */
	    mqual = 0;
	    cigar_tmp = "*";
	    isize = 0;
	} else {
	    mqual = r->mqual;
	    cigar_tmp = cigar;
	}

	fprintf(fp, "%.*s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%.*s\t%.*s",
		tname_len, tname,
		flag,
		c->name,
		pos+offset,
		mqual,
		cigar_tmp,
		mate_ref,
		iend,
		isize,
		len, S,
		len, Q);

	if (s->parent_type == GT_Library) {
	    if (last_lrec != s->parent_rec) {
		last_lrec  = s->parent_rec;
		lib = cache_search(io, GT_Library, s->parent_rec);
	    }
	    if (lib->name)
		fprintf(fp, "\tRG:Z:%s", lib->name);
	    else
		fprintf(fp, "\tRG:Z:rg#%d", lib->rec);
	}

	if (s->aux_len)
	    fprintf(fp, "\t%s", sam_aux_stringify(s->sam_aux, s->aux_len));

	fprintf(fp, "\n");

	if (s != sorig)
	    free(s);
    }
    contig_iter_del(ci);
    cache_decr(io, c);

    if (Q) free(Q);
    if (S) free(S);
    if (cigar) free(cigar);

    return 0;
}

/* #define FASTQ_COMPLEMENTED */
static int export_contig_fastq(GapIO *io, FILE *fp,
			       int crec, int start, int end) {
    contig_iterator *ci = contig_iter_new(io, crec, 0, CITER_FIRST,
					  start, end);
    rangec_t *r;
    int qalloc = 0;
    char *q = NULL;
    char *name;
    int name_len;
    
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	int len = s->len < 0 ? -s->len : s->len;
	int i;

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
	}
	for (i = 0; i < len; i++) {
	    int v = '!' + s->conf[i];
	    if (v < '!') v = '!';
	    if (v > 255) v = 255;
	    q[i] = v;
	}

	if (s->name_len) {
	    name_len = s->name_len;
	    name = s->name;
	} else {
	    name = false_name(io, s, 1, &name_len);
	}

#ifdef FASTQ_COMPLEMENTED
	fprintf(fp, "@%.*s %d\n%.*s\n+\n%.*s\n",
		name_len, name, s != origs, len, s->seq, len, q);

	if (s != origs)
	    free(s);
#else
	fprintf(fp, "@%.*s\n%.*s\n+\n%.*s\n",
		name_len, name, len, s->seq, len, q);
#endif
    }
    contig_iter_del(ci);

    if (q)
	free(q);

    return 0;
}

static int export_contig_fasta(GapIO *io, FILE *fp,
			       int crec, int start, int end) {
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

/*
 * A FIFO queue. We stack up sequences in here until we have sufficient
 * data to punt the entire thing out including all the annotations.
 */
typedef struct fifo {
    struct fifo *next;
    rangec_t r; /* Make this a generic payload? */
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

static void fifo_queue_push(fifo_queue_t *f, rangec_t *r) {
    fifo_t *i = malloc(sizeof(*i));

    i->r    = *r;
    i->next = NULL;

    if (!f->head)
	f->head = i;
    if (f->tail)
	f->tail->next = i;
    f->tail = i;
}

static fifo_t *fifo_queue_head(fifo_queue_t *f) {
    return f->head;
}

static fifo_t *fifo_queue_pop(fifo_queue_t *f) {
    fifo_t *i;

    i = f->head;

    if (f->head)
	f->head = f->head->next;

    if (f->head == NULL)
	f->tail = NULL;

    return i;
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
    /* Best guess at template name */
    if ((cp = strchr(name, '.')) ||
	(cp = strchr(name, '/')))
	fprintf(fp, "TN=%.*s\n", (int)(cp-name), name);
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
}

static int export_contig_baf(GapIO *io, FILE *fp,
			     int crec, int start, int end) {
    contig_iterator *ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
						  start, end,
						  GRANGE_FLAG_ISANY);
    rangec_t *r;
    int qalloc = 0;
    char *q = NULL, *S = NULL, *cp;
    contig_t *c;
    char *name;
    int name_len;
    fifo_queue_t *fq = fifo_queue_create(), *tq = fifo_queue_create();
    fifo_t *fi;
    int last_start = 0;

    /* Contig record */
    c = (contig_t *)cache_search(io, GT_Contig, crec);
    fprintf(fp, "CO=%s\nLN=%d\n\n", c->name, c->end-c->start+1);

    /* Seq/tag records */
    while (r = contig_iter_next(io, ci)) {
	/* Add new items to fifo */
	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    fifo_queue_push(fq, r);
	    last_start = r->start;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    if (r->pair_rec == crec) {
		/* Contig tag */
		anno_ele_t *a = cache_search(io, GT_AnnoEle, r->rec);
		char type[5];

		fprintf(fp, "AN=%s\nAT=C\n", type2str(r->mqual, type));
		fprintf(fp, "LO=@%d\n", r->start);
		if (r->start != r->end)
		    fprintf(fp, "LL=%d\n", r->end - r->start + 1);
		fprintf(fp, "AN=%s\n", type2str(r->mqual, type));

		if (a->comment && *a->comment) {
		    char *escaped = escape_C_string(a->comment);
		    fprintf(fp, "TX=%s\n", escaped);
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
	    sprintf(name+strlen(name), ".%d", s->rec);
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
	    s->trace_name_len ? s->trace_name_len : strlen(name),
	    s->trace_name_len ? s->trace_name     : name);

    //fprintf(fp, "Insert_size %d %d\n", fixme, fixme);

    /* Best guess at template name */
    if ((cp = strchr(name, '.')) ||
	(cp = strchr(name, '/')))
	fprintf(fp, "Template %.*s\n", (int)(cp-name), name);

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
	    fputc(s->seq[i+j] != '*' ? s->seq[i+j] : '-', fp);
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
}

static int export_contig_caf(GapIO *io, FILE *fp,
			     int crec, int start, int end) {
    contig_iterator *ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
						  start, end,
						  GRANGE_FLAG_ISANY);
    rangec_t *r;
    int qalloc = 0;
    char *q = NULL, *S = NULL, *cp;
    contig_t *c;
    char *name;
    int name_len;
    fifo_queue_t *fq = fifo_queue_create(), *tq = fifo_queue_create();
    fifo_t *fi;
    int last_start = 0;
    dstring_t *ds = dstring_create(NULL);
    int first_base, last_base, len;
    char *cons;
    float *qual;
    int i, wrap;

    /* Seq/tag records */
    while (r = contig_iter_next(io, ci)) {
	/* Add new items to fifo */
	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    fifo_queue_push(fq, r);
	    last_start = r->start;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    if (r->pair_rec == crec) {
		/* Contig tag */
		/* FIXME: to do still */
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
    fprintf(fp, "\n");

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
	int q = qual[i]+.5;
	wrap += fprintf(fp, "%d ", q);
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
			     int crec, int start, int end) {
    contig_iterator *ci;
    rangec_t *r;
    int qalloc = 0;
    char *q = NULL, *S = NULL;
    contig_t *c;
    int nreads = 0;
    int i, j, len, nBS, last, first_base, last_base;
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
    len = last_base - first_base + 1;

    fprintf(fp, "\nCO %s %d %d %d U\n",
	    c->name, len, nreads, nBS);


    /* Contig sequence */
    cons = (char *)xmalloc(len);
    qual = (float *)xmalloc(len * sizeof(*qual));
    calculate_consensus_simple(io, crec, first_base, last_base, cons, qual);
    for (i = 0; i < len; i += 50) {
	fprintf(fp, "%.*s\n", len-i > 50 ? 50 : len-i, &cons[i]);
    }

    /* Contig Quality */
    fprintf(fp, "\nBQ\n");
    for (i = j = 0; i < len; i++) {
	int q;
	if (cons[i] == '*')
	    continue;

	q = qual[i]+0.5;
	if (q < 0) q = 0;
	if (q > 97) q = 97;
	fprintf(fp, " %d", q);

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
	int i;

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

	fprintf(fp, "DS CHROMAT_FILE: %s PHD_FILE: %s.phd.1 TIME: Thu Jan 01 00:00:01 GMT 1970\n", s->name, s->name);

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
			  char *fn, int fixmates) {
    int i;
    FILE *fp;
    
    if (NULL == (fp = fopen(fn, "w"))) {
	perror(fn);
	return -1;
    }

    /* Header */
    switch (format) {
    case FORMAT_ACE:
	export_header_ace(io, fp, cc, cv);
	break;

    case FORMAT_SAM:
	export_header_sam(io, fp, cc, cv);
	break;
    }

    /* Per contig */
    for (i = 0; i < cc; i++) {
	switch (format) {
	case FORMAT_SAM:
	    export_contig_sam(io, fp, cv[i].contig, cv[i].start, cv[i].end,
			      fixmates);
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

    fclose(fp);

    return 0;
}

/* ------------------------------------------------------------------------
 * export_tags implementation
 */

static int export_header_tags_gff(GapIO *io, FILE *fp,
				  int cc, contig_list_t *cv) {
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
			   int crec, int start, int end,
			   int consensus, int unpadded) {
    contig_iterator *ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST,
						  start, end,
						  GRANGE_FLAG_ISANNO);
    rangec_t *r;
    contig_t *c;
    char *con = NULL;
    int *map = NULL;

    c = cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    /* Generate a padded to unpadded mapping table */
    if (unpadded) {
	int i, np;
	if (NULL == (con = malloc(c->end - c->start + 2))) {
	    return -1;
	}
	if (NULL == (map = malloc((c->end - c->start + 2) *
				  sizeof(int)))) {
	    free(con);
	    return -1;
	}
	calculate_consensus_simple(io, crec, c->start, c->end, con, NULL);

	for (np = 0, i = c->start; i <= c->end; i++) {
	    map[i-c->start] = i - np;
	    if (con[i-c->start] == '*')
		np++;
	}
    }

    /* Export */
    while (r = contig_iter_next(io, ci)) {
	anno_ele_t *a = cache_search(io, GT_AnnoEle, r->rec);
	char type[5], type2[5];
	char *name, *escaped;
	int st, en;

	cache_incr(io, a);

	if (!consensus && (r->flags & GRANGE_FLAG_TAG_SEQ)) {
	    int seq_start, seq_end;
	    seq_t *s;

	    sequence_get_position(io, r->pair_rec, NULL,
				  &seq_start, &seq_end, NULL);

	    st = r->start - seq_start;
	    en = r->end - seq_start;

	    s = cache_search(io, GT_Seq, r->pair_rec);
	    name = s->name;
	} else {
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

	if (a->comment && *a->comment) {
	    escaped = escape_C_string(a->comment);

	    fprintf(fp, "%s\tgap5\t%s\t%d\t%d\t.\t.\t.\tAnno \"%s\"\n",
		    name, type2str(r->mqual, type2),
		    st, en,
		    escaped);
	    free(escaped);
	} else {
	    fprintf(fp, "%s\tgap5\t%s\t%d\t%d\t.\t.\t.\n",
		    name, type2str(r->mqual, type2),
		    st, en);
	}

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
