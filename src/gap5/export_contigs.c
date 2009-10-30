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

/* Sequence formats */
#define FORMAT_FASTA 0
#define FORMAT_FASTQ 1
#define FORMAT_SAM   2
#define FORMAT_BAF   3
#define FORMAT_ACE   4

/* Annotation formats */
#define FORMAT_GFF   5


static int export_contigs(GapIO *io, int cc, contig_list_t *cv, int format,
			  char *fn);
static int export_tags(GapIO *io, int cc, contig_list_t *cv, int format,
		       int consensus, int unpadded, char *fn);

/* ------------------------------------------------------------------------
 * Tcl interfaces
 */
typedef struct {
    GapIO *io;
    char *inlist;
    char *format;
    char *outfile;
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
    else
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    res = export_contigs(args.io, rargc, rargv, format_code, args.outfile);

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

    return 0;
}

static int export_contig_sam(GapIO *io, FILE *fp,
			     int crec, int start, int end) {
    contig_iterator *ci = contig_iter_new(io, crec, 0, CITER_FIRST,
					  start, end);
    rangec_t *r;
    contig_t *c;
    int qalloc = 0;
    char *Q = NULL, *S = NULL;
    int cg_alloc = 0;
    char *cigar = NULL;
    int insert_sizes = 0; /* set to 1 to also set MPOS/ISIZE */
    int offset;

    c = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    /* Sam can only have coordinates from 1 onwards, so shift if needed */
    offset = c->start <= 0
	? 1-c->start
	: 0;

    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	seq_t *sorig = s;
	char *tname, *cp;
	int len = s->len < 0 ? -s->len : s->len, olen = len;
	int i, j, flag, tname_len, pos;
	int left, right, op, oplen, cglen;
	int iend, isize;
	char *mate_ref;

	if (s->len < 0) {
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
	if (r->pair_rec) flag |= 0x03;
	if (sorig->len < 0)  flag |= 0x10;
	if (r->pair_rec) {
	    if ((r->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_REV)
		flag |= 0x40;
	    else
		flag |= 0x80;
	}

	pos = r->start;

	/* Generate cigar string. */
	cp = cigar;
	left  = s->left-1;
	right = s->right;
	pos += left;
	
	op = 'S'; oplen = 0;
	cglen = 0;
	for (i = j = 0; i < olen; i++) {
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

	assert(cglen == len);

	/* Compute insert sizes */
	if (insert_sizes && r->pair_rec) {
	    int other_c, other_st, other_en, other_dir;

	    sequence_get_position(io, r->pair_rec, &other_c, &other_st,
				  &other_en, &other_dir);
	    
	    iend = other_st < other_en ? other_st : other_en;

	    if (crec == other_c) {
		int l, r;
		mate_ref = "=";
		l = (sorig->len >= 0) ? pos : pos - sorig->len - 1;
		r = other_dir ? other_st : other_en;
		isize = r-l+1;
	    } else {
		contig_t *oc = cache_search(io, GT_Contig, other_c);
		mate_ref = oc->name;
		isize = 0;
	    }
	} else {
	    iend = 0;
	    isize = 0;
	    mate_ref = r->pair_rec ? "=" : "*";
	}

	fprintf(fp, "%.*s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%.*s\t%.*s\n",
		tname_len, tname,
		flag,
		c->name,
		pos+offset,
		r->mqual,
		cigar,
		mate_ref,
		iend,
		isize,
		len, S,
		len, Q);

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
	if (s->len < 0) {
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

    printf("POP SEQ %d @ %d,%d\n", fi->r.rec, fi->r.start, fi->r.end);

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
    fprintf(fp, "DR=%d\n", s->len >= 0 ? 1 : -1);
    fprintf(fp,"AP=%d\n", 
	    s->len >= 0
	    ? fi->r.start + s->left-1
	    : fi->r.start + (-s->len - s->right));
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

	printf("pop TAG %d @ %d,%d\n", ti->r.rec, ti->r.start, ti->r.end);

	/* Tag is for this seq. */
	a = cache_search(io, GT_AnnoEle, ti->r.rec);
	fprintf(fp, "AN=%s\n", type2str(ti->r.mqual, type));
	if (s->len >= 0) {
	    fprintf(fp, "LO=%d\n", ti->r.start - (fi->r.start-1));
	} else {
	    fprintf(fp, "LO=%d\n", -s->len+1 - (ti->r.end - (fi->r.start-1)));
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
	    printf("push SEQ %d @ %d,%d\n", r->rec, r->start, r->end);
	    last_start = r->start;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    /* Technically this should be a list and not a fifo */
	    fifo_queue_push(tq, r);
	    printf("push TAG %d @ %d,%d\n", r->rec, r->start, r->end);
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

    ci = contig_iter_new(io, crec, 0, CITER_FIRST, c->start, c->end);
    last = 1;
    nBS = 0;
    first_base = 999999;
    last_base = -999999;
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	nreads++;
	if (r->start + s->right-1 > last) {
	    nBS++;
	    last = r->start + s->right-1;
	}
	if (first_base > r->start + s->left-1)
	    first_base = r->start + s->left-1;
	if (last_base < r->start + s->right-1)
	    last_base = r->start + s->right-1;
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
		"UC"[s->len < 0], r->start - first_base);
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

	if (s->len < 0) {
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
			  char *fn) {
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
	    export_contig_sam(io, fp, cv[i].contig, cv[i].start, cv[i].end);
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
