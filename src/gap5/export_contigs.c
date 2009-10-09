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

#define FORMAT_FASTA 0
#define FORMAT_FASTQ 1
#define FORMAT_SAM   2
#define FORMAT_BAF   3
#define FORMAT_ACE   4

static int export_contigs(GapIO *io, int cc, contig_list_t *cv, int format,
			  char *fn);

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
	{"-format",	ARG_STR, 1, "fastq",  offsetof(ec_arg, format)},
	{"-outfile",    ARG_STR, 1, "out.aln",offsetof(ec_arg, outfile)},
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

    c = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

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
		pos,
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

#ifdef FASTQ_COMPLEMENTED
	fprintf(fp, "@%.*s %d\n%.*s\n+\n%.*s\n",
		s->name_len, s->name, s != origs, len, s->seq, len, q);

	if (s != origs)
	    free(s);
#else
	fprintf(fp, "@%.*s\n%.*s\n+\n%.*s\n",
		s->name_len, s->name, len, s->seq, len, q);
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
    
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	int len = s->len < 0 ? -s->len : s->len;

	fprintf(fp, ">%.*s\n%.*s\n", s->name_len, s->name, len, s->seq);
    }
    contig_iter_del(ci);

    return 0;
}

static int export_contig_baf(GapIO *io, FILE *fp,
			     int crec, int start, int end) {
    contig_iterator *ci = contig_iter_new(io, crec, 0, CITER_FIRST,
					  start, end);
    rangec_t *r;
    int qalloc = 0;
    char *q = NULL, *S = NULL, *cp;
    contig_t *c;

    /* Contig record */
    c = (contig_t *)cache_search(io, GT_Contig, crec);
    fprintf(fp, "CO=%s\nLN=%d\n\n", c->name, c->end-c->start+1);

    /* Seq records */
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	int len = s->len < 0 ? -s->len : s->len;
	int i;

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

	fprintf(fp, "RD=%.*s\n", s->name_len, s->name);
	fprintf(fp, "DR=%d\n", s->len >= 0 ? 1 : -1);
	fprintf(fp,"AP=%d\n", 
		s->len >= 0
		? r->start + s->left-1
		: r->start + (-s->len - s->right));
	fprintf(fp, "QL=%d\n", s->left);
	fprintf(fp, "QR=%d\n", s->right);
	fprintf(fp, "PR=%d\n",
		(s->flags & SEQ_END_MASK) == SEQ_END_FWD ? 0 : 1);
	/* Best guess at template name */
	if ((cp = strchr(s->name, '.')) ||
	    (cp = strchr(s->name, '/')))
	    fprintf(fp, "TN=%.*s\n", (int)(cp-s->name), s->name);
	if (s->trace_name_len)
	    fprintf(fp,"TR=%.*s\n", s->trace_name_len, s->trace_name);
	if (s->alignment_len)
	    fprintf(fp, "AL=%.*s\n", s->alignment_len, s->alignment);
	if (r->mqual != 255)
	    fprintf(fp, "MQ=%d\n", r->mqual);
	fprintf(fp, "SQ=%.*s\n", len, S);
	fprintf(fp, "FQ=%.*s\n\n", len, q);
    }
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
	}
    }

    fclose(fp);

    return 0;
}
