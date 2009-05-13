#include <tcl.h>

#include <tg_gio.h>
#include "export_contigs.h"
#include "gap_cli_arg.h"
#include "newgap_structs.h"
#include "list_proc.h"
#include "xalloc.h"
#include "dna_utils.h"

#define FORMAT_FASTA 0
#define FORMAT_FASTQ 1
#define FORMAT_SAM   2
#define FORMAT_BAF   3

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

    c = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, c);
    
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	char *cp;
	int len = s->len < 0 ? -s->len : s->len;
	int i, flag;

	if (len > qalloc) {
	    qalloc = len;
	    Q = realloc(Q, qalloc);
	    S = realloc(S, qalloc);
	}

	if (NULL == (cp = strchr(s->name, '/')))
	    cp = s->name + s->name_len;

	for (i = 0; i < len; i++) {
	    int v = '!' + s->conf[i];
	    if (v < '!') v = '!';
	    if (v > 255) v = 255;
	    Q[i] = v;

	    if (s->seq[i] == '-')
		S[i] = 'n';
	    else if (s->seq[i] == '*')
		S[i] = '-';
	    else
		S[i] = s->seq[i];
	}

	if (s->len < 0) {
	    complement_seq(S, len);
	    reverse_dna(Q, len);
	}

	flag = 0x02;
	if (r->pair_rec) flag |= 0x01;
	if (s->len < 0)  flag |= 0x10;
	if ((r->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_REV)
	    flag |= 0x40;
	else
	    flag |= 0x80;

	fprintf(fp, "%.*s\t%d\t%s\t%d\t%d\t%dM\t*\t0\t0\t%.*s\t%.*s\n",
		(int)(cp - s->name), s->name,
		flag,
		c->name,
		r->start,
		r->mqual,
		len,
		len, S,
		len, Q);
    }
    contig_iter_del(ci);
    cache_decr(io, c);

    if (Q) free(Q);
    if (S) free(S);

    return 0;
}

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

	fprintf(fp, "@%.*s\n%.*s\n+\n%.*s\n",
		s->name_len, s->name, len, s->seq, len, q);
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
	if (cp = strchr(s->name, '.'))
	    fprintf(fp, "TN=%.*s\n", (int)(cp-s->name), s->name);
	fprintf(fp,"TR=%.*s\n",
		s->trace_name_len ? s->trace_name_len : s->name_len,
		s->trace_name_len ? s->trace_name     : s->name);
	if (s->alignment_len)
	    fprintf(fp, "AL=%.*s\n", s->alignment_len, s->alignment);
	fprintf(fp, "SQ=%.*s\n", len, S);
	fprintf(fp, "FQ=%.*s\n\n", len, q);
    }
    contig_iter_del(ci);

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
	}
    }

    fclose(fp);

    return 0;
}
