#ifndef _HAPLO_H
#define _HAPLO_H

#include "IO.h"
#include "dstring.h"
#include "haplo.h"

typedef struct {
    int tmplate;
    char base;
    int1 conf;
    double tscore;
} seq_base_t;

typedef struct {
    int pos;
    double score;
    seq_base_t *seqs;
    int nseqs;
} snp_t;

dstring_t *haplo_snps(GapIO *io, int contig, int start, int end,
		      int discrep_cutoff, int snp_cutoff,
		      int min_base_qual, int two_alleles);

void free_snps(snp_t *snps, int nsnps);

dstring_t *haplo_split(GapIO *io, snp_t *snp, int nsnps, int verbosity,
		       double min_score, int two_pass, int fast_mode,
		       double correlation_offset, int max_sets);

int calc_template_consensus(GapIO *io, int contig,
			    int start, int end,
			    int *templates, int ntemplates,
			    char **cons, float **qual);

int calc_template_depth(GapIO *io, int contig, int start, int end,
			int *tdepth);

#endif
