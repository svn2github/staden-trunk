#ifndef _PRIMLIB_H
#define _PRIMLIB_H

#include "primer3.h"

typedef struct {
    double min_tm;	/* Minimum, maximum and optimum temperature */
    double max_tm;
    double opt_tm;
    double min_gc;	/* GC Content - as a percentage */
    double max_gc;
    double opt_gc;
    double min_len;	/* Length, in bases */
    double max_len;
    double opt_len;
    double max_end_stability;
    double salt_conc;
    double dna_conc;
    double self_any;
    double self_end;
    double gc_clamp;
    double max_poly_x;
} primlib_args;

typedef struct {
    primer_args p3args;
    primer_state *p3state;

    /* public */
    int nprimers;
    primer_rec *primers;
} primlib_state;

void primlib_set_args(primlib_state *state, primlib_args *args);

primlib_state *primlib_create(void);

int primlib_choose(primlib_state *state, char *seq);

void primlib_destroy(primlib_state *state);

primlib_args *primlib_str2args(char *str);

#endif /* _PRIMLIB_H */
