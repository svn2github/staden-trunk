#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "gap_globals.h"
#include "primlib.h"

/*
 * Copies values from primlib_args into the primlib_state structure.
 * Only values in primlib_args that are non zero are copied over,
 * with everything else being left as the defaults.
 *
 * All arguments are double, but some will be converted to int or short
 * during the copy. If you wish to set a parameter to be a zero value pass
 * it over as 0.0001 (to avoid the non-zero check above).
 */
void primlib_set_args(primlib_state *state, primlib_args *args) {
    if (!state || !args)
	return;

    if (args->min_tm)
	state->p3args.min_tm = args->min_tm;
    if (args->max_tm)
	state->p3args.max_tm = args->max_tm;
    if (args->opt_tm)
	state->p3args.opt_tm = args->opt_tm;

    if (args->min_gc)
	state->p3args.min_gc = args->min_gc;
    if (args->max_gc)
	state->p3args.max_gc = args->max_gc;
    if (args->opt_gc)
	state->p3args.opt_gc_content = args->opt_gc;

    if (args->min_len)
	state->p3args.primer_min_size = (int)args->min_len;
    if (args->max_len)
	state->p3args.primer_max_size = (int)args->max_len;
    if (args->opt_len)
	state->p3args.primer_opt_size = (int)args->opt_len;

    if (args->max_end_stability)
	state->p3args.max_end_stability = args->max_end_stability;

    if (args->salt_conc)
	state->p3args.salt_conc = args->salt_conc;
    if (args->dna_conc)
	state->p3args.dna_conc = args->dna_conc;

    if (args->self_any)
	state->p3args.self_any = (short)args->self_any * 100;
    if (args->self_end)
	state->p3args.self_end = (short)args->self_end * 100;

    if (args->gc_clamp)
	state->p3args.gc_clamp = args->gc_clamp;

    if (args->max_poly_x)
	state->p3args.max_poly_x = args->max_poly_x;

    if (args->num_return)
	state->p3args.num_return = args->num_return;
}

primlib_state *primlib_create(void)
{
    primlib_state *state = malloc(sizeof(*state));
    if (!state)
	return NULL;

    memset(state, 0, sizeof(*state));

    set_default_global_primer_args(&state->p3args);
    /*state->p3args.primer_task = pick_hyb_probe_only;*/
    state->p3args.primer_task = pick_left_only;
    state->p3args.liberal_base = 1; /* converts - to N */
    state->p3state = primer3_create();

    return state;
}

int
primlib_choose(primlib_state *state, char *seq)
{
    seq_args sa;

    if (!state)
	return -1;

    /* Initialise seq_args structure */
    memset(&sa, 0, sizeof(sa));
    sa.start_codon_pos = PR_DEFAULT_START_CODON_POS;
    sa.sequence = seq;
    sa.incl_l = strlen(seq);
    sa.incl_s = state->p3args.first_base_index;

    memset(&state->p3args.glob_err, 0, sizeof(state->p3args.glob_err));
    if (0 != primer3_choose(state->p3state, &state->p3args, &sa)) {
	if (sa.error.data || state->p3args.glob_err.data) {
	    printf("primer3 failed: ");
	    if (sa.error.data)
		printf("'%s' ", sa.error.data);	
	    if (state->p3args.glob_err.data)
		printf("'%s'", state->p3args.glob_err.data);
	    printf("\n");
	}

	state->nprimers = 0;
	return -1;
    }

    /* Convert primer3 results to primlib data structures */
    state->nprimers = state->p3state->n_f;
    state->primers = state->p3state->f;

    return 0;
}

int
primlib_choose_pcr(primlib_state *state, char *seq, int target, int tlen)
{
    seq_args sa;

    if (!state)
	return -1;

    /* Initialise seq_args structure */
    memset(&sa, 0, sizeof(sa));
    sa.start_codon_pos = PR_DEFAULT_START_CODON_POS;
    sa.sequence = seq;
    sa.incl_l = strlen(seq);
    sa.incl_s = state->p3args.first_base_index;

    sa.tar[0][0] = target;
    sa.tar[0][1] = tlen;
    sa.num_targets = 1;

    memset(&state->p3args.glob_err, 0, sizeof(state->p3args.glob_err));
    if (0 != primer3_choose(state->p3state, &state->p3args, &sa)) {
	if (sa.error.data || state->p3args.glob_err.data) {
	    printf("primer3 failed: ");
	    if (sa.error.data)
		printf("'%s' ", sa.error.data);	
	    if (state->p3args.glob_err.data)
		printf("'%s'", state->p3args.glob_err.data);
	    printf("\n");
	}

	state->nprimers = 0;
	return -1;
    }

    /* Convert primer3 results to primlib data structures */
    state->nprimers = state->p3state->n_f;
    state->primers = state->p3state->f;
    state->npairs = state->p3state->best_pairs.num_pairs;
    state->pairs = state->p3state->best_pairs.pairs;

    return 0;
}

void primlib_destroy(primlib_state *state)
{
    if (!state)
	return;

    primer3_destroy(state->p3state);
    free(state);
}

static void primlib_set_arg_by_str(primlib_args *args,
				   char *name, int name_len,
				   char *value, int value_len) {
    char tmpbuf[256];

    /* nul-terminate value by copying into a tmp buffer */
    strncpy(tmpbuf, value, value_len < 255 ? value_len : 255);
    tmpbuf[value_len < 255 ? value_len : 255] = 0;

    if (strncmp(name, "min_tm", name_len) == 0) {
	args->min_tm = atof(tmpbuf);
    } else if (strncmp(name, "max_tm", name_len) == 0) {
	args->max_tm = atof(tmpbuf);
    } else if (strncmp(name, "opt_tm", name_len) == 0) {
	args->opt_tm = atof(tmpbuf);
    } else if (strncmp(name, "min_gc", name_len) == 0) {
	args->min_gc = atof(tmpbuf);
    } else if (strncmp(name, "max_gc", name_len) == 0) {
	args->max_gc = atof(tmpbuf);
    } else if (strncmp(name, "opt_gc", name_len) == 0) {
	args->opt_gc = atof(tmpbuf);
    } else if (strncmp(name, "min_len", name_len) == 0) {
	args->min_len = atof(tmpbuf);
    } else if (strncmp(name, "max_len", name_len) == 0) {
	args->max_len = atof(tmpbuf);
    } else if (strncmp(name, "opt_len", name_len) == 0) {
	args->opt_len = atof(tmpbuf);
    } else if (strncmp(name, "max_end_stability", name_len) == 0) {
	args->max_end_stability = atof(tmpbuf);
    } else if (strncmp(name, "salt_conc", name_len) == 0) {
	args->salt_conc = atof(tmpbuf);
    } else if (strncmp(name, "self_any", name_len) == 0) {
	args->self_any = atof(tmpbuf);
    } else if (strncmp(name, "self_end", name_len) == 0) {
	args->self_end = atof(tmpbuf);
    } else if (strncmp(name, "gc_clamp", name_len) == 0) {
	args->gc_clamp = atoi(tmpbuf);
    } else if (strncmp(name, "max_poly_x", name_len) == 0) {
	args->max_poly_x = atoi(tmpbuf);
    } else if (strncmp(name, "num_return", name_len) == 0) {
	args->num_return = atof(tmpbuf);
    } else {
	fprintf(stderr, "Unknown keyword '%.*s'\n",
		name_len, name);
    }
}

/*
 * Converts from 'string value string value ...' syntax to a primlib_args
 * structure suitable for passing into primlib_set_args. The various
 * string keywords follow the same names as the variables used within the
 * primlib_args structure.
 *
 * Returns:
 *    A pointer to a malloc()ed primlib_args structure. The caller should
 *    free this.
 *    NULL on failure.
 */
primlib_args *primlib_str2args(char *str) {
    enum {name_start, name_end, value_start, value_end} state;
    char *cp, *name = NULL, *value = NULL;
    int name_len = 0, value_len = 0;
    primlib_args *args;

    if (NULL == (args = calloc(1, sizeof(*args))))
	return NULL;
    
    state = name_start;
    cp = str;
    do {
	switch (state) {
	case name_start:
	    if (!isspace(*cp)) {
		state = name_end;
		name = cp;
	    }
	    break;

	case name_end:
	    if (isspace(*cp)) {
		name_len = cp-name;
		state = value_start;
	    }
	    break;

	case value_start:
	    if (!isspace(*cp)) {
		state = value_end;
		value = cp;
	    }
	    break;

	case value_end:
	    if (isspace(*cp) || *cp == 0) {
		value_len = cp-value;
		primlib_set_arg_by_str(args, name, name_len, value, value_len);
		state = name_start;
	    }
	    break;
	}
    } while (*cp++);

    return args;
}
