#ifndef _FINISH_H
#define _FINISH_H

/*
 * We have a choice of using either OSP or Primer3. Defaults to Primer3 unless
 * USE_OSP is defined.
 */
/* #define USE_OSP */

#include <tcl.h>
#include "IO.h"
#include "vseqs.h"
#include "template.h"
#include "finish_hash.h"
#include "primlib.h"

#define CLASS_STRAND_TOP	1
#define CLASS_STRAND_BOTTOM	2
#define CLASS_SEQ_DEPTH_GT	3
#define CLASS_TEMP_DEPTH_GT	4
#define CLASS_CONFIDENCE_GT	5
#define CLASS_CHEMISTRY		6
#define CLASS_CONTIG_LEFT_END	7
#define CLASS_CONTIG_RIGHT_END	8
#define CLASS_LOW_COMPLEXITY	9
#define CLASS_SEQ_DEPTH_GE	10
#define CLASS_TEMP_DEPTH_GE	11
#define CLASS_CONFIDENCE_GE	12
#define CLASS_POLY_A		13
#define CLASS_POLY_C		14
#define CLASS_POLY_G		15
#define CLASS_POLY_T		16
#define CLASS_POLY_K		17
#define CLASS_POLY_M		18
#define CLASS_POLY_R		19
#define CLASS_POLY_S		20
#define CLASS_POLY_W		21
#define CLASS_POLY_Y		22


typedef struct {
    int bit;		/* Bit number */
    int type;		/* CLASS_* defines above */
    int arg;		/* Optional paremeter */
} con_bits_t;

typedef struct {
    /* General purpose options */
    int    use_avg_insert;	/* Average insert size vs minimum size? */
    double mandatory_ratio;	/* If >ratio mandatory problems left, reject */
    int    prob_mandatory;	/* Mask for probs =>must solve these */
    double max_score_drop;	/* If score/max_score < max_drop then reject */
    double min_template_score;	/* Minimum score for template */
    double min_score;		/* Absolute minimum score for experiment */
    int	   find_dup_templates;	/* Try to identify duplicate templates */
    int    dust_level;		/* Dust low-complex filter threshold */
    int    min_extension;	/* Minimum contig extension to be worthwhile */
    int    svec_as_cvec;	/* Treat SVEC as clone ends; useful for ESTs */

    /* Options for resequencing */
    int	   reseq_length;	/* Expected length of normal resequencing */
    int    reseq_nsolutions;	/* Number of re-sequences to pick */

    /* Options for long gel readings */
    int    long_length;		/* Expected length of long reading */
    int    long_nsolutions;	/* Number of long reads to pick */

    /* Options for primer walking */
    int    pwalk_search_dist;	/* Extra dist +/- temp for primer mismatches */
    double pwalk_max_match;	/* Uniqueness %age. Above this => reject */
    int    pwalk_osp_score;	/* Primers must have score <= this figure */
    int    pwalk_noligos;	/* Number of best OSP primers to explore */
    int    pwalk_ntemplates;	/* Number of templates to use for each oligo */
    int    pwalk_offset1;	/* Location offset from problem base to find */
    int    pwalk_offset2;	/*   primers within. offset1 > offset2 */
    int    pwalk_length;	/* Expected reading length */
    int    pwalk_nsolutions;	/* Number of primers to use */
    int    pwalk_seq_gap;	/* Gap between primer and base calling */
    int    pwalk_consistent_only; /* Only use known consistent templates */
    double pwalk_max_err;	/* Maximum prob of error in oligo sequence */
    int    pwalk_min_qual;	/* Min consensus qual allowed in a primer */
    double pwalk_max_err2;	/* As max_err/min_qual, but weaker */
    int    pwalk_min_qual2;	/*    constraints for contig ends */
    int    pwalk_end_dist;	/* How close to end before max_err2 is used? */
    int    pwalk_use_template;	/* Max number of times to use a single temp. */
    double pwalk_use_template_score; /* Cost factor for over use, 0=>reject */
    double pwalk_dup_template_cost; /* Cost multiplier for dup templates */
    int	   pwalk_tag_cons;	/* Tags on consensus(1) or reads(0) */
    int	   pwalk_prob_mask;	/* If prob[i] & mask then don't use in oligo */
    char  *pwalk_tag_type;	/* Type of tag to create; defaults to PRIM */
    primlib_args primer_args;	/* Temperature, GC %, length, conc, etc */
    
    /* Debug options */
    int    debug[10];		/* Debug levels. */
} finish_options_t;

/* Only used temporarily during command line argument parsing */
typedef struct {
    char *contig;
    char *ccontigs;
    char *eseq;
    char *skip_template_file;
    char *avail_template_file;
    char *external_seq_file;
    char *pscores;
    char *mscores;
    char *output_file;
} finish_args_t;

typedef struct {
    finish_options_t   opts;		/* Options */
    finish_args_t      args;		/* Used in CLI parsing */
    GapIO	      *io;		/* GapIO handle */
    int		       contig;		/* Contig number */
    int		       cvec_left;	/* Cloning vector at left end */
    int		       cvec_right;	/* Cloning vector at right end */
    int		       start;		/* Start pos in contig (1..N) */
    int		       end;		/* End pos in contig (1..N) */
    int		       length;		/* Length (end-start+1) */
    vcontig_t	      *vc;		/* Virtual contig structure */
    char 	      *cons;		/* Masked Consensus */
    char	      *filtered;	/* cons with low-complexity removed */
    float	      *qual;		/* Quality */
    float	      *orig_qual;	/* Quality (before virtual seqs) */
    int		       alloc_len;	/* Allocated length of cons/qual */
    con_bits_t	      *classify;	/* Classification mechanism */
    int		       nclassify;	/* Number of items in 'classify' */
    unsigned int      *base_bits;	/* Base classification bits */
    unsigned int      *prob_bits;	/* Problem bits, during run */
    unsigned int      *orig_prob_bits;	/* Problem bits, at start of run */
    unsigned int      *solution_bits;	/* Solution bits */
    unsigned int      *tag_mask;	/* Mask this problem (skipped tag) */
    template_c	     **tarr;		/* Template info array (template.c) */
    int		      *template_dup;	/* Holds duplicate templates */
    char	      *prob_script;	/* Tcl function to find problems */
    char	      *solu_script;	/* Tcl function to find solutions */
    int                left_extent;	/* Max left posn in vc, <= 0 */
    int                right_extent;	/* Max right posn in vc, >= len-1 */
    int		      *template_used;	/* How many times we use each temp. */ 
    int		      *template_skip;	/* Skip these templates */
    int		       fake_searched;	/* Have we looked for fake templates?*/
    int		       count[10];	/* Count of each experiment type */
    float	       cost[10];	/* Cost of each experiment type */
    float 	       pscore[32];	/* Score for each problem type */
    float 	       mscore[32];	/* Score added for each mand. prob. */
    Tcl_DString        tag_list;	/* Tcl list of tags to add to db */
    char 	     **skip_tags;	/* Tag types to mask out solutions */
    int		       nskip_tags;	/* Number of elements in skip_tags */
    char	      *external_seq;	/* External sequence (eg vector) */
    int		       external_seq_len;/* Length of external_seq */
    Hash	      *external_seq_h;	/* Hash of external_seq */
    char	      *external_seq_rev;/* Reverse complement of external_seq*/
    char	      *all_cons;	/* Complete consensus sequence (cat) */
    int	      	       all_cons_len;	/* Length of all_cons */
    Hash	      *all_cons_h;	/* Hash of all_cons */

    Tcl_Command	       command_token;	/* From Tcl_CreateObjCommand() */

    FILE	      *out_fp;		/* Where to send the experiment info */
} finish_t;

typedef struct {
    char primer[100];	/* DNA sequence, complemented if appropriate */
    double score;	/* Combined primer score */
    double osp_score;	/* Score from OSP */
    double p_err;	/* Chance of sequencing error in primer */
    double secbind;	/* Secondary binding score (high is bad) */
    int start;		/* Left end of primer in consensus */
    int end;		/* Right end of primer in consensus */
    int dir;		/* Direction of primer: 1 == fwd, 2 == rev */
    int homopolymer;	/* Max length of homopolymer found in primer */
} experiment_walk_t;

typedef struct {
    int dummy;
} experiment_long_t;

typedef struct experiments {
    GReadings r;
    double score;	/* Experiment score */
    double cost;
    int expt_id;
    int group_id;
    int group_num;
    int type;
    int nsolutions;	/* How many solutions of this type we desire */
    double t_score;	/* Template score */
    int t_dir;		/* Template direction */
    void (*log_func)(FILE *fp, finish_t *f, struct experiments *e, int contig,
		     Tcl_DString *tagds);
    union {
	experiment_walk_t e_walk;
	experiment_long_t e_long;
    } data;
} experiments_t;

#define EXPERIMENT_LONG 1
#define EXPERIMENT_VPWALK 2
#define EXPERIMENT_RESEQ 3
#define EXPERIMENT_CPWALK 4
#define EXPERIMENT_REVERSE 4

#define FIN_DEBUG_REGEXP 5
#define FIN_DEBUG_DUST 6
#define FIN_DEBUG_SCORE 7
#define FIN_DEBUG 0 /* general */

int Finish_Init(Tcl_Interp *interp);

/*
 * finishing_rules
 *
 * Calls the Tcl finishing_rules function on a series of base classification
 * bit patterns to produce a series of problem bit patterns.
 */
unsigned int *finishing_rules(Tcl_Interp *interp,
			      finish_t *fin,
			      int mask_offset,
			      char *script,
			      unsigned int *classbits,
			      int len);

/*
 * finishing_solutions
 *
 * Calls the Tcl tcl_find_solutions function on a series of base and problem
 * bit patterns to produce a series of solution bit patterns.
 */
unsigned int *finishing_solutions(Tcl_Interp *interp,
				  char *solu_script,
				  unsigned int *classbits,
				  unsigned int *probbits,
				  int len);

#endif /* _FINISH_H */
