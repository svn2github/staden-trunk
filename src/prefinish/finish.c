#define FINISH_VERSION "1.21"

#include <tcl.h>
#include <limits.h>
#include <ctype.h>

#include "IO.h"
#include "gap_cli_arg.h"
#include "gap_globals.h"
#include "qual.h"
#include "misc.h"
#include "finish.h"
#include "finish_main.h"
#include "finish_utils.h"
#include "finish_long.h"
#include "finish_walk.h"
#include "finish_filter.h"
#include "vseqs.h"
#include "template.h"
#include "list_proc.h"
#include "tagUtils.h"
#include "consen.h"
#include "dna_utils.h"
#include "dust.h"
#include "sequence_formats.h"

/*#define NO_TCL_REGEXP
  #define SYSV_REGEX*/
#include "reg_exp.h"
#define REGEXP_TEMPLATES

/*
 * ---------------------------------------------------------------------------
 * Local prototypes.
 * ---------------------------------------------------------------------------
 */
static int tcl_finish_cmd(ClientData clientData, Tcl_Interp *interp,
			  int objc, Tcl_Obj *CONST objv[]);

static int tcl_finish_obj_cmd(ClientData clientData, Tcl_Interp *interp,
			      int objc, Tcl_Obj *CONST objv[]);

static int tcl_finish_configure(finish_t *fin, Tcl_Interp *interp,
				int objc, Tcl_Obj *CONST objv[]);

static int tcl_classify_bases(finish_t *fin, Tcl_Interp *interp,
			      int objc, Tcl_Obj *CONST objv[]);

static int tcl_find_problems(finish_t *fin, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]);

static int tcl_dump_problems(finish_t *fin, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]);

static int tcl_implement_solutions(finish_t *fin, Tcl_Interp *interp,
				   int objc, Tcl_Obj *CONST objv[]);

static void tcl_finish_delete(ClientData clientdata);

static con_bits_t *parse_con_bits(Tcl_Interp *interp, char *str, int *nbits);

/*
 * ---------------------------------------------------------------------------
 * Finish (de)allocation
 * ---------------------------------------------------------------------------
 */

/*
 * Allocates a new finish structure and initialises it to blank.
 */
static finish_t *finish_new(void) {
    finish_t *fin;
    int i;

    if (NULL == (fin = (finish_t *)xmalloc(sizeof(*fin))))
	return NULL;

    /* Initialise finish_t structure */
    fin->opts.use_avg_insert	= 0;
    fin->opts.mandatory_ratio   = 0.4;
    fin->opts.prob_mandatory    = 0;
    fin->opts.max_score_drop    = 0.1;
    fin->opts.min_template_score= 0.25;
    fin->opts.min_score		= 0.6;
    fin->opts.find_dup_templates= 1;
    fin->opts.dust_level	= 18;
    fin->opts.min_extension	= 50;
    fin->opts.reseq_length	= 400;
    fin->opts.reseq_nsolutions	= 4;
    fin->opts.long_length       = 700;
    fin->opts.long_nsolutions	= 4;
    fin->opts.pwalk_search_dist = 500;
    fin->opts.pwalk_max_match   = 8;
    fin->opts.pwalk_osp_score   = 16;
    fin->opts.pwalk_noligos     = 2;
    fin->opts.pwalk_ntemplates  = 4;
    fin->opts.pwalk_offset1     = 80;
    fin->opts.pwalk_offset2     = 10;
    fin->opts.pwalk_length      = 400;
    fin->opts.pwalk_nsolutions  = 400;
    fin->opts.pwalk_seq_gap     = 20;
    fin->opts.pwalk_consistent_only = 0;
    fin->opts.pwalk_max_err     = 0.02;
    fin->opts.pwalk_min_qual    = 25;
    fin->opts.pwalk_max_err2    = 0.1;
    fin->opts.pwalk_min_qual2   = 15;
    fin->opts.pwalk_end_dist    = 700;
    fin->opts.pwalk_use_template = 1;
    fin->opts.pwalk_use_template_score = 0.1;
    fin->opts.pwalk_dup_template_cost = 1.1;
    fin->opts.pwalk_tag_cons	= 0;
    fin->opts.pwalk_prob_mask	= 0;
    fin->opts.pwalk_tag_type    = NULL;
    memset(&fin->opts.primer_args, 0, sizeof(fin->opts.primer_args));
    fin->io = NULL;
    fin->contig = 0;
    fin->cvec_left = 0;
    fin->cvec_right = 0;
    fin->start = 0;
    fin->length = 0;
    fin->end = 0;
    fin->vc = NULL;
    fin->cons = NULL;
    fin->qual = NULL;
    fin->filtered = NULL;
    fin->orig_qual = NULL;
    fin->alloc_len = 0;
    fin->classify = NULL;
    fin->nclassify = 0;
    fin->base_bits = NULL;
    fin->prob_bits = NULL;
    fin->tag_mask = NULL;
    fin->orig_prob_bits = NULL;
    fin->solution_bits = NULL;
    fin->tarr = NULL;
    fin->template_dup = NULL;
    fin->prob_script = NULL;
    fin->solu_script = NULL;
    fin->left_extent = 0;
    fin->right_extent = 0;
    fin->template_used = NULL;
    fin->template_skip = NULL;
    fin->skip_tags = NULL;
    fin->nskip_tags = 0;
    fin->external_seq = NULL;
    fin->external_seq_rev = NULL;
    fin->external_seq_h = NULL;
    fin->all_cons = NULL;
    fin->all_cons_len = 0;
    fin->all_cons_h = NULL;
    fin->command_token = NULL;
    Tcl_DStringInit(&fin->tag_list);
    for (i = 0; i < 10; i++)
	fin->count[i] = 0;
    fin->cost[EXPERIMENT_RESEQ]  = 1.0;
    fin->cost[EXPERIMENT_LONG]   = 2.0;
    fin->cost[EXPERIMENT_VPWALK] = 3.0;
    fin->cost[EXPERIMENT_CPWALK] = 8.0;
    fin->cost[EXPERIMENT_REVERSE] = 1.0;
    for (i = 0; i < 32; i++) {
	fin->pscore[i] = 1.0; 
	fin->mscore[i] = 1.0;
    }
    memset(fin->opts.debug, 0, 10*sizeof(*fin->opts.debug));
    fin->out_fp = NULL;

    return fin;
}

/*
 * Just deallocates the bits of a finish object that relate to the specific
 * contig being analysed.
 */
static void finish_del_contig(finish_t *fin) {
    if (fin->cons) {
	xfree(fin->cons);
	fin->cons = NULL;
    }

    if (fin->qual) {
	xfree(fin->qual);
	fin->qual = NULL;
    }

    if (fin->filtered) {
	xfree(fin->filtered);
	fin->filtered = NULL;
    }

    if (fin->orig_qual) {
	xfree(fin->orig_qual);
	fin->orig_qual = NULL;
    }

    if (fin->vc) {
	del_vcontig(fin->vc);
	fin->vc = NULL;
    }

    if (fin->classify) {
	xfree(fin->classify);
	fin->classify = NULL;
    }

    if (fin->base_bits) {
	xfree(fin->base_bits);
	fin->base_bits = NULL;
    }

    if (fin->prob_bits) {
	xfree(fin->prob_bits);
	fin->prob_bits = NULL;
    }

    if (fin->tag_mask) {
	xfree(fin->tag_mask);
	fin->tag_mask = NULL;
    }

    if (fin->orig_prob_bits) {
	xfree(fin->orig_prob_bits);
	fin->orig_prob_bits = NULL;
    }

    if (fin->solution_bits) {
	xfree(fin->solution_bits);
	fin->solution_bits = NULL;
    }

    if (fin->tarr) {
	uninit_template_checks(fin->io, fin->tarr);
	fin->tarr = NULL;
    }

    if (fin->template_dup) {
	xfree(fin->template_dup);
	fin->template_dup = NULL;
    }
}

/*
 * Destroys and frees a finish_t structure.
 */
static void finish_del(finish_t *fin) {
    if (!fin)
	return;

    finish_del_contig(fin);

    if (fin->prob_script) {
	xfree(fin->prob_script);
	fin->prob_script = NULL;
    }

    if (fin->solu_script) {
	xfree(fin->solu_script);
	fin->solu_script = NULL;
    }

    if (fin->template_used) {
	xfree(fin->template_used);
	fin->template_used = NULL;
    }

    if (fin->template_skip) {
	xfree(fin->template_skip);
	fin->template_skip = NULL;
    }

    Tcl_DStringFree(&fin->tag_list);

    if (fin->external_seq) {
	xfree(fin->external_seq);
    }

    if (fin->external_seq_rev) {
	xfree(fin->external_seq_rev);
    }

    if (fin->external_seq_h) {
	free_hash8n(fin->external_seq_h);
	fin->external_seq_h = NULL;
    }

    if (fin->all_cons) {
	xfree(fin->all_cons);
	fin->all_cons = NULL;
    }

    if (fin->all_cons_h) {
	free_hash8n(fin->all_cons_h);
	fin->all_cons_h = NULL;
    }

    xfree(fin);
}

/*
 * ---------------------------------------------------------------------------
 * Tcl initialisation and top-level registered functions
 * ---------------------------------------------------------------------------
 */

/*
 * This is called when the library is dynamically linked in with the calling
 * program. Use it to initialise any tables and to register the necessary
 * commands.
 */
int Prefinish_Init(Tcl_Interp *interp) {
    if (NULL == Tcl_CreateObjCommand(interp, "finish", tcl_finish_cmd,
				     (ClientData)NULL,
				     (Tcl_CmdDeleteProc *) NULL))
	return TCL_ERROR;

    return TCL_OK;
}

int Prefinish_SafeInit(Tcl_Interp *interp) {
    return Prefinish_Init(interp);
}

/*
 * tcl_finish_cmd
 *
 * Initialises a new the finish object and registers a Tcl command associated
 * to that instance of finish.
 */
/*ARGSUSED*/
static int tcl_finish_cmd(ClientData clientData, Tcl_Interp *interp,
			  int objc, Tcl_Obj *CONST objv[]) {
    finish_t *fin;

    if (objc < 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 Tcl_GetString(objv[0]), " objName ?options?\"",
			 (char *)NULL);
	return TCL_ERROR;
    }

    /* Allocate obj and initialise 'fin' */
    if (NULL == (fin = finish_new()))
	return TCL_ERROR;

    printf("Initialising prefinish version %s\n", FINISH_VERSION);

    /* Add a new command to handle finish methods */
    fin->command_token =
	Tcl_CreateObjCommand(interp, Tcl_GetString(objv[1]),
			     tcl_finish_obj_cmd, (ClientData)fin,
			     (Tcl_CmdDeleteProc *)tcl_finish_delete);

    /* Configure the new object */
    if (objc >= 2)
	tcl_finish_configure(fin, interp, objc - 1, objv + 1);

    /* Return result */
    Tcl_SetObjResult(interp, objv[1]);

    return TCL_OK;
}

/*
 * tcl_finish_obj_cmd
 *
 * This is the command registered for a specific finish object.
 * It Implements the finish methods for that object.
 */
static int tcl_finish_obj_cmd(ClientData clientData, Tcl_Interp *interp,
			      int objc, Tcl_Obj *CONST objv[]) {
    finish_t *fin = (finish_t *)clientData;
    int index;
    static char *finishCmds[] = {
	"configure",		"classify",		"find_problems",
	"implement_solutions",	"delete",		"dump_problems",
	(char *)NULL
    };
    enum finishCmds {
	FIN_CONFIGURE,		FIN_CLASSIFY,		FIN_FIND_PROBLEMS,
	FIN_IMPLEMENT_SOLUTIONS,FIN_DELETE,		FIN_DUMP_PROBLEMS
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "command ?option val? ...");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], finishCmds, "command", 0,
			    &index) != TCL_OK) {
	return TCL_ERROR;
    }

    switch ((enum finishCmds)index) {
    case FIN_CONFIGURE:
	if (fin->opts.debug[FIN_DEBUG])
	    puts("Finish configure");
	return tcl_finish_configure(fin, interp, objc - 1, objv + 1);
	
    case FIN_CLASSIFY:
	if (fin->opts.debug[FIN_DEBUG])
	    puts("Finish classify");
	return tcl_classify_bases(fin, interp, objc - 1, objv + 1);
	
    case FIN_FIND_PROBLEMS:
	if (fin->opts.debug[FIN_DEBUG])
	    puts("Finish find problems");
	return tcl_find_problems(fin, interp, objc - 1, objv + 1);
	
    case FIN_IMPLEMENT_SOLUTIONS:
	if (fin->opts.debug[FIN_DEBUG])
	    puts("Finish implement solutions");
	return tcl_implement_solutions(fin, interp, objc - 1, objv + 1);

    case FIN_DUMP_PROBLEMS:
	if (fin->opts.debug[FIN_DEBUG])
	    puts("Finish dump problems");
	return tcl_dump_problems(fin, interp, objc - 1, objv + 1);

    case FIN_DELETE:
	if (objc != 2) {
	    Tcl_WrongNumArgs(interp, 2, objv, "finish delete");
	    return TCL_ERROR;
	}
	Tcl_DeleteCommandFromToken(interp, fin->command_token);
	break;
    }

    return TCL_OK;
}

/*
 * Configures the finish_t.template_skip[] array.
 * skip_template_file is the name of a file of template names.
 * 'skip' is false if this is for -available_templates, otherwise it is true
 * for the -skip_templates option.
 * Set skip_template_file to "" to clear any existing setup.
 */
static int configure_skip_templates(finish_t *fin,
				    Tcl_Interp *interp,
				    char *skip_template_file,
				    int skip) {
    FILE *fp;
    char *r_exp = NULL;
    int r_exp_len = 0, r_exp_alen = 0;
    char line[1024];
    int r_exp_elements;

    if (NULL == (fp = fopen(skip_template_file, "r"))) {
	verror(ERR_WARN, "finish_init", "Could not open file '%s'",
	       skip_template_file);
	return TCL_OK; /* Not strictly correct, but it aids our interface */
    }
    
    /* Initialise array of templates to reject automatically */
    if (*skip_template_file == 0 && fin->template_skip) {
	xfree(fin->template_skip);
	fin->template_skip = NULL;
    }
    if (!fin->template_skip) {
	fin->template_skip = (int *)xcalloc(Ntemplates(fin->io)+1,
					    sizeof(int));
	if (!fin->template_skip)
	    return TCL_ERROR;
    }

    /* Picking available templates? If so initialiase everything to reject */
    if (!skip) {
	int i;
	for (i = 1; i <= Ntemplates(fin->io); i++) {
	    fin->template_skip[i] = 1;
	}
    }
    
#ifdef REGEXP_TEMPLATES
    do {
	r_exp_elements = 0;
	r_exp_len = 0;
	if (r_exp)
	    *r_exp = 0;
	/*
	 * Combine lines together in exp|exp|exp|... format for up to
	 * 80 expressions in a row. More than this and we run the risk of
	 * the regexp compiler from breaking.
	 */
#define REGBLOCKS 80
	while (fgets(line, 1024, fp) && ++r_exp_elements < REGBLOCKS) {
	    char *cp;
	    int len;
	    
	    if ((cp = strchr(line, '\n')))
		*cp = 0;

	    if (*line == 0)
		continue;
	    
	    /* (Re)allocate and initialise r_exp, our regular expression */
	    len = strlen(line);
	    r_exp_len += len + 3;
	    if (r_exp_len > r_exp_alen) {
		char *r_exp_new;
		r_exp_alen = r_exp_len * 2 + 1;
		
		if (NULL == (r_exp_new = xrealloc(r_exp, r_exp_alen))) {
		    verror(ERR_WARN, "finish_init",
			   "Not enough memory to build regexp");
		    xfree(r_exp);
		    return TCL_ERROR;
		}
		if (!r_exp) {
		    *r_exp_new = 0;
		}
		r_exp = r_exp_new;
	    }
	    
	    /* Add "|line" to regexp */
	    if (*r_exp) {
		char *ptr = r_exp + strlen(r_exp);
		sprintf(ptr, "|%s", line);
	    } else {
		sprintf(r_exp, "%s", line);
	    }
	}
    
	/* Compiler our up-to-80-element expression */
	if (r_exp) {
	    char *comp;
	    if (fin->opts.debug[FIN_DEBUG_REGEXP] > 0)
		printf("Compiling regexp...\n");
	    if (NULL == (comp = REGCMP(interp, r_exp))) {
		if (fin->opts.debug[FIN_DEBUG_REGEXP] > 0)
		    printf("Failed!\n");
		verror(ERR_WARN, "finish_init",
		       "Could not compile regexp '%s'", r_exp);
	    } else {
		int i;
		if (fin->opts.debug[FIN_DEBUG_REGEXP] > 0)
		    printf("Done\n");
		/* Iterate around all templates matching this expression */
		for (i = 1; i <= Ntemplates(fin->io); i++) {
		    char *tname;

		    if (fin->template_skip[i] == skip)
			/* Already matched in a previous block.. */
			continue;

		    if (!(tname = get_template_name(fin->io, i)))
			continue;

		    if (REGEX(interp, tname, comp)) {
			if (fin->opts.debug[FIN_DEBUG_REGEXP] > 0)
			    printf("%sing template '%s'\n",
				   skip ? "Skipp" : "Us", tname);
			fin->template_skip[i] = skip;
		    }
		}
	    
		REGFREE(interp, comp);
	    }
	}
    } while (r_exp_elements == REGBLOCKS);
    xfree(r_exp);

    if (fin->opts.debug[FIN_DEBUG_REGEXP] > 0)
	printf("Regexp matching done\n");

#else /* ifdef REGEXP_TEMPLATES */
    while (fgets(line, 1024, fp)) {
	char *cp;
	int len;
	    
	if ((cp = strchr(line, '\n')))
	    *cp = 0;

	if (*line == 0)
	    continue;
	    
	{
	    int tnum;
	    if (!(tnum = template_name_to_number(fin->io, line))) {
		verror(ERR_WARN, "finish_init",
		       "Template '%s' not found in database", line);
	    } else {
		fin->template_skip[tnum] = skip;
	    }
	}
    }
#endif
    fclose(fp);

    return TCL_OK;
}

#if 0
static void finish_filter(finish_t *fin) {
    int i, clen;

    if (fin->opts.debug[FIN_DEBUG_DUST])
	puts("Filtering using dust...");

    clen = io_clength(fin->io, fin->contig);
    if (NULL == (fin->filtered = (char *)xmalloc(clen)))
	return;

    /* Filter using Dust */
    memcpy(fin->filtered, fin->cons, clen);
    set_dust_level(fin->opts.dust_level);
    dust(clen, fin->filtered);

    /* Look for low complexity data with 32 of the end, if so extend to ends */
    for (i = 0; i < clen && i < 32; i++) {
	if (fin->filtered[i] == 'N') {
	    for (i = 0; i < 32 && i < clen; i++)
		fin->filtered[i] = 'N';
	    break;
	}
    }

    for (i = 0; clen-1-i >= 0 && i < 32; i++) {
	if (fin->filtered[clen-1-i] == 'N') {
	    for (i = 0; clen-1-i >= 0 && i < 32; i++)
		fin->filtered[clen-1-i] = 'N';
	    break;
	}
    }

    if (fin->opts.debug[FIN_DEBUG_DUST] > 1)
	printf("%.*s\n", clen, fin->filtered);
}
#endif

/*
 * Called from tcl_finish_configure to read a sequence from a file and
 * returns the result. Handles formats supported by seq_utils get_seq().
 *
 * Returns a malloced string containing the sequence on success. This should
 *         be freed by the caller.
 * Returns NULL on failure.
 */
static char *read_external_seq_file(char *filename) {
    char *seq = NULL;
    int seq_len, r;

    r = get_seq(&seq, 0 /* unused! */, &seq_len, filename,
		NULL /* entryname */);

    if (r != 0 || seq == NULL)
	return NULL;

    /* realloc and null-term as I'm not sure if get_seq() does this for us */
    seq = (char *)xrealloc(seq, seq_len+1);
    if (!seq)
	return NULL;
    seq[seq_len] = 0;

    return seq;
}


/*
 * Called when "$fin configure ..." is called. Also called with the additonal
 * command line arguments after creating a new finish object.
 */
static int tcl_finish_configure(finish_t *fin, Tcl_Interp *interp,
				int objc, Tcl_Obj *CONST objv[]) {
    int i;
    Hidden_params hparams;
    int maxseq = gap4_global_get_maxseq();
    char *eseq = NULL; /* external sequence */
    int eseq_alloced = 0;

    /* The mapping of the argument strings to our structure above */
    /* The 'def' section of cli_args is not used */
    cli_args conf[] = {
	{"-io",       	       ARG_IO,    1, 0,
	        offsetof(finish_t, io)},
	{"-contig",   	       ARG_STR,   1, NULL,
	        offsetof(finish_t, args.contig)},
	{"-check_contigs",     ARG_STR,   1, NULL,
	        offsetof(finish_t, args.ccontigs)},
	{"-external_seq",      ARG_STR,   1, NULL,
	        offsetof(finish_t, args.eseq)},
	{"-skip_template_file",ARG_STR,   1, NULL,
	        offsetof(finish_t, args.skip_template_file)},
	{"-avail_template_file",ARG_STR,   1, NULL,
	        offsetof(finish_t, args.avail_template_file)},
	{"-external_seq_file",ARG_STR,   1, NULL,
	        offsetof(finish_t, args.external_seq_file)},
	{"-output_file",ARG_STR,   1, NULL,
	        offsetof(finish_t, args.output_file)},
	{"-use_avg_insert",    ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.use_avg_insert)},
	{"-mandatory_ratio",   ARG_DBL, 1, NULL,
	 	offsetof(finish_t, opts.mandatory_ratio)},
	{"-prob_mandatory",    ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.prob_mandatory)},
	{"-max_score_drop",    ARG_DBL, 1, NULL,
	        offsetof(finish_t, opts.max_score_drop)},
	{"-min_template_score",    ARG_DBL, 1, NULL,
	        offsetof(finish_t, opts.min_template_score)},
	{"-min_score",    ARG_DBL, 1, NULL,
	        offsetof(finish_t, opts.min_score)},
	{"-find_dup_templates",    ARG_INT, 1, NULL,
	        offsetof(finish_t, opts.find_dup_templates)},
	{"-dust_level",       ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.dust_level)},
	{"-min_extension",       ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.min_extension)},
	{"-reseq_length",      ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.reseq_length)},
	{"-reseq_nsolutions",  ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.reseq_nsolutions)},
	{"-long_length",       ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.long_length)},
	{"-long_nsolutions",   ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.long_nsolutions)},
	{"-pwalk_search_dist", ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_search_dist)},
	{"-pwalk_max_match",   ARG_DBL, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_max_match)},
	{"-pwalk_osp_score",   ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_osp_score)},
	{"-pwalk_noligos",     ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_noligos)},
	{"-pwalk_ntemplates",  ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_ntemplates)},
	{"-pwalk_offset1",     ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_offset1)},
	{"-pwalk_offset2",     ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_offset2)},
	{"-pwalk_length",      ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_length)},
	{"-pwalk_nsolutions",  ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_nsolutions)},
	{"-pwalk_seq_gap",     ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_seq_gap)},
	{"-pwalk_consistent_only", ARG_INT,   1, NULL,
	 	offsetof(finish_t, opts.pwalk_consistent_only)},
	{"-pwalk_max_err",     ARG_DBL, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_max_err)},
	{"-pwalk_min_qual",    ARG_INT, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_min_qual)},
	{"-pwalk_max_err2",     ARG_DBL, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_max_err2)},	
	{"-pwalk_min_qual2",    ARG_INT, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_min_qual2)},
	{"-pwalk_end_dist",    ARG_INT, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_end_dist)},
	{"-pwalk_use_template",ARG_INT, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_use_template)},
	{"-pwalk_use_template_score",  ARG_DBL, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_use_template_score)},
	{"-pwalk_dup_template_cost",  ARG_DBL, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_dup_template_cost)},
	{"-pwalk_tag_cons",    ARG_INT, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_tag_cons)},
	{"-pwalk_prob_mask",    ARG_INT, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_prob_mask)},
	{"-pwalk_tag_type",    ARG_STR, 1, NULL,
	 	offsetof(finish_t, opts.pwalk_tag_type)},
	{"-reseq_cost", ARG_FLOAT, 1, NULL,
		 offsetof(finish_t, cost[EXPERIMENT_RESEQ])},
	{"-long_cost",  ARG_FLOAT, 1, NULL,
		 offsetof(finish_t, cost[EXPERIMENT_LONG])},
	{"-vpwalk_cost", ARG_FLOAT, 1, NULL,
		 offsetof(finish_t, cost[EXPERIMENT_VPWALK])},
	{"-cpwalk_cost", ARG_FLOAT, 1, NULL,
		 offsetof(finish_t, cost[EXPERIMENT_CPWALK])},
	{"-reverse_cost", ARG_FLOAT, 1, NULL,
		 offsetof(finish_t, cost[EXPERIMENT_REVERSE])},
	{"-pscores", ARG_STR, 1, NULL, offsetof(finish_t, args.pscores)},
	{"-mscores", ARG_STR, 1, NULL, offsetof(finish_t, args.mscores)},
	{"-primer_min_tm", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.min_tm)},
	{"-primer_max_tm", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.max_tm)},
	{"-primer_opt_tm", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.opt_tm)},
	{"-primer_min_gc", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.min_gc)},
	{"-primer_max_gc", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.max_gc)},
	{"-primer_opt_gc", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.opt_gc)},
	{"-primer_min_len", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.min_len)},
	{"-primer_max_len", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.max_len)},
	{"-primer_opt_len", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.opt_len)},
	{"-primer_max_end_stability", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.max_end_stability)},
	{"-primer_salt_conc", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.salt_conc)},
	{"-primer_dna_conc", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.dna_conc)},
	{"-primer_self_any", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.self_any)},
	{"-primer_self_end", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.self_end)},
	{"-primer_gc_clamp", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.gc_clamp)},
	{"-primer_max_poly_x", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.max_poly_x)},
	{"-primer_max_end_stability", ARG_DBL, 1, NULL,
	         offsetof(finish_t, opts.primer_args.max_end_stability)},
	{"-debug0", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[0])},
	{"-debug1", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[1])},
	{"-debug2", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[2])},
	{"-debug3", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[3])},
	{"-debug4", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[4])},
	{"-debug5", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[5])},
	{"-debug6", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[6])},
	{"-debug7", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[7])},
	{"-debug8", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[8])},
	{"-debug9", ARG_INT, 1, NULL,
	 	 offsetof(finish_t, opts.debug[9])},
	{NULL,      0,       0, NULL, 0}
    };

    /*
     * First things first, add a header to the output window. This shows the
     * date and function name.
     */
    if (fin->opts.debug[FIN_DEBUG])
	vfuncheader("finish_init");

    /* Parse the arguments */
    memset(&fin->args, 0, sizeof(fin->args));
    if (-1 == gap_parse_obj_config(conf, fin, objc, (Tcl_Obj **)objv)) {
	return TCL_ERROR;
    }
    
    /* Contigs to check primers against - hash them */
    if (fin->args.ccontigs && *fin->args.ccontigs) {
	int num_check_contigs;
	contig_list_t *check_contigs = NULL;
	Contig_parms *clist;

	if (fin->opts.debug[FIN_DEBUG])
	    printf("Hashing consensus seqs\n");

	if (!fin->io)
	    return TCL_ERROR;

	active_list_contigs(fin->io, fin->args.ccontigs, &num_check_contigs,
			    &check_contigs);
	clist = get_contig_list(0, fin->io, num_check_contigs, check_contigs);
	
	/* Deallocate old data if required */
	if (fin->all_cons) {
	    xfree(fin->all_cons);
	    fin->all_cons = NULL;
	}
	fin->all_cons_len = 0;
	if (fin->all_cons_h) {
	    free_hash8n(fin->all_cons_h);
	    fin->all_cons_h = NULL;
	}

	/* Compute the entire consensus, used for chromosomal primer match */
	if (NULL == (fin->all_cons = (char *)xmalloc(maxseq))) {
	    return TCL_ERROR;
	}
	hparams.do_it = 0;
	if (make_consensus(ADDTITLE | NORMALCONSENSUS, fin->io, fin->all_cons,
			   NULL /*qual*/, clist, num_check_contigs,
			   &fin->all_cons_len, max_gel_len(fin->io), maxseq,
			   hparams, gap4_global_get_consensus_cutoff() )) {
	    xfree(clist);
	    return TCL_ERROR;
	}

	depad_seq(fin->all_cons, &fin->all_cons_len, NULL);

	/* hash it */
	if (init_hash8n(fin->all_cons_len, FIN_MAXPRIMERLEN,
			4 /* word_length */,
			0 /* max_matches - unused */,
			0 /* min_match - unused */,
			1 /* job */,
			&fin->all_cons_h)) {
	    verror(ERR_WARN, "finish_init", "Failed to hash consenus");
	    fin->all_cons_h = NULL;
	    xfree(fin->all_cons);
	    fin->all_cons = NULL;
	}
	
	fin->all_cons_h->seq1 = fin->all_cons;
	fin->all_cons_h->seq1_len = fin->all_cons_len;

	hash_seqn(fin->all_cons_h, 1);
	store_hashn(fin->all_cons_h);

	if (check_contigs)
	    xfree(check_contigs);
	xfree(clist);
    }

    /* External_seq_file takes precedence over external_seq argument */
    if (fin->args.external_seq_file) {
	eseq = read_external_seq_file(fin->args.external_seq_file);
	eseq_alloced = 1;
    } else if (fin->args.eseq && *fin->args.eseq) {
	eseq = fin->args.eseq;
	eseq_alloced = 0;
    } else {
	eseq = NULL;
	eseq_alloced = 0;
    }

    /* External sequence to check primers against - hash this too */
    if (eseq) {
	int i;

	/* Deallocate old data if required */
	if (fin->external_seq) {
	    xfree(fin->external_seq);
	    fin->external_seq = NULL;
	}
	if (fin->external_seq_rev) {
	    xfree(fin->external_seq_rev);
	    fin->external_seq_rev = NULL;
	}
	fin->external_seq_len = 0;
	if (fin->external_seq_h) {
	    free_hash8n(fin->external_seq_h);
	    fin->external_seq_h = NULL;
	}

	/* Copy and depad sequence */
	fin->external_seq_len = strlen(eseq);
	fin->external_seq = (char *)xmalloc(fin->external_seq_len+1);
	strcpy(fin->external_seq, eseq);
	depad_seq(fin->external_seq, &fin->external_seq_len, NULL);

	/* Uppercase it */
	for (i = 0; i < fin->external_seq_len; i++)
	    fin->external_seq[i] = toupper(fin->external_seq[i]);

#if !defined(USE_OSP)
	/* For primer3 only, reverse and complement it */
	fin->external_seq_rev = strdup(fin->external_seq);
	complement_seq(fin->external_seq_rev, fin->external_seq_len);
#endif

	/* Hash it */
	if (init_hash8n(fin->external_seq_len, FIN_MAXPRIMERLEN,
			4 /* word_length */,
			0 /* max_matches - unused */,
			0 /* min_match - unused */,
			1 /* job */,
			&fin->external_seq_h)) {
	    verror(ERR_WARN, "finish_init", "Failed to hash external_seq");
	    fin->external_seq_h = NULL;
	    xfree(fin->external_seq);
	    fin->external_seq = NULL;
	}
	
	fin->external_seq_h->seq1 = fin->external_seq;
	fin->external_seq_h->seq1_len = fin->external_seq_len;

	hash_seqn(fin->external_seq_h, 1);
	store_hashn(fin->external_seq_h);
    }

    if (eseq_alloced && eseq) {
	xfree(eseq);
    }

    if (fin->args.avail_template_file) {
	if (TCL_OK != configure_skip_templates(fin, interp,
					       fin->args.avail_template_file,
					       0))
	    return TCL_ERROR;
    }

    if (fin->args.skip_template_file) {
	if (TCL_OK != configure_skip_templates(fin, interp,
					       fin->args.skip_template_file,
					       1))
	    return TCL_ERROR;
    }

    if (fin->args.output_file) {
	if (fin->out_fp)
	    fclose(fin->out_fp);
	if (NULL == (fin->out_fp = fopen(fin->args.output_file, "w"))) {
	    verror(ERR_WARN, "finish_init",
		   "Could not output to '%s'", fin->args.output_file);
	    return TCL_ERROR;
	}
    }

    /* Which contig are we working on? */
    if (fin->args.contig) {
	int num_contigs;
	contig_list_t *contigs = NULL;

	if (!fin->io)
	    return TCL_ERROR;

	/* Deallocate old structures */
	finish_del_contig(fin);

	/* Find which contig? */
	active_list_contigs(fin->io, fin->args.contig, &num_contigs, &contigs);
	if (num_contigs == 0)
	    return TCL_ERROR;
    
	fin->contig = contigs[0].contig;
	fin->start = contigs[0].start;
	fin->end = contigs[0].end;
	fin->length = fin->end - fin->start + 1;
	fin->vc = new_vcontig(fin->io, contigs[0].contig);

	/*
	 * Compute consensus - over entire contig as we need to evaluate
	 * sequences which span only part of our 'start -> end' region.
	 */
	if (NULL == (fin->cons = (char *)xmalloc(io_clength(fin->io,
							    fin->contig))))
	    return TCL_ERROR;

	if (NULL == (fin->qual = (float *)xmalloc(io_clength(fin->io,
							     fin->contig) *
						  sizeof(float))))
	    return TCL_ERROR;

	if (NULL == (fin->orig_qual = (float *)xmalloc(io_clength(fin->io,
								  fin->contig)
						       * sizeof(float))))
	    return TCL_ERROR;

	fin->alloc_len = io_clength(fin->io, fin->contig);
	fin->left_extent = 0;
	fin->right_extent = io_clength(fin->io, fin->contig)-1;

	calc_consensus(fin->contig, 1, io_clength(fin->io, fin->contig),
		       CON_SUM, fin->cons, NULL, fin->qual, NULL,
		       gap4_global_get_consensus_cutoff(),
		       gap4_global_get_quality_cutoff(),
		       database_info, (void *)(fin->io));
	fin->vc->cons = fin->cons;
	memcpy(fin->orig_qual, fin->qual, io_clength(fin->io, fin->contig)
	       * sizeof(fin->orig_qual[0]));

	/*
	 * Low complexity filtering using the Dust algorithm.
	 */
	finish_filter(fin);

	/* Initialise our count of how many times each template is used */
	if (fin->template_used)
	    xfree(fin->template_used);
	fin->template_used = (int *)xcalloc(Ntemplates(fin->io)+1,
					    sizeof(int));
	if (!fin->template_used)
	    return TCL_ERROR;

	/* Check for cloning (cosmid, bac, etc) vector at contig ends */
	find_cloning_vector(fin->io, fin->contig,
			    &fin->cvec_left, &fin->cvec_right);

	if (contigs)
	    xfree(contigs);
    }

    /* problem scores */
    if (fin->args.pscores) {
	char **scores;
	int nscores;
	SplitList(fin->args.pscores, &nscores, &scores);

	for (i = 0; i < nscores; i++) {
	    fin->pscore[i] = atof(scores[i]);
	}
	Tcl_Free((char *)scores);
    }

    /* mandatory scores */
    if (fin->args.mscores) {
	char **scores;
	int nscores;
	SplitList(fin->args.mscores, &nscores, &scores);

	for (i = 0; i < nscores; i++) {
	    fin->mscore[i] = atof(scores[i]);
	}
	Tcl_Free((char *)scores);
    }

    return TCL_OK;
}

static int tcl_classify_bases(finish_t *fin, Tcl_Interp *interp,
			      int objc, Tcl_Obj *CONST objv[]) {
    con_bits_t *con_bits;
    int bitsc;

    /* A structure definition to store the arguments in */
    typedef struct {
	char *bits;
    } cp_args;

    /* The mapping of the argument strings to our structure above */
    cp_args args;
    cli_args a[] = {
	{"-bits",     ARG_STR, 1, NULL, offsetof(cp_args, bits)},
	{NULL,      0,       0, NULL, 0}
    };

    /*
     * First things first, add a header to the output window. This shows the
     * date and function name.
     */
    vfuncheader("classify_bases");

    /* Parse the arguments */
    if (-1 == gap_parse_obj_args(a, &args, objc, (Tcl_Obj **)objv)) {
	return TCL_ERROR;
    }

    /* Parse and build bits structure */
    con_bits = parse_con_bits(interp, args.bits, &bitsc);

    fin->classify = con_bits;
    fin->nclassify = bitsc;

    fin->base_bits = classify_bases(fin, fin->start, fin->end, NULL,
				    database_info, (void *)(fin->io));

    if (NULL == fin->base_bits)
	return TCL_ERROR;

    return TCL_OK;
}


/*
 * Given a tag 'curr_tag' on a reading 'r', we convert the tag corrdinates
 * into absolute contig coordinates and return these in start_p and end_p
 * accorindgly.
 */
static void consensus_tag_pos(finish_t *fin,
			      GAnnotations *curr_tag,
			      int rnum,
			      int *start_p,
			      int *end_p) {
    GReadings r;
    int start, end;

    gel_read(fin->io, rnum, r);

    start = r.position - r.start +
	(r.sense
	 ? r.length - (curr_tag->length + curr_tag->position - 1)
	 : curr_tag->position - 1);
    end = start + curr_tag->length - 1;
    if (end > r.position + r.sequence_length - 1)
	end = r.position + r.sequence_length - 1;

    if (start < 1)
	start = 1;
    if (end > fin->length)
	end = fin->length;

    if (start_p)
	*start_p = start;
    if (end_p)
	*end_p = end;
}


static int tcl_find_problems(finish_t *fin, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    Tcl_Obj *prob_objv[2];
    Tcl_Obj *solu_objv[3];
    int i;
    int init_orig_prob;
    GAnnotations *curr_tag;
    int tag_start, tag_end;
    int do_ctags;

    typedef struct {
	Tcl_Obj *prob;
	Tcl_Obj *solu;
	char *tag_list;
    } fp_args;

    fp_args args;
    cli_args a[] = {
	{"-problem_command",  ARG_OBJ, 1, NULL, offsetof(fp_args, prob)},
	{"-solution_command", ARG_OBJ, 1, NULL, offsetof(fp_args, solu)},
	{"-tag_types",	      ARG_STR, 1, "",   offsetof(fp_args, tag_list)},
	{NULL,      0,       0, NULL, 0}
    };

    vfuncheader("find_problems");

    /* Parse the arguments */
    if (-1 == gap_parse_obj_args(a, &args, objc, (Tcl_Obj **)objv)) {
	return TCL_ERROR;
    }

    if (fin->prob_bits)
	xfree(fin->prob_bits);
    fin->prob_bits = (unsigned int *)xmalloc(fin->length*sizeof(unsigned int));
    if (fin->orig_prob_bits) {
	/* Already initialised - don't change! */
	init_orig_prob = 0;
    } else {
	fin->orig_prob_bits =
	    (unsigned int *)xmalloc(fin->length*sizeof(unsigned int));
	init_orig_prob = 1;
    }
    if (!fin->prob_bits || !fin->orig_prob_bits)
	return TCL_ERROR;

    if (fin->solution_bits)
	xfree(fin->solution_bits);
    fin->solution_bits = (unsigned int *)xmalloc(fin->length *
						 sizeof(unsigned int));
    if (!fin->solution_bits)
	return TCL_ERROR;

    if (fin->prob_script)
	xfree(fin->prob_script);
    fin->prob_script = strdup(Tcl_GetString(args.prob));

    if (fin->solu_script)
	xfree(fin->solu_script);
    fin->solu_script = strdup(Tcl_GetString(args.solu));

    if (fin->tag_mask) {
	xfree(fin->tag_mask);
	fin->tag_mask = NULL;
    }

    prob_objv[0] = args.prob;
    prob_objv[1] = Tcl_NewIntObj(0);
    Tcl_IncrRefCount(prob_objv[0]);
    Tcl_IncrRefCount(prob_objv[1]);

    solu_objv[0] = args.solu;
    solu_objv[1] = prob_objv[1];
    Tcl_IncrRefCount(solu_objv[0]);

    SplitList(args.tag_list, &fin->nskip_tags, &fin->skip_tags);
    if (fin->nskip_tags) {
	curr_tag = vtagget(fin->io, -fin->contig, fin->nskip_tags,
			   fin->skip_tags);
	fin->tag_mask = (unsigned int *)xcalloc(fin->length,
						sizeof(unsigned int));
    }
    if (!curr_tag) {
	do_ctags = 0;
    } else {
	do_ctags = 1;
	tag_start = curr_tag->position;
	tag_end = curr_tag->position + curr_tag->length-1;
    }

    for (i = 0; i < fin->length; i++) {
	Tcl_Obj *pobj;

	/* Tagged regions, by definition, have no problems */
	if (fin->nskip_tags && do_ctags) {
	    while (i+1 > tag_end) {
		/* Use next tag */
		curr_tag = vtagget(fin->io, 0, fin->nskip_tags,
				   fin->skip_tags);

		if (curr_tag) {
		    char type[5];
		    tag_start = curr_tag->position;
		    if (curr_tag->position + curr_tag->length-1 > tag_end)
			tag_end = curr_tag->position + curr_tag->length-1;
		    if (fin->opts.debug[FIN_DEBUG])
			printf("Skipping problems in tagged cons "
			       "%d-%d (type %s)\n",
			       tag_start, tag_end,
			       type2str(curr_tag->type, type));
		} else {
		    do_ctags = 0;
		    break;
		}
	    }

	    if (i+1 >= tag_start && i+1 <= tag_end) {
		fin->tag_mask[i] = 1;
		fin->prob_bits[i] = 0;
		fin->orig_prob_bits[i] = 0;
		fin->solution_bits[i] = 0;
		continue;
	    }
	}

	/* Call problem_command */
	Tcl_SetIntObj(prob_objv[1], fin->base_bits[i]);
	if (TCL_OK != Tcl_EvalObjv(interp, 2, prob_objv, TCL_EVAL_GLOBAL)) {
	    fprintf(stderr, "Eval failed '%s'\n", Tcl_GetStringResult(interp));
	    continue;
	}

	pobj = Tcl_GetObjResult(interp);
	Tcl_IncrRefCount(pobj);

	/* Store result in fin->prob_bits[i] */
	if (TCL_OK != Tcl_GetIntFromObj(interp, pobj,
					(int *)&fin->prob_bits[i])) {
	    fprintf(stderr, "Tcl_GetIntFromObj: %s\n",
		    Tcl_GetStringResult(interp));
	    Tcl_DecrRefCount(pobj);
	    continue;
	}
	if (init_orig_prob)
	    fin->orig_prob_bits[i] = fin->prob_bits[i];

	/* Optimise - if no problems then do not call solution command */
	if (!fin->prob_bits[i]) {
	    fin->solution_bits[i] = 0;
	    Tcl_DecrRefCount(pobj);
	    continue;
	}

	/* Call solution_command */
	solu_objv[2] = pobj;
	if (TCL_OK != Tcl_EvalObjv(interp, 3, solu_objv, TCL_EVAL_GLOBAL)) {
	    fprintf(stderr, "Eval failed '%s'\n", Tcl_GetStringResult(interp));
	    Tcl_DecrRefCount(pobj);
	    continue;
	}

	/* Store result in fin->solution_bits[i] */
	Tcl_DecrRefCount(pobj);
	
	pobj = Tcl_GetObjResult(interp);
	Tcl_IncrRefCount(pobj);
	if (TCL_OK != Tcl_GetIntFromObj(interp, pobj,
					(int *)&fin->solution_bits[i])) {
	    fprintf(stderr, "Tcl_GetIntFromObj: %s\n",
		    Tcl_GetStringResult(interp));
	    Tcl_DecrRefCount(pobj);
	    continue;
	}
	
	Tcl_DecrRefCount(pobj);

	if (fin->opts.debug[FIN_DEBUG] > 2)
	    printf("%d: bits 0x%x, probs 0x%x, soln 0x%x\n",
		   i,
		   fin->base_bits[i],
		   fin->prob_bits[i],
		   fin->solution_bits[i]);
    }

    /* Skip tagged regions in sequences too.. */
    if (fin->nskip_tags) {
	int rnum;

	for (rnum = io_clnbr(fin->io, fin->contig);
	     rnum;
	     rnum = io_rnbr(fin->io, rnum)) {
	    curr_tag = vtagget(fin->io, rnum, fin->nskip_tags,
			       fin->skip_tags);

	    if (!curr_tag)
		continue;

	    while (curr_tag && curr_tag != (GAnnotations *)-1) {
		char type[5];
		int start, end;
		int i;

		consensus_tag_pos(fin, curr_tag, rnum, &start, &end);

		if (fin->opts.debug[FIN_DEBUG] > 1 && start <= end)
		    printf("Skipping problems in tagged seq #%d "
			   "%d-%d (type %s)\n",
			   rnum, start, end, type2str(curr_tag->type, type));
		for (i = start; i <= end; i++) {
		    fin->tag_mask[i-1] = 1;
		    fin->prob_bits[i-1] = 0;
		    fin->solution_bits[i-1] = 0;
		}

		curr_tag = vtagget(fin->io, 0, fin->nskip_tags,
				   fin->skip_tags);
	    }
	}
    }

    Tcl_DecrRefCount(prob_objv[0]);
    Tcl_DecrRefCount(prob_objv[1]);
    Tcl_DecrRefCount(solu_objv[0]);
    ckfree((char *)fin->skip_tags);

    return TCL_OK;
}


/*
 * tcl_dump_problems() - part of the finishing software.
 *
 * Returns a Tcl list (object) of the problem array
 */
static int tcl_dump_problems(finish_t *fin, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    Tcl_Obj *lobj, *obj;
    int i;

    /* Create a list object */
    lobj = Tcl_NewListObj(0, NULL);
    if (!lobj)
	return TCL_ERROR;
    Tcl_IncrRefCount(lobj);

    /* Create integer objects for each problem and add them to the list */
    for (i = 0; i < fin->length; i++) {
	obj = Tcl_NewIntObj(fin->prob_bits[i]);
	if (!obj) {
	    Tcl_DecrRefCount(lobj);
	    return TCL_ERROR;
	}

	Tcl_ListObjAppendElement(interp, lobj, obj);
    }

    /* Set the result and return */
    Tcl_SetObjResult(interp, lobj);
    Tcl_DecrRefCount(lobj);

    return TCL_OK;
}


/*
 * Replace A,C,G,T with d,e,f,i to mask.
 * Note that when used on the consensus this requires the liberal_base
 * flag in primer3 to be defined in order for it to treat these as Ns
 * instead of errors.
 */
static void tag_mask(char *sequence, int start, int length) {
    int i;

    for (i = 0; i < length; i++) {
	int pos = start-1 + i;
	switch(sequence[pos]) {
	case 'a':
	case 'A':
	    sequence[pos] = 'd';
	    break;
	case 'c':
	case 'C':
	    sequence[pos] = 'e';
	    break;
	case 'g':
	case 'G':
	    sequence[pos] = 'f';
	    break;
	case 't':
	case 'T':
	    sequence[pos] = 'i';
	    break;
	default:
	    sequence[pos] = '-';
	}
    }
}


/*
 * tcl_implement_solutions() - part of the finishing software.
 *
 * Given a the base classification, problem and solution bit-patterns,
 * this finds possible implementations of the solutions.
 */
static int tcl_implement_solutions(finish_t *fin, Tcl_Interp *interp,
				   int objc, Tcl_Obj *CONST objv[]) {
    GAnnotations *curr_tag;

    typedef struct {
	char *prob;
	char *solu;
	char *tag_list;
    } is_args;

    is_args args;
    cli_args a[] = {
	{"-problem_command",  ARG_STR, 1, "",   offsetof(is_args, prob)},
	{"-solution_command", ARG_STR, 1, "",   offsetof(is_args, solu)},
	{"-tag_types",	      ARG_STR, 1, "",   offsetof(is_args, tag_list)},
	{NULL,      0,       0, NULL, 0}
    };
 
    vfuncheader("implement_solutions");

    /* Parse the arguments */
    if (-1 == gap_parse_obj_args(a, &args, objc, (Tcl_Obj **)objv)) {
	return TCL_ERROR;
    }

    if (*args.prob) {
	if (fin->prob_script)
	    xfree(fin->prob_script);
	fin->prob_script = strdup(args.prob);
    }

    if (*args.solu) {
	if (fin->solu_script)
	    xfree(fin->solu_script);
	fin->solu_script = strdup(args.solu);
    }

    SplitList(args.tag_list, &fin->nskip_tags, &fin->skip_tags);
    
    if (fin->nskip_tags) {
	int rnum;

	/* Consensus tags */
	curr_tag = vtagget(fin->io, -fin->contig,
			   fin->nskip_tags, fin->skip_tags);
	while (curr_tag && curr_tag != (GAnnotations *)-1) {
	    /* Exclude from consensus => don't choose oligos here */
	    char type[5];

	    if (fin->opts.debug[FIN_DEBUG] > 1)
		printf("Excluding %d-%d (cons, type %s)\n", 
		       curr_tag->position,
		       curr_tag->position + curr_tag->length-1,
		       type2str(curr_tag->type, type));

	    tag_mask(fin->cons, curr_tag->position, curr_tag->length);

	    curr_tag = vtagget(fin->io, 0, fin->nskip_tags, fin->skip_tags);
	}

	/* Reading tags */
	for (rnum = io_clnbr(fin->io, fin->contig);
	     rnum;
	     rnum = io_rnbr(fin->io, rnum)) {

	    curr_tag = vtagget(fin->io, rnum, fin->nskip_tags, fin->skip_tags);

	    if (!curr_tag)
		continue;

	    while (curr_tag && curr_tag != (GAnnotations *)-1) {
		/* Exclude from consensus => don't choose oligos here */
		char type[5];
		int start, end;

		consensus_tag_pos(fin, curr_tag, rnum, &start, &end);

		if (fin->opts.debug[FIN_DEBUG] > 1 && start <= end)
		    printf("Excluding %d-%d (seq %d, type %s)\n", 
			   start, end,
			   rnum, type2str(curr_tag->type, type));
		
		tag_mask(fin->cons, start, end-start+1 /* length */);

		curr_tag = vtagget(fin->io, 0,
				   fin->nskip_tags, fin->skip_tags);
	    }
	}
    }

    implement_solutions(interp, fin);
    ckfree((char *)fin->skip_tags);

    return TCL_OK;
}

/*
 * Callback from Tcl when a finish object command is destroyed.
 */
static void tcl_finish_delete(ClientData clientdata) {
    finish_t *fin = (finish_t *)clientdata;

    if (fin->opts.debug[FIN_DEBUG])
	puts("Deleting finish object");

    puts("");
    printf("Total number of long reads:         %d\n",
	   fin->count[EXPERIMENT_LONG]);
    printf("Total number of resequences:        %d\n",
	   fin->count[EXPERIMENT_RESEQ]);
    printf("Total number of vector walks:       %d\n",
	   fin->count[EXPERIMENT_VPWALK]);
    printf("Total number of chromosomal walks:  %d\n",
	   fin->count[EXPERIMENT_CPWALK]);
    printf("Total number of reverse sequences:  %d\n",
	   fin->count[EXPERIMENT_REVERSE]);
    puts("");

    finish_del(fin);
}

/*
 * ---------------------------------------------------------------------------
 * Tcl utility functions
 * ---------------------------------------------------------------------------
 */

/*
 * Parses the bit classification string.
 */
static con_bits_t *parse_con_bits(Tcl_Interp *interp, char *str, int *nbits) {
    con_bits_t *con_bits;
    int bitsc;
    char **bitsv;
    int i;

    /* Parse and build bits structure */
    Tcl_SplitList(interp, str, &bitsc, &bitsv);
    con_bits = (con_bits_t *)xmalloc(bitsc * sizeof(*con_bits));

    for (i = 0; i < bitsc; i++) {
	int argc;
	char **argv;
	if (-1 == Tcl_SplitList(interp, bitsv[i], &argc, &argv))
	    continue;
	if (argc < 2)
	    continue;
	con_bits[i].bit = atoi(argv[0]);
	if (strcmp(argv[1], "strand_top") == 0) {
	    con_bits[i].type = CLASS_STRAND_TOP;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "strand_bottom") == 0) {
	    con_bits[i].type = CLASS_STRAND_BOTTOM;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "sequence_depth_gt") == 0) {
	    con_bits[i].type = CLASS_SEQ_DEPTH_GT;
	    con_bits[i].arg = argc > 2 ? atoi(argv[2]) : 2;

	} else if (strcmp(argv[1], "sequence_depth_ge") == 0) {
	    con_bits[i].type = CLASS_SEQ_DEPTH_GE;
	    con_bits[i].arg = argc > 2 ? atoi(argv[2]) : 3;

	} else if (strcmp(argv[1], "template_depth_gt") == 0) {
	    con_bits[i].type = CLASS_TEMP_DEPTH_GT;
	    con_bits[i].arg = argc > 2 ? atoi(argv[2]) : 1;

	} else if (strcmp(argv[1], "template_depth_ge") == 0) {
	    con_bits[i].type = CLASS_TEMP_DEPTH_GE;
	    con_bits[i].arg = argc > 2 ? atoi(argv[2]) : 2;

	} else if (strcmp(argv[1], "confidence_gt") == 0) {
	    con_bits[i].type = CLASS_CONFIDENCE_GT;
	    con_bits[i].arg = argc > 2 ? atoi(argv[2]) : 14;

	} else if (strcmp(argv[1], "confidence_ge") == 0) {
	    con_bits[i].type = CLASS_CONFIDENCE_GE;
	    con_bits[i].arg = argc > 2 ? atoi(argv[2]) : 15;

	} else if (strcmp(argv[1], "chemistry") == 0) {
	    con_bits[i].type = CLASS_CHEMISTRY;
	    con_bits[i].arg = argc > 2 ? atoi(argv[2]) : 17;

	} else if (strcmp(argv[1], "contig_left_end") == 0) {
	    con_bits[i].type = CLASS_CONTIG_LEFT_END;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "contig_right_end") == 0) {
	    con_bits[i].type = CLASS_CONTIG_RIGHT_END;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "low_complexity") == 0) {
	    con_bits[i].type = CLASS_LOW_COMPLEXITY;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_A") == 0) {
	    con_bits[i].type = CLASS_POLY_A;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_C") == 0) {
	    con_bits[i].type = CLASS_POLY_C;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_G") == 0) {
	    con_bits[i].type = CLASS_POLY_G;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_T") == 0) {
	    con_bits[i].type = CLASS_POLY_T;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_K") == 0) {
	    con_bits[i].type = CLASS_POLY_K;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_M") == 0) {
	    con_bits[i].type = CLASS_POLY_M;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_R") == 0) {
	    con_bits[i].type = CLASS_POLY_R;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_S") == 0) {
	    con_bits[i].type = CLASS_POLY_S;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_W") == 0) {
	    con_bits[i].type = CLASS_POLY_Y;
	    con_bits[i].arg = 0;

	} else if (strcmp(argv[1], "poly_Y") == 0) {
	    con_bits[i].type = CLASS_POLY_Y;
	    con_bits[i].arg = 0;

	} else {
	    verror(ERR_WARN, "classify_bases", "Unknown class type '%s'",
		   argv[1]);
	}
	Tcl_Free((char *)argv);
    }
    Tcl_Free((char *)bitsv);

    *nbits = bitsc;
    return con_bits;
}


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
			      int len)
{
    int i;
    Tcl_Obj *objv[2];
    unsigned int *probs;

    if (!script)
	return NULL;

    probs = (unsigned int *)xmalloc(len * sizeof(*probs));
    if (NULL == probs)
	return NULL;

    objv[0] = Tcl_NewStringObj(script, -1);
    objv[1] = Tcl_NewIntObj(0);
    Tcl_IncrRefCount(objv[0]);
    Tcl_IncrRefCount(objv[1]);

    for (i = 0; i < len; i++) {
	if (fin->tag_mask &&
	    mask_offset + i < fin->length &&
	    fin->tag_mask[mask_offset+i]) {
	    probs[i] = 0;
	} else {
	    Tcl_SetIntObj(objv[1], classbits[i]);
	    Tcl_EvalObjv(interp, 2, objv, 0);
	    Tcl_GetIntFromObj(interp, Tcl_GetObjResult(interp),
			      (int *)&probs[i]);
	}
    }

    Tcl_DecrRefCount(objv[0]);
    Tcl_DecrRefCount(objv[1]);

    return probs;
}

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
				  int len)
{
    int i;
    Tcl_Obj *objv[3];
    unsigned int *soln;

    soln = (unsigned int *)xmalloc(len * sizeof(*soln));
    if (NULL == soln)
	return NULL;

    objv[0] = Tcl_NewStringObj(solu_script, -1);
    objv[1] = Tcl_NewIntObj(0);
    objv[2] = Tcl_NewIntObj(1);
    Tcl_IncrRefCount(objv[0]);
    Tcl_IncrRefCount(objv[1]);
    Tcl_IncrRefCount(objv[2]);

    for (i = 0; i < len; i++) {
	Tcl_SetIntObj(objv[1], classbits[i]);
	Tcl_SetIntObj(objv[2], probbits[i]);
	Tcl_EvalObjv(interp, 3, objv, 0);
	Tcl_GetIntFromObj(interp, Tcl_GetObjResult(interp), (int *)&soln[i]);
    }

    Tcl_DecrRefCount(objv[0]);
    Tcl_DecrRefCount(objv[1]);
    Tcl_DecrRefCount(objv[2]);

    return soln;
}

