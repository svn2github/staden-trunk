#include <tcl.h>
#include <tk.h>
#include <string.h>
#include <ctype.h>

#include "seq_raster.h"
#include "sip_globals.h"
#include "sip_structs.h"
#include "cli_arg.h"
#include "seq_results.h"
#include "sip_results.h"
#include "text_output.h"
#include "align.h"
#include "sip_align.h"
#include "dna_utils.h"
#include "sip_hash.h"
#include "sip_find_identity.h"
#include "sequence_formats.h"
#include "probs.h"
#include "readpam.h"
#include "sip_similar_spans.h"
#include "sip_quick_scan.h"
#include "tcl_utils.h"
#include "sip_sendto.h"
#include "sip_sim.h"
#include "sequtils.h"
#include "tkRaster.h"

/*
 * display the expected and observed scores
 */
int tcl_sip_find_probs(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{    
    find_prob_arg args;
    char *seq1;
    char *seq2;
    int seq1_len, seq2_len;
    int seq1_type, seq2_type;
    int seq1_num, seq2_num;

    cli_args a[] = {
	{"-win_len",  ARG_INT, 1, NULL, offsetof(find_prob_arg, win_len)},
	{"-seq_id_h", ARG_INT, 1, NULL, offsetof(find_prob_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL, offsetof(find_prob_arg, seq_id_v)},
	{"-start_h",  ARG_INT, 1, "1",  offsetof(find_prob_arg, start_h)},
	{"-end_h",    ARG_INT, 1, NULL, offsetof(find_prob_arg, end_h)},
	{"-start_v",  ARG_INT, 1, "1",  offsetof(find_prob_arg, start_v)},
	{"-end_v",    ARG_INT, 1, NULL, offsetof(find_prob_arg, end_v)},
	{"-type_h",   ARG_INT, 1, "-1", offsetof(find_prob_arg, type_h)},
	{"-type_v",   ARG_INT, 1, "-1", offsetof(find_prob_arg, type_v)},
	{"-use_av_comp", ARG_INT, 1, "0", offsetof(find_prob_arg, use_av_comp)},
	{NULL,       0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* get first and second sequence saved using extract_sequence */
    seq1_num = GetSeqNum(args.seq_id_h);
    seq2_num = GetSeqNum(args.seq_id_v);

    if (seq1_num == -1) {
	verror(ERR_WARN, "find probabilities", "horizontal sequence undefined");
	return TCL_OK;
    } else if (seq2_num == -1) {
	verror(ERR_WARN, "find probabilities", "vertical sequence undefined");
	return TCL_OK;
    }

    seq1 = GetSeqSequence(seq1_num);
    seq1_len = GetSeqLength(seq1_num);
    seq2 = GetSeqSequence(seq2_num);
    seq2_len = GetSeqLength(seq2_num);

    if (args.start_h < 1)
	args.start_h = 1;

    if (args.end_h > seq1_len)
	args.end_h = seq1_len;

    if (args.start_v < 1)
	args.start_v = 1;

    if (args.end_v > seq2_len)
	args.end_v = seq2_len;

    if (args.type_h == -1) 
	seq1_type = GetSeqType(seq1_num);
    else 
	seq1_type = args.type_h;

    if (args.type_v == -1) 
	seq2_type = GetSeqType(seq2_num);
    else
	seq2_type = args.type_v;

    if (args.use_av_comp) {
	seq1_type = PROTEIN;
	seq2_type = PROTEIN;
    }

    if (seq1_type != seq2_type) {
	verror(ERR_WARN, "find score", "sequences must both be either DNA or protein");
	return TCL_OK;
    } else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
        set_score_matrix(get_matrix_file(PROTEIN));
    } else if (seq1_type == DNA) {
	set_char_set(DNA);
        set_score_matrix(get_matrix_file(DNA));
    }

    FindProbs(seq1, seq2, args.start_h, args.end_h, args.start_v, args.end_v, 
	      args.win_len, seq1_type, args.use_av_comp);
   
    return TCL_OK;
}

/*
 * finds the expected score to give num_matches for a given window length
 */
int tcl_sip_find_score(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{    
    find_score_arg args;
    char *seq1;
    char *seq2;
    int seq1_len, seq2_len;
    int seq1_type, seq2_type;
    int seq1_num, seq2_num;
    int score;

    cli_args a[] = {
	{"-win_len",     ARG_INT, 1, NULL, offsetof(find_score_arg, win_len)},
	{"-num_matches", ARG_INT, 1, NULL, 
	     offsetof(find_score_arg, num_matches)},
	{"-seq_id_h", ARG_INT, 1, NULL, offsetof(find_score_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL, offsetof(find_score_arg, seq_id_v)},
	{"-start_h",  ARG_INT, 1, "1",  offsetof(find_score_arg, start_h)},
	{"-end_h",    ARG_INT, 1, NULL, offsetof(find_score_arg, end_h)},
	{"-start_v",  ARG_INT, 1, "1",  offsetof(find_score_arg, start_v)},
	{"-end_v",    ARG_INT, 1, NULL, offsetof(find_score_arg, end_v)},
	{"-use_av_comp", ARG_INT, 1, "0", offsetof(find_score_arg, use_av_comp)},
	{NULL,       0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq1_num = GetSeqNum(args.seq_id_h);
    seq2_num = GetSeqNum(args.seq_id_v);

    seq1 = GetSeqSequence(seq1_num);
    seq2 = GetSeqSequence(seq2_num);
    seq1_type = GetSeqType(seq1_num);
    seq2_type = GetSeqType(seq2_num);
    seq1_len = GetSeqLength(seq1_num);
    seq2_len = GetSeqLength(seq2_num);

    if (args.start_h < 1)
	args.start_h = 1;

    if (args.end_h > seq1_len)
	args.end_h = seq1_len;

    if (args.start_v < 1)
	args.start_v = 1;

    if (args.end_v > seq2_len)
	args.end_v = seq2_len;

    seq1_len = args.end_h - args.start_h + 1;
    seq2_len = args.end_v - args.start_v + 1;

    if (args.use_av_comp) {
	seq1_len = seq1_len / 3;
	seq2_len = seq2_len / 3;
	seq1_type = PROTEIN;
	seq2_type = PROTEIN;
    }

    if (seq1_type != seq2_type) {
	verror(ERR_WARN, "find score", "sequences must both be either DNA or protein");
	return TCL_OK;
    } else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
        set_score_matrix(get_matrix_file(PROTEIN));
    } else if (seq1_type == DNA) {
	set_char_set(DNA);
        set_score_matrix(get_matrix_file(DNA));
    }

    score = FindScore(seq1_len, seq2_len, args.win_len, args.num_matches);
    vTcl_SetResult(interp, "%d", score);
   
    return TCL_OK;
}

/* HACK should go somewhere else!!!
 * returns the position in the sequence of the first base ie ignores leading
 * padding chars
 */
int 
find_first_base(char *seq,
		int seq_length,
		char pad)
{
    int i = 0;

    while ((seq[i] == pad) && (i < seq_length)) {
	i++;
    }
    return i;
}

int sip_list(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    list_arg args;
    seq_result *result;
    int seq_num;
    int num_id;
    char **result_id;
    int i;

    cli_args a[] = {
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(list_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(list_arg, result_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return -1;

    if (Tcl_SplitList(interp, args.result_id, &num_id, &result_id) != TCL_OK) {
	return -1;
    }

    seq_num = GetSeqNum(args.seq_id);;
    for (i = 0; i < num_id; i++) {
	result = result_data(atoi(result_id[i]), seq_num);
	result->txt_func(result);
    }
    Tcl_Free((char *)result_id);
    return 0;
}

int sip_similar_spans_create(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc, 
			     char *argv[]) 
{
    similar_spans_arg args;
    int id;

    cli_args a[] = {
	{"-seq_id_h", ARG_INT, 1, NULL,  offsetof(similar_spans_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL,  offsetof(similar_spans_arg, seq_id_v)},
	{"-win_len",  ARG_INT, 1, "15", offsetof(similar_spans_arg, win_len)},
	{"-min_match",ARG_INT, 1, "13", offsetof(similar_spans_arg,min_match)},
	{"-start_h",  ARG_INT, 1, "1",  offsetof(similar_spans_arg, start_h)},
	{"-end_h",    ARG_INT, 1, "-1", offsetof(similar_spans_arg, end_h)},
	{"-start_v",  ARG_INT, 1, "1",  offsetof(similar_spans_arg, start_v)},
	{"-end_v",    ARG_INT, 1, "-1", offsetof(similar_spans_arg, end_v)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;
 
    if (-1 == init_sip_similar_spans_create(interp, args.seq_id_h,
					    args.seq_id_v, args.start_h,
					    args.end_h, args.start_v,
					    args.end_v, args.win_len,
					    args.min_match, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    
    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int sip_similar_spans_plot(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_v)},
	{"-result_id", ARG_INT, 1, NULL, offsetof(plot_arg, result_id)},
	{"-window",    ARG_STR, 1, NULL, offsetof(plot_arg, raster)},
	{"-raster_id", ARG_INT, 1, NULL, offsetof(plot_arg, raster_id)},
	{"-fill",      ARG_STR, 1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT, 1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	printf("failure in sip_similar_spans_plot\n");
	return TCL_ERROR;
    }
    
    if (-1 == init_sip_similar_spans_plot(interp, args.seq_id_h, args.seq_id_v,
					  args.result_id,
					  args.raster, args.raster_id,
					  args.colour, args.line_width)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    return TCL_OK;
}

int tcl_sip_similar_spans(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      sip_similar_spans_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      sip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      sip_similar_spans_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int sip_global_align_create(ClientData clientData, 
			    Tcl_Interp *interp, 
			    int argc, 
			    char *argv[]) 
{
    align_seqs_arg args;
    int id;

    cli_args a[] = {
	{"-seq_id_h", ARG_INT, 1, NULL, offsetof(align_seqs_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL, offsetof(align_seqs_arg, seq_id_v)},
	{"-start_h",  ARG_INT, 1, "1",  offsetof(align_seqs_arg, start_h)},
	{"-end_h",    ARG_INT, 1, "-1", offsetof(align_seqs_arg, end_h)},
	{"-start_v",  ARG_INT, 1, "1",  offsetof(align_seqs_arg, start_v)},
	{"-end_v",    ARG_INT, 1, "-1", offsetof(align_seqs_arg, end_v)},
	{"-match",    ARG_INT, 1, NULL, offsetof(align_seqs_arg, match)},
	{"-mismatch", ARG_INT, 1, NULL, offsetof(align_seqs_arg, mismatch)},
	{"-start_gap",ARG_INT, 1, NULL, offsetof(align_seqs_arg, start_gap)},
	{"-cont_gap", ARG_INT, 1, NULL, offsetof(align_seqs_arg, cont_gap)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "Global alignment", "failure to parse input\n");
	return TCL_OK;
    }

    if (-1 == init_sip_global_align_create(interp, args.seq_id_h,
					   args.seq_id_v, args.start_h,
					   args.end_h, args.start_v,
					   args.end_v, args.match, 
					   args.mismatch, args.start_gap,
					   args.cont_gap, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    
    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int sip_global_align_plot(ClientData clientData, 
			    Tcl_Interp *interp, 
			    int argc, 
			    char *argv[]) 
{
    plot_arg args;
    
    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_v)},
	{"-result_id", ARG_INT, 1, NULL, offsetof(plot_arg, result_id)},
	{"-window",    ARG_STR, 1, NULL, offsetof(plot_arg, raster)},
	{"-raster_id", ARG_INT, 1, NULL, offsetof(plot_arg, raster_id)},
	{"-fill",      ARG_STR, 1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT, 1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	printf("failure in sip_global_align_plot\n");
	return TCL_ERROR;
    }
    
    if (-1 == init_sip_global_align_plot(interp, args.seq_id_h, args.seq_id_v,
					 args.result_id, args.raster, 
					 args.raster_id, args.colour, 
					 args.line_width)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    return TCL_OK;
}

int tcl_sip_global_align(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      sip_global_align_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      sip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      sip_global_align_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int sip_matching_words_create(ClientData clientData, 
			      Tcl_Interp *interp, 
			      int argc, 
			      char *argv[]) 
{
    identity_arg args;
    int id;

    cli_args a[] = {
	{"-seq_id_h", ARG_INT, 1, NULL, offsetof(identity_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL, offsetof(identity_arg, seq_id_v)},
	{"-start_h",  ARG_INT, 1, "1", offsetof(identity_arg, start_h)},
	{"-end_h",    ARG_INT, 1, "-1", offsetof(identity_arg, end_h)},
	{"-start_v",  ARG_INT, 1, "1", offsetof(identity_arg, start_v)},
	{"-end_v",    ARG_INT, 1, "-1", offsetof(identity_arg, end_v)},
	{"-word_len", ARG_INT, 1, NULL, offsetof(identity_arg, word_len)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "Find matching words", "failure to parse input\n");
	return TCL_OK;
    }
 
    if (-1 == init_sip_matching_words_create(interp, args.seq_id_h,
					     args.seq_id_v, args.start_h,
					     args.end_h, args.start_v,
					     args.end_v, args.word_len,
					     &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    
    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int sip_matching_words_plot(ClientData clientData, 
			    Tcl_Interp *interp, 
			    int argc, 
			    char *argv[]) 
{
    plot_arg args;
    
    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_v)},
	{"-result_id", ARG_INT, 1, NULL, offsetof(plot_arg, result_id)},
	{"-window",    ARG_STR, 1, NULL, offsetof(plot_arg, raster)},
	{"-raster_id", ARG_INT, 1, NULL, offsetof(plot_arg, raster_id)},
	{"-fill",      ARG_STR, 1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT, 1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	printf("failure in sip_matching_words_plot\n");
	return TCL_ERROR;
    }
    
    if (-1 == init_sip_matching_words_plot(interp, args.seq_id_h, args.seq_id_v,
					   args.result_id,
					   args.raster, args.raster_id,
					   args.colour, args.line_width)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    return TCL_OK;
}

int tcl_sip_matching_words(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      sip_matching_words_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      sip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      sip_matching_words_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int sip_best_diagonals_create(ClientData clientData, 
			      Tcl_Interp *interp, 
			      int argc, 
			      char *argv[]) 
{
    quick_scan_arg args;
    int id;

    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(quick_scan_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, NULL, offsetof(quick_scan_arg, seq_id_v)},
	{"-start_h",   ARG_INT, 1, "1", offsetof(quick_scan_arg, start_h)},
	{"-end_h",     ARG_INT, 1, "-1", offsetof(quick_scan_arg, end_h)},
	{"-start_v",   ARG_INT, 1, "1", offsetof(quick_scan_arg, start_v)},
	{"-end_v",     ARG_INT, 1, "-1", offsetof(quick_scan_arg, end_v)},
	{"-win_len",   ARG_INT, 1, "15", offsetof(quick_scan_arg,win_len)},
	{"-min_match", ARG_INT, 1, "13", offsetof(quick_scan_arg, min_match)},
	{"-word_len",  ARG_INT, 1, NULL, offsetof(quick_scan_arg, word_len)},
	{"-sd",        ARG_FLOAT, 1, NULL, offsetof(quick_scan_arg, min_sd)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "Find best diagonals", "failure to parse input\n");
	return TCL_OK;
    }
    if (-1 == init_sip_best_diagonals_create(interp, args.seq_id_h,
					     args.seq_id_v, args.start_h,
					     args.end_h, args.start_v,
					     args.end_v, args.win_len,
					     args.min_match, args.word_len,
					     args.min_sd, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    
    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int sip_best_diagonals_plot(ClientData clientData, 
			    Tcl_Interp *interp, 
			    int argc, 
			    char *argv[]) 
{
    plot_arg args;
    
    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_v)},
	{"-result_id", ARG_INT, 1, NULL, offsetof(plot_arg, result_id)},
	{"-window",    ARG_STR, 1, NULL, offsetof(plot_arg, raster)},
	{"-raster_id", ARG_INT, 1, NULL, offsetof(plot_arg, raster_id)},
	{"-fill",      ARG_STR, 1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT, 1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	printf("failure in sip_matching_words_plot\n");
	return TCL_ERROR;
    }
    
    if (-1 == init_sip_best_diagonals_plot(interp, args.seq_id_h, args.seq_id_v,
					   args.result_id,
					   args.raster, args.raster_id,
					   args.colour, args.line_width)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    return TCL_OK;
}

int tcl_sip_best_diagonals(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      sip_best_diagonals_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      sip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      sip_best_diagonals_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

/*
 * set the default DNA or protein score matrix file
 */
int
SetScoreMatrix(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{
    set_score_matrix_arg args;
    cli_args a[] = {
	{"-file",  ARG_STR, 1, NULL, offsetof(set_score_matrix_arg, file)},
	{"-type",  ARG_INT, 1, NULL, offsetof(set_score_matrix_arg, type)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
#ifdef DEBUG
    printf("SET file %s type %d \n", args.file, args.type);
#endif

    vfuncheader("Change score matrix");

    if (strcmp(args.file, "<identity>") == 0) {
	if (-1 == set_matrix_file(NULL, args.type)) {
	    verror(ERR_WARN, "set score matrix", "unable to set the identity"
		   "score matrix");
	} else {
	    vmessage("Current dna score matrix file is %s\n", args.file);
	}

    } else {
	if (-1 == set_matrix_file(args.file, args.type)) {
	    verror(ERR_WARN, "set score matrix", "unable to set the score"
		   "matrix %s", args.file);
	} else {
	    vmessage("Current protein score matrix file is %s\n", args.file);
	}
    }
    return TCL_OK;
}

/* 
 * return the current score matrix for type 'type' to tcl
 */
int
GetScoreMatrix(ClientData clientData, 
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[]) 
{
    get_score_matrix_arg args;
    char *matrix_file;

    cli_args a[] = {
	{"-type",  ARG_INT, 1, NULL, offsetof(get_score_matrix_arg, type)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    matrix_file = get_matrix_name(args.type);

#ifdef DEBUG
    printf("GET file %s type %d \n", matrix_file, args.type);
#endif
    if (matrix_file == NULL) {
	vTcl_SetResult(interp, "<identity>");	    
    } else {
	vTcl_SetResult(interp, "%s", get_matrix_name(args.type));	    
    }
    return TCL_OK;
}

/*
 * set global variable determining whether to replot temporary results
 */
int
SetReplotTemp(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    set_replot_temp(atoi(argv[1]));

    return TCL_OK;
}

/*
 * set global variable determining whether to replot temporary results
 */
int
GetReplotTemp(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    int temp;

    temp = get_replot_temp();
    vTcl_SetResult(interp, "%d", temp);	    
    return TCL_OK;
}

/*
 * set global variable of max num of matches that can be saved
 */
int
SetMaxMatches(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    set_max_matches(atoi(argv[1]));
    return TCL_OK;
}

/*
 * return global variable of max num of matches that can be saved
 */
int
GetMaxMatches(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    vTcl_SetResult(interp, "%d", get_max_matches());
    return TCL_OK;
}

/*
 * set global variable of default num of matches that can be saved
 */
int
SetDefMatches(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    set_def_matches(atoi(argv[1]));
    return TCL_OK;
}

/*
 * return global variable of default num of matches that can be saved
 */
int
GetDefMatches(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    vTcl_SetResult(interp, "%d", get_def_matches());
    return TCL_OK;
}

int
SipFreeAll(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    /* SipFreeResults(); */
  
    return TCL_OK;
}

int
SetRemoveDup(ClientData clientData, 
	   Tcl_Interp *interp, 
	   int argc, 
	   char *argv[]) 
{
    /* 1 = yes, 0 = no */
    set_remove_dup(atoi(argv[1]));
    return TCL_OK;
}

int
GetRemoveDup(ClientData clientData, 
	   Tcl_Interp *interp, 
	   int argc, 
	   char *argv[]) 
{
    /* 1 = yes, 0 = no */

    vTcl_SetResult(interp, "%d", get_remove_dup());
    return TCL_OK;
}

/* 
 * Local alignment using the Smith-Waterman algorithm and the sim.c code of
 * Huang.
 * There are 2 modes of operation, 
 * 1) finding the best num_align alignments 
 * 2) finding all alignments above a certain score. In this case, we need to 
 * set num_align to be a large number and then stop the search when then
 * score drops below score_align
 */
int sip_local_align_create(ClientData clientData, 
			      Tcl_Interp *interp, 
			      int argc, 
			      char *argv[]) 
{
    sim_arg args;
    int id;

    cli_args a[] = {
	{"-seq_id_h",     ARG_INT, 1, NULL,  offsetof(sim_arg, seq_id_h)},
	{"-seq_id_v",     ARG_INT, 1, NULL,  offsetof(sim_arg, seq_id_v)},
	{"-start_h",      ARG_INT, 1, "1",   offsetof(sim_arg, start_h)},
	{"-end_h",        ARG_INT, 1, NULL,  offsetof(sim_arg, end_h)},
	{"-start_v",      ARG_INT, 1, "1",   offsetof(sim_arg, start_v)},
	{"-end_v",        ARG_INT, 1, NULL,  offsetof(sim_arg, end_v)},
	{"-num_alignments",ARG_INT, 1, "1",   offsetof(sim_arg, num_align)},
	{"-score_alignment",ARG_FLOAT, 1, "-1",offsetof(sim_arg, score_align)},
	{"-match",        ARG_FLOAT, 1, "1",   offsetof(sim_arg, match)},
	{"-transition",   ARG_FLOAT, 1, "-1",  offsetof(sim_arg, transition)},
	{"-transversion", ARG_FLOAT, 1, "-1", offsetof(sim_arg, transversion)},
	{"-start_gap",    ARG_FLOAT, 1, "6.0", offsetof(sim_arg, start_gap)},
	{"-cont_gap",     ARG_FLOAT, 1, "0.2", offsetof(sim_arg, cont_gap)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "Local alignment", "failure to parse input\n");
	return TCL_OK;
    }
    if (-1 == init_sip_local_align_create(interp, args.seq_id_h,
					  args.seq_id_v, args.start_h,
					  args.end_h, args.start_v,
					  args.end_v, args.num_align,
					  args.score_align, args.match,
					  args.transition, args.transversion,
					  args.start_gap, args.cont_gap, &id)){
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    
    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int sip_local_align_plot(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    plot_arg args;
    
    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_v)},
	{"-result_id", ARG_INT, 1, NULL, offsetof(plot_arg, result_id)},
	{"-window",    ARG_STR, 1, NULL, offsetof(plot_arg, raster)},
	{"-raster_id", ARG_INT, 1, NULL, offsetof(plot_arg, raster_id)},
	{"-fill",      ARG_STR, 1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT, 1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "Local alignment", "failure to parse input\n");
	return TCL_ERROR;
    }
    
    if (-1 == init_sip_local_align_plot(interp, args.seq_id_h, args.seq_id_v,
					args.result_id,
					args.raster, args.raster_id,
					args.colour, args.line_width)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    return TCL_OK;
}

int tcl_sip_local_align(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, 
			char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      sip_local_align_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      sip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      sip_local_align_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

/* set up tcl commands which call C procedures */
/*****************************************************************************/
/*                                 Sip_Init                                  */
/*****************************************************************************/
int
Sip_Init(Tcl_Interp *interp) {
    sip_init_globals(interp);

    Tcl_CreateCommand(interp, "sip_similar_spans", tcl_sip_similar_spans, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sip_find_probs", tcl_sip_find_probs, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sip_find_score", tcl_sip_find_score, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sip_global_align", tcl_sip_global_align, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sip_matching_words", tcl_sip_matching_words, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sip_best_diagonals", tcl_sip_best_diagonals, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "set_score_matrix", SetScoreMatrix, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_score_matrix", GetScoreMatrix, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "set_replot_temp", SetReplotTemp, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_replot_temp", GetReplotTemp, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "set_max_matches", SetMaxMatches, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_max_matches", GetMaxMatches, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "set_def_matches", SetDefMatches, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_def_matches", GetDefMatches, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sip_free_all", SipFreeAll, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "set_remove_dup", SetRemoveDup, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_remove_dup", GetRemoveDup, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sip_local_align", tcl_sip_local_align, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    {
        char *s, c[10];
	/*
	 * Set packages(name). This is done to prevent subsequent reloading
	 * of this library (for efficiency reasons). The only reason that this
	 * is necessary is that currently gap4 dynamically links with some
	 * libraries at link time. When they're all at run time this won't
	 * be necessary.
	 */
	if (s = Tcl_GetVar2(interp, "packages", "sip", TCL_GLOBAL_ONLY))
	    sprintf(c, "%d", atoi(s)|2);
	else
	    strcpy(c, "2");
	Tcl_SetVar2(interp, "packages", "sip", c, TCL_GLOBAL_ONLY);
    }

    return TCL_OK;
}
