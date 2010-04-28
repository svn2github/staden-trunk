#include <staden_config.h>

#include <tcl.h>
#include <tk.h>
#include <string.h>
#include <math.h>

/* added for parsing features */
#include <stdio.h>
#include "codon_content.h"
/* added for parsing features */

#include "misc.h"
#include "xalloc.h"
#include "nip_structs.h"
#include "cli_arg.h"
#include "text_output.h"
#include "tkSeqed.h"
#include "tkSheet_common.h"
#include "seq_raster.h"
#include "tkRaster.h"
#include "seq_results.h"
#include "sequence_formats.h" /* DNA PROTEIN */
#include "dna_utils.h"
#include "base_comp.h"
#include "seqed.h"
#include "nip_base_comp.h"
#include "ruler_tick.h"
#include "nip_globals.h"
#include "nip_results.h"
#include "nip_raster.h"
#include "nip_gene_search.h"
#include "renz_utils.h"
#include "trna_search.h"
#include "nip_stop_codon.h"
#include "nip_trna_search.h"
#include "splice_search.h"
#include "nip_splice_search.h"
#include "tcl_utils.h"
#include "canvas_box.h"
#include "nip_canvas_box.h"
#include "genetic_code.h"
#include "seqed_translate.h"
#include "nip_restriction_enzymes.h"
#include "sequtils.h"
#include "splice_search.h"
#include "nip_string_search.h"
#include "seqed_write.h"
#include "edge.h"
#include "nip_sendto.h"
#include "open_reading_frames.h"
#include "dinuc_freqs.h"
#include "nip_wtmatrix_search.h"

static char *get_subseq(char *dnaseq, Featcds **key_index, int k);
static int get_transl_table(Featcds **key_index);
static int init_genetic_code_ft(Tcl_Interp *interp, int idx);
static char *TranslateSubseq(char *subseq, int rf);

int nip_list(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    s_codon_arg args;
    seq_result *result;
    int seq_num;
    int num_id;
    char **result_id;
    int i;

    cli_args a[] = {
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(s_codon_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(s_codon_arg, result_id)},
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

int nip_base_comp_create(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    pbc_arg args;
    int id;

    cli_args a[] = {
	{"-win_len",   ARG_INT, 1, "31", offsetof(pbc_arg, win_len)},
	{"-a",         ARG_INT, 1, "1",  offsetof(pbc_arg, a)},
	{"-c",         ARG_INT, 1, "0",  offsetof(pbc_arg, c)},
	{"-g",         ARG_INT, 1, "0",  offsetof(pbc_arg, g)},
	{"-t",         ARG_INT, 1, "1",  offsetof(pbc_arg, t)},
	{"-start",     ARG_INT, 1, "1",  offsetof(pbc_arg, start)},
	{"-end",       ARG_INT, 1, "-1", offsetof(pbc_arg, end)},
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(pbc_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;
 
    if (-1 == init_nip_base_comp_create(interp, args.seq_id, args.start, 
					args.end, args.win_len, args.a,
					args.c, args.g, args.t, &id)) {
      vTcl_SetResult(interp, "%d", -1);
      return TCL_OK;
    }

    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int nip_base_comp_plot(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "nip base composition", "failure to parse input\n");
	return TCL_ERROR;
    }
    
    if (-1 == init_nip_base_comp_plot(interp, args.seq_id, 
				      atoi(args.result_id),
				      args.raster, atoi(args.raster_id),
				      args.colour, args.line_width)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    return TCL_OK;
}

int tcl_nip_base_comp(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      nip_base_comp_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_base_comp_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int CountBaseComp(ClientData clientData,
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[])
{
    set_range_arg args;
    char *seq;
    int seq_len, seq_num;
    Tcl_DString input_params;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(set_range_arg, seq_id)},
	{"-start",  ARG_INT, 1, "1",  offsetof(set_range_arg, start)},
	{"-end",    ARG_INT, 1, "-1", offsetof(set_range_arg, end)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("sequence composition");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq_num = GetSeqNum(args.seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = seq_len;
    }
    seq_len = args.end - args.start + 1;

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n",
		       GetSeqName(seq_num), args.start, args.end);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    sequence_info(GetSeqName(seq_num), GetSeqSequence(seq_num),
		  args.start, args.end, GetSeqStructure(seq_num), 
		  GetSeqType(seq_num));
    
    return TCL_OK;
}

int nip_codon_pref_create(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{
    gene_arg args;
    int id[3];

    cli_args a[] = {
	{"-codon_table", ARG_STR, 1, NULL, offsetof(gene_arg, codon_table)},
	{"-win_len",     ARG_INT, 1, "0",  offsetof(gene_arg, win_len)},
	{"-start",       ARG_INT, 1, "1",  offsetof(gene_arg, start)},
	{"-end",         ARG_INT, 1, "-1", offsetof(gene_arg, end)},
	{"-option",      ARG_INT, 1, NULL, offsetof(gene_arg, option)},
	{"-seq_id",      ARG_INT, 1, NULL, offsetof(gene_arg, seq_id)},
	{NULL,           0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	return TCL_ERROR;
    }
    if (-1 == init_nip_codon_pref_create(interp, args.seq_id, args.start, 
					 args.end,
					 args.codon_table, args.win_len,
					 args.option, id)) {
      vTcl_SetResult(interp, "%d %d %d", -1, -1, -1);
      return TCL_OK;
    }

    vTcl_SetResult(interp, "%d %d %d", id[0], id[1], id[2]);
    return TCL_OK;
}

int nip_gene_search_plot(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    plot_arg args;
    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "nip gene search plot", "failure to parse input\n");
	return TCL_ERROR;
    }

    if (-1 == init_nip_gene_search_plot(interp, args.raster, args.raster_id, 
					args.seq_id, args.result_id,
					args.colour, args.line_width))
      return TCL_ERROR;
    return TCL_OK;
}

int tcl_nip_codon_pref(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
	nip_codon_pref_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
	nip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
	nip_gene_search_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int nip_author_test_create(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{
    author_arg args;
    int id[3];

    cli_args a[] = {
	{"-codon_table", ARG_STR, 1, NULL, offsetof(author_arg, codon_table)},
	{"-error",       ARG_DOUBLE, 1, "0.1", offsetof(author_arg, error)},
	{"-start",       ARG_INT, 1, "1",  offsetof(author_arg, start)},
	{"-end",         ARG_INT, 1, "-1", offsetof(author_arg, end)},
	{"-seq_id",      ARG_INT, 1, NULL, offsetof(author_arg, seq_id)},
	{NULL,           0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	return TCL_ERROR;
    }
    if (-1 == init_nip_author_test_create(interp, args.seq_id, args.start, 
					  args.end, args.codon_table, 
					  args.error, id)) {
      vTcl_SetResult(interp, "%d %d %d", -1, -1, -1);
      return TCL_OK;
    }

    vTcl_SetResult(interp, "%d %d %d", id[0], id[1], id[2]);
    return TCL_OK;
    
}

int tcl_nip_author_test(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, 
			char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      nip_author_test_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_gene_search_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int nip_base_bias_create(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    gene_arg args;
    int id;

    cli_args a[] = {
	{"-win_len",    ARG_INT, 1, "0",  offsetof(gene_arg, win_len)},
	{"-start",      ARG_INT, 1, "1",  offsetof(gene_arg, start)},
	{"-end",        ARG_INT, 1, "-1", offsetof(gene_arg, end)},
	{"-seq_id",     ARG_INT, 1, NULL, offsetof(gene_arg, seq_id)},
	{NULL,          0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	return TCL_ERROR;
    }
    if (-1 == init_nip_base_bias_create(interp, args.seq_id, args.start, 
					args.end, args.win_len, &id)) {
      vTcl_SetResult(interp, "%d", -1);
      return TCL_OK;
    }

    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
    
}
int nip_base_bias_plot(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{
    plot_arg args;
    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	return TCL_ERROR;
    }
    if (-1 == init_nip_base_bias_plot(interp, args.raster, 
				      args.raster_id, 
				      args.seq_id, args.result_id,
				      args.colour, args.line_width))
	return TCL_ERROR;
    return TCL_OK;
}

int tcl_nip_base_bias(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      nip_base_bias_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_base_bias_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int ValidCodonTable(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[])
{
    FILE *in_file;
    int res;
    codon_arg args;
    double codon_usage_table[4][4][4];

    cli_args a[] = {
	{"-codon_table", ARG_STR, 1, NULL, offsetof(codon_arg, codon_table)},
	{NULL,          0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
 
    in_file = fopen(args.codon_table,"r");
    if (!in_file) {
	vTcl_SetResult(interp, "%d", 0);
	return TCL_OK;
    }

    res = read_cod_table(in_file, codon_usage_table);
    fclose ( in_file );
    if (!res) {
	vTcl_SetResult(interp, "%d", 0);
	return TCL_OK;
    }
    vTcl_SetResult(interp, "%d", 1);
    return TCL_OK;
}

int nip_trna_search_create(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{
    trna_arg args;
    int id;

    cli_args a[] = {
	{"-start",     ARG_INT, 1, "1",  offsetof(trna_arg, start)},
	{"-end",       ARG_INT, 1, "-1", offsetof(trna_arg, end)},
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(trna_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;

    if (-1 == init_nip_trna_search_create(interp, args.seq_id, args.start,
					  args.end, &id)) {
      vTcl_SetResult(interp, "%d", -1);
      return TCL_OK;
    }

    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int nip_trna_search_plot(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{"-tick_ht",   ARG_FLOAT, 1, "20", offsetof(plot_arg, tick_ht)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "nip trna search", "unable to parse input\n");
	return TCL_ERROR;
    }

    if (-1 == init_nip_trna_search_plot(interp, args.seq_id, 
					atoi(args.result_id),
					args.raster, atoi(args.raster_id),
					args.colour, args.line_width,
					args.tick_ht)) {
      return TCL_ERROR;
    }
    return TCL_OK;
}

int tcl_nip_trna_search(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, 
			char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      nip_trna_search_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_trna_search_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int
nip_stop_codons_create(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{
    s_codon_arg args;
    int id[3];

    cli_args a[] = {
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(s_codon_arg, seq_id)},
	{"-start",      ARG_INT,   1, "1",  offsetof(s_codon_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(s_codon_arg, end)},
	{"-strand",     ARG_STR,   1, "+",  offsetof(s_codon_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;

    init_nip_stop_codons_create(args.seq_id, args.start, args.end,
				args.strand, id);
    vTcl_SetResult(interp, "%d %d %d", id[0], id[1], id[2]);
    return TCL_OK;
}

int
nip_stop_codons_plot(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{"-tick_ht",   ARG_FLOAT, 1, "20", offsetof(plot_arg, tick_ht)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "nip_stop_codons_plot", "failed to parse input\n");
	return TCL_ERROR;
    }

    if (-1 == init_nip_stop_codons_plot(interp, args.raster, args.raster_id,
					args.seq_id, args.result_id,
					args.colour, args.line_width,
					args.tick_ht))
      return TCL_ERROR;

    return TCL_OK;
}

int tcl_nip_stop_codons(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, 
			char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
      nip_stop_codons_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_stop_codons_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int
nip_start_codons_create(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, 
			char *argv[]) 
{
    s_codon_arg args;
    int id[3];

    cli_args a[] = {
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(s_codon_arg, seq_id)},
	{"-start",      ARG_INT,   1, "1",  offsetof(s_codon_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(s_codon_arg, end)},
	{"-strand",     ARG_STR,   1, "+",  offsetof(s_codon_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;

    init_nip_start_codons_create(args.seq_id, args.start, args.end,
				 args.strand, id);

    vTcl_SetResult(interp, "%d %d %d", id[0], id[1], id[2]);
    return TCL_OK;
}

int tcl_nip_start_codons(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{

    if (strcmp(argv[1], "create") == 0) {
      nip_start_codons_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_stop_codons_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int
nip_splice_search_create(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    splice_arg args;
    int id[3];

    cli_args a[] = {
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(splice_arg, seq_id)},
	{"-start",     ARG_INT,   1, "1",  offsetof(splice_arg, start)},
	{"-end",       ARG_INT,   1, "-1", offsetof(splice_arg, end)},
	{"-donor",     ARG_STR,   1, NULL, offsetof(splice_arg, donor)},
	{"-acceptor",  ARG_STR,   1, NULL, offsetof(splice_arg, acceptor)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;

    if (-1 == init_splice_search_create(args.seq_id, args.start, args.end,
					args.donor, args.acceptor, id)) {
	vTcl_SetResult(interp, "%d %d %d", -1, -1, -1);
	return TCL_OK;
    }

    vTcl_SetResult(interp, "%d %d %d", id[0], id[1], id[2]);
    return TCL_OK;
}

int
nip_splice_search_plot(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{"-tick_ht",   ARG_FLOAT, 1, "20", offsetof(plot_arg, tick_ht)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "nip splice search", "failure to parse input\n");
	return TCL_ERROR;
    }
    if (-1 == (init_splice_search_plot(interp, args.raster, 
				       atoi(args.raster_id), 
				       args.result_id, args.seq_id, 
				       args.colour, args.line_width, 
				       args.tick_ht)))
      return TCL_ERROR;

    return TCL_OK;
}

int tcl_splice_search(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{

    if (strcmp(argv[1], "create") == 0) {
      nip_splice_search_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv); 
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_splice_search_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int nip_string_search_create(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc, 
			     char *argv[]) 
{
    string_arg args;
    int id;

    cli_args a[] = {
	{"-strand",     ARG_STR,   1, "+",  offsetof(string_arg, strand)},
	{"-min_pmatch", ARG_FLOAT, 1, "75.",offsetof(string_arg, match)},
	{"-string",     ARG_STR,   1, NULL, offsetof(string_arg, string)},
	{"-use_iub",    ARG_INT,   1, "1", offsetof(string_arg, use_iub_code)},
	{"-start",      ARG_INT,   1, "1",  offsetof(string_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(string_arg, end)},
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(string_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;

    if (-1 == init_nip_string_search_create(args.strand, args.match,
					    args.string, args.use_iub_code, 
					    args.start, args.end, args.seq_id, 
					    &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }

    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;

}

int nip_string_search_plot(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{"-tick_ht",   ARG_FLOAT, 1, "20", offsetof(plot_arg, tick_ht)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "nip string search", "unable to parse input\n");
	return TCL_ERROR;
    }

    if (-1 == init_nip_string_search_plot(interp, args.raster, 
					  atoi(args.raster_id),
					  atoi(args.result_id), args.seq_id,
					  args.colour, args.line_width,
					  args.tick_ht))
      return TCL_ERROR;
    return TCL_OK;
}

int tcl_nip_string_search(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{

    if (strcmp(argv[1], "create") == 0) {
      nip_string_search_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv); 
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_string_search_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int nip_wtmatrix_search_create(ClientData clientData, 
			       Tcl_Interp *interp, 
			       int argc, 
			       char *argv[]) 
{
    wtmatrix_arg args;
    int id;

    cli_args a[] = {
	{"-start",      ARG_INT,   1, "1",  offsetof(wtmatrix_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(wtmatrix_arg, end)},
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(wtmatrix_arg, seq_id)},
	{"-wt_matrix",  ARG_STR,   1, NULL, offsetof(wtmatrix_arg, wt_matrix)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return TCL_ERROR;

    if (-1 == init_nip_wtmatrix_search_create(args.start, args.end, 
					      args.seq_id, args.wt_matrix, 
					      &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int nip_wtmatrix_search_plot(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc, 
			     char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-window",    ARG_STR,   1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_STR,   1, NULL, offsetof(plot_arg, raster_id)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(plot_arg, result_id)},
	{"-fill",      ARG_STR,   1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT,   1, "1",  offsetof(plot_arg, line_width)},
	{"-tick_ht",   ARG_FLOAT, 1, "20", offsetof(plot_arg, tick_ht)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "nip weight matrix search", "failure to parse input\n");
	return TCL_ERROR;
    }

    if (-1 == init_nip_wtmatrix_search_plot(interp, args.seq_id, 
					    atoi(args.result_id),
					    args.raster, atoi(args.raster_id), 
					    args.colour,
					    args.line_width, args.tick_ht))
      return TCL_ERROR;

    return TCL_OK;
}


int tcl_nip_wtmatrix_search(ClientData clientData, 
			    Tcl_Interp *interp, 
			    int argc, 
			    char *argv[]) 
{

    if (strcmp(argv[1], "create") == 0) {
      nip_wtmatrix_search_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
      nip_list(clientData, interp, argc, argv); 
    } else if (strcmp(argv[1], "plot") == 0) {
      nip_wtmatrix_search_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int CodonUsage(ClientData clientData, 
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[]) 
{
    codon_usage_arg args;
    char *seq;
    int seq_len;
    int i;
    int seq_num;
    double codon_table[4][4][4];
    double codon_table_t[4][4][4];
    double codon_table_c[4][4][4];
    int num_ranges;
    int num_range;
    FILE *fp = NULL;
    char *sequence = NULL;
    char **ranges = NULL;
    char **srange = NULL;
    range *c_range = NULL;
    int   retval   = TCL_OK;
    Tcl_DString input_params;

    cli_args a[] = {
	{"-seq_id",  ARG_INT,   1, NULL, offsetof(codon_usage_arg, seq_id)},
	{"-outfile", ARG_STR,   1, "", offsetof(codon_usage_arg, filename)},
	{"-concatenate", ARG_INT,1, "0",offsetof(codon_usage_arg, concat)},
	{"-format",  ARG_INT,   1, "1",  offsetof(codon_usage_arg, format)},
	{"-strand",  ARG_STR,   1, "+",  offsetof(codon_usage_arg, strand)},
	{"-table",   ARG_STR,   1, "",   offsetof(codon_usage_arg, table)},
	{"-range",   ARG_STR,   1, "",   offsetof(codon_usage_arg, range)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vfuncheader("codon usage");

    set_char_set(DNA);

    if (Tcl_SplitList(interp, args.range, &num_ranges, &ranges) != TCL_OK) {
	return TCL_ERROR;
    }

    if (NULL == (c_range = (range *)xmalloc(num_ranges * sizeof(range))))
    {
	retval = TCL_ERROR;
	goto cleanup;
    }
    
    for (i = 0; i < num_ranges; i++) {
	if (Tcl_SplitList(interp, ranges[i], &num_range, &srange) != TCL_OK)
	{
	    retval = TCL_ERROR;
	    goto cleanup;
	}
	if (num_range != 2) {
	    verror(ERR_WARN, "CodonUsage", 
		   "incorrect number of values in range\n");
	    goto cleanup;
	}
	c_range[i].start = atoi(srange[0]);
	c_range[i].end = atoi(srange[1]);
    }

    seq_num = GetSeqNum(args.seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);

    {
	char c_totals[16];
	char strand[8];

	if (args.format == 1) {
	    strcpy(c_totals, "observed counts");
	} else {
	    strcpy(c_totals, "percentage");
	} 

	if (strcmp(args.strand, "+") == 0) {
	    strcpy(strand, "forward");
	} else {
	    strcpy(strand, "reverse");
	}
	vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
			   "save table as %s\nuse existing codon table %s\n"
			   "show codon totals as %s\nstrand %s\n",
			   GetSeqName(seq_num), c_range[0].start, 
			   c_range[0].end,
			   args.filename, args.table, c_totals, strand);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    /* initialise table */
    init_codon_table(codon_table);

    if (args.concat) {
      /* if concatenating, init table_t and read table_c */
	init_codon_table(codon_table_t);

	if (NULL == (fp = fopen(args.table, "r"))) {
	    verror(ERR_WARN, "codon usage", "Unable to open codon usage table %s", args.table);
	    goto cleanup;            
	}
	/* read in table for output later */
	if (0 == read_cod_table(fp, codon_table_c)) {
	  verror(ERR_WARN, "codon usage", "Unable to read codon usage table %s\n", args.table);
	  goto cleanup;
	}
	fclose(fp);
        fp = NULL;
    } else {
      /* not concatenating */
      if (strcmp(args.table, "") == 0) {
	  init_codon_table(codon_table_t);
      } else {
	/* read in table from file */
	if (NULL == (fp = fopen(args.table, "r"))) {
	    verror(ERR_WARN, "codon usage", "Unable to open codon usage table %s", args.table);
	    goto cleanup;
	}
	if (0 == read_cod_table(fp, codon_table_t)) {
	        verror(ERR_WARN, "codon usage", "Unable to read codon usage table %s\n", args.table);
		goto cleanup;
	}
	fclose(fp);
        fp = NULL;
      }
    }

    sequence = strdup(seq);
    /* complement */
    if (strcmp(args.strand, "-") == 0) {
	complement_seq(sequence, seq_len);
    }

    vmessage("Sequence %s\n", GetSeqName(seq_num));

    for (i = 0; i < num_ranges; i++) {
	vmessage("Range from %d to %d \n", c_range[i].start, c_range[i].end);
	calc_codon_usage(&sequence[c_range[i].start-1], 
			 c_range[i].end - c_range[i].start + 1, 
			 codon_table);
	if (args.format == 2) {
	    codon_table_percent(codon_table);
	}
	/* comment out for now because ranges aren't supported yet */
#ifdef FIXME
	write_screen_cod_table(codon_table);
#endif
	calc_codon_usage(&sequence[c_range[i].start-1], 
			 c_range[i].end - c_range[i].start + 1, 
			 codon_table_t);
	init_codon_table(codon_table);
    }
    
    if (args.format == 2) {
	codon_table_percent(codon_table_t);
    }

    if (strcmp(args.filename, "") != 0) {
      if (NULL == (fp = fopen(args.filename, "w"))) {
	verror(ERR_WARN, "codon usage", "Unable to save codon usage table");
	goto cleanup;
      }

      /* write to file */
      if (args.concat) {
	  write_cod_table(fp, codon_table_c);
      } 
      write_cod_table(fp, codon_table_t);
      fclose(fp);
      fp = NULL;
    }

    /* write to text window */
    if (args.concat) {
        write_screen_cod_table(codon_table_c);
    } 
    write_screen_cod_table(codon_table_t);


cleanup:
    if(fp)       fclose(fp);
    if(sequence) xfree(sequence);
    if(c_range)  xfree(c_range);
    if(ranges)   Tcl_Free( (char*) ranges );
    if(srange)   Tcl_Free( (char*) srange );
    return retval;
}

int NipPlotRenz(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    nip_renz_arg args;
    cursor_s cursor;
    tick_s *tick;
    ruler_s *ruler;
    int id;
    out_canvas *output;
    Tcl_DString input_params;

    cli_args a[] = {
	{"-file",	ARG_STR, 1, NULL, offsetof(nip_renz_arg, filename)},
	{"-frame",	ARG_STR, 1, NULL, offsetof(nip_renz_arg, frame)},
	{"-win_names",	ARG_STR, 1, NULL, offsetof(nip_renz_arg, win_name)},
	{"-window",	ARG_STR, 1, NULL, offsetof(nip_renz_arg, plot)},
	{"-win_ruler",	ARG_STR, 1, NULL, offsetof(nip_renz_arg, win_ruler)},
	{"-enzymes",	ARG_STR, 1, NULL, offsetof(nip_renz_arg, inlist)},
	{"-num_enzymes",ARG_INT, 1, NULL, offsetof(nip_renz_arg, num_items)},
	{"-text_offset",ARG_INT, 1, NULL, offsetof(nip_renz_arg, text_offset)},
	{"-text_fill",  ARG_STR, 1, NULL, offsetof(nip_renz_arg, text_fill)},
	{"-tick_height",ARG_INT, 1, "-1", offsetof(nip_renz_arg, tick_ht)},
	{"-tick_width", ARG_INT, 1, "-1", offsetof(nip_renz_arg, tick_wd)},
	{"-tick_fill",  ARG_STR, 1,   "", offsetof(nip_renz_arg, tick_fill)},
	{"-cursor_width",ARG_INT,1, "-1", offsetof(nip_renz_arg, cursor_wd)},
	{"-cursor_fill", ARG_STR, 1,  "", offsetof(nip_renz_arg, cursor_fill)},
	{"-yoffset",	 ARG_INT, 1, NULL, offsetof(nip_renz_arg, yoffset)},
	{"-seq_id",	 ARG_INT, 1, NULL, offsetof(nip_renz_arg, seq_id)},
	{"-start",	 ARG_INT, 1, "1",  offsetof(nip_renz_arg, start)},
	{"-end",	 ARG_INT, 1, "-1", offsetof(nip_renz_arg, end)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (output = (out_canvas *)xmalloc(sizeof(out_canvas))))
	return TCL_OK;

    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = GetSeqLength(GetSeqNum(args.seq_id));
    }

    vfuncheader("restriction enzyme plot");
    set_char_set(DNA);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "enzymes: %s\n",
		       GetSeqName(GetSeqNum(args.seq_id)), args.start, 
		       args.end,
		       args.inlist);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
    
    cursor = cursor_struct(interp, tk_utils_defs, "R_ENZ", args.cursor_wd, 
			   args.cursor_fill);

    tick = tick_struct(interp, tk_utils_defs, "R_ENZ", args.tick_wd, args.tick_ht, 
		       args.tick_fill);

    ruler = ruler_struct(interp, tk_utils_defs, "R_ENZ", 0);

    ruler->start = args.start;
    ruler->end = args.end;
    strcpy(ruler->window, args.win_ruler);

    output->interp = interp;
    id = nip_renz_reg(interp, args.seq_id, output, args.filename, args.frame, 
		      args.win_name, args.plot, args.inlist, args.num_items,
		      args.start, args.end, args.text_offset, args.text_fill,
		      tick, args.yoffset, ruler, cursor);

    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
    
}

int NipTranslate(ClientData clientData,
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[])
{
    nip_trans_arg args;
    char *seq;
    int seq_len, prot_len;
    int i, j, k;
    int frame;
    int width;
    int seq_num;
    Tcl_DString input_params;
    int transl_table_number;
    int num_selcds, num_cds;
    char *protein_id;
    char buffer[1024];
    char *line     = NULL;
    char *sequence = NULL;
    char *prot_seq = NULL;
    char **selcds  = NULL;  
    char *subseq   = NULL;
    int retval     = TCL_ERROR;


    /* added for parsing feature tables */

    cli_args a[] = {
	{"-seq_id",      ARG_INT, 1, NULL, offsetof(nip_trans_arg, seq_id)},
	{"-start",       ARG_INT, 1, "1",  offsetof(nip_trans_arg, start)},
	{"-end",         ARG_INT, 1, "-1", offsetof(nip_trans_arg, end)},
	{"-line_length", ARG_INT, 1, NULL, offsetof(nip_trans_arg, line_length)},
	{"-size",     ARG_INT, 1, NULL, offsetof(nip_trans_arg, size)},
	{"-feat",     ARG_INT, 1, NULL, offsetof(nip_trans_arg, feat)},
	{"-selcds",     ARG_STR, 1, NULL, offsetof(nip_trans_arg, selcds)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    vfuncheader("translation");

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
   
    seq_num = GetSeqNum(args.seq_id);
    seq     = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    if (args.feat==2) { /* use selected range */
    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = seq_len;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "line length %d display as %d letter use %d(1 for feature table and 2 for entry box)\n",
		       GetSeqName(seq_num), args.start, args.end,
		       args.line_length, args.size, args.feat);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
    
    if (NULL == (sequence = (char *)xmalloc((seq_len+5) * sizeof(char))))
	goto end;
    
    if (NULL == (line = (char *)xmalloc((args.line_length+4) * sizeof(char))))
	goto end;
    
    sequence[0] = '-';
    sequence[1] = '-';
    strcpy(&(sequence[2]), seq);
    sequence[seq_len+2] = '-';
    sequence[seq_len+3] = '-';
    sequence[seq_len+4] = '\0';

    frame = 1;
    for (j = args.start; j < args.end; j += args.line_length) {
	if (args.end - j + 1 < args.line_length) {
	    width = args.end - j + 1;
	} else {
	    width = args.line_length;
	}
	for (k = 1; k < 4; k++) {
	    line[0] = ' ';
	    seqed_write_translation(&sequence[j-1], k, args.size, j, width, 1, &line[1]);
	    vmessage("%s\n", line);
	}
	line[0] = ' ';
	seqed_write_sequence(&sequence[j+1], j+1, width, &line[1]);
	vmessage("%s\n", line);

	seqed_write_ruler(j, width, &line[1]);
	vmessage("%s\n", line);
	
	seqed_write_complement(&sequence[j+1], j+1, width, &line[1]);
	vmessage("%s\n", line);

	for (k = 4; k < 7; k++) {
	    line[0] = ' ';
	    seqed_write_translation(&sequence[j-1], k, args.size, j, width, 1, &line[1]);
	    vmessage("%s\n", line);
	}
	vmessage("\n");
    }
    retval = TCL_OK;
    goto end;
    }
    if (args.feat==1){ /* use feature table */
      
      if ( GetSeqKeyIndex(seq_num) == NULL){
	verror(ERR_WARN, "Translation", "Error in translation\n");
	goto end;
      }
      transl_table_number = get_transl_table(GetSeqKeyIndex(seq_num));

      if( init_genetic_code_ft(interp, transl_table_number) )
	goto end;

      set_dna_lookup();  

      /* create selecing CDS array */
      if (Tcl_SplitList(interp, args.selcds, &num_selcds, &selcds) != TCL_OK)
	goto end;

      num_cds = GetSeqKeyIndex(seq_num)[0]->id;
      for(i = 1; i <= num_cds; i++){
	for(j = 0; j < num_selcds; j++){
	  if(!strcmp(selcds[j], GetSeqCdsExpr(seq_num, i))){  
	    subseq = get_subseq(seq, GetSeqKeyIndex(seq_num),i);
	    prot_seq = TranslateSubseq(subseq, 0);
	    if (!prot_seq) {
		xfree(subseq);
		goto end;
	    }
            protein_id = GetSeqProteinId(seq_num, i);
            if( protein_id ) {
		strcpy( buffer, protein_id );                
		buffer[strlen(buffer)-1] = 0;
		vmessage( ">%s\n", &buffer[13] );
	    } else {
	        vmessage(">UNKNOWN\n");
	    }
	    prot_len = strlen(prot_seq);
	    for (k=0; k<prot_len; k += args.line_length) {
	      vmessage( "%.*s\n", args.line_length, &prot_seq[k] );
	    }
	    xfree(subseq);
	    xfree(prot_seq);
	  }
	}
      }
      retval = TCL_OK;
    } 


end:
    if(line)     xfree(line);
    if(sequence) xfree(sequence);
    if(selcds)   Tcl_Free((char*)selcds);
    return retval;
}

int TranslateORFToFeatureTable(ClientData clientData,
			       Tcl_Interp *interp, 
			       int argc, 
			       char *argv[])
{
    trans_ft_arg args;
    char *seq;
    int seq_len, seq_num;
    Tcl_DString input_params;

    cli_args a[] = {
	{"-seq_id",  ARG_INT, 1, NULL, offsetof(trans_ft_arg, seq_id)},
	{"-start",   ARG_INT, 1, "1",  offsetof(trans_ft_arg, start)},
	{"-end",     ARG_INT, 1, "-1", offsetof(trans_ft_arg, end)},
	{"-min_orf", ARG_INT, 1, "30", offsetof(trans_ft_arg, min_orf)},
	{"-strand",  ARG_STR, 1, "+",  offsetof(trans_ft_arg, strand)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("Translate open reading frames to protein: write as feature table");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq_num = GetSeqNum(args.seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = seq_len;
    }

   /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    { 
        char strand[8];
	if (strcmp(args.strand, "+") == 0) {
	    strcpy(strand, "forward");
	} else if (strcmp(args.strand, "-") == 0) {
            strcpy(strand, "reverse");
	} else {
	    strcpy(strand, "both");
	}
	vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "strand %s minimum ORF in codons %d\n",
		       GetSeqName(seq_num), args.start, args.end,
		       strand, args.min_orf);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
  
    vmessage("Sequence %s\n", GetSeqName(seq_num));
    if (strcmp(args.strand, "+") == 0 || strcmp(args.strand, "both") == 0) {
	write_screen_open_frames_f_ft(seq, seq_len, args.start, args.end,
				     args.min_orf);
    }

    if (strcmp(args.strand, "-") == 0 || strcmp(args.strand, "both") == 0) {
	write_screen_open_frames_r_ft(seq, seq_len, args.start, args.end,
				      args.min_orf);
    }
    
    return TCL_OK;
}

int TranslateORFToFasta(ClientData clientData,
			Tcl_Interp *interp, 
			int argc, 
			char *argv[])
{
    trans_fasta_arg args;
    char *seq;
    int seq_len, seq_num;
    FILE *fp;
    Tcl_DString input_params;

    cli_args a[] = {
	{"-seq_id",   ARG_INT, 1, NULL, offsetof(trans_fasta_arg, seq_id)},
	{"-start",    ARG_INT, 1, "1",  offsetof(trans_fasta_arg, start)},
	{"-end",      ARG_INT, 1, "-1", offsetof(trans_fasta_arg, end)},
	{"-min_orf",  ARG_INT, 1, "30", offsetof(trans_fasta_arg, min_orf)},
	{"-strand",   ARG_STR, 1, "+",  offsetof(trans_fasta_arg, strand)},
	{"-filename", ARG_STR, 1, NULL, offsetof(trans_fasta_arg, filename)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("Translate open reading frames to protein: write as fasta file");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq_num = GetSeqNum(args.seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = seq_len;
    }

    if (NULL == (fp = fopen(args.filename, "w"))) {
	verror(ERR_WARN, "Translate open reading frames to protein", 
	       "Unable to open file %s\n", args.filename);
	return TCL_OK;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    { 
        char strand[8];
	if (strcmp(args.strand, "+") == 0) {
	    strcpy(strand, "forward");
	} else if (strcmp(args.strand, "-") == 0) {
            strcpy(strand, "reverse");
	} else {
	    strcpy(strand, "both");
	}
	vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "strand %s minimum ORF in codons %d fasta filename %s\n",
		       GetSeqName(seq_num), args.start, args.end,
		       strand, args.min_orf, args.filename);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    vmessage("Sequence %s\n", GetSeqName(seq_num));
    if (strcmp(args.strand, "+") == 0 || strcmp(args.strand, "both") == 0) {
	write_screen_open_frames_f(seq, seq_len, args.start, args.end,
				   args.min_orf);
	write_open_frames_f(fp, seq, seq_len, args.start, args.end,
			    args.min_orf);
    }

    if (strcmp(args.strand, "-") == 0 || strcmp(args.strand, "both") == 0) {
	write_screen_open_frames_r(seq, seq_len, args.start, args.end,
				   args.min_orf);
	write_open_frames_r(fp, seq, seq_len, args.start, args.end,
			    args.min_orf);
    }
    fclose(fp);
    return TCL_OK;
}

int tcl_load_genetic_code(ClientData clientData, Tcl_Interp *interp,
			  int argc, char *argv[])
{
    read_enz_arg args;
    FILE *fp;

    cli_args a[] = {
	{"-filename",	ARG_STR, 1, NULL, offsetof(read_enz_arg, filename)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (fp = fopen(args.filename, "r"))) {
	Tcl_SetResult(interp, "unable to open file\n", TCL_STATIC);
	return TCL_ERROR;
    }

    if (0 == read_global_genetic_code(fp)) {
	verror(ERR_WARN, "load_genetic_code",
	       "Could not read genetic code. Using standard code.");
	init_genetic_code();
	vTcl_SetResult(interp, "%d", -1);
    } else {
	vTcl_SetResult(interp, "%d", 0);
    }
    fclose(fp);
    

    return TCL_OK;
}

int CountDinucFreq(ClientData clientData,
		   Tcl_Interp *interp, 
		   int argc, 
		   char *argv[])
{
    set_range_arg args;
    double obs_freqs[5][5];
    double exp_freqs[5][5];
    int seq_num;
    char *seq;
    char base[]="ACGT";
    int i,j;
    Tcl_DString input_params;
    

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(set_range_arg, seq_id)},
	{"-start",  ARG_INT, 1, "1",  offsetof(set_range_arg, start)},
	{"-end",    ARG_INT, 1, "-1", offsetof(set_range_arg, end)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("count dinucleotide frequencies");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq_num = GetSeqNum(args.seq_id);
    seq = GetSeqSequence(seq_num);
    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = GetSeqLength(seq_num);
    }

   /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n",
		       GetSeqName(seq_num), args.start, args.end);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);


    calc_dinuc_freqs (seq, args.start, args.end, obs_freqs);
    calc_expected_dinuc_freqs (seq, args.start, args.end, exp_freqs);

    vmessage("Sequence %s\n", GetSeqName(seq_num));

    vmessage("        A                C                G                T\n");
    vmessage("     Obs    Expected  Obs    Expected  Obs    Expected  Obs    Expected\n");
    for (i = 0; i < 4; i++ ) {
	vmessage(" %c", base[i]);
	for (j = 0; j < 4; j++ ) {
	    vmessage("  %7.2f %7.2f", obs_freqs[i][j], exp_freqs[i][j]);
	}
	vmessage("\n");
    }

    return TCL_OK;
}


int
NipScrollCanvas(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    nip_scroll_arg args;
    seq_result *result;
    seq_reg_info info;
    renz_res *data;

    cli_args a[] = {
	{"-id",             ARG_INT, 1, NULL, offsetof(nip_scroll_arg, id)},
	{"-xscrollcommand", ARG_STR, 1, "", offsetof(nip_scroll_arg, xscroll)},
	{"-yscrollcommand", ARG_STR, 1, "", offsetof(nip_scroll_arg, yscroll)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    if (strcmp(args.xscroll, "") != 0) {
	canvasScrollX(interp, data->re_win, data->win_list, 
		      data->num_wins,
		      data->world->visible, data->canvas, args.xscroll);
    }
    if (strcmp(args.yscroll, "") != 0) {
	canvasScrollY(interp, data->re_win, data->win_list, 
		      data->num_wins,
		      data->world->visible, data->canvas, args.yscroll);
    }    
    return TCL_OK;
}

int 
NipZoomCanvas(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    box *zoom;
    nip_zoom_arg args;
    seq_result *result;
    out_canvas *output;
    seq_reg_info info;
    renz_res *data;

    cli_args a[] = {
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(nip_zoom_arg, seq_id)},
	{"-id",        ARG_INT, 1, NULL, offsetof(nip_zoom_arg, id)},
	{"-x1",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, x1)},
	{"-y1",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, y1)},
	{"-x2",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, x2)},
	{"-y2",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, y2)},
	{"-direction", ARG_STR, 1, "b", offsetof(nip_zoom_arg, scroll)},
	{NULL,	0,	0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;
    output = result->output;

    /* zoom back */
    if (args.x1 == -1 && args.y1 == -1 && args.x2 == -1 && args.y2 == -1) {
	canvasZoomback (interp, data->canvas, data->re_win, data->world,
			data->win_list, data->num_wins, &data->zoom);
    } else {

	/* zoom */
	if (NULL == (zoom = (box *)xmalloc(sizeof(box))))
	    return TCL_OK;
	zoom->x1 = args.x1;
	zoom->y1 = args.y1;
	zoom->x2 = args.x2;
	zoom->y2 = args.y2;

	canvasZoom(interp, data->canvas, data->re_win, data->world,
		   data->win_list, data->num_wins, &data->zoom, zoom, 
		   args.scroll[0]);
	xfree(zoom);
    }

    /* redraw ruler ticks by redrawing ruler */
    
    draw_single_ruler(interp, data->ruler, data->canvas, data->ruler->start,
		      data->ruler->end, 1);
    scaleSingleCanvas(interp, data->world, data->canvas, data->ruler->window, 
		      'x', "all");
		      
    /* redraw cursor */
    nip_canvas_cursor_refresh(interp, args.seq_id, output->cursor, 
			      output->cursor, data->canvas, data->win_list, 
			      data->num_wins, result->id, 
			      &output->cursor_visible, data->world, 1);

    return TCL_OK;
}

int NipCanvasCursorX(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{
    nip_cursor_arg args;
    seq_result *result;
    seq_reg_info info;
    renz_res *data;
    double wx, wy;
    char *label;
    
    cli_args a[] = {
	{"-id", ARG_INT, 1, NULL, offsetof(nip_cursor_arg, id)},
	{"-x",  ARG_INT, 1, NULL, offsetof(nip_cursor_arg, cx)},
	{NULL,	0,	0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    CanvasToWorld(data->canvas, args.cx, 0, &wx, &wy);
    label = get_default_astring(interp, tk_utils_defs, w("R_ENZ.CURSOR"));

    canvasCursorX(interp, data->canvas, data->frame, label, 
		  data->cursor.colour, data->cursor.width, args.cx, wx, 
		  data->win_list, data->num_wins);
    
    xfree(label);
    return TCL_OK;
}
int NipResizeCanvas(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{
    nip_resize_arg args;
    seq_result *result;
    seq_reg_info info;
    renz_res *data;
    
    cli_args a[] = {
	{"-id", ARG_INT, 1, NULL, offsetof(nip_resize_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    resizeCanvas(interp, data->re_win, data->win_list, data->num_wins,
		 data->world->visible, data->world->total, data->canvas);
    
    return TCL_OK;
}

/*
 * return the name of an enzyme given a file of names and an index into that
 * file
 */
int 
NipGetREnzName(ClientData clientData,
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[])
{
    nip_enz_name_arg args;
    seq_result *result;
    seq_reg_info info;
    renz_res *data;

    cli_args a[] = {
	{"-id",	    ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, enzyme)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    vTcl_SetResult(interp, "%s", data->r_enzyme[args.enzyme].name);

    return TCL_OK;
}

/*
 * return the name of an enzyme given a file of names and an index into that
 * file
 */
int 
NipGetREnzInfo(ClientData clientData,
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[])
{
    nip_enz_name_arg args;
    seq_reg_generic gen;

    cli_args a[] = {
	{"-id",	    ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, enzyme)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = SEQ_GENERIC;
    gen.task = TASK_NIP_RENZ_INFO;
    gen.data = (void *)&args.enzyme;

    vfuncgroup(5, "restriction enzymes");
    seq_result_notify(args.id, (seq_reg_data *)&gen, 0); 

    return TCL_OK;
}

/*
 * print out restriction enzyme info. 
 * Option 0 = PrintEnzymeByEnzyme
 * Option 1 = OrderOnPosition
 */
int 
NipREnzInfo(ClientData clientData,
	    Tcl_Interp *interp, 
	    int argc, 
	    char *argv[])
{
    nip_enz_info_arg args;
    seq_result *result;
    int seq_num;
    renz_res *data;

    cli_args a[] = {
	{"-result_id", ARG_INT, 1, NULL, offsetof(nip_enz_info_arg, result_id)},
	{"-option",    ARG_INT, 1, NULL, offsetof(nip_enz_info_arg, 
						  print_opt)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    result = seq_id_to_result(args.result_id);
    seq_num = GetSeqNum(result->seq_id[0]);
    data = result->data;

    nip_renz_info(seq_num, data, data->ruler->start, args.print_opt);

    
    return TCL_OK;
}

int NipCanvasToWorld(ClientData clientData,
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[])
{
    nip_world_arg args;
    double wx, wy;
    renz_res *r;
    seq_result *result;
    seq_reg_info info;

    cli_args a[] = {
	{"-id", ARG_INT, 1, NULL, offsetof(nip_world_arg, id)},
	{"-x",  ARG_INT, 1, NULL, offsetof(nip_world_arg, cx)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    r = result->data;
   
    CanvasToWorld(r->canvas, args.cx, 0, &wx, &wy);
    vTcl_SetResult(interp, "%d", (int)wx);
    return TCL_OK;
}


int
FreeNip(ClientData clientData, 
	Tcl_Interp *interp, 
	int argc, 
	char *argv[]) 
{

    Tcl_DeleteInterp(interp);
    return TCL_OK;
}

/* set up tcl commands which call C procedures */
/*****************************************************************************/
/*                                 Nip_Init                                  */
/*****************************************************************************/
int
NipCmds_Init(Tcl_Interp *interp) {

    nip_init_globals(interp); 
    Tcl_CreateCommand(interp, "nip_base_comp", tcl_nip_base_comp, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);

    /* no longer used */
#if 0
    Tcl_CreateCommand(interp, "base_pref", BasePref, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
#endif
    Tcl_CreateCommand(interp, "nip_codon_pref", tcl_nip_codon_pref, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_author_test", tcl_nip_author_test, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_base_bias", tcl_nip_base_bias, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "valid_codon_table", ValidCodonTable, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_trna_search", tcl_nip_trna_search, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_stop_codons", tcl_nip_stop_codons, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_start_codons", tcl_nip_start_codons, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "splice_search", tcl_splice_search, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_string_search", tcl_nip_string_search, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_wtmatrix_search", tcl_nip_wtmatrix_search, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "codon_usage", CodonUsage, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_plot_renz", NipPlotRenz, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_get_renz_name", NipGetREnzName, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_get_renz_info", NipGetREnzInfo, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_renz_info", NipREnzInfo, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_translate", NipTranslate, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "count_base_comp", CountBaseComp, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "translate_orf_to_feature_table", 
		      TranslateORFToFeatureTable, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "translate_orf_to_fasta", 
		      TranslateORFToFasta, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "load_genetic_code",
		      tcl_load_genetic_code, (ClientData)NULL, 
		      (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "count_dinuc_freq", CountDinucFreq, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);


    Tcl_CreateCommand(interp, "nip_scroll_canvas", NipScrollCanvas, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_zoom_canvas", NipZoomCanvas, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_canvas_cursor_x", NipCanvasCursorX, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_resize_canvas", NipResizeCanvas, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nip_canvas_to_world", NipCanvasToWorld, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "free_nip", FreeNip, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
}

/* added for parsing features */
static int get_transl_table(Featcds **key_index){

    int i,k;
    char *tt;
    char *tmp;
    tt = (char *)malloc((20) * sizeof(char));
    tmp = (char *)malloc((20) * sizeof(char));
    for(k=1; k <= key_index[0]->id; k++){  
	for(i=0; i < number_quas; i++){
	    if((key_index[0][k].qualifier[i] != NULL 
                 && !strncmp(&key_index[0][k].qualifier[i][0], 
                 "/transl_table=", 14))){
		tt=strchr(key_index[0][k].qualifier[i],'=');
		strcpy(tmp,&tt[1]);
		if(tmp != NULL) return(atoi(tmp));
		else return 1;
	    } 
	}/*for(k=1*/
    }
    free(tmp);
    free(tt);
    return 1;
}

static int init_genetic_code_ft(Tcl_Interp *interp, int idx ) {
    FILE *fp;
    char infile[1024];
    char *dir = get_default_string(interp, nip_defs, w("GENETIC_CODE_DIR"));
    sprintf( infile, "%s/%s", dir, genetic_code_ft[idx] );
    fp = fopen( infile, "r" );
    if( fp )
    {
	read_global_genetic_code (fp);
	fclose(fp);
	return 0;
    }
    else
    {
        verror( ERR_WARN, "Translation", "Unable to open genetic code file %s.\n", infile );
        return 1;
    }
}


/*
 * Translates subseq and returns an allocated string.
 * Caller is expected to deallocate returned buffer with xfree.
 *
 * Returns NULL on failure.
 */
static char *TranslateSubseq(char *subseq, int rf)
{
    int i;
    char *prot_seq;
    int cnt = 0;
    int length = strlen(subseq);
  
    if (NULL == (prot_seq = (char *)xmalloc(((length/3)+3) * sizeof(char))))
	return NULL;

    for (i = rf; i < length-2; i+=3) {
	prot_seq[cnt++] = codon_to_acid1(&subseq[i]);
    }
   
    prot_seq[cnt-1] = '\0';
    return (prot_seq);
}

/*
 * Allocates and returns a subsequence based on CDS feature number 'k'.
 * Caller is expected to deallocate returned buffer with xfree.
 */
static char *get_subseq(char *dnaseq, Featcds **key_index, int k){

    char *type_loca, *type_range;
    char *range_seq = NULL;
    char *sub_seq = NULL;
    int start_pos, end_pos;
    BasePos *current;
 
    if (NULL == (range_seq=(char*)xmalloc((strlen(dnaseq))*sizeof(char))))
	return NULL;
    if (NULL == (sub_seq=(char*)xmalloc((strlen(dnaseq))*sizeof(char))))
	return NULL;
    strcpy(sub_seq, "");

    type_loca = key_index[0][k].type_loca;
    for(current=key_index[0][k].loca; current!=NULL; current=current->next){
	type_range=current->type_range;
	start_pos=current->start_pos;
	end_pos=current->end_pos;
	/* extract corresponding sub_sequence and to form coding sequence */ 
	strncpy(range_seq, &dnaseq[start_pos-1], end_pos - start_pos+1);   	
	range_seq[end_pos - start_pos+1] = 0;   
	if(!strcmp(type_range,"c"))
	    (void) complement_seq( range_seq, strlen(range_seq));
	strcat(sub_seq,range_seq);
    }    
    if(!strcmp(type_loca,"c") || !strcmp(type_loca,"cj"))
	(void) complement_seq(sub_seq, strlen(sub_seq));      
    xfree(range_seq);
    return (sub_seq);
}

/* added for parsing features */
