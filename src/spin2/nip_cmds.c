#include <tcl.h>
#include <tk.h>
#include <string.h>
#include <math.h>

/* added for parsing features */
#include <stdio.h>
#include "codon_content.h"
/* added for parsing features */

#include "seq_reg.h"
#include "misc.h"
#include "xalloc.h"
#include "nip_structs.h"
#include "cli_arg.h"
#include "text_output.h"
#include "tkSeqed.h"
#include "tkSheet_common.h"
#include "seq_results.h"
#include "sequence_formats.h" /* DNA PROTEIN */
#include "dna_utils.h"
#include "base_comp.h"
#include "seqed.h"
#include "nip_base_comp.h"
#include "ruler_tick.h"
#include "nip_globals.h"
#include "nip_results.h"
#include "nip_gene_search.h"
#include "renz_utils.h"
#include "trna_search.h"
#include "nip_stop_codon.h"
#include "nip_trna_search.h"
#include "splice_search.h"
#include "nip_splice_search.h"
#include "tcl_utils.h"
#include "canvas_box.h"
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
#include "licence.h"
#include "tkCanvGraph.h"
#include "tclCanvGraph.h"
#include "spin_cli_arg.h"
#include "ft_viewer.h"
#include "spin_globals.h"
#include "spin_structs.h"

static char *get_subseq(char *dnaseq, Featcds **key_index, int k);
static int get_transl_table(Featcds **key_index);
static int init_genetic_code_ft(Tcl_Interp *interp, int idx);
static char *TranslateSubseq(char *subseq, int rf);

int nip_list(ClientData clientData,
	     Tcl_Interp *interp,
	     int objc,
	     Tcl_Obj *CONST objv[])
{
    nip_list_arg args;
    seq_result *result;
    int seq_num;
    int num_id;
    char **result_id;
    int i;

    cli_args a[] = {
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(nip_list_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(nip_list_arg, result_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1]))
	return TCL_ERROR;

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
			 int objc,
			 Tcl_Obj *CONST objv[])
{
    pbc_arg args;
    int id;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph;
    Tcl_Obj *id_obj;

    cli_args a[] = {
	{"-win_len",   ARG_INT, 1, "31", offsetof(pbc_arg, win_len)},
	{"-a",         ARG_INT, 1, "1",  offsetof(pbc_arg, a)},
	{"-c",         ARG_INT, 1, "0",  offsetof(pbc_arg, c)},
	{"-g",         ARG_INT, 1, "0",  offsetof(pbc_arg, g)},
	{"-t",         ARG_INT, 1, "1",  offsetof(pbc_arg, t)},
	{"-start",     ARG_INT, 1, "1",  offsetof(pbc_arg, start)},
	{"-end",       ARG_INT, 1, "-1", offsetof(pbc_arg, end)},
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(pbc_arg, seq_id)},
	{"-strand",    ARG_INT, 1, "1", offsetof(pbc_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1]))
	return TCL_ERROR;

    if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);

    if (-1 == init_nip_base_comp_create(interp, args.seq_id, args.start,
					args.end, args.win_len, args.a,
					args.c, args.g, args.t, args.strand,
					&graph, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }

    id_obj = Tcl_NewIntObj(id);
    Tcl_ListObjAppendElement(interp, list, id_obj);
    Tcl_ListObjAppendElement(interp, list, graph);

    Tcl_SetObjResult(interp, list);

    /* FIXME - don't understand Tcl_DecrRefCount - if I do this, if frees
     * everything in the list which I don't want but presumably I should
     * do this somewhere!
     */
    /* Tcl_DecrRefCount(list); */

    return TCL_OK;
}

int nip_base_comp_plot(ClientData clientData,
		       Tcl_Interp *interp,
		       int objc,
		       Tcl_Obj *CONST objv[])
{
    plot_arg args;

    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
	{"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-width", ARG_INT, 1, "-1", offsetof(plot_arg, line_width)},
	{"-fill", ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-orientation", ARG_INT, 1, "-1", offsetof(plot_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "nip base composition", "failure to parse input\n");
	return TCL_ERROR;
    }

    if (args.line_width == -1) {
	args.line_width = get_default_int(interp, nip_defs,
					   "NIP.PBC.L_WIDTH");
    }

    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, nip_defs, "NIP.PBC.COLOUR");
    }

    if (args.orientation == -1) {
	args.orientation = HORIZONTAL;
    }

    init_nip_base_comp_plot(interp, args.seq_id, args.result_id,
			    args.element, args.container, args.results,
			    args.container_id, args.element_id,
			    args.element_type, args.line_width, args.colour,
			    args.orientation);
    return TCL_OK;
}

int tcl_nip_base_comp(ClientData clientData,
		      Tcl_Interp *interp,
		      int objc,
		      Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_base_comp_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_base_comp_plot(clientData, interp, objc, objv);
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

    vfuncheader("base composition");
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
			  int objc,
			  Tcl_Obj *CONST objv[])
{
    gene_arg args;
    int id[3];
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph[3];
    Tcl_Obj *id_obj;
    int i;

    cli_args a[] = {
	{"-codon_table", ARG_STR, 1, NULL, offsetof(gene_arg, codon_table)},
	{"-win_len",     ARG_INT, 1, "0",  offsetof(gene_arg, win_len)},
	{"-start",       ARG_INT, 1, "1",  offsetof(gene_arg, start)},
	{"-end",         ARG_INT, 1, "-1", offsetof(gene_arg, end)},
	{"-option",      ARG_INT, 1, NULL, offsetof(gene_arg, option)},
	{"-seq_id",      ARG_INT, 1, NULL, offsetof(gene_arg, seq_id)},
	{"-strand",      ARG_INT, 1, "1", offsetof(gene_arg, strand)},
	{NULL,           0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	return TCL_ERROR;
    }
    if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);
    if (-1 == init_nip_codon_pref_create(interp, args.seq_id, args.start,
					 args.end, args.codon_table,
					 args.win_len, args.option,
					 args.strand, graph, id)) {
	vTcl_SetResult(interp, "%d %d %d", -1, -1, -1);
	return TCL_OK;
    }

    for (i = 0; i < 3 ; i++) {
	id_obj = Tcl_NewIntObj(id[i]);
	Tcl_ListObjAppendElement(interp, list, id_obj);
	Tcl_ListObjAppendElement(interp, list, graph[i]);
    }
    Tcl_SetObjResult(interp, list);
    /* Tcl_DecrRefCount(list); */

    /* vTcl_SetResult(interp, "%d %d %d", id[0], id[1], id[2]); */
    return TCL_OK;
}

int nip_gene_search_plot(ClientData clientData,
			 Tcl_Interp *interp,
			 int objc,
			 Tcl_Obj *CONST objv[])
{
    plot_arg args;

    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
	{"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-width", ARG_INT, 1, "-1", offsetof(plot_arg, line_width)},
	{"-fill", ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-orientation", ARG_INT, 1, "-1", offsetof(plot_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "gene search plot", "failure to parse input\n");
	return TCL_ERROR;
    }

    if (args.line_width == -1) {
	args.line_width = get_default_int(interp, nip_defs,
					   "NIP.PGS.L_WIDTH");
    }

    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, nip_defs, "NIP.CODONPREF.COLOUR.F1");
    }

    if (args.orientation == -1) {
	args.orientation = HORIZONTAL;
    }

    init_nip_gene_search_plot(interp, args.seq_id, args.result_id,
			      args.element, args.container, args.results,
			      args.container_id, args.element_id,
			      args.element_type, args.line_width, args.colour,
			      args.orientation);

    return TCL_OK;
}

int tcl_nip_codon_pref(ClientData clientData,
		       Tcl_Interp *interp,
		       int objc,
		       Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_codon_pref_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_gene_search_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

int nip_author_test_create(ClientData clientData,
			  Tcl_Interp *interp,
			  int objc,
			  Tcl_Obj *CONST objv[])
{
    author_arg args;
    int id[3];
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph[3];
    Tcl_Obj *id_obj;
    int i;

    cli_args a[] = {
	{"-codon_table", ARG_STR, 1, NULL, offsetof(author_arg, codon_table)},
	{"-error",       ARG_DOUBLE, 1, "0.1", offsetof(author_arg, error)},
	{"-start",       ARG_INT, 1, "1",  offsetof(author_arg, start)},
	{"-end",         ARG_INT, 1, "-1", offsetof(author_arg, end)},
	{"-seq_id",      ARG_INT, 1, NULL, offsetof(author_arg, seq_id)},
	{"-strand",      ARG_INT, 1, "1", offsetof(author_arg, strand)},
	{NULL,           0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	return TCL_ERROR;
    }

    if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);

    if (-1 == init_nip_author_test_create(interp, args.seq_id, args.start,
					  args.end, args.codon_table,
					  args.error, args.strand, graph,
					  id)) {
      vTcl_SetResult(interp, "%d %d %d", -1, -1, -1);
      return TCL_OK;
    }

    for (i = 0; i < 3 ; i++) {
	id_obj = Tcl_NewIntObj(id[i]);
	Tcl_ListObjAppendElement(interp, list, id_obj);
	Tcl_ListObjAppendElement(interp, list, graph[i]);
    }
    Tcl_SetObjResult(interp, list);

    return TCL_OK;

}

int tcl_nip_author_test(ClientData clientData,
			Tcl_Interp *interp,
			int objc,
			Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_author_test_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_gene_search_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

int nip_base_bias_create(ClientData clientData,
			 Tcl_Interp *interp,
			 int objc,
			 Tcl_Obj *CONST objv[])
{
    gene_arg args;
    int id;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph;
    Tcl_Obj *id_obj;

    cli_args a[] = {
	{"-win_len",    ARG_INT, 1, "0",  offsetof(gene_arg, win_len)},
	{"-start",      ARG_INT, 1, "1",  offsetof(gene_arg, start)},
	{"-end",        ARG_INT, 1, "-1", offsetof(gene_arg, end)},
	{"-seq_id",     ARG_INT, 1, NULL, offsetof(gene_arg, seq_id)},
	{"-strand",     ARG_INT, 1, "1", offsetof(gene_arg, strand)},
	{NULL,          0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	return TCL_ERROR;
    }
    if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);
    if (-1 == init_nip_base_bias_create(interp, args.seq_id, args.start,
					args.end, args.win_len,
					args.strand, &graph, &id)) {
      vTcl_SetResult(interp, "%d", -1);
      return TCL_OK;
    }

    id_obj = Tcl_NewIntObj(id);
    Tcl_ListObjAppendElement(interp, list, id_obj);
    Tcl_ListObjAppendElement(interp, list, graph);

    Tcl_SetObjResult(interp, list);

    return TCL_OK;

}

int tcl_nip_base_bias(ClientData clientData,
		      Tcl_Interp *interp,
		      int objc,
		      Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);
    if (strcmp(cmd, "create") == 0) {
	nip_base_bias_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_gene_search_plot(clientData, interp, objc, objv);
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
			   int objc,
			   Tcl_Obj *CONST objv[])
{
    trna_arg args;
    int id;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph;
    Tcl_Obj *id_obj;

    cli_args a[] = {
	{"-start",     ARG_INT, 1, "1",  offsetof(trna_arg, start)},
	{"-end",       ARG_INT, 1, "-1", offsetof(trna_arg, end)},
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(trna_arg, seq_id)},
	{"-strand",    ARG_INT, 1, "1", offsetof(trna_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	return TCL_ERROR;
    }

    if (-1 == init_nip_trna_search_create(interp, args.seq_id, args.start,
					  args.end, args.strand, &graph, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }

     if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);
    id_obj = Tcl_NewIntObj(id);
    Tcl_ListObjAppendElement(interp, list, id_obj);
    Tcl_ListObjAppendElement(interp, list, graph);

    Tcl_SetObjResult(interp, list);
    return TCL_OK;
}

int nip_trna_search_plot(ClientData clientData,
			 Tcl_Interp *interp,
			 int objc,
			 Tcl_Obj *CONST objv[])
{
    plot_arg args;
    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
	{"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-width", ARG_INT, 1, "-1", offsetof(plot_arg, line_width)},
	{"-fill", ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-tick_ht",   ARG_FLOAT, 1, "-1", offsetof(plot_arg, tick_ht)},
	{"-orientation", ARG_INT, 1, "-1", offsetof(plot_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "nip trna search", "unable to parse input\n");
	return TCL_ERROR;
    }

    if (args.line_width == -1) {
	args.line_width = get_default_int(interp, nip_defs,
					   "NIP.TRNA.L_WIDTH");
    }

    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, nip_defs, "NIP.TRNA.COLOUR");
    }

    if (args.orientation == -1) {
	args.orientation = HORIZONTAL;
    }

    if (-1 == init_nip_trna_search_plot(interp, args.seq_id, args.result_id,
					args.element, args.container,
					args.results,
					args.container_id, args.element_id,
					args.element_type, args.line_width,
					args.colour,
					args.tick_ht, args.orientation)) {
      return TCL_ERROR;
    }
    return TCL_OK;
}

int tcl_nip_trna_search(ClientData clientData,
			Tcl_Interp *interp,
			int objc,
			Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_trna_search_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_trna_search_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

int
nip_stop_codons_create(ClientData clientData,
		       Tcl_Interp *interp,
		       int objc,
		       Tcl_Obj *CONST objv[])
{
    s_codon_arg args;
    int id[3];
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph[3];
    Tcl_Obj *id_obj;
    int i;

    cli_args a[] = {
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(s_codon_arg, seq_id)},
	{"-start",      ARG_INT,   1, "1",  offsetof(s_codon_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(s_codon_arg, end)},
	{"-strand",     ARG_INT,   1, "1",  offsetof(s_codon_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1]))
	return TCL_ERROR;

     if (!list)
	return TCL_ERROR;

     Tcl_IncrRefCount(list);

     if (-1 == init_nip_stop_codons_create(args.seq_id, args.start, args.end,
					   args.strand, graph, id)) {
	 vTcl_SetResult(interp, "%d %d %d", -1, -1, -1);
	 return TCL_OK;

     }

     for (i = 0; i < 3; i++) {
	 if (id[i] != -1) {
	     id_obj = Tcl_NewIntObj(id[i]);
	     Tcl_ListObjAppendElement(interp, list, id_obj);
	     Tcl_ListObjAppendElement(interp, list, graph[i]);
	 }
     }
    Tcl_SetObjResult(interp, list);
    return TCL_OK;
}

int
nip_stop_codons_plot(ClientData clientData,
		     Tcl_Interp *interp,
		     int objc,
		     Tcl_Obj *CONST objv[])
{
    plot_arg args;

    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
	{"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-width", ARG_INT, 1, "-1", offsetof(plot_arg, line_width)},
	{"-fill", ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-tick_ht",   ARG_FLOAT, 1, "-1", offsetof(plot_arg, tick_ht)},
	{"-orientation", ARG_INT, 1, "-1", offsetof(plot_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "nip_stop_codons_plot", "failed to parse input\n");
	return TCL_ERROR;
    }

    if (args.line_width == -1) {
	args.line_width = get_default_int(interp, nip_defs,
					   "NIP.STOP_CODON.L_WIDTH");
    }

    if (args.tick_ht == -1) {
	args.tick_ht = get_default_double(interp, nip_defs,
					  "NIP.STOP_CODON.TICK_HT");
    }
    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, nip_defs, "NIP.STOP_CODON.COLOUR.F1");
    }
    if (args.orientation == -1) {
	args.orientation = HORIZONTAL;
    }

    if (-1 == init_nip_stop_codons_plot(interp, args.seq_id,
					args.result_id,
					args.element, args.container,
					args.results,
					args.container_id, args.element_id,
					args.element_type, args.line_width,
					args.colour,
					args.tick_ht, args.orientation)) {
	return TCL_ERROR;
    }

    return TCL_OK;
}

int tcl_nip_stop_codons(ClientData clientData,
			Tcl_Interp *interp,
			int objc,
			Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_stop_codons_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_stop_codons_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

int
nip_start_codons_create(ClientData clientData,
			Tcl_Interp *interp,
			int objc,
			Tcl_Obj *CONST objv[])
{
    s_codon_arg args;
    int id[3];
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph[3];
    Tcl_Obj *id_obj;
    int i;

    cli_args a[] = {
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(s_codon_arg, seq_id)},
	{"-start",      ARG_INT,   1, "1",  offsetof(s_codon_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(s_codon_arg, end)},
	{"-strand",     ARG_INT,   1, "1",  offsetof(s_codon_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1]))
	return TCL_ERROR;

     if (!list)
	return TCL_ERROR;

     Tcl_IncrRefCount(list);

    if (-1 == init_nip_start_codons_create(args.seq_id, args.start, args.end,
					   args.strand, graph, id)) {
	vTcl_SetResult(interp, "%d %d %d", -1, -1, -1);
	return TCL_OK;
    }

     for (i = 0; i < 3 ; i++) { 
	 if (id[i] != -1) {
	     id_obj = Tcl_NewIntObj(id[i]);
	     Tcl_ListObjAppendElement(interp, list, id_obj);
	     Tcl_ListObjAppendElement(interp, list, graph[i]);
	 }
     }
    Tcl_SetObjResult(interp, list);
    return TCL_OK;
}

int tcl_nip_start_codons(ClientData clientData,
			 Tcl_Interp *interp,
			 int objc,
			 Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_start_codons_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_stop_codons_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

int
nip_splice_search_create(ClientData clientData,
			 Tcl_Interp *interp,
			 int objc,
			 Tcl_Obj *CONST objv[])
{
    splice_arg args;
    int id[6];
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph[6];
    Tcl_Obj *id_obj;
    int i;

    cli_args a[] = {
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(splice_arg, seq_id)},
	{"-start",     ARG_INT,   1, "1",  offsetof(splice_arg, start)},
	{"-end",       ARG_INT,   1, "-1", offsetof(splice_arg, end)},
	{"-donor",     ARG_STR,   1, NULL, offsetof(splice_arg, donor)},
	{"-acceptor",  ARG_STR,   1, NULL, offsetof(splice_arg, acceptor)},
	{"-strand",    ARG_INT,   1, "1", offsetof(splice_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	printf("failure to parse args nip_splice_search_create\n");
	return TCL_ERROR;
    }

    if (!list) {
	printf("failed to create list\n");
	return TCL_ERROR;
    }
     Tcl_IncrRefCount(list);

     if (-1 == init_splice_search_create(args.seq_id, args.start, args.end,
					 args.donor, args.acceptor,
					 args.strand, graph, id)) {
	return TCL_OK;
    }

     for (i = 0; i < 6; i++) {
	 if (id[i] != -1) {
	     id_obj = Tcl_NewIntObj(id[i]);
	     Tcl_ListObjAppendElement(interp, list, id_obj);
	     Tcl_ListObjAppendElement(interp, list, graph[i]);
	 }
     }
    Tcl_SetObjResult(interp, list);
    return TCL_OK;
}

int
nip_splice_search_plot(ClientData clientData,
		       Tcl_Interp *interp,
		       int objc,
		       Tcl_Obj *CONST objv[])
{
    plot_arg args;

    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
	{"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-width", ARG_INT, 1, "-1", offsetof(plot_arg, line_width)},
	{"-fill", ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-tick_ht",   ARG_FLOAT, 1, "-1", offsetof(plot_arg, tick_ht)},
	{"-orientation", ARG_INT, 1, "-1", offsetof(plot_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "nip splice search", "failure to parse input\n");
	return TCL_ERROR;
    }
    if (args.line_width == -1) {
	args.line_width = get_default_int(interp, nip_defs,
					   "NIP.SPLICE.L_WIDTH");
    }

    if (args.tick_ht == -1) {
	args.tick_ht = get_default_double(interp, nip_defs,
					  "NIP.SPLICE.TICK_HT");
    }
    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, nip_defs, "NIP.SPLICE.COLOUR.F1");
    }

    if (args.orientation == -1) {
	args.orientation = HORIZONTAL;
    }

    if (-1 == (init_splice_search_plot(interp, args.seq_id, args.result_id,
					args.element, args.container,
					args.results,
					args.container_id, args.element_id,
					args.element_type, args.line_width,
					args.colour,
					args.tick_ht, args.orientation)))
      return TCL_ERROR;

    return TCL_OK;
}

int tcl_splice_search(ClientData clientData,
		     Tcl_Interp *interp,
		     int objc,
		     Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_splice_search_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_splice_search_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

int nip_string_search_create(ClientData clientData,
			     Tcl_Interp *interp,
			     int objc,
			     Tcl_Obj *CONST objv[])
{
    string_arg args;
    int id;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph;
    Tcl_Obj *id_obj;

    cli_args a[] = {
	{"-strand",     ARG_INT,   1, "1",  offsetof(string_arg, strand)},
	{"-min_pmatch", ARG_FLOAT, 1, "75.",offsetof(string_arg, match)},
	{"-string",     ARG_STR,   1, NULL, offsetof(string_arg, string)},
	{"-use_iub",    ARG_INT,   1, "1", offsetof(string_arg, use_iub_code)},
	{"-start",      ARG_INT,   1, "1",  offsetof(string_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(string_arg, end)},
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(string_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1]))
	return TCL_ERROR;

    if (-1 == init_nip_string_search_create(args.strand, args.match,
					    args.string, args.use_iub_code,
					    args.start, args.end, args.seq_id,
					    &graph, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }

     if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);
    id_obj = Tcl_NewIntObj(id);
    Tcl_ListObjAppendElement(interp, list, id_obj);
    Tcl_ListObjAppendElement(interp, list, graph);

    Tcl_SetObjResult(interp, list);
    /* Tcl_DecrRefCount(list); */

    return TCL_OK;

}

int nip_string_search_plot(ClientData clientData,
			   Tcl_Interp *interp,
			   int objc,
			   Tcl_Obj *CONST objv[])
{
    plot_arg args;

    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
	{"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-width", ARG_INT, 1, "-1", offsetof(plot_arg, line_width)},
	{"-fill", ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-tick_ht",   ARG_FLOAT, 1, "-1", offsetof(plot_arg, tick_ht)},
	{"-orientation", ARG_INT, 1, "-1", offsetof(plot_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "nip string search", "unable to parse input\n");
	return TCL_ERROR;
    }

    if (args.line_width == -1) {
	args.line_width = get_default_int(interp, nip_defs,
					   "NIP.STRING_SEARCH.L_WIDTH");
    }

    if (args.tick_ht == -1) {
	args.tick_ht = get_default_double(interp, nip_defs,
					  "NIP.STRING_SEARCH.TICK_HT");
    }
    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, nip_defs, "NIP.STRING_SEATCH.COLOUR");
    }

    if (args.orientation == -1) {
	args.orientation = HORIZONTAL;
    }

    if (-1 == init_nip_string_search_plot(interp, args.seq_id,
					  args.result_id,
					  args.element, args.container,
					  args.results,
					  args.container_id, args.element_id,
					  args.element_type, args.line_width,
					  args.colour,
					  args.tick_ht, args.orientation)) {
	return TCL_ERROR;
    }
    return TCL_OK;
}

int tcl_nip_string_search(ClientData clientData,
			  Tcl_Interp *interp,
			  int objc,
			  Tcl_Obj *CONST objv[])
{
   char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_string_search_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_string_search_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

int nip_wtmatrix_search_create(ClientData clientData,
			       Tcl_Interp *interp,
			       int objc,
			       Tcl_Obj *CONST objv[])
{
    wtmatrix_arg args;
    int id;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph;
    Tcl_Obj *id_obj;

    cli_args a[] = {
	{"-start",      ARG_INT,   1, "1",  offsetof(wtmatrix_arg, start)},
	{"-end",        ARG_INT,   1, "-1", offsetof(wtmatrix_arg, end)},
	{"-seq_id",     ARG_INT,   1, NULL, offsetof(wtmatrix_arg, seq_id)},
	{"-wt_matrix",  ARG_STR,   1, NULL, offsetof(wtmatrix_arg, wt_matrix)},
	{"-strand",     ARG_INT,   1, "1",   offsetof(wtmatrix_arg, strand)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	printf("failure to parse args nip_wtmatrix_search_create\n");
	return TCL_ERROR;
    }

    if (!list) {
	printf("failed to create list\n");
	return TCL_ERROR;
    }
    Tcl_IncrRefCount(list);

    if (-1 == init_nip_wtmatrix_search_create(args.start, args.end,
					      args.seq_id, args.wt_matrix,
					      args.strand, &graph, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    id_obj = Tcl_NewIntObj(id);
    Tcl_ListObjAppendElement(interp, list, id_obj);
    Tcl_ListObjAppendElement(interp, list, graph);

    Tcl_SetObjResult(interp, list);
    return TCL_OK;
}

int nip_wtmatrix_search_plot(ClientData clientData,
			     Tcl_Interp *interp,
			     int objc,
			     Tcl_Obj *CONST objv[])
{
    plot_arg args;

    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
	{"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-width",     ARG_INT, 1, "-1", offsetof(plot_arg, line_width)},
	{"-fill",      ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-tick_ht",   ARG_FLOAT, 1, "-1", offsetof(plot_arg, tick_ht)},
	{"-orientation", ARG_INT, 1, "-1", offsetof(plot_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "nip splice search", "failure to parse input\n");
	return TCL_ERROR;
    }

    if (args.line_width == -1) {
	args.line_width = get_default_int(interp, nip_defs,
					   "NIP.WTMATRIX_SEARCH.L_WIDTH");
    }

    if (args.tick_ht == -1) {
	args.tick_ht = get_default_double(interp, nip_defs,
					  "NIP.WTMATRIX_SEARCH.TICK_HT");
    }
    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, nip_defs, "NIP.WTMATRIX_SEARCH.COLOUR");
    }

    if (args.orientation == -1) {
	args.orientation = HORIZONTAL;
    }

    if (-1 == init_nip_wtmatrix_search_plot(interp, args.seq_id,
					    args.result_id,
					    args.element, args.container,
					    args.results,
					    args.container_id,
					    args.element_id,
					    args.element_type,
					    args.line_width,
					    args.colour,
					    args.tick_ht,
					    args.orientation))
	return TCL_ERROR;

    return TCL_OK;
}


int tcl_nip_wtmatrix_search(ClientData clientData,
			    Tcl_Interp *interp,
			    int objc,
			    Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	nip_wtmatrix_search_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	nip_wtmatrix_search_plot(clientData, interp, objc, objv);
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
    FILE *fp = NULL;
    char *sequence = NULL;
    int num_ranges;
    char **ranges;
    int num_range;
    char **srange;
    range *c_range;
    int   retval   = TCL_OK;
    Tcl_DString input_params;

    cli_args a[] = {
	{"-seq_id",  ARG_INT,   1, NULL, offsetof(codon_usage_arg, seq_id)},
	{"-outfile", ARG_STR,   1, "", offsetof(codon_usage_arg, filename)},
	{"-concatenate", ARG_INT,1, "0",offsetof(codon_usage_arg, concat)},
	{"-format",  ARG_INT,   1, "1",  offsetof(codon_usage_arg, format)},
	{"-strand",  ARG_INT,   1, "1",  offsetof(codon_usage_arg, strand)},
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

	if (args.strand & TOP_S) {
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
    if (args.strand & BOTTOM_S) {
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
    char *subseq   = NULL;
    char *prot_seq = NULL;
    char **selcds  = NULL;
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
	{"-strand",  ARG_INT, 1, "1",  offsetof(trans_ft_arg, strand)},
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
	if (args.strand & TOP_S && args.strand & BOTTOM_S) {
	    strcpy(strand, "both");
	} else if (args.strand & TOP_S) {
	    strcpy(strand, "forward");
	} else {
            strcpy(strand, "reverse");
	}
	vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "strand %s minimum ORF in codons %d\n",
		       GetSeqName(seq_num), args.start, args.end,
		       strand, args.min_orf);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    vmessage("Sequence %s\n", GetSeqName(seq_num));
    if (args.strand & TOP_S) {
	write_screen_open_frames_f_ft(seq, seq_len, args.start, args.end,
				     args.min_orf);
    }

    if (args.strand & BOTTOM_S) {
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
	{"-strand",   ARG_INT, 1, "1",  offsetof(trans_fasta_arg, strand)},
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
	if (args.strand & TOP_S && args.strand & BOTTOM_S) {
	    strcpy(strand, "both");
	} else if (args.strand & TOP_S) {
	    strcpy(strand, "forward");
	} else {
            strcpy(strand, "reverse");
	}
	vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "strand %s minimum ORF in codons %d fasta filename %s\n",
		       GetSeqName(seq_num), args.start, args.end,
		       strand, args.min_orf, args.filename);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    vmessage("Sequence %s\n", GetSeqName(seq_num));
    if (args.strand & TOP_S) {
	write_screen_open_frames_f(seq, seq_len, args.start, args.end,
				   args.min_orf);
	write_open_frames_f(fp, seq, seq_len, args.start, args.end,
			    args.min_orf);
    }

    if (args.strand & BOTTOM_S) {
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

int
FreeNip(ClientData clientData,
	Tcl_Interp *interp,
	int argc,
	char *argv[])
{

    Tcl_DeleteInterp(interp);
    return TCL_OK;
}


typedef struct {
    int start;
    int end;
    int seq_id;
    int strand;
    int orientation;
} ft_viewer_arg;

int ft_viewer_create(ClientData clientData,
		     Tcl_Interp *interp,
		     int objc,
		     Tcl_Obj *CONST objv[])
{
    ft_viewer_arg args;
    int id;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *graph;
    Tcl_Obj *id_obj;

    cli_args a[] = {
	{"-start",     ARG_INT, 1, "1",  offsetof(ft_viewer_arg, start)},
	{"-end",       ARG_INT, 1, "-1", offsetof(ft_viewer_arg, end)},
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(ft_viewer_arg, seq_id)},
	{"-strand",    ARG_INT, 1, "1", offsetof(ft_viewer_arg, strand)},
	{"-orientation",ARG_INT,   1, "1", offsetof(ft_viewer_arg, orientation)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1]))
	return TCL_ERROR;

    if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);

    if (-1 == init_ft_viewer_create(interp, args.seq_id, args.start,
				    args.end, args.strand, args.orientation,
				    &graph, &id)) {
	vTcl_SetResult(interp, "%d", -1);
	return TCL_OK;
    }
    id_obj = Tcl_NewIntObj(id);
    Tcl_ListObjAppendElement(interp, list, id_obj);
    Tcl_ListObjAppendElement(interp, list, graph);

    Tcl_SetObjResult(interp, list);

    return TCL_OK;
}

int ft_viewer_plot(ClientData clientData,
		   Tcl_Interp *interp,
		   int objc,
		   Tcl_Obj *CONST objv[])
{
    plot_arg args;

    cli_args a[] = {
	{"-element",    ARG_STR,   1, NULL, offsetof(plot_arg, element)},
	{"-container", ARG_STR,   1, NULL, offsetof(plot_arg, container)},
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(plot_arg, seq_id)},
	{"-result_id", ARG_INT,   1, NULL, offsetof(plot_arg, result_id)},
        {"-results",   ARG_OBJ,   1, NULL, offsetof(plot_arg, results)},
	{"-container_id", ARG_INT, 1, NULL, offsetof(plot_arg, container_id)},
	{"-element_id", ARG_INT, 1, NULL, offsetof(plot_arg, element_id)},
	{"-element_type", ARG_STR, 1, NULL, offsetof(plot_arg, element_type)},
	{"-tick_ht", ARG_FLOAT, 1, "-1", offsetof(plot_arg, tick_ht)},
	{"-fill", ARG_STR, 1, "", offsetof(plot_arg, colour)},
	{"-offset", ARG_FLOAT, 1, "-1", offsetof(plot_arg, offset)},
	{"-orientation",    ARG_INT,   1, "1", offsetof(plot_arg, orientation)},
	{"-origin",    ARG_DOUBLE,   1, "-1.0", offsetof(plot_arg, origin)},
	{"-display_mode", ARG_INT, 1, "-1", offsetof(plot_arg, display_mode)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == spin_parse_obj_args(a, &args, objc-1, &objv[1])) {
	verror(ERR_WARN, "ft viewer plot", "failure to parse input\n");
	return TCL_ERROR;
    }

    if (args.orientation == CIRCLE) {
	if (args.offset == -1) {
	    args.offset = get_default_int(interp, spin_defs,
					  "FT.CIRCLE.OFFSET");
	}
	if (args.tick_ht == -1) {
	    args.tick_ht = get_default_int(interp, spin_defs,
					   "FT.CIRCLE.HEIGHT");
	}

    } else {
	if (args.offset == -1) {
	    args.offset = get_default_int(interp, spin_defs,
					  "FT.SINGLE.OFFSET");
	}
	if (args.tick_ht == -1) {
	    args.tick_ht = get_default_int(interp, spin_defs,
					   "FT.SINGLE.HEIGHT");
	}
    }
    if (args.origin == -1.0) {
	args.origin = get_default_double(interp, spin_defs,
					 "FT.CIRCLE.ORIGIN");
    }

    if (strcmp(args.colour, "") == 0) {
	args.colour = get_default_string(interp, spin_defs, "FT.COLOUR");
    }

    init_ft_viewer_plot(interp, args.seq_id, args.result_id,
			args.element, args.container, args.results,
			args.container_id, args.element_id, args.colour,
			args.orientation, args.tick_ht,	args.offset,
			args.element_type, args.origin, args.display_mode);
    return TCL_OK;
}

int tcl_ft_viewer(ClientData clientData,
		  Tcl_Interp *interp,
		  int objc,
		  Tcl_Obj *CONST objv[])
{
    char *cmd;

    cmd = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(cmd, "create") == 0) {
	ft_viewer_create(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "list") == 0) {
	nip_list(clientData, interp, objc, objv);
    } else if (strcmp(cmd, "plot") == 0) {
	ft_viewer_plot(clientData, interp, objc, objv);
    }
    return TCL_OK;
}

typedef struct {
    int id;
    int display_mode;
    float offset;
    char *strand;
} ft_update_arg;

int
tcl_ft_viewer_update(ClientData clientData,
		     Tcl_Interp *interp,
		     int argc,
		     char *argv[])
{
    ft_update_arg args;
    seq_result *s_result;
    seq_reg_info info;
#if 0
    plot_data *result;
    feats *extra_data;
#endif

    cli_args a[] = {
	{"-id", ARG_INT, 1, NULL, offsetof(ft_update_arg, id)},
	{"-display_mode", ARG_INT, 1, "-1", offsetof(ft_update_arg, display_mode)},
	{"-offset", ARG_FLOAT, 1, "1", offsetof(ft_update_arg, offset)},
	{"-strand", ARG_INT, 1, "3", offsetof(ft_update_arg, strand)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (args.display_mode != -1) {
	s_result = seq_id_to_result(args.id);

#if 0
	data = s_result->data;
	result = find_plot_data(e, s_result->id);
	(feats *)extra_data = result->extra_data;


	if (args.strand & TOP_S) ||
	    (strcmp(args.strand, "both") == 0)) {
	    data->forward.display = args.display_mode;
	    data->forward.offset = args.offset;
	}
	if ((strcmp(args.strand, "-") == 0) ||
	    (strcmp(args.strand, "both") == 0)) {
	    data->reverse.display = args.display_mode;
	    data->reverse.offset = args.offset;
	}
#endif
    }
    info.job = SEQ_PLOT;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    return TCL_OK;
}

/* set up tcl commands which call C procedures */
/*****************************************************************************/
/*                                 Nip_Init                                  */
/*****************************************************************************/
int
NipCmds_Init(Tcl_Interp *interp) {

    nip_init_globals(interp);
    Tcl_CreateObjCommand(interp, "nip_base_comp", tcl_nip_base_comp,
			 (ClientData) NULL,
			 (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_codon_pref", tcl_nip_codon_pref,
			 (ClientData) NULL,
			 (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_author_test", tcl_nip_author_test,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_base_bias", tcl_nip_base_bias,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "valid_codon_table", ValidCodonTable,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_trna_search", tcl_nip_trna_search,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_stop_codons", tcl_nip_stop_codons,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_start_codons", tcl_nip_start_codons,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "splice_search", tcl_splice_search,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_string_search", tcl_nip_string_search,
			 (ClientData) NULL,
			 (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "nip_wtmatrix_search",
	tcl_nip_wtmatrix_search,
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
    Tcl_CreateObjCommand(interp, "ft_viewer", tcl_ft_viewer,
			 (ClientData) NULL,
			 (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ft_viewer_update", tcl_ft_viewer_update,
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
