#include <tcl.h>

#include "IO.h"
#include "gap_cli_arg.h"
#include "list_proc.h"
#include "haplo.h"
#include "xalloc.h"

typedef struct {
    GapIO *io;
    char *inlist;
    int discrep_cutoff;
    int snp_cutoff;
    int min_base_qual;
    int two_alleles;
} snp_arg;

static int haplo_snp_cmd(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[]) {
    int rargc, i;
    contig_list_t *rargv;
    snp_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(snp_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(snp_arg, inlist)},
	{"-discrep_cutoff",
	                ARG_INT, 1, "40",  offsetof(snp_arg, discrep_cutoff)},
	{"-min_base_qual",
	                ARG_INT, 1, "0",   offsetof(snp_arg, min_base_qual)},
	{"-snp_cutoff", ARG_INT, 1, "40",  offsetof(snp_arg, snp_cutoff)},
	{"-two_alleles",ARG_INT, 1, "0",   offsetof(snp_arg, two_alleles)},
	{NULL,	    0,	     0, NULL, 0}
    };
    dstring_t *ds;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    ds = dstring_create(NULL);
    for (i = 0; i < rargc; i++) {
	dstring_t *ds2;
	ds2 = haplo_snps(args.io,
			 rargv[i].contig, rargv[i].start, rargv[i].end,
			 args.discrep_cutoff, args.snp_cutoff,
			 args.min_base_qual, args.two_alleles);

	if (ds2) {
	    dstring_appendf(ds, "{%s} ", dstring_str(ds2));
	}

	dstring_destroy(ds2);
    }

    Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);

    dstring_destroy(ds);
    xfree(rargv);
    return TCL_OK;
}


typedef struct {
    GapIO *io;
    Tcl_Obj *list;
    int verbosity;
    double minscore;
    int twopass;
    int fast;
    double correlation_offset;
    int max_sets;
} iolist_arg;


/**
 * This is the Tcl interface between the Tcl haplo_split C code.
 * It converts a snp Tcl list back into appropriate data structures to
 * start re-grouping templates by.
 */
static int haplo_split_cmd(Tcl_Interp *interp,
			   int objc,
			   Tcl_Obj *CONST objv[]) {
    iolist_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,     1, NULL, offsetof(iolist_arg, io)},
	{"-snps",	ARG_OBJ,    1, NULL, offsetof(iolist_arg, list)},
	{"-verbosity",	ARG_INT,    1, "0",  offsetof(iolist_arg, verbosity)},
	{"-minscore",   ARG_DOUBLE, 1, "0",  offsetof(iolist_arg, minscore)},
	{"-twopass",    ARG_INT,    1, "0",  offsetof(iolist_arg, twopass)},
	{"-fast",       ARG_INT,    1, "0",  offsetof(iolist_arg, fast)},
	{"-correlation_offset",
	                ARG_DOUBLE, 1, ".9", offsetof(iolist_arg,
						      correlation_offset)},
	{"-max_sets",   ARG_INT,    1, "99", offsetof(iolist_arg, max_sets)},
	{NULL,	    0,	     0, NULL, 0}
    };
    int robjc, i, j;
    Tcl_Obj **robjv;
    snp_t *snp = NULL;
    int nsnps = 0;
    dstring_t *ds;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    /* Split the list */
    Tcl_ListObjGetElements(interp, args.list, &robjc, &robjv);
    
    /* Allocate SNP memory */
    nsnps = robjc;
    if (!nsnps) {
	Tcl_ResetResult(interp);
	return TCL_OK;
    }

    if (NULL == (snp = (snp_t *)xmalloc(nsnps * sizeof(*snp))))
	return TCL_ERROR;

    /* Convert TCL lists into snp_t/seq_base_t structs */
    for (i = 0; i < robjc; i++) {
	int nele;
	Tcl_Obj *posobj, *scoreobj;
	int pos;
	double score;

	/* Get the position and length */
	if (TCL_OK != Tcl_ListObjLength(interp, robjv[i], &nele))
	    return TCL_ERROR;

	if (TCL_OK != Tcl_ListObjIndex(interp, robjv[i], 0, &posobj) ||
	    TCL_OK != Tcl_GetIntFromObj(interp, posobj, &pos))
	    return TCL_ERROR;

	if (TCL_OK != Tcl_ListObjIndex(interp, robjv[i], 1, &scoreobj) ||
	    TCL_OK != Tcl_GetDoubleFromObj(interp, scoreobj, &score))
	    return TCL_ERROR;

	snp[i].pos = pos;
	snp[i].score = score;
	snp[i].seqs = NULL;
	snp[i].nseqs = 0;

	/* Loop around the template sequences */
	for (j = 2; j < nele; j++) {
	    int n4;
	    Tcl_Obj *seqobj, *obj;
	    int tnum, conf;
	    double tscore;
	    char *seq;
	    int index;

	    if (TCL_OK != Tcl_ListObjIndex(interp, robjv[i], j, &seqobj))
		return TCL_ERROR;
	    
	    if (TCL_OK != Tcl_ListObjLength(interp, seqobj, &n4) ||
		n4 != 4)
		return TCL_ERROR;

	    if (TCL_OK != Tcl_ListObjIndex(interp, seqobj, 0, &obj) ||
		TCL_OK != Tcl_GetIntFromObj(interp, obj, &tnum))
		return TCL_ERROR;

	    if (TCL_OK != Tcl_ListObjIndex(interp, seqobj, 1, &obj) ||
		TCL_OK != Tcl_GetDoubleFromObj(interp, obj, &tscore))
		return TCL_ERROR;

	    if (TCL_OK != Tcl_ListObjIndex(interp, seqobj, 2, &obj) ||
		NULL == (seq = Tcl_GetString(obj)))
		return TCL_ERROR;

	    if (TCL_OK != Tcl_ListObjIndex(interp, seqobj, 3, &obj) ||
		TCL_OK != Tcl_GetIntFromObj(interp, obj, &conf))
		return TCL_ERROR;

	    index = snp[i].nseqs++;
	    snp[i].seqs = xrealloc(snp[i].seqs,
				   snp[i].nseqs * sizeof(seq_base_t));

	    snp[i].seqs[index].tmplate = tnum;
	    snp[i].seqs[index].tscore = tscore;
	    snp[i].seqs[index].base = *seq;
	    snp[i].seqs[index].conf = conf;
	}
    }

    ds = haplo_split(args.io, snp, nsnps, args.verbosity, args.minscore,
		     args.twopass, args.fast, args.correlation_offset,
		     args.max_sets);
    Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);

    dstring_destroy(ds);
    free_snps(snp, nsnps);

    return TCL_OK;
}


#if 0
typedef struct {
    GapIO *io;
    char *list;
} iolist_arg;

static int haplo_split_cmd(Tcl_Interp *interp,
			   int objc,
			   Tcl_Obj *CONST objv[]) {
    iolist_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(iolist_arg, io)},
	{"-snps",	ARG_STR, 1, NULL,  offsetof(iolist_arg, list)},
	{NULL,	    0,	     0, NULL, 0}
    };
    int count;
    char *cp;
    enum state_e {
	SNP_START, SNP_END, SNP_POS, SNP_TEMPLATE, SNP_EXIT
    };
    enum state_e state = SNP_START;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    count = 0;
    for (cp = args.list;*cp && state != SNP_EXIT;) {
	switch(state) {
	case SNP_START:
	    while (*cp && isspace(*cp))
		cp++;

	    if (*cp++ != '{')
		goto error;

	    state = SNP_START;
	    /* flow through */

	case SNP_POS: {
	    int pos, index;
	    double score;

	    if (2 != sscanf(cp, "%d %lf%n", &pos, &score, &index))
		goto error;
	    cp += index;

	    state = SNP_TEMPLATE;
	    /* flow through */
	}

	case SNP_TEMPLATE: {
	    int tnum, conf, index;
	    double tscore;
	    char base;

	    while (*cp && isspace(*cp))
		cp++;

	    if (4 != sscanf(cp, "{%d %lf %c %d}%n",
			    &tnum, &tscore, &base, &conf, &index))
		goto error;
	    cp += index;

	    while (*cp && isspace(*cp))
		cp++;

	    if (*cp && *cp == '}')
		state = SNP_END;

	    break;
	}

	case SNP_END:
	    puts("");

	    cp++; /* } */
	    while (*cp && isspace(*cp))
		cp++;

	    state = *cp ? SNP_START : SNP_EXIT;
	    break;
	}
    }

    return TCL_OK;

 error:
    fprintf(stderr, "Unexpected symbol in snp list\n");
    return TCL_ERROR;
}
#endif


typedef struct {
    GapIO *io;
    char *contig;
    char *tlist;
} cons_arg;

static int haplo_consensus_cmd(Tcl_Interp *interp,
			       int objc,
			       Tcl_Obj *CONST objv[]) {
    int i, cnum, clen;
    cons_arg args;
    char **rargv;
    int rargc;
    int *templates;
    char *cons;
    unsigned char *qual_str;
    float *qual;
    Tcl_Obj *res, *obj[2];

    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(cons_arg, io)},
	{"-templates",	ARG_STR, 1, NULL,  offsetof(cons_arg, tlist)},
	{"-contig",	ARG_STR, 1, NULL,  offsetof(cons_arg, contig)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (-1 == (cnum = get_contig_num(args.io, args.contig, GGN_ID)))
	return TCL_ERROR;

    if (TCL_OK != Tcl_SplitList(interp, args.tlist, &rargc, &rargv))
	return TCL_ERROR;

    if (rargc) {
	templates = (int *)xmalloc(rargc * sizeof(int *));
	for (i = 0; i < rargc; i++) {
	    templates[i] = template_name_to_number(args.io, rargv[i]);
	}
    } else {
	templates = NULL;
    }
    if (rargv)
	Tcl_Free((char *)rargv);

    clen = io_clength(args.io, cnum);
    calc_template_consensus(args.io, cnum, 1, clen, templates, rargc,
			    &cons, &qual);
    qual_str = (unsigned char *)xmalloc(clen);
    for (i = 0; i < clen; i++) {
	qual_str[i] = (unsigned char)qual[i];
    }
    obj[0] = Tcl_NewStringObj(cons, clen);
    obj[1] = Tcl_NewStringObj((char *)qual_str, clen);
    res = Tcl_NewListObj(2, obj);
    Tcl_SetObjResult(interp, res);

    xfree(cons);
    xfree(qual);
    xfree(qual_str);
    if (templates)
	xfree(templates);

    return TCL_OK;
}


typedef struct {
    GapIO *io;
    char *inlist;
} tdepth_arg;

static int haplo_tdepth_cmd(Tcl_Interp *interp,
			    int objc,
			    Tcl_Obj *CONST objv[]) {
    int rargc, i, j, start, end, last, depth;
    contig_list_t *rargv;
    tdepth_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(tdepth_arg, io)},
	{"-contig",	ARG_STR, 1, NULL,  offsetof(tdepth_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };
    int *tdepth;
    dstring_t *ds;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);
    tdepth = (int *)xcalloc(io_clength(args.io, rargv[0].contig), sizeof(int));
    depth = calc_template_depth(args.io, rargv[0].contig, rargv[0].start,
				rargv[0].end, tdepth);

    /* Convert the depth array into a series of (start, end, depth) ranges */
    ds = dstring_create(NULL);
    dstring_appendf(ds, "%d ", depth);
    start = rargv[0].start;
    end = rargv[0].end;
    last = tdepth[0];
    for (j = 0, i = 1; i < end-start+1; i++) {
	if (tdepth[i] != last) {
	    dstring_appendf(ds, "%d %d %d ", j+start, i+start, last);
	    j = i;
	    last = tdepth[i];
	}
    }
    dstring_appendf(ds, "%d %d %d ", j+start, i+start, last);

    Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);

    xfree(rargv);
    xfree(tdepth);
    dstring_destroy(ds);

    return TCL_OK;
}

/**
 * Returns a list of SNP locations and scores. The format is a Tcl
 * list of lists, with each sub-list containing the SNP location, the
 * score and a string of base calls present at that coordinate.
 */
static int tcl_haplo_cmd(ClientData clientData, Tcl_Interp *interp,
			 int objc, Tcl_Obj *CONST objv[]) {
    int opt;
    static char *subcmds[] = {
	"snps",         "split",        "consensus",    "tdepth",
	(char *)NULL
    };
    enum subcmds {
	HAPLO_SNPS,     HAPLO_SPLIT,    HAPLO_CONS,     HAPLO_TDEPTH
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "options ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], subcmds, "option", 0, &opt)
	!= TCL_OK) {
	return TCL_ERROR;
    }

    switch ((enum subcmds)opt) {
    case HAPLO_SNPS:
	return haplo_snp_cmd(interp, objc-1, objv+1);

    case HAPLO_SPLIT:
	return haplo_split_cmd(interp, objc-1, objv+1);

    case HAPLO_CONS:
	return haplo_consensus_cmd(interp, objc-1, objv+1);

    case HAPLO_TDEPTH:
	return haplo_tdepth_cmd(interp, objc-1, objv+1);
    }

    return TCL_OK;
}


int Haplo_Init(Tcl_Interp *interp) {
    if (NULL == Tcl_CreateObjCommand(interp, "haplo", tcl_haplo_cmd,
                                     (ClientData)NULL,
                                     (Tcl_CmdDeleteProc *) NULL))
        return TCL_ERROR;

    return TCL_OK;
}

int Haplo_SafeInit(Tcl_Interp *interp) {
    return Haplo_Init(interp);
}
