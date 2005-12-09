#include <tk.h>

#include "io_utils.h"
#include "fort.h"
#include "IO.h"
#include "gap_globals.h"
#include "fij.h"
#include "io-reg.h"
#include "contig_selector.h"
#include "newgap_cmds.h" /* GetInterp() */
#include "contigEditor.h"
#include "misc.h"
#include "text_output.h"
#include "complement.h"
#include "consen.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"

static int counter;
static int counter_max;
static mobj_fij *global_match;

void *fij_obj_func(int job, void *jdata, obj_fij *obj,
		      mobj_fij *fij) {
    static char buf[80];
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(fij->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(fij->io, cs_id, 0);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	if (io_rdonly(fij->io) && ((obj->c1 > 0 && obj->c2 < 0) ||
				   (obj->c1 < 0 && obj->c2 > 0))) {
	    return "Information\0Hide\0IGNORE\0"
		"IGNORE\0SEPARATOR\0Remove\0";
	} else {
	    return "Information\0Hide\0Invoke join editor *\0"
		"Invoke contig editors\0SEPARATOR\0Remove\0";
	}

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */

	    start_message();
	    vmessage("FIJ match\n");
	    vmessage("    From contig %s(#%d) at %d\n",
		     get_contig_name(fij->io, ABS(obj->c1)),
		     io_clnbr(fij->io, ABS(obj->c1)), obj->pos1);
	    vmessage("    With contig %s(#%d) at %d\n",
		     get_contig_name(fij->io, ABS(obj->c2)),
		     io_clnbr(fij->io, ABS(obj->c2)), obj->pos2);
	    vmessage("    Length %d, score %d, mismatch %2.2f%%\n\n",
		     obj->length, obj->score, ((float)obj->percent)/10000);
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, (obj_match *)obj,
		     (mobj_repeat *)fij, csplot_hash);
	    break;

	case -2: /* default */
        case 2: /* Invoke join editor */ {
	    int cnum[2], llino[2], pos[2];

	    obj->flags |= OBJ_FLAG_VISITED;
	    fij->current = obj - (obj_fij *)fij->match;
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(fij), NULL);

	    cnum[0] = abs(obj->c1);
	    cnum[1] = abs(obj->c2);

	    /* Complement a contig if needed */
	    if ((obj->c1 > 0) != (obj->c2 > 0)) {
		if (cnum[0] == cnum[1]) {
		    verror(ERR_WARN, "join_editor",
			   "cannot display the same contig in two "
			   "different orientations");
		    break;
		}
		if (io_rdonly(fij->io)) {
		    bell();
		    break;
		}

		if (io_clength(fij->io, ABS(obj->c1))
		    < io_clength(fij->io, ABS(obj->c2))) {
		    if (-1 == complement_contig(fij->io, ABS(obj->c1)))
			if (-1 == complement_contig(fij->io, ABS(obj->c2)))
			    return NULL;
		} else {
		    if (-1 == complement_contig(fij->io, ABS(obj->c2)))
			if (-1 == complement_contig(fij->io, ABS(obj->c1)))
			    return NULL;
		}
	    }

	    /*
	     * NB: obj->pos1 now may not be the same value as when this
	     * function was entered, due to the complementing!
	     */
	    pos[0] = obj->pos1;
	    pos[1] = obj->pos2;

	    llino[0] = io_clnbr(fij->io, cnum[0]);
	    llino[1] = io_clnbr(fij->io, cnum[1]);

	    join_contig(GetInterp(), fij->io, cnum, llino, pos,
			consensus_cutoff, quality_cutoff);
	    break;
	}

	case 3: /* Invoke contig editors */ {
	    int cnum, llino, pos, reveal;

	    cnum  = ABS(obj->c1);
	    llino = io_clnbr(fij->io, cnum);
	    pos   = obj->pos1;
	    reveal= (obj->pos1 <= 0 ||
		     obj->pos2 <= 0 ||
		     obj->pos1 >= io_clength(fij->io, ABS(obj->c1)) ||
		     obj->pos2 >= io_clength(fij->io, ABS(obj->c2))) ? 1 : 0;

	    edit_contig(GetInterp(), fij->io, cnum, llino, pos,
			consensus_cutoff, quality_cutoff, reveal, NULL);

	    cnum  = ABS(obj->c2);
	    llino = io_clnbr(fij->io, cnum);
	    pos   = obj->pos2;

	    edit_contig(GetInterp(), fij->io, cnum, llino, pos,
			consensus_cutoff, quality_cutoff, reveal, NULL);
	    break;
	}

	case 4: /* Remove */
	    obj_remove(GetInterp(), cs->window, (obj_match *)obj,
		       (mobj_repeat *)fij, csplot_hash);
	    break;

        }
        break;

    case OBJ_GET_BRIEF:
	sprintf(buf,
		"FIJ: %c#%d@%d with %c#%d@%d, len %d, score %d, mis %2.2f%%",
		obj->c1 > 0 ? '+' : '-',
		io_clnbr(fij->io, ABS(obj->c1)), obj->pos1,
		obj->c2 > 0 ? '+' : '-',
		io_clnbr(fij->io, ABS(obj->c2)), obj->pos2,
		obj->length, obj->score, ((float)obj->percent)/10000);
	return buf;
    }

    return NULL;
}

static int sort_func(const void *p1, const void *p2) {
    obj_fij *m1 = (obj_fij *)p1, *m2 = (obj_fij *)p2;
    return m2->score - m1->score;
}

/*
 * Match callback.
 * 'obj' is a match contained within the 'fij' list.
 */
void fij_callback(GapIO *io, int contig, void *fdata, reg_data *jdata) {
    mobj_fij *r = (mobj_fij *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id, 0);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Find Internal Joins");
	break;


    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, (mobj_repeat *)r,
			csplot_hash, cs->window);
	break;


    case REG_COMPLEMENT:

	csmatch_complement(io, contig, (mobj_repeat *)r,
			   csplot_hash, cs->window);
	break;


    case REG_GET_OPS:

	if (r->all_hidden)
	    jdata->get_ops.ops = "PLACEHOLDER\0PLACEHOLDER\0Information\0"
		"PLACEHOLDER\0Hide all\0Reveal all\0Sort Matches\0"
		    "SEPARATOR\0Remove\0";
	else
	    jdata->get_ops.ops = "Use for 'Next'\0Reset 'Next'\0Information\0"
		"Configure\0Hide all\0Reveal all\0Sort Matches\0"
		    "SEPARATOR\0Remove\0";
	break;


    case REG_INVOKE_OP:

	switch (jdata->invoke_op.op) {
	case 0: /* Next */
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(r), NULL);
	    break;
	case 1: /* Reset Next */
	    csmatch_reset_next((mobj_repeat *)r);
	    break;
	case 2: /* Information */
	    csmatch_info((mobj_repeat *)r, "Find Internal Joins");
	    break;
	case 3: /* Configure */
	    csmatch_configure(io, cs->window, (mobj_repeat *)r);
	    break;
	case 4: /* Hide all */
	    csmatch_hide(GetInterp(), cs->window, (mobj_repeat *)r,
			 csplot_hash);
	    break;
	case 5: /* Reveal all */
	    csmatch_reveal(GetInterp(), cs->window, (mobj_repeat *)r,
			   csplot_hash);
	    break;
	case 6: /* Sort */
	    qsort(r->match, r->num_match, sizeof(obj_fij), sort_func);
	    csmatch_reset_hash(csplot_hash, (mobj_repeat *)r);
	    r->current = -1;
	    break;
	case 7: /* Remove */
	    csmatch_remove(io, cs->window,
			   (mobj_repeat *)r,
			   csplot_hash);
	    break;
	}
	break;


    case REG_PARAMS:

	jdata->params.string = r->params;
	break;


    case REG_NUMBER_CHANGE:

	csmatch_renumber(io, contig, jdata->number.number,
			 (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_ORDER:
#ifdef DEBUG
	printf("find internal joins REG_ORDER\n");
#endif
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs->window,
		       (mobj_repeat *)r,
		       csplot_hash);
	break;

    case REG_DELETE:

	csmatch_contig_delete(io, (mobj_repeat *)r, contig,
			      cs->window, csplot_hash);
	break;

    case REG_LENGTH:
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;
    }
}

/*
 * Main find internal joins algorithm entry.
 *
 * Mode should be one of COMPARE_ALL (0) or COMPARE_SINGLE (1)
 * When mode is COMPARE_SINGLE, the first element of contig_array should
 * be the single contig to compare (against all other elements in
 * contig_array). Otherwise, all elements of contig_array are compared
 * against all other elements.
 */
int
fij(GapIO *io,
    int mask,
    int compare_mode, /* could be sorted out to save work below*/
    int min_overlap,
    double max_percent_mismatch,
    int word_len,
    double max_prob,
    int min_match,
    int band,
    int window_len,
    int max_unknown,
    double min_conf,
    int use_conf,
    int use_hidden,
    int max_alignment,
    int num_contigs,
    contig_list_t *contig_array)
{
    char *consensus;
    mobj_fij *FIJMatch;
    int i, id;
    char *val;
    int ret;
    int max_consensus;
    Contig_parms *contig_list;
    int database_size, number_of_contigs;
    int task_mask;
    int consensus_length, max_read_length;
    static char buf[80];
    Hidden_params p;

    /* FIXME: get these from elsewhere */

    int gap_open, gap_extend;

    gap_open = gopenval; /* note we were using 8 and 6, not 8 and 7 */
    gap_extend = gextendval;

    p.min = p.max = p.verbose = 0;
    p.do_it = use_hidden;
    p.use_conf = use_conf;
    p.test_mode = 0;
    p.start = 0;
    p.lwin1 = 0;
    p.lcnt1 = 0;
    p.rwin1 = window_len;
    p.rcnt1 = max_unknown;
    p.qual_val = min_conf;
    p.window_len = window_len;
    p.gap_open = gopenval;
    p.gap_extend = gextendval;
    p.band = 30; /*FIXME: hardwired 30 bases band for aligning hidden data */
    max_consensus = maxseq;  /* global! */

    max_read_length = find_max_gel_len(io, 0, 0);

    if ((consensus = (char *)xmalloc(max_consensus * sizeof(char)))==NULL){
	return(-1);
    }

    if ((FIJMatch = (mobj_fij *)xmalloc(sizeof(mobj_fij)))==NULL){
	xfree ( consensus );
	return -1;
    }

    counter_max = 14;
    if (NULL == (FIJMatch->match = (obj_fij *)xmalloc(counter_max *
						      sizeof(obj_fij)))) {
	xfree ( consensus );
	xfree(FIJMatch);
	return -1;
    }
    database_size = io_dbsize(io);
    number_of_contigs = num_contigs;
    if ( ! (contig_list = get_contig_list ( database_size, io,
			   number_of_contigs, contig_array ))) {
	xfree ( consensus );
	xfree(FIJMatch->match);
	xfree(FIJMatch);
	return -5;
    }

    global_match = FIJMatch;
    counter = 0;

    task_mask = ADDTITLE | NORMALCONSENSUS;
    if ( mask == 2 ) task_mask |= MARKING;
    if ( mask == 3 ) task_mask |= MASKING;
    if ( p.do_it ) task_mask |= ADDHIDDENDATA;

    consensus_length = 0;
    if ( ret = make_consensus ( task_mask, io, consensus, NULL,
			  contig_list, number_of_contigs,
			  &consensus_length,
			  max_read_length,
			  max_consensus,
			  p,
			  consensus_cutoff ) ) {

	xfree ( consensus );
	xfree(FIJMatch);
	xfree(FIJMatch->match);
	xfree(contig_list);
	return -1;
    }

    if ( ret = do_it_fij ( consensus, consensus_length,
			   word_len,
			   min_overlap, max_percent_mismatch, compare_mode,
			   band, gap_open, gap_extend, max_prob, min_match,
			   max_alignment, contig_list, number_of_contigs ) ) {

	xfree(FIJMatch->match);
	xfree(FIJMatch);
	xfree(contig_list);
	xfree ( consensus );
	return -1;
    }

    if (0 == counter) {
	vmessage("No joins found \n");
	xfree(FIJMatch->match);
	xfree(FIJMatch);
	xfree(contig_list);
	xfree ( consensus );
	return 0;
    }

    sprintf(buf, " Number of potential joins found   %d", counter);
    vmessage("%s\n", buf);

    FIJMatch->num_match = counter;
    FIJMatch->io = io;
    strcpy(FIJMatch->tagname, CPtr2Tcl(FIJMatch));

    val = get_default_string(GetInterp(), gap_defs, "FIJ.COLOUR");
    strcpy(FIJMatch->colour, val);

    FIJMatch->linewidth = get_default_int(GetInterp(), gap_defs,
					  "FIJ.LINEWIDTH");

    FIJMatch->params = (char *)xmalloc(100);
    if (FIJMatch->params)
	sprintf(FIJMatch->params, "Unknown at present");
    FIJMatch->all_hidden = 0;
    FIJMatch->current = -1;
    FIJMatch->reg_func = fij_callback;
    FIJMatch->match_type = REG_TYPE_FIJ;

    for (i = 0; i < counter; i++){
	obj_fij *match = &FIJMatch->match[i];

	if (match->c1 < 0) {
	    match->c1 =  rnumtocnum(io, ABS(match->c1)) * -1;
	    match->pos1 = io_clength(io, -match->c1) - match->pos1 + 1;
	} else {
	    match->c1 =  rnumtocnum(io, ABS(match->c1));
	}
	if (match->c2 < 0) {
	    match->c2 =  rnumtocnum(io, ABS(match->c2)) * -1;
	    match->pos2 = io_clength(io, -match->c2) - match->pos2 + 1;
	} else {
	    match->c2 =  rnumtocnum(io, ABS(match->c2));
	}
    }

    /* Sort matches */
    qsort(FIJMatch->match, FIJMatch->num_match, sizeof(obj_fij), sort_func);

    PlotRepeats(io, (mobj_repeat *)FIJMatch);
    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(FIJMatch), NULL);

    /*
     * Register find internal joins with each of the contigs used.
     * Currently we assume that this is all.
     */
    if (counter) {
	id = register_id();
	for (i = 1; i <= NumContigs(io); i++) {
	    contig_register(io, i, fij_callback, (void *)FIJMatch, id,
			    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			    REG_NUMBER_CHANGE | REG_ORDER, REG_TYPE_FIJ);
	}
    }

    xfree(contig_list);
    xfree ( consensus );

    return 0;
} /* end fij */

/* store hits of find internal joins */
void
buffij(int c1,
	int pos1,
	int c2,
	int pos2,
	int len,
        int score,
        double percent)
{
    global_match->match[counter].func = fij_obj_func;
    global_match->match[counter].data = global_match;

    global_match->match[counter].c1 = c1;
    global_match->match[counter].pos1 = pos2;
    global_match->match[counter].c2 = c2;
    global_match->match[counter].pos2 = pos1;
    global_match->match[counter].length = len;
    global_match->match[counter].score = score;
    global_match->match[counter].percent = 10000 * percent;
    global_match->match[counter].flags = 0;

    if (++counter >= counter_max) {
	counter_max *= 2;
	global_match->match = (obj_fij *)xrealloc(global_match->match,
						    counter_max *
						    sizeof(obj_fij));
    }

}
