#include "tg_gio.h"
#include "gap_globals.h"
#include "misc.h"
#include "cs-object.h"
#include "contig_selector.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "find_repeats.h"
#include "consen.h"
#include "gap_hash.h"
#include "gap4_compat.h"
#include "editor_view.h"
#include "tk-io-reg.h"

/*
 * Match callback.
 * 'obj' is a match contained within the 'repeat' list.
 */
void *repeat_obj_func(int job, void *jdata, obj_match *obj,
		      mobj_repeat *repeat) {
    static char buf[80];
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(repeat->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(repeat->io, cs_id);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	if (io_rdonly(repeat->io) && ((obj->c1 > 0 && obj->c2 < 0) ||
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
	    vmessage("Repeat match (%s)\n",
		     ((obj->c1 > 0) == (obj->c2 > 0)) ? "direct" : "inverted");
	    vmessage("    From contig %s(#%"PRIrec") at %d\n",
		     get_contig_name(repeat->io, ABS(obj->c1)),
		     io_clnbr(repeat->io, ABS(obj->c1)), obj->pos1);
	    vmessage("    With contig %s(#%"PRIrec") at %d\n",
		     get_contig_name(repeat->io, ABS(obj->c2)),
		     io_clnbr(repeat->io, ABS(obj->c2)), obj->pos2);
	    vmessage("    Length %d\n\n", obj->length);
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, obj,
		     (mobj_repeat *)repeat, csplot_hash);
	    break;

	case -2: /* default */
	case 2: /* Invoke join editor */ {
	    tg_rec cnum[2], llino[2];
	    int pos[2];

	    obj->flags |= OBJ_FLAG_VISITED;
	    repeat->current = obj - repeat->match;
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(repeat), NULL);

	    cnum[0] = ABS(obj->c1);
	    cnum[1] = ABS(obj->c2);

	    /* Complement a contig if needed */
	    if ((obj->c1 > 0) != (obj->c2 > 0)) {
		if (cnum[0] == cnum[1]) {
		    verror(ERR_WARN, "join_editor",
			   "cannot display the same contig in two "
			   "different orientations");
		    break;
		}
		if (io_rdonly(repeat->io)) {
		    bell();
		    break;
		}

		if (io_clength(repeat->io, ABS(obj->c1)) <
		    io_clength(repeat->io, ABS(obj->c2))) {
		    if (-1 == complement_contig(repeat->io, ABS(obj->c1)))
			if (-1 == complement_contig(repeat->io, ABS(obj->c2)))
			    return NULL;
		} else {
		    if (-1 == complement_contig(repeat->io, ABS(obj->c2)))
			if (-1 == complement_contig(repeat->io, ABS(obj->c1)))
			    return NULL;
		}
	    }

	    /*
	     * NB: obj->pos1 now may not be the same value as when this
	     * function was entered, due to the complementing!
	     */
	    pos[0] = obj->pos1;
	    pos[1] = obj->pos2;

	    llino[0] = io_clnbr(repeat->io, cnum[0]);
	    llino[1] = io_clnbr(repeat->io, cnum[1]);

	    join_contig(repeat->io, cnum, llino, pos);
	    break;
	}

	case 3: /* Invoke contig editors */ {
	    tg_rec cnum, llino;
	    int pos;

	    cnum  = ABS(obj->c1);
	    llino = io_clnbr(repeat->io, cnum);
	    pos   = obj->pos1;
	    edit_contig(repeat->io, cnum, llino, pos);
    
	    cnum  = ABS(obj->c2);
	    llino = io_clnbr(repeat->io, cnum);
	    pos   = obj->pos2;
	    edit_contig(repeat->io, cnum, llino, pos);
	    break;
	}

	case 4: /* Remove */
	    obj_remove(GetInterp(), cs->window, obj,
		     (mobj_repeat *)repeat, csplot_hash);
	    break;

	}
	break;

    case OBJ_GET_BRIEF:
	sprintf(buf, "Repeat: %c#%"PRIrec"@%d with %c#%"PRIrec"@%d, len %d",
		obj->c1 > 0 ? '+' : '-',
		io_clnbr(repeat->io, ABS(obj->c1)), obj->pos1,
		obj->c2 > 0 ? '+' : '-',
		io_clnbr(repeat->io, ABS(obj->c2)), obj->pos2,
		obj->length);
	return buf;
    }

    return NULL;
}

static int sort_func(const void *p1, const void *p2) {
    obj_match *m1 = (obj_match *)p1, *m2 = (obj_match *)p2;
    return m2->score - m1->score;
}

void repeat_callback(GapIO *io, tg_rec contig, void *fdata, reg_data *jdata) {
    mobj_repeat *r = (mobj_repeat *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Repeat search");
	break;


    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, r, csplot_hash, cs->window);
	break;


    case REG_COMPLEMENT:

	csmatch_complement(io, contig, r, csplot_hash, cs->window);
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
	    csmatch_info((mobj_repeat *)r, "Find Repeats");
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
	    qsort(r->match, r->num_match, sizeof(obj_match), sort_func);
	    csmatch_reset_hash(csplot_hash, (mobj_repeat *)r);
	    r->current = -1;
	    break;
	case 7: /* Remove */
	    csmatch_remove(io, cs->window, (mobj_repeat *)r,
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
	printf("find repeats REG_ORDER\n");
#endif
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs->window, (mobj_repeat *)r, csplot_hash);
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

void
plot_rpt(GapIO *io,
	 int nres,
	 tg_rec c1[],
	 int pos1[],
	 tg_rec c2[],
	 int pos2[],
	 int len[])
{
    int i, id;
    mobj_repeat *repeat;
    obj_match *matches = NULL;
    char *val;

    /* If nres is zero - do nothing */
    if (0 == nres)
	return;

    if (NULL == (repeat = (mobj_repeat *)xmalloc(sizeof(mobj_repeat)))) {
	f_proc_return();
    }

    if (NULL == (matches = (obj_match *)xmalloc(nres * sizeof(obj_match)))) {
	xfree(repeat);
	f_proc_return();
    }

    repeat->num_match = nres;
    repeat->match = matches;
    repeat->io = io;
    strcpy(repeat->tagname, CPtr2Tcl(repeat));

    val = get_default_string(GetInterp(), gap5_defs,"FINDREP.COLOUR");
    strcpy(repeat->colour, val);

    repeat->linewidth = get_default_int(GetInterp(), gap5_defs,
					"FINDREP.LINEWIDTH");

    repeat->params = (char *)xmalloc(100);
    if (repeat->params)
	sprintf(repeat->params, "Unknown at present");
    repeat->all_hidden = 0;
    repeat->current = -1;
    repeat->reg_func = repeat_callback;
    repeat->match_type = REG_TYPE_REPEAT;

    /* Create and plot our match array */
    for (i= 0; i < nres; i++){
	matches[i].func = repeat_obj_func;
	matches[i].data = repeat;
        matches[i].c1 = rnumtocnum(io, ABS(c1[i])) * (c1[i] < 1 ? -1 : 1);
        matches[i].pos1 = pos1[i];
        matches[i].c2 = rnumtocnum(io, ABS(c2[i])) * (c2[i] < 1 ? -1 : 1);
        matches[i].pos2 = pos2[i];
        matches[i].length = len[i];
        matches[i].score = len[i];
	matches[i].flags = 0;
    }

    /* Sort matches */
    qsort(repeat->match, repeat->num_match, sizeof(obj_match), sort_func);

    PlotRepeats(io, repeat);
    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(repeat), NULL);

    /*
     * Register the repeat search with each of the contigs used.
     * Currently we assume that this is all.
     */
    id = register_id();
    contig_register(io, 0, repeat_callback, (void *)repeat, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
		    REG_NUMBER_CHANGE | REG_ORDER, REG_TYPE_REPEAT);
    update_results(io);
}


int repeat_search (
	           int mode,		/* 1=f,2=r,3=b */
		   int min_match,	/* the minimum match length */
		   int **seq1_match,	/* positions of matches in seq1 */
		   int **seq2_match,	/* positions of matches in seq2 */
		   int **len_match,	/* length of matches */
		   int max_matches,	/* maximum number of matches */
		   char *seq1,		/* seq1 */
		   int seq1_len,	/* size of seq1 and its hash array */
		   int *num_f_matches,
		   int *num_r_matches
		   );


int
find_repeats(GapIO *io,
             int mode,
             int min_match,
	     int mask,
	     float percd, /* del */
	     int num_contigs,
	     contig_list_t *contig_array,
	     char *out_name)
{
    tg_rec *crec1, *crec2;
    int *pos1, *pos2, *len; /* position in crec1/2 and joint length */
    char *consensus;
    int max_read_length, database_size, number_of_contigs;
    int consensus_length, ret, task_mask;
    int max_matches;
    int i,j;
    int num_f_matches, num_r_matches;
    Contig_parms *contig_list;
    Hidden_params p;

    p.min = p.max = p.verbose = p.use_conf = p.qual_val = p.window_len =0;
    p.test_mode = 0;
    p.start = 0;
    p.lwin1 = 0;
    p.lcnt1 = 0;
    p.rwin1 = 0;
    p.rcnt1 = 0;
    p.do_it = 0;
    p.gap_open = 12;
    p.gap_extend = 4;

    max_matches = 10000; /* FIXME: make this adjustable */
    consensus = NULL;
    contig_list = NULL;
    crec1 = crec2 = NULL;
    pos1 = pos2 = len = NULL;

    if ((pos1 = (int *)xmalloc(max_matches * sizeof(int)))==NULL){
	goto bail_out;
    }
    if ((pos2 = (int *)xmalloc(max_matches * sizeof(int)))==NULL){
	goto bail_out;
    }
    if ((len = (int *)xmalloc(max_matches * sizeof(int)))==NULL){
	goto bail_out;
    }

    max_read_length = find_max_gel_len(io, 0, 0);

    database_size = io_dbsize(io);
    number_of_contigs = num_contigs;
    if ( ! (contig_list = get_contig_list ( database_size, io,
			   number_of_contigs, contig_array ))) {
	goto bail_out;
    }

    task_mask = ADDTITLE | NORMALCONSENSUS;
    if ( mask == 3 ) task_mask |= MASKING;

    consensus_length = 0;

/*    printf("TASK mask %d mask %d mode %d\n",task_mask,mask,mode); */
    if ( ret = make_consensus ( task_mask, io, &consensus, NULL,
			  contig_list, number_of_contigs,
			  &consensus_length,
			  max_read_length,
			  p,
			  consensus_cutoff ) ) {
	goto bail_out;
    }

    ret =  repeat_search ( mode, min_match, &pos2, &pos1, &len, max_matches,
			   consensus, consensus_length,
			   &num_f_matches, &num_r_matches );

    if( ret < 0 ) {
	goto bail_out;
    }

    /* get the output arrays filled in to fit with old routines!!*/
    if ((crec1 = (tg_rec *)xmalloc((num_f_matches + num_r_matches + 1) * sizeof(tg_rec)))==NULL){
	goto bail_out;
    }

    if ((crec2 = (tg_rec *)xmalloc((num_f_matches + num_r_matches + 1) * sizeof(tg_rec)))==NULL){
	goto bail_out;
    }

    for ( i=0;i<num_f_matches; i++ ) {
	if ( -1 != (j = (contig_listel_from_con_pos ( contig_list, number_of_contigs,
					pos1[i] )))) {
	    crec1[i] = contig_list[j].contig_left_gel;
	    pos1[i] += contig_list[j].contig_start
		     - contig_list[j].contig_start_offset - 1;
	}
	else {
	    goto bail_out;
	}
	if ( -1 != (j = (contig_listel_from_con_pos ( contig_list, number_of_contigs,
					pos2[i] )))) {
	    crec2[i] = contig_list[j].contig_left_gel;
	    pos2[i] += contig_list[j].contig_start
		     - contig_list[j].contig_start_offset - 1;
	}
	else {
	    goto bail_out;
	}
    }

    for ( i=num_f_matches; i<num_f_matches+num_r_matches;i++ ) {
	if ( -1 != (j = (contig_listel_from_con_pos ( contig_list, number_of_contigs,
					pos1[i] )))) {
	    crec1[i] = -contig_list[j].contig_left_gel;
	    pos1[i] += contig_list[j].contig_start
		     - contig_list[j].contig_start_offset - 1;;
	}
	else {
	    goto bail_out;
	}
	if ( -1 != (j = (contig_listel_from_con_pos ( contig_list, number_of_contigs,
					pos2[i] )))) {
	    crec2[i] = contig_list[j].contig_left_gel;
	    pos2[i] += contig_list[j].contig_start
		     - contig_list[j].contig_start_offset - 1;
	}
	else {
	    goto bail_out;
	}
    }

/* end new */

    flush2t(io);

    /* FIXME write_tags needs arguments changing from f_int to int */
#if 0
    if (out_name) {
	write_tags(io, out_name, num_f_matches+num_r_matches, crec1, pos1, crec2, pos2, len);
    }
#endif

    plot_rpt(io, num_f_matches+num_r_matches, crec1, pos1, crec2, pos2, len);

	if ( crec1 ) xfree(crec1);
	if ( pos1 )  xfree(pos1);
	if ( crec2 ) xfree(crec2);
	if ( pos2 )  xfree(pos2);
	if ( len )   xfree(len);
	if ( consensus )   xfree(consensus);
	if ( contig_list ) xfree(contig_list);

	return 0;
 bail_out:

	if ( crec1 ) xfree(crec1);
	if ( pos1 )  xfree(pos1);
	if ( crec2 ) xfree(crec2);
	if ( pos2 )  xfree(pos2);
	if ( len )   xfree(len);
	if ( consensus )   xfree(consensus);
	if ( contig_list ) xfree(contig_list);

	return -1;
} /* end FindRepeats */
