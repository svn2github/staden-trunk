#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "io-reg.h"
#include "canvas_box.h"
#include "gap_canvas_box.h"
#include "gap_cli_arg.h"
#include "restriction_enzymes.h"
#include "restriction_enzyme_map.h"
#include "getfile.h"
#include "list_proc.h"
#include "newgap_cmds.h"
#include "misc.h"                                        /* need for strdup */
#include "qual.h"
#include "dna_utils.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void renz_callback(GapIO *io, int contig,
			      void *fdata,
			      reg_data *jdata);
void display_renz(Tcl_Interp *interp, GapIO *io, obj_renz *r);

/*
 *---------------------------------------------------------------------------
 * Registration
 *---------------------------------------------------------------------------
 */
/*
 * Provide information about a set of matches - the meta object.
 */
int
renz_info(char *window,
	  int contig,
	  R_Match *match,
	  int num_match,
	  int sequence_type,
	  R_Enz *r_enzyme,
	  int num_enzymes,
	  char *name, 
	  GapIO *io,
	  int i,
	  int lreg,
	  int rreg,
	  int print_opt) 
{
    int sequence_len;
    char *sequence;
    
    vfuncheader("%s result list", name);
    
    vmessage("Contig %s (#%d) \n", get_contig_name(io, contig),
	     io_clnbr(io, contig));

    vmessage("Number of enzymes = %d\n", num_enzymes);
    vmessage("Number of matches = %d\n", num_match);

    /* generate a consensus sequence of the contig */
    sequence_len = rreg - lreg + 1;

    if (NULL == (sequence = (char *)xmalloc((sequence_len + 1) * 
					   sizeof(char ))))
	return 0;

    calc_consensus(contig, lreg, rreg, CON_SUM, sequence, 
		   NULL, NULL, NULL, consensus_cutoff, quality_cutoff, 
		   database_info, io);
    depad_seq(sequence, &sequence_len, NULL);

    if (print_opt == 0) {
	start_message();
	if (0 == PrintEnzymeByEnzyme(r_enzyme, match, 
				     num_match, 
				     num_enzymes, sequence, 
				     sequence_len, 
				     sequence_type, lreg, 1)) {
	    vmessage("No enzyme cut sites found\n");
	}
	end_message(window);
    } else {
	start_message();
	if (0 == OrderOnPosition(r_enzyme, match, num_match,sequence,
				 sequence_len, sequence_type, lreg)) {
	    vmessage("No enzyme cut sites found\n");
	}
	end_message(window);
    }
    xfree(sequence);
    return 1;
}


/*
 * Callback routines for the restriction enzyme plot in the template display
 */

static int update_obj_t_renz(GapIO *io, obj_t_renz *r, int contig) {
    int i, newl;

    for (i = 0; i < r->num_contigs; i++) {
	if (r->c_match[i].contig == contig) {
	    break;
	}
    }
    newl = ABS(io_clength(io, r->c_match[i].contig));

    r->start = 1;
    r->end = newl;
    r->c_match[i].length = newl;
    return 0;
}


/*
 * Callback routines for the stand alone restriction enzyme plot
 */

static void
update_obj_renz(GapIO *io, obj_renz *r, int contig) {
    int newl;

    newl = ABS(io_clength(io, r->c_match.contig));

    /* if need to replot, lose original lreg to rreg info */
    r->start = 1;
    r->end = newl;
    r->c_match.length = newl;

    /* reset ruler from and to */
    r->ruler->start = 1;
    r->ruler->end = newl;
}

static void gap_renz_shutdown(GapIO *io, obj_renz *r) {
    char cmd[1024];

    contig_deregister(io, r->c_match.contig, renz_callback, (void *)r);

    delete_contig_cursor(io, r->c_match.contig, r->cursor->id, 0);

    renz_shutdown(r->r_enzyme, r->num_enzymes, r->c_match.match, r->canvas, 
		  r->world, r->zoom);

    sprintf(cmd, "DeleteREnzPlot %s %s\n", r->frame, r->window);
    if (TCL_ERROR == Tcl_Eval(r->interp, cmd)) {
	printf("%s\n", Tcl_GetStringResult(r->interp));
    }

    if (r->xhair.colour) free(r->xhair.colour);
    if (r->tick->colour) free(r->tick->colour);
    if (r->text_colour) free(r->text_colour);
    free_win_list(r->win_list, r->num_wins);
    free_ruler_struct(r->ruler);
    xfree(r->tick);

    xfree(r);
}

static void
renz_renumber(GapIO *io,
	      obj_renz *r, 
	      int old_contig,
	      int new_contig)
{
    if (r->c_match.contig == old_contig) {
	r->c_match.contig = r->c_match.contig > 0 ? new_contig : -new_contig;
    }
}

int
renz_replot(Tcl_Interp *interp,
	    GapIO *io,
	    obj_renz *r)
{
    char *sequence;
    int sequence_len;
    R_Match *match;
    int total_matches;
    int *depad_to_pad;
    int i;
    int start, end;

    start = MAX(1, r->start - r->max_overlap);
    end = MIN(ABS(io_clength(io, r->c_match.contig)), r->end + r->max_overlap);

    /* generate a consensus sequence of the contig */
    sequence_len = end - start + 1;
    if (NULL == (sequence = (char *)xmalloc((sequence_len + 1) * 
					    sizeof(char ))))
	return 0;
    if (NULL == (depad_to_pad = (int *)xmalloc((sequence_len + 1) *
					       sizeof(int))))
	return 0;

    /* Compute pad-stripped consensus */
    calc_consensus(r->c_match.contig, start, end, CON_SUM, sequence, 
		   NULL, NULL, NULL, consensus_cutoff, quality_cutoff, 
		   database_info, io);
    depad_seq(sequence, &sequence_len, depad_to_pad);
    
    if (r->c_match.match)
	xfree(r->c_match.match);
	
    if (NULL == (match = (R_Match*)xcalloc(MAXMATCHES, sizeof(R_Match))))
	return 0;

    FindMatches(r->r_enzyme, r->num_enzymes, sequence, sequence_len,
		r->sequence_type, &match, &total_matches);

    /* Reposition matches back into the padded data */
    for (i = 0; i < total_matches; i++) {
	if (match[i].cut_pos >= sequence_len) {
	    match[i].padded_cut_pos = depad_to_pad[sequence_len-1] +
		match[i].cut_pos - sequence_len + 1;
	} else if (match[i].cut_pos >= 0) {
	    match[i].padded_cut_pos = depad_to_pad[match[i].cut_pos];
	} else {
	    match[i].padded_cut_pos = match[i].cut_pos;
	}
	match[i].cut_pos = match[i].cut_pos - (r->start - start);
	match[i].padded_cut_pos = match[i].padded_cut_pos - (r->start - start);
    }

    r->c_match.match = match;
    r->c_match.num_match = total_matches;

    xfree(sequence);
    xfree(depad_to_pad);
    /* display_renz(interp, io, r); */
    plot_renz_matches(interp, r->window, r->names_win, r->text_offset,
		      r->text_colour, r->yoffset, r->num_enzymes,
		      r->r_enzyme, r->ruler,
		      io_clength(io, r->c_match.contig),
		      r->c_match.num_match, r->c_match.match, r->tick,
		      r->frame, r->world, r->canvas, r->win_list,
		      r->num_wins, &r->zoom);
    return 1;
}

/*
 * Callback for the stand alone restriction enzyme plot
 */
static void 
renz_callback(GapIO *io, 
	      int contig,
	      void *fdata,
	      reg_data *jdata) 
{
    obj_renz *r = (obj_renz *)fdata;

    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Restriction enzymes");
	    return;
	}
    case REG_JOIN_TO:
#ifdef DEBUG
	printf("renz JOIN TO contig %d join %d\n", contig, jdata->join.contig);
#endif
	r->c_match.contig = jdata->join.contig;
	/* NOTE - flow through to update_obj_renz() */

    case REG_COMPLEMENT:
    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("renz COMPLEMENT LENGTH contig %d\n", contig);
#endif
	    update_obj_renz(io, r, contig);
	    renz_replot(r->interp, io, r);
	    return;
	}
    case REG_QUIT:
    case REG_DELETE:
	{
	    gap_renz_shutdown(io, r);
	    return;
	}

    case REG_GET_OPS:
	jdata->get_ops.ops = 
	    "Output enzyme by enzyme\0"
	    "Output ordered on position\0"
		"SEPARATOR\0"
	    "Remove\0";

	return;

    case REG_INVOKE_OP:
	switch (jdata->invoke_op.op) {

	case 0: /* Information print enzyme by enzyme */
	    renz_info(r->window, r->c_match.contig, r->c_match.match,
		      r->c_match.num_match, r->sequence_type,
		      r->r_enzyme, r->num_enzymes, 
		      "Restriction Enzymes", io, 0, 
		      r->start, r->end, 0);
	    break;
	case 1: /* Information order */
	    renz_info(r->window, r->c_match.contig, r->c_match.match,
		      r->c_match.num_match, r->sequence_type,
		      r->r_enzyme, r->num_enzymes, 
		      "Restriction Enzymes", io, 0, 
		      r->start, r->end, 1);
	    break;
	case 2: /* Remove */
	    gap_renz_shutdown(io, r);
	    break;
	}
	break;
    case REG_PARAMS:
	/* jdata->params.string = r->params; */
	return;

    case REG_CURSOR_NOTIFY:
	canvas_cursor_refresh(r->interp, io, r->c_match.contig,
			      jdata->cursor_notify.cursor, r->cursor,
			      r->canvas, r->win_list, r->num_wins, r->id,
			      0, &r->cursor_visible, r->world, 1);
	return;
	

    case REG_NUMBER_CHANGE:
#ifdef DEBUG
	printf("renz NUMBER_CHANGE %d\n", jdata->number.number);
#endif
	renz_renumber(io, r, contig, jdata->number.number); 
	renz_replot(r->interp, io, r);
	return; 
    case REG_GENERIC:
	switch (jdata->generic.task) {
	case TASK_CANVAS_SCROLLX: 
	    {
		char *scroll = (char *)jdata->generic.data;
		
		canvasScrollX(r->interp, r->window, r->win_list, r->num_wins,
			      r->world->visible, r->canvas, scroll);
		break;
	    }
	case TASK_CANVAS_SCROLLY:
	    {
		char *scroll = (char *)jdata->generic.data;
		
		canvasScrollY(r->interp, r->window, r->win_list, r->num_wins,
			      r->world->visible, r->canvas, scroll);

		break;
	    }
	case TASK_CANVAS_ZOOMBACK:
	    {
		canvasZoomback(r->interp, r->canvas, r->window, r->world,
			       r->win_list, r->num_wins, &r->zoom);

		/* redraw ruler ticks by redrawing ruler */
		
		draw_single_ruler(r->interp, r->ruler, r->canvas, 
				  r->ruler->start, r->ruler->end, 1);
		scaleSingleCanvas(r->interp, r->world, r->canvas, 
				  r->ruler->window, 'x', "all");
		
		canvas_cursor_refresh(r->interp, io, r->c_match.contig,
				      r->cursor, r->cursor, r->canvas, 
				      r->win_list, r->num_wins, r->id, 0, 
				      &r->cursor_visible, r->world, 0);
		
		break;
	    }
	case TASK_CANVAS_ZOOM:
	    {
		s_zoom *szoom = (s_zoom *)jdata->generic.data;
		canvasZoom(r->interp, r->canvas, r->window, r->world,
			   r->win_list, r->num_wins, &r->zoom, szoom->zoom, 
			   szoom->scroll);

		/* redraw ruler ticks by redrawing ruler */
		
		draw_single_ruler(r->interp, r->ruler, r->canvas,  
				  r->ruler->start, r->ruler->end, 1);
		scaleSingleCanvas(r->interp, r->world, r->canvas, 
				  r->ruler->window, 'x', "all");
		
		canvas_cursor_refresh(r->interp, io, r->c_match.contig,
				      r->cursor, r->cursor, r->canvas, 
				      r->win_list, r->num_wins, r->id, 0, 
				      &r->cursor_visible, r->world, 0);

		break;
	    }
	case TASK_CANVAS_CURSOR_X:
	    {    
		char *label;
		int *cx = (int *)jdata->generic.data;
		double wx, wy;

		CanvasToWorld(r->canvas, *cx, 0, &wx, &wy);
		label = get_default_string(r->interp, gap5_defs,
					   "R_ENZ.CURSOR");

		canvasCursorX(r->interp, r->canvas, r->frame, label, 
			     r->xhair.colour, r->xhair.width, *cx, wx, 
			      r->win_list, r->num_wins);

		break;
	    }
	case TASK_CANVAS_RESIZE:
	    {
		char scroll_args[20];
		/* resize restriction enzyme display window */
		resizeCanvas(r->interp, r->window, r->win_list, r->num_wins,
			      r->world->visible, r->world->total, r->canvas);
		sprintf(scroll_args, "scroll 0 units");
		canvasScrollX(r->interp, r->window, r->win_list, r->num_wins,
			      r->world->visible, r->canvas, scroll_args);
		break;
	    }
	case TASK_CANVAS_WORLD:
	    {
		task_world_t *tw = (task_world_t *)jdata->generic.data;
		int cx = tw->canvasx;
		double wx, wy;

		CanvasToWorld(r->canvas, cx, 0, &wx, &wy);
		tw->basex = wx;
		break;
	    }
	case TASK_RENZ_INFO: 
	    {
		/* info on all cuts by a specific restriction enzyme */
		int *enz_name = (int *)jdata->generic.data;
		int i;
	        R_Match *tmp_match;
		int cnt = 0;
		int sequence_len;
		int sequence_type = 0;
		char *sequence;

#ifdef DEBUG
		printf("TASK_RENZ_INFO name %d\n", *enz_name);
#endif
		/* find the sequence */
		sequence_len = r->end - r->start + 1;
		if (NULL == (sequence = (char *)xmalloc((sequence_len + 1) * 
							sizeof(char ))))
		    return;
		
		calc_consensus(r->c_match.contig, r->start, r->end, CON_SUM, 
			       sequence, NULL, NULL, NULL, consensus_cutoff, 
			       quality_cutoff, database_info, io);
		
		if (NULL == (tmp_match = (R_Match *)malloc(r->c_match.num_match * 
							   sizeof(R_Match))))
		    return;
		
		/* 
		 * make tmp match array containing only those matches for the 
		 * selected enzyme
		 */
		for (i = 0; i < r->c_match.num_match; i++) {
		    if (*enz_name == r->c_match.match[i].enz_name) {
			tmp_match[cnt++] = r->c_match.match[i];
		    }
		}
		start_message();
		PrintEnzymeByEnzyme(r->r_enzyme, tmp_match, cnt, 
				    r->num_enzymes, sequence, sequence_len, 
				    sequence_type, r->start, 0);
		end_message(r->window);
		xfree(sequence);
		xfree(tmp_match);
		break;
	    }
	}
    }
}


/*****************************************************************************/
/*                  stand alone restriction enzyme display                   */
/*****************************************************************************/
/*
 * plot the r enzyme names and cut sites in the stand alone plot
 */
void
display_renz(Tcl_Interp *interp,
	     GapIO *io,
	     obj_renz *r)
{
    char cmd[1024];
    int item;
    int contig_len;
    int offset;
    int y_offset;
    int t_offset;
    int i;
    int cut_site;

    sprintf(cmd, "%s delete all", r->window);
    Tcl_Eval(interp, cmd);

    sprintf(cmd, "%s delete all", r->names_win);
    Tcl_Eval(interp, cmd);

    contig_len = ABS(io_clength(io, r->c_match.contig));

/*    sequence_len = r->rreg - r->lreg + 1; */

    t_offset = r->text_offset;
    y_offset = r->yoffset;

    /* for each selected enzyme, display its name and cut sites */
    for (item = 0; item < r->num_enzymes; item++) {

	sprintf(cmd,"%s create text 10 %d -text %s -anchor w -fill %s "
		"-tag {S re_%d}",
		r->names_win, t_offset, r->r_enzyme[item].name, 
		r->text_colour, item);
	Tcl_Eval(interp, cmd);
	sprintf(cmd, "%s create line %d %d %d %d -tag contig -fill %s",
		r->window, 0, y_offset, contig_len, y_offset, r->ruler->colour);
	Tcl_Eval(interp, cmd);

        offset = r->yoffset;
	for (i = 0; i < r->c_match.num_match; i++) {
	    if (item == r->c_match.match[i].enz_name) {
	
		offset = (item * r->tick->ht) + r->yoffset;
		cut_site = r->start - 1 + r->c_match.match[i].padded_cut_pos;

		PlotStickMap(interp, r->window, cut_site, cut_site, 0, offset, 
			     r->tick->ht, r->tick->line_width, r->tick->colour, 
			     item, 1, contig_len);
	    }
	}
	y_offset += r->tick->ht;
	t_offset += r->tick->ht;
    }

    /* draw final line */
    sprintf(cmd, "%s create line %d %d %d %d -tag contig -fill %s",
	    r->window, 0, y_offset, contig_len, y_offset, r->ruler->colour);
    Tcl_Eval(interp, cmd);

    /* One extra for the blank line at the bottom */
    y_offset += r->tick->ht;

    /* add back selection rectangles */
    if (TCL_ERROR == Tcl_VarEval(interp, "ReSelectRect ", r->frame, " ",
				 r->names_win, NULL))
	printf("display_renz: %s\n", Tcl_GetStringResult(interp));

    r->world->total->x1 = 1;
    r->world->total->x2 = contig_len;
    r->world->total->y1 = 1;
    r->world->total->y2 = y_offset;

    memcpy(r->world->visible, r->world->total, sizeof(d_box));
    r->world->visible->y2 = r->canvas->height;

    SetCanvasCoords(interp, r->world->visible->x1, r->world->visible->y1, 
		    r->world->visible->x2, r->world->visible->y2, r->canvas);

    draw_single_ruler(interp, r->ruler, r->canvas, r->ruler->start, 
		      r->ruler->end, 1);

    scaleCanvas(interp, r->win_list, r->num_wins, "all", r->world->visible, 
		r->canvas);
    scrollRegion(interp, r->win_list, r->num_wins, r->world->total, r->canvas);

    /* remove all current zooming info */
    freeZoom(&r->zoom);

    /* add first zoom */
    pushZoom(&r->zoom, r->world->visible);
}


/*
 * stand alone restriction enzyme display
 * register data and do first plot
 */
int
renz_reg(Tcl_Interp *interp,
	 GapIO *io,
	 char *filename,
	 char *frame,
	 char *names_win,
	 char *re_win,
	 char *inlist,
	 int num_items,
	 int contig_num,
	 int lreg,
	 int rreg,
	 int text_offset,
	 char *text_fill,
	 tick_s *tick,
	 int yoffset,
	 ruler_s *ruler,
	 cursor_s xhair) 
{
    int num_enzymes;
    R_Enz *r_enzyme;
    int sequence_type = 0;                      /* circle if 1; linear if 0 */
    obj_renz *r;
    int id;

    if (NULL == (r = (obj_renz *)xmalloc(sizeof(obj_renz)))) {
	return 0;
    }
    id = register_id();
    r->id = id;

    strcpy(r->window, re_win);
    strcpy(r->frame, frame);
    strcpy(r->names_win, names_win);
    
    r->tick = tick;
    r->ruler = ruler;
    r->xhair = xhair;
    r->interp = interp;

    /* create list of windows in the restriction enzyme display */
    if (NULL == (r->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;

    r->num_wins = 0;
    addWindow(r->win_list, &r->num_wins, r->window, 'b', r->id);
    addWindow(r->win_list, &r->num_wins, r->ruler->window, 'x', r->id);
    addWindow(r->win_list, &r->num_wins, r->names_win, 'y', r->id);

    if (NULL == (r->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (r->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (r->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (r->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    initCanvas(r->interp, r->canvas, r->window);
    createZoom(&r->zoom);

    /*
     * open the file of enzyme names and parse the selected enzymes into the
     * r_enzyme array of size num_enzymes
     */
    open_renz_file(filename, inlist, num_items, &r_enzyme, &num_enzymes);

    r->max_overlap = find_max_cut_dist(r_enzyme, num_enzymes);

    r->c_match.contig = contig_num;
    r->start = lreg;
    r->end = rreg;
    r->c_match.match = NULL;

    r->r_enzyme = r_enzyme;
    r->num_enzymes = num_enzymes;
    r->tick = tick;
    r->yoffset = yoffset;
    r->sequence_type = sequence_type;
    r->text_offset = text_offset;
    r->text_colour = strdup(text_fill);
    r->cursor = create_contig_cursor(io, contig_num, 0, id);
    r->cursor_visible = 0;

    renz_replot(interp, io, r);

    contig_register(io, contig_num, renz_callback, (void *)r, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS
		    | REG_CURSOR_NOTIFY | REG_NUMBER_CHANGE | REG_GENERIC,
		    REG_TYPE_RESTRICTION);
    canvas_cursor_refresh(r->interp, io, r->c_match.contig,
			  r->cursor, r->cursor,
			  r->canvas, r->win_list, r->num_wins, r->id,
			  0, &r->cursor_visible, r->world, 1);

    return id;
}

/*
 * create tags from stand alone restriction enzyme display
 */
int
Create_REnz_Tags(GapIO *io,
		 int contig_num,
		 obj_renz *r,
		 char *enz_list,
		 char **enz_ids,
		 int num_ids)
{
#if 0
    int i,j;
    int seq_len;
    int start;
    int length;
    char comments[1024];
    char *seq;
    char num[4];
    reg_anno ra;
    char *item_str;
    int item;
    int cnt = 0;
    int n_tags = 0;
    
    /* check on status of contig selector - only allow addition of tags if
     * the contig editor is not being used
     */
    if (contig_lock_write(io, contig_num) == -1) {
	verror(ERR_WARN, "create restriction enzyme tags", "Contig is busy");
	return -1;
    }

    /* set inlist to be the active list */
    if (-1 == set_active_list(enz_list)) {
	return -1;
    }

    /* get the first item on the inlist */
    item_str = get_active_list_item();
    item = atoi(item_str); /* convert to integer */

    while (item_str) {

	for (i = 0; i < r->c_match.num_match; i++) {

	    /* found a match for the enzyme */
	    /* note that r->c_match contains an entry for each cut site */
	    if (r->c_match.match[i].enz_name == item) {

		strcpy(comments, r->r_enzyme[item].name);
		/* number of recognition sequences */
		for (j = 0; j < r->r_enzyme[item].num_seq; j++) {
		    
		    if (r->c_match.match[i].enz_seq == j) {
#ifdef DEBUG
			printf("cut_pos %d cut_site %d seq %s\n", 
			       r->c_match.match[i].padded_cut_pos,
			       r->r_enzyme[item].cut_site[j],
			       r->r_enzyme[item].seq[j]);
#endif
			start = (r->start - 1 +
				 r->c_match.match[i].padded_cut_pos) - 
			    r->r_enzyme[item].cut_site[j];
			seq_len = strlen(r->r_enzyme[item].seq[j]);
			/* 
			   if (r->r_enzyme[item].cut_site[j] < 0) { 
			   length = seq_len -  r->r_enzyme[item].cut_site[j];
			   start = r->c_match[i].padded_cut_pos;
			   } else if(r->r_enzyme[item].cut_site[j]>= seq_len) {
			   length = r->r_enzyme[item].cut_site[j];
			   } else {
			   length = seq_len;
			   }
			   */
			length = seq_len;
			seq = AddCutSites(r->r_enzyme[item].seq[j],  
					  r->r_enzyme[item].cut_site[j]);
			strcat(comments, "\n");
			strcat(comments, seq);
			strcat(comments, "\t");
			sprintf(num, "%d", r->r_enzyme[item].cut_site[j]);
			strcat(comments, num);
			strcat(comments, "\n");
			insert_NEW_tag(io, (tag_id)-contig_num, start, length, 
				       enz_ids[cnt], comments, 2);
			n_tags++;
		    }
		}
	    }
	}
	item_str = get_active_list_item();
	cnt++;
	if (!item_str)
	    break;
	item = atoi(item_str); /* convert to integer */
    }
    ra.job = REG_ANNO;
    contig_notify(io, contig_num, (reg_data *)&ra);
    return (n_tags);
#else
    return 0;
#endif
}

/*
 *---------------------------------------------------------------------------
 * Tcl command interfaces.
 *---------------------------------------------------------------------------
 */

typedef struct {
    char *filename;
} read_enz_arg;

typedef struct {
    GapIO *io;
    int id;
    char *filename;
    char *frame;
    char *plot;
    char *inlist;
    int num_items;
    char *contig;
    int tick_ht;
    int tick_wd;
    char *tick_fill;
    int yoffset;
} single_enz_arg;

typedef struct {
    GapIO *io;
    char *filename;
    char *frame;
    char *win_name;
    char *plot;
    char *win_ruler;
    char *inlist;
    int num_items;
    char *contigs;
    int text_offset;
    char *text_fill;
    int tick_ht;
    int tick_wd;
    char *tick_fill;
    int cursor_wd;
    char *cursor_fill;
    int yoffset;
} renz_arg;

typedef struct {
    int item;
    GapIO *io;
    int id;
    tg_rec contig;
} enz_name_arg;

typedef struct {
    int enzyme;
    GapIO *io;
    int id;
    tg_rec contig;
} enz_info_arg;

typedef struct {
    GapIO *io;
    char *contigs;
    int id;
    char *enz_list;
    char *id_list;
} enz_tag_arg;

typedef struct {
    GapIO *io;
    int id;
    int cx;
    int cy;
    tg_rec cnum;
} t_cursor_arg;

int tcl_read_enz_file(ClientData clientData,
		      Tcl_Interp *interp,
		      int objc,
		      Tcl_Obj *CONST objv[])
{
    read_enz_arg args;
    int num_enzymes;
    int i;
    char **names;

    cli_args a[] = {
	{"-file",   ARG_STR, 1, NULL, offsetof(read_enz_arg, filename)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    /* open_r_enz_file(args.filename, &r_enzyme, &num_enzymes); */
    if (0 == r_enz_file_names(args.filename, &names, &num_enzymes))
	return TCL_OK;

    for (i =0; i < num_enzymes; i++) {
	Tcl_AppendElement(interp, names[i]);
	xfree(names[i]);
    }

    if (num_enzymes)
	xfree(names);

    return TCL_OK;
}

/*
 * finds and plots cut sites generated by selected enzymes provided in 'inlist'
 * from a file (filename) on a contig consensus sequence given by 'contig_num'
 * The names of each enzyme are display in canvas 'win_name' and the cuts in
 * canvas 'plot'
 */
int
PlotREnz(ClientData clientData,
	 Tcl_Interp *interp,
	 int objc,
	 Tcl_Obj *CONST objv[])
{
    renz_arg args;
    contig_list_t *contigs;
    int num_contigs;
    ruler_s *ruler;
    int id;
    tick_s *tick;
    cursor_s cursor;

    cli_args a[] = {
	{"-io",		 ARG_IO,  1, NULL, offsetof(renz_arg, io)},
	{"-file",	 ARG_STR, 1, NULL, offsetof(renz_arg, filename)},
	{"-frame",	 ARG_STR, 1, NULL, offsetof(renz_arg, frame)},
	{"-win_names",	 ARG_STR, 1, NULL, offsetof(renz_arg, win_name)},
	{"-window",	 ARG_STR, 1, NULL, offsetof(renz_arg, plot)},
	{"-win_ruler",	 ARG_STR, 1, NULL, offsetof(renz_arg, win_ruler)},
	{"-enzymes",	 ARG_STR, 1, NULL, offsetof(renz_arg, inlist)},
	{"-num_enzymes", ARG_INT, 1, NULL, offsetof(renz_arg, num_items)},
	{"-contigs",	 ARG_STR, 1, NULL, offsetof(renz_arg, contigs)},
	{"-text_offset", ARG_INT, 1, NULL, offsetof(renz_arg, text_offset)},
	{"-text_fill",   ARG_STR, 1, NULL, offsetof(renz_arg, text_fill)},
	{"-tick_height", ARG_INT, 1, "-1", offsetof(renz_arg, tick_ht)},
	{"-tick_width",  ARG_INT, 1, "-1", offsetof(renz_arg, tick_wd)},
	{"-tick_fill",   ARG_STR, 1,   "", offsetof(renz_arg, tick_fill)},
	{"-cursor_width",ARG_INT,1, "-1", offsetof(renz_arg, cursor_wd)},
	{"-cursor_fill", ARG_STR, 1,  "", offsetof(renz_arg, cursor_fill)},
	{"-yoffset",	 ARG_INT, 1, NULL, offsetof(renz_arg, yoffset)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncgroup(5, "restriction enzymes");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contigs);

    if (num_contigs != 1) {
	printf("ONLY DEAL WITH SINGLE CONTIG \n");
    }

    cursor = cursor_struct(interp, tk_utils_defs, "R_ENZ", args.cursor_wd,
			   args.cursor_fill);

    tick = tick_struct(interp, tk_utils_defs, "R_ENZ", args.tick_wd,
		       args.tick_ht, args.tick_fill);

    ruler = ruler_struct(interp, tk_utils_defs, "R_ENZ", 0);

    ruler->start = contigs[0].start;
    ruler->end = contigs[0].end;
    strcpy(ruler->window, args.win_ruler);

    id = renz_reg(interp, args.io, args.filename, args.frame,
		  args.win_name, args.plot, args.inlist, args.num_items,
		  contigs[0].contig, contigs[0].start, contigs[0].end,
		  args.text_offset, args.text_fill, tick, args.yoffset,
		  ruler, cursor);
    vTcl_SetResult(interp, "%d", id);
    xfree(contigs);
    return TCL_OK;
}

/*
 * return the name of an enzyme given a file of names and an index into that
 * file
 */
int
GetREnzName(ClientData clientData,
	    Tcl_Interp *interp,
	    int objc,
	    Tcl_Obj *CONST objv[])
{
    enz_name_arg args;
    obj_t_renz *r;

    cli_args a[] = {
	{"-io",	    ARG_IO, 1, NULL, offsetof(enz_name_arg, io)},
	{"-id",	    ARG_INT, 1, NULL, offsetof(enz_name_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(enz_name_arg, item)},
	{"-cnum",   ARG_REC, 1, NULL, offsetof(enz_name_arg, contig)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    r = result_data(args.io, args.id, args.contig);
    /* printf("name %s \n", r->r_enzyme[args.item].name); */
    if (r) {
	vTcl_SetResult(interp, "%s", r->r_enzyme[args.item].name);
    } else {
	vTcl_SetResult(interp, "No renz plot for id %d, contig %d\n",
		       args.id, args.contig);
	return TCL_ERROR;
    }

    return TCL_OK;
}

/*
 * return info on an enzyme given a contig_register id and an item index into
 * the r_enzyme array
 */
int
GetREnzInfo(ClientData clientData,
	    Tcl_Interp *interp,
	    int objc,
	    Tcl_Obj *CONST objv[])
{
    enz_info_arg args;
    reg_generic gen;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(enz_info_arg, io)},
	{"-id",	    ARG_INT, 1, NULL, offsetof(enz_info_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(enz_info_arg, enzyme)},
	{"-cnum",   ARG_REC, 1, NULL, offsetof(enz_info_arg, contig)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_RENZ_INFO;
    gen.data = (void *)&args.enzyme;

    vfuncgroup(5, "restriction enzymes");
    result_notify(args.io, args.id, (reg_data *)&gen, args.contig);

    return TCL_OK;
}

int
CreateREnzTags(ClientData clientData,
	       Tcl_Interp *interp,
	       int objc,
	       Tcl_Obj *CONST objv[])
{
    contig_list_t *contigs;
    int num_contigs;
    enz_tag_arg args;
    obj_renz *r;
    int result;
    char **enz_ids = NULL;
    int num_ids;

    cli_args a[] = {
	{"-io",	      ARG_IO,  1, NULL, offsetof(enz_tag_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(enz_tag_arg, id)},
	{"-contigs",  ARG_STR, 1, NULL, offsetof(enz_tag_arg, contigs)},
	{"-enum",     ARG_STR, 1, NULL, offsetof(enz_tag_arg, enz_list)},
	{"-tag_types",ARG_STR, 1, NULL, offsetof(enz_tag_arg, id_list)},
	{NULL,	      0,       0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contigs);
    if (num_contigs == 0 || !contigs) {
	if (contigs) xfree(contigs);
	return TCL_OK;
    }
    if (num_contigs != 1) {
	printf("ERROR: only supported for single contig. Processing first"
	       " contig only\n");
    }

    r = result_data(args.io, args.id, contigs[0].contig);

    /* create reading name array */
    if (Tcl_SplitList(interp, args.id_list, &num_ids, &enz_ids)
	!= TCL_OK) {
	return TCL_ERROR;
    }

    result = Create_REnz_Tags(args.io, contigs[0].contig, r, args.enz_list,
			      enz_ids, num_ids);
    vTcl_SetResult(interp, "%d", result);
    xfree(contigs);
    Tcl_Free((char *)enz_ids);
    return TCL_OK;
}

int Canvas_To_World(ClientData clientData,
		    Tcl_Interp *interp,
		    int objc,
		    Tcl_Obj *CONST objv[])
{
    t_cursor_arg args;
    reg_generic gen;
    double wx;
    task_world_t tw;

    cli_args a[] = {
	{"-io",	  ARG_IO,  1, NULL, offsetof(t_cursor_arg, io)},
	{"-id",   ARG_INT, 1, NULL, offsetof(t_cursor_arg, id)},
	{"-cnum", ARG_REC, 1, "0",  offsetof(t_cursor_arg, cnum)},
	{"-x",    ARG_INT, 1, NULL, offsetof(t_cursor_arg, cx)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_WORLD;
    tw.canvasx = args.cx;
    tw.cnum = args.cnum;
    gen.data = (void *)&tw;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    wx = ((task_world_t *)gen.data)->basex;
    vTcl_SetResult(interp, "%d", (int)wx);
    return TCL_OK;
}

int REnz_Init(Tcl_Interp *interp) {
    Tcl_CreateObjCommand(interp, "read_enz_file", tcl_read_enz_file,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "plot_renz", PlotREnz,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "get_r_enz_name", GetREnzName,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "get_r_enz_info", GetREnzInfo,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "create_renz_tags", CreateREnzTags,
			 (ClientData) NULL, NULL);

    Tcl_CreateObjCommand(interp, "canvas_to_world", Canvas_To_World,
			 (ClientData)NULL, NULL);
    return TCL_OK;
}
