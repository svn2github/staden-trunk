#include <tk.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#include "IO.h"
#include "io-reg.h"
#include "canvas_box.h"
#include "gap_canvas_box.h"
#include "newgap_cmds.h"
#include "gap_globals.h"
#include "misc.h"
#include "confidence_graph.h"
#include "tclXkeylist.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "consen.h"
#include "qual.h" 
#include "consistency_display.h"
#include "text_output.h"
#include "ruler_display.h"

/*
 * Maintains confidence graph with the contig registration scheme.
 * The graph is drawn as part of a consistency display
 */

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void confidence_callback(GapIO *io, int contig, void *fdata,
				reg_data *jdata);


/*
 *---------------------------------------------------------------------------
 * Internal functions
 *---------------------------------------------------------------------------
 */


/*
 */
void plot_confidence(Tcl_Interp *interp,
		     float *qual,
		     int seq_len,
		     char *c_win,
		     GapIO *io,
		     int offset,
		     int width,
		     char *colour,
		     float min,
		     float max)
{
    int i;
    /*
     * "%d %f %d %f" (see below) is at most 10+1+>8+1+10+1+>8+1 chars => 40.
     * (>8 as 8 for 0.xxxx, more for 10.xxxx. Guess at max 100.xxx = 44 chars)
     * Therefore, cmd needs to be at least 100*44 long.
     * Set to 10000 for paranoia!
     */
    char cmd[10000], *tmpp;
    int start, end, len;

    start = 0; end = 100;
    if (seq_len < 100)
        end = seq_len-1;

    if (strcmp(get_default_string(interp, gap_defs, 
				  "CONFIDENCE_GRAPH.PLOT_TYPE"), "line")==0) { 
	while (start  < seq_len-1) {
	    len = sprintf(cmd, "%s create line ", c_win);
	    tmpp = cmd + len;

	    for (i = start; i < end; i++) {
		len = sprintf(tmpp, "%d %f %d %f ",
			      i+offset, max - qual[i] + min, 
			      i+1+offset, max - qual[i+1] + min);
		tmpp += len;
	    }

	    sprintf(tmpp, "-fill %s -width %d", colour, width);
	    
	    start += 100;
	    end += 100;
	    if (end > seq_len-1)
		end = seq_len-1;
	    
	    Tcl_Eval(interp, cmd);
	}
    } else {
	for (i = 0; i < seq_len-1; i++) {
	    sprintf(cmd, "%s create line %d %f %d %f -fill %s -width %d -capstyle round", 
		    c_win, i+offset, max - qual[i] + min, i+1+offset, 
		    max - qual[i] + min, colour, width);
	    
	    Tcl_Eval(interp, cmd);
	}
    }

#if 0
    for (i = 0; i < seq_len-1; i++) {
      /*
        sprintf(cmd, "%s create line %d %f %d %f -fill %s -width %d", 
		c_win, i+offset, max - qual[i] + min, i+1+offset, 
		max - qual[i+1] + min, colour, width);
      */
        sprintf(cmd, "%s create line %d %f %d %f -fill %s -width %d -capstyle round", 
		c_win, i+offset, max - qual[i] + min, i+1+offset, 
		max - qual[i] + min, colour, width);

#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);
    }
#endif
}

int calc_confidence(GapIO *io,
		    int contig,
		    int start,
		    int end,
		    int mode,
		    float *qual, /* out */
		    float *min,  /* out */
		    float *max)  /* out */
{
    char *con;
    int i;

    con = (char *)xmalloc((end - start + 1) * sizeof(*con));

    if (!con)
        return -1;

    calc_consensus(contig, start, end,
		   CON_SUM, con, NULL,
		   mode == CONFIDENCE_GRAPH_DISCREP ? NULL : qual,
		   mode == CONFIDENCE_GRAPH_DISCREP ? qual : NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);
    
    for (i = 0; i < end - start + 1; i++) {
        if (qual[i] > *max)
	    *max = qual[i];
	if (qual[i] < *min)
	    *min = qual[i];

#ifdef DEBUG    
        printf("%d %f\n", i, qual[i]);
#endif
    }
    /* always set minimum to 0.0 */
    *min = 0.0;
    xfree(con);
    return 0;
}

/*
 * plot the vertical ruler
 */
void plot_confidence_ruler(Tcl_Interp *interp,
			   obj_confidence_graph *conf,
			   CanvasPtr *canvas,
			   WorldPtr *world)
{

    display_ruler_v(interp, canvas, conf->ruler, conf->t_min, conf->t_max);
    scaleSingleCanvas(interp, world, canvas, conf->ruler->window, 'y', "all");  
}

/*
 * confidence graph
 */
void display_confidence_graph(GapIO *io, obj_confidence_graph *conf) {
    char cmd[1024];
    int i;
    obj_consistency_disp *c;
    int length;
    int win_num;

    c = result_data(io, conf->cons_id, 0);

    sprintf(cmd, "%s delete all", conf->c_win);
    Tcl_Eval(c->interp, cmd);

    win_num = get_consistency_win_num(c, conf->id);
    for (i = 0; i < c->num_contigs; i++) {
        if (c->num_contigs > 1) {
	    length = ABS(io_clength(io, c->contigs[i]));
	} else {
	    length = c->end - c->start + 1;
	}
        plot_confidence(c->interp, conf->qual[i], length, conf->c_win, io, 
			c->start + c->contig_offset[c->contigs[i]].offset, 
			conf->linewidth, conf->colour, conf->t_min, 
			conf->t_max);
    }

    plot_confidence_ruler(c->interp, conf, c->win_list[win_num]->canvas,
			  c->win_list[win_num]->world);

    scaleCanvas(c->interp, &c->win_list[win_num], 1, "all", 
		c->win_list[win_num]->world->visible, 
		c->win_list[win_num]->canvas); 

    scrollRegion(c->interp, &c->win_list[win_num], 1, 
		 c->win_list[win_num]->world->total, 
		 c->win_list[win_num]->canvas);

    consistency_update_cursors(io, c, 0);
}

/*
 * Removes the confidence graph (and unplots etc).
 */
static void confidence_shutdown(GapIO *io, obj_confidence_graph *conf) {
    char cmd[1024];
    int i;
    obj_consistency_disp *c;
    int win_num;

    if (NULL != (c = result_data(io, conf->cons_id, 0))) {
	win_num = get_consistency_win_num(c, conf->id);
	delete_consistency_window(c, win_num);
    }
    for (i = 0; i < c->num_contigs; i++) {
      contig_deregister(io, c->contigs[i], confidence_callback, 
			(void *)conf);
    }
    sprintf(cmd, "DeleteConfidenceGraph %d %s %s %d\n", *handle_io(io), 
	    conf->frame, conf->c_win, conf->cons_id);
    if (TCL_ERROR == Tcl_Eval(c->interp, cmd)) 
      printf("confidence_shutdown: %s\n", c->interp->result);

    if (conf->qual) {
        for (i = 0; i < c->num_contigs; i++) {
	  xfree(conf->qual[i]);
	}
	xfree(conf->qual);
    }

    if (conf->min)
	xfree(conf->min);
    if (conf->max)
	xfree(conf->max);

    free_ruler_struct(conf->ruler);

    xfree(conf);
    if (c->num_wins == 0)
        consistency_shutdown(io, c);
    
}

/*
 * Updates an obj_confidence_graph structure.
 * This includes recalculating the length, reallocing the memory used,
 * and recalculating the consensus sequence.
 *
 * Returns 0 for success, -1 for failure.
 */
static int update_obj_confidence_graph(GapIO *io, 
				       obj_confidence_graph *conf, 
				       int contig) {
    int i, ret, newl;
    obj_consistency_disp *c;
    int start, end;
    int win_num;

#ifdef DEBUG
    printf("update_obj_confidence_graph\n");
#endif
    c = result_data(io, conf->cons_id, 0);
    
    for (i = 0; i < c->num_contigs; i++) {
	if (c->contigs[i] == contig) {
	    break;
	}
    }

    /* 
     * if single contig, start & end are a range otherwise set start & end
     * to length of contig
    */
    if (c->num_contigs == 1) {
        start = c->start;
	end = c->end;
    } else {
        start = 1;
	end = ABS(io_clength(io, c->contigs[i]));
    }
    newl = end - start + 1;
    if (conf->qual[i]) {
      if (NULL == (conf->qual[i] = 
		   (float *)xrealloc(conf->qual[i], newl*sizeof(float))))
	return -1;
    }

    /* re-initialise */
    conf->max[i] = -FLT_MIN;
    conf->min[i] = FLT_MAX;

    conf->t_max = -FLT_MIN;
    conf->t_min = FLT_MAX;

    ret = calc_confidence(io, c->contigs[i], start, end, conf->mode,
			  conf->qual[i], &conf->min[i], &conf->max[i]);

    /* find total min and max over all contigs */
    for (i = 0; i < c->num_contigs; i++) {
#ifdef DEBUG
        printf("max[%d]=%f\n", i, conf->max[i]);
#endif
        if (conf->t_max < conf->max[i])
	    conf->t_max = conf->max[i];
    }

    conf->t_min = 0;

    /* need to reset y coords */
    win_num = get_consistency_win_num(c, conf->id);
    c->win_list[win_num]->world->total->y1 = conf->t_min;
    c->win_list[win_num]->world->total->y2 = conf->t_max;
    c->win_list[win_num]->world->visible->y1 = conf->t_min;
    c->win_list[win_num]->world->visible->y2 = conf->t_max;
    
    return ret;
}

#if 0
static void confidence_list(GapIO *io, 
			    obj_consistency_disp *c,
			    obj_confidence_graph *conf,
			    int contig)
{
    int i, j;
    int length;

    for (i = 0; i < c->num_contigs; i++) {
        if (c->num_contigs > 1) {
	    length = ABS(io_clength(io, c->contigs[i]));
	} else {
	    length = c->end - c->start + 1;
	}
	vmessage("Contig %s (#%d)\n", 
		 get_contig_name(io, contig), io_clnbr(io, contig));
	for (j = 0; j < length; j++) {
	  vmessage("%d %f\n", j + c->start, conf->qual[i][j]);
	}
	vmessage("\n");
    }
}
#endif

/* scroll the confidence ruler */
static void confidenceScrollY(Tcl_Interp *interp,
			      obj_confidence_graph *conf,
			      CanvasPtr *canvas,
			      d_box *visible,
			      char *scroll_args)
{
    char cmd[1024];
    double wx;
  
    sprintf(cmd, "eval %s yview %s", conf->ruler->window, scroll_args);
    Tcl_Eval(interp, cmd);

    /* find new top edge of canvas in canvasy coords */
    Tcl_VarEval(interp, conf->ruler->window, " canvasy 0", NULL);
    canvas->y = atoi(interp->result);

    /* find new top and bottom edges of canvas in world coords */
    CanvasToWorld(canvas, 0, canvas->y, &wx, &visible->y1);
    CanvasToWorld(canvas, 0, canvas->y + canvas->height, &wx, &visible->y2);
    SetCanvasCoords(interp, visible->x1, visible->y1, 
		    visible->x2, visible->y2, canvas);
}

static void confidence_join(GapIO *io,
			    obj_consistency_disp *c,
			    obj_confidence_graph *conf,
			    int index)
{
    int length;

    length = c->num_contigs - index - 1;

    if (conf->qual[index])
	xfree(conf->qual[index]);
    memmove(&(conf->qual)[index], &(conf->qual)[index+1], 
	    length * sizeof(float*));
    memmove(&(conf->min)[index], &(conf->min)[index+1], 
	    length * sizeof(int));
    memmove(&(conf->max)[index], &(conf->max)[index+1], 
	    length * sizeof(int));
}

/*
 * The callback from the contig registration scheme.
 */
static void confidence_callback(GapIO *io, int contig, void *fdata,
				reg_data *jdata) {
    obj_confidence_graph *conf = (obj_confidence_graph *)fdata;
    obj_consistency_disp *c;
    int win_num;

    c = result_data(io, conf->cons_id, 0);

    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, 
		    (conf->mode == CONFIDENCE_GRAPH_DISCREP)
		    ? "Discrepancy graph"
		    : "Confidence graph");
	    return;
	}
    case REG_COMPLEMENT:
    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("confidence COMPLEMENT LENGTH contig %d\n", contig);
#endif
	    if (update_obj_confidence_graph(io, conf, contig) == 0) {
		display_confidence_graph(io, conf);
		return;
	    }
	    /* else
	     *    flow through to REG_QUIT
	     */
        }
    case REG_QUIT:
    case REG_DELETE:
	{
#ifdef DEBUG
	    printf("confidence QUIT DELETE \n");
#endif
	    confidence_shutdown(io, conf);
	    return;
	}

    case REG_GET_OPS:
	{
	    jdata->get_ops.ops =
#if 0
		"List\0"
		    "SEPARATOR\0"
#endif
		"Remove\0";
	    return;
	}
    case REG_INVOKE_OP:
	{
	    switch (jdata->invoke_op.op) {
#if 0
	    case 0: /* List Quality */
	        start_message();
		confidence_list(io, c, conf, contig);
		end_message();
		break;
#endif
	    case 0: /* Remove */
		confidence_shutdown(io, conf);
		break;
	    }
	    return;
	}
    case REG_PARAMS:
	{
	    puts("REG_PARAMS");
	    return;
	}
    case REG_GENERIC:
	switch (jdata->generic.task) {
	case TASK_CONS_ID:
	  jdata->generic.data = (void *)conf->cons_id;
	case TASK_DISPLAY_TICKS:
	    {
		/* int *ticks = (int *)jdata->generic.data; */
		char cmd[100];
		int win_num;
#ifdef DEBUG
		printf("Confidence display_ticks\n");
#endif
		/* get confidence window */
		win_num = get_consistency_win_num(c, conf->id);

		sprintf(cmd, "%s delete all", conf->ruler->window);
		Tcl_Eval(c->interp, cmd);
		plot_confidence_ruler(c->interp, conf, 
				      c->win_list[win_num]->canvas,
				      c->win_list[win_num]->world);
		break;
	    }
	case TASK_CANVAS_SCROLLY: 
	    {
		char *scroll = (char *)jdata->generic.data;

		win_num = get_consistency_win_num(c, conf->id);

		/* scroll the ruler */
		confidenceScrollY(c->interp, conf, 
				  c->win_list[win_num]->canvas,
				  c->win_list[win_num]->world->visible,
				  scroll);
		
		/* scroll the confidence plot */
		consistencyScrollY(c->interp, conf->c_win, 
				   c->win_list[win_num]->scroll,
				   c->win_list[win_num]->world->visible, 
				   c->win_list[win_num]->canvas, scroll);
		
		break;
	    }
	case TASK_CANVAS_REDRAW:
#ifdef DEBUG
	    printf("confidence TASK_CANVAS_REDRAW\n");
#endif
	    display_confidence_graph(io, conf);
	    break;
	case TASK_CANVAS_RESIZE:
	  {
	      char cmd[1024];
	      int win_num;
	      /* resize confidence graph window */
#ifdef DEBUG
	    printf("Confidence TASK_CANVAS_RESIZE\n");
#endif
	    win_num = get_consistency_win_num(c, conf->id);
	    sprintf(cmd, "%s delete all", conf->ruler->window);
	    Tcl_Eval(c->interp, cmd);
	    plot_confidence_ruler(c->interp, conf, 
				  c->win_list[win_num]->canvas, 
				  c->win_list[win_num]->world);
	  }
	    break;

	case TASK_CANVAS_WORLD:
	    {
		task_world_t *tw = (task_world_t *)jdata->generic.data;
		int cx = tw->canvasx;
		double wx, wy;
#ifdef DEBUG
		printf("CONFIDENCE TASK_CANVAS_WORLD\n");
#endif
		win_num = get_consistency_win_num(c, conf->id);
		CanvasToWorld(c->win_list[win_num]->canvas, cx, 0, &wx, &wy);
		tw->basex = wx - c->contig_offset[tw->cnum].offset;
		break;
	    }
	case TASK_CONS_WORLD:
	    {
		d_box *cw = (d_box *)jdata->generic.data;

		win_num = get_consistency_win_num(c, conf->id);
		cw->x1 = c->win_list[win_num]->world->total->x1;
		cw->y1 = c->win_list[win_num]->world->total->y1;
		cw->x2 = c->win_list[win_num]->world->total->x2;
		cw->y2 = c->win_list[win_num]->world->total->y2;
		
		break;
	    }
	case TASK_CONS_JOIN:
	    {
#ifdef DEBUG
	      printf("confidence JOIN_TO contig %d join %d\n", 
		     contig, jdata->join.contig); 
#endif
	      confidence_join(io, c, conf, (int)jdata->generic.data);
	      break;
	    }
	case TASK_CONS_CURSOR_DELETE:
	    {
	      Tcl_VarEval(c->interp, conf->c_win, " delete cursor_y", NULL);
	      break;
	    }

	case TASK_CANVAS_CURSOR_Y:
	    {    
		char *label;
		int *cy = (int *)jdata->generic.data;
		double wx, wy;
		int win_num;

		win_num = get_consistency_win_num(c, conf->id);
		CanvasToWorld(c->win_list[win_num]->canvas, 0, *cy, &wx, &wy);
		wy = conf->t_max - wy + conf->t_min;

		label = get_default_string(c->interp, gap_defs,
					   "CONSISTENCY_DISPLAY.CURSORY");
		canvasCursorY(c->interp, c->win_list[win_num]->canvas, 
			      c->frame, label, 
			      c->xhair.colour, c->xhair.width, *cy, wy, 
			      c->win_list, c->num_wins);
		break;
	    }
	}
    }
}


/*
 * Registers and initialises a quality buffer for a particular contig.
 */
int confidence_graph_reg(GapIO *io, 
			 Tcl_Interp *interp,
			 char *frame,
			 char *conf_win,
			 int cons_id,
			 ruler_s *ruler,
			 int mode) 
{
    obj_confidence_graph *conf;
    obj_consistency_disp *c;
    int id;
    int i;
    int start, end, length;
    char *val;

    c = result_data(io, cons_id, 0);

    /* first check that there are enough available windows */
    if (c->num_wins >= MAX_NUM_WINS)
	return -1;

    if (NULL == (conf = (obj_confidence_graph *)xmalloc(sizeof(obj_confidence_graph))))
	return -1;

    if (NULL == (conf->qual = (float **)xmalloc(c->num_contigs*sizeof(float*))))
	return -1;

    if (NULL == (conf->min = (float *)xmalloc(c->num_contigs*sizeof(float))))
	return -1;
    if (NULL == (conf->max = (float *)xmalloc(c->num_contigs*sizeof(float))))
	return -1;

    id = register_id();
    
    conf->id = id;
    conf->cons_id = cons_id;
    strcpy(conf->c_win, conf_win);
    strcpy(conf->frame, frame);
    conf->linewidth = get_default_int(interp, gap_defs,
				      "CONFIDENCE_GRAPH.LINEWIDTH");
    val = get_default_string(interp, gap_defs, "CONFIDENCE_GRAPH.COLOUR"); 
    strcpy(conf->colour, val);

    conf->t_max = -FLT_MIN;
    conf->t_min = FLT_MAX;
    conf->ruler = ruler;
    conf->mode = mode;
    
    for (i = 0; i < c->num_contigs; i++) {

        if (c->num_contigs > 1) {
	    length = ABS(io_clength(io, c->contigs[i]));
	    start = 1;
	    end = length;
	} else {
	    length = c->end - c->start + 1;
	    start = c->start;
	    end = c->end;
	}

	if (NULL==(conf->qual[i]=(float *)xmalloc(length * sizeof(float)))) {
	    return -1;
	}

	conf->max[i] = -FLT_MIN;
	conf->min[i] = FLT_MAX;

        calc_confidence(io, c->contigs[i], start, end, conf->mode,
			conf->qual[i], &conf->min[i], &conf->max[i]);
	/* find total min and max over all contigs */
        if (conf->t_max < conf->max[i]) {
#ifdef DEBUG
	  printf("max[%d]=%f\n", i, conf->max[i]);
#endif
	  conf->t_max = conf->max[i];
	}
	conf->t_min = 0;
    }

    add_consistency_window(io, c, conf_win, 'b', id, c->orig_total->x1,
			   conf->t_min, c->orig_total->x2, conf->t_max);

    display_confidence_graph(io, conf);

    for (i = 0; i < c->num_contigs; i++) {
        contig_register(io, c->contigs[i], confidence_callback, 
			(void *)conf, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_CURSOR_NOTIFY | REG_GENERIC, 
			REG_TYPE_CONFIDENCE);
    }
#ifdef DEBUG
    printf("CONFIDENCE ID %d\n", id);
#endif
    return id;
}
