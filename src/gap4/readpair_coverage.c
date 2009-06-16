#include <tk.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#include "io-reg.h"
#include "template_display.h"
#include "IO.h"
#include "canvas_box.h"
#include "gap_canvas_box.h"
#include "newgap_cmds.h"
#include "gap_globals.h"
#include "misc.h"
#include "reading_coverage.h"
#include "readpair_coverage.h"
#include "tclXkeylist.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "consen.h"
#include "qual.h" 
#include "consistency_display.h"
#include "text_output.h"
#include "ruler_display.h"
#include "template.h"

/*
 * Maintains readpair_coverage with the contig registration scheme.
 * The histogram is drawn as part of a consistency display
 */

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void readpair_coverage_callback(GapIO *io, int contig, void *fdata,
				      reg_data *jdata);


/*
 *---------------------------------------------------------------------------
 * Internal functions
 *---------------------------------------------------------------------------
 */


/*
 */
void plot_readpair_coverage(Tcl_Interp *interp,
			    int *histogram,
			    int seq_len,
			    char *c_win,
			    GapIO *io,
			    int offset,
			    int width,
			    char *colour,
			    int min,
			    int max)
{
    int i, start;
    char cmd[1024];
    
    start = 1;

    for (i = 2; i <= seq_len; i++) {
#ifdef DEBUG
      printf("h[%d]=%d\n", i, histogram[i]);
#endif
      if (histogram[i] != histogram[i-1]) {
	/* horizonal line */
	sprintf(cmd, "%s create line %d %d %d %d -fill %s -width %d", 
		  c_win, start+offset-1, max - histogram[start] + min, 
		  i-1+offset-1, max - histogram[i-1] + min, colour, width);

#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);
	/* vertical line */
	  sprintf(cmd, "%s create line %d %d %d %d -fill %s -width %d", 
		  c_win, i-1+offset-1, max - histogram[i-1] + min, 
		  i+offset-1, max - histogram[i] + min, colour, width);

#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);
	  start = i;
      }
    }
    /* last horizonal line */
    sprintf(cmd, "%s create line %d %d %d %d -fill %s -width %d", 
	    c_win, start+offset-1, max - histogram[start] + min, 
	    i-1+offset-1, max - histogram[i-1] + min, colour, width);
    
#ifdef DEBUG
    printf("%s\n", cmd);
#endif
    Tcl_Eval(interp, cmd);
}


int calc_readpair_coverage(GapIO *io,
			   int contig,
			   int start,
			   int end,
			   int *histogram, /* out */
			   int *min,    /* out */
			   int *max)  /* out */
{
    int i, j;
    int x1, x2, tx1, tx2;
    template_c *t;
    template_c **tarr;
    int status;
    int ntemplates = Ntemplates(io);
    item_t *item;
    gel_cont_t *gc;
    int tmp;

    tx1 = INT_MAX;
    tx2 = INT_MIN;

    if (ntemplates == 0)
	return -1;

    /* Initialise templates */
    if (NULL == (tarr = init_template_checks(io, 1, &contig, 1)))
	return -1;

    check_all_templates(io, tarr); 

    /* for each template */
    for (i = 1; i <= ntemplates; i++){
	t = tarr[i];

	/* if the template is in the contig */
	if (t) {
	    for (item = head(t->gel_cont); item; item=item->next) {
	      gc = (gel_cont_t *)(item->data);

	      if (gc->contig != contig)
		continue;

	      tmp = t->consistency;
	      get_template_positions(io, t, gc->contig);
	      t->consistency |= tmp;
    
	      status = getStatus(t);
	      if (status == FORW_REV_READINGS) {
		x1 = MIN(MIN(t->start, t->end), t->min);
		x2 = MAX(MAX(t->start, t->end), t->max);
#ifdef DEBUG
		printf("FORW_REV_READINGS %d x1 %d x2 %d\n", i, x1, x2);
#endif
		for (j = x1; j <= x2; j++) {
		  if (j >= start && j <= end) {
		    histogram[j-start+1]++;
		    
		    if (histogram[j-start] > *max)
		      *max = histogram[j-start];
		    if (histogram[j-start] < *min)
		      *min = histogram[j-start];
		  }
		}
	      }
	      break;
	    }
	}
    }

    if (tarr) 
      uninit_template_checks(io, tarr);

    /* always set minimum to 0 */
    *min = 0;
    return 0;
}

/*
 * plot the vertical ruler
 */
void plot_readpair_coverage_ruler(Tcl_Interp *interp,
				  obj_reading_coverage *rcov,
				  CanvasPtr *canvas,
				  WorldPtr *world)
{

    display_ruler_v(interp, canvas, rcov->ruler, rcov->t_min, rcov->t_max);
    scaleSingleCanvas(interp, world, canvas, rcov->ruler->window, 'y', "all");  
}

/*
 * readpair_coverage histogram
 */
void display_readpair_coverage(GapIO *io, obj_reading_coverage *rcov) {
    char cmd[1024];
    int i;
    obj_consistency_disp *c;
    int length;
    int win_num;

    c = result_data(io, rcov->cons_id, 0);

    sprintf(cmd, "%s delete all", rcov->c_win);
    Tcl_Eval(c->interp, cmd);

    win_num = get_consistency_win_num(c, rcov->id);
    for (i = 0; i < c->num_contigs; i++) {
        if (c->num_contigs > 1) {
	    length = ABS(io_clength(io, c->contigs[i]));
	} else {
	    length = c->end - c->start + 1;
	}
        plot_readpair_coverage(c->interp, rcov->histogram1[i], length, 
			       rcov->c_win, io, 
			       c->start+c->contig_offset[c->contigs[i]].offset, 
			       rcov->linewidth, rcov->colour1, rcov->t_min, 
			       rcov->t_max);
    }

    plot_readpair_coverage_ruler(c->interp, rcov, c->win_list[win_num]->canvas,
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
 * Removes the readpair_coverage (and unplots etc).
 */
static void readpair_coverage_shutdown(GapIO *io, obj_reading_coverage *rcov) {
    char cmd[1024];
    int i;
    obj_consistency_disp *c;
    int win_num;

    if (NULL != (c = result_data(io, rcov->cons_id, 0))) {
	win_num = get_consistency_win_num(c, rcov->id);
	if (win_num != -1) {
	  delete_consistency_window(c, win_num);
	}
    }
    for (i = 0; i < c->num_contigs; i++) {
      contig_deregister(io, c->contigs[i], readpair_coverage_callback, 
			(void *)rcov);
    }
    sprintf(cmd, "DeleteReadPairCoverage %d %s %s %d\n", *handle_io(io), 
	    rcov->frame, rcov->c_win, rcov->cons_id);
    if (TCL_ERROR == Tcl_Eval(c->interp, cmd)) 
      printf("readpair_coverage_shutdown: %s\n", Tcl_GetStringResult(c->interp));

    if (rcov->histogram1) {
        for (i = 0; i < c->num_contigs; i++) {
	  xfree(rcov->histogram1[i]);
	}
	xfree(rcov->histogram1);
    }

    if (rcov->min)
        xfree(rcov->min);

    if (rcov->max)
        xfree(rcov->max);

    free_ruler_struct(rcov->ruler);

    xfree(rcov);
    if (c->num_wins == 0)
        consistency_shutdown(io, c);
}

/*
 * Updates an obj_reading_coverage structure.
 * This includes recalculating the length, reallocing the memory used,
 * and recalculating the consensus sequence.
 *
 * Returns 0 for success, -1 for failure.
 */
static int update_obj_readpair_coverage(GapIO *io, 
					obj_reading_coverage *rcov, 
					int contig) {
    int i, j, ret, newl;
    obj_consistency_disp *c;
    int start, end;
    int win_num;

#ifdef DEBUG
    printf("UPDATE_OBJ_READPAIR_COVERAGE\n");
#endif
    c = result_data(io, rcov->cons_id, 0);
    
    /* find contig index */
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

    if (rcov->histogram1[i]) {
      if (NULL == (rcov->histogram1[i] = 
		   (int *)xrealloc(rcov->histogram1[i], (newl+1)*sizeof(int))))
	return -1;
    }

    /* re-initialise */
    rcov->max[i] = INT_MIN;
    rcov->min[i] = INT_MAX;

    rcov->t_max = INT_MIN;
    rcov->t_min = INT_MAX;

    for (j = 0; j <= newl; j++) {
      rcov->histogram1[i][j] = 0;
    }

    ret = calc_readpair_coverage(io, c->contigs[i], start, end, 
				 rcov->histogram1[i], &rcov->min[i], 
				 &rcov->max[i]);

    /* find total min and max over all contigs */
    for (i = 0; i < c->num_contigs; i++) {
        if (rcov->t_max < rcov->max[i])
	    rcov->t_max = rcov->max[i];
    }

    rcov->t_min = 0;

    /* need to reset y coords */
    win_num = get_consistency_win_num(c, rcov->id);
    c->win_list[win_num]->world->total->y1 = rcov->t_min;
    c->win_list[win_num]->world->total->y2 = rcov->t_max;
    c->win_list[win_num]->world->visible->y1 = rcov->t_min;
    c->win_list[win_num]->world->visible->y2 = rcov->t_max;
    
    return ret;
}

#if 0
static void readpair_coverage_list(GapIO *io, 
				   obj_consistency_disp *c,
				   obj_reading_coverage *rcov,
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
	  vmessage("%d %d\n", c->start+j, rcov->histogram1[i][j]);
	}
	vmessage("\n");
    }
}
#endif

/* scroll the readpair_coverage ruler */
static void readpair_coverageScrollY(Tcl_Interp *interp,
				     obj_reading_coverage *rcov,
				     CanvasPtr *canvas,
				     d_box *visible,
				     char *scroll_args)
{
    char cmd[1024];
    double wx;
  
    sprintf(cmd, "eval %s yview %s", rcov->ruler->window, scroll_args);
    Tcl_Eval(interp, cmd);

    /* find new top edge of canvas in canvasy coords */
    Tcl_VarEval(interp, rcov->ruler->window, " canvasy 0", NULL);
    canvas->y = atoi(Tcl_GetStringResult(interp));

    /* find new top and bottom edges of canvas in world coords */
    CanvasToWorld(canvas, 0, canvas->y, &wx, &visible->y1);
    CanvasToWorld(canvas, 0, canvas->y + canvas->height, &wx, &visible->y2);
    SetCanvasCoords(interp, visible->x1, visible->y1, 
		    visible->x2, visible->y2, canvas);
}

static void readpair_coverage_join(GapIO *io,
				   obj_consistency_disp *c,
				   obj_reading_coverage *rcov, 
				   int index)
{
    int length;

    length = c->num_contigs - index - 1;

    memmove(&(rcov->histogram1)[index], &(rcov->histogram1)[index+1], 
	    length * sizeof(int*));
    memmove(&(rcov->min)[index], &(rcov->min)[index+1], 
	    length * sizeof(int));
    memmove(&(rcov->max)[index], &(rcov->max)[index+1], 
	    length * sizeof(int));
}

/*
 * The callback from the contig registration scheme.
 */
static void readpair_coverage_callback(GapIO *io, int contig, void *fdata,
				       reg_data *jdata) {
    obj_reading_coverage *rcov = (obj_reading_coverage *)fdata;
    obj_consistency_disp *c;
    int win_num;

    c = result_data(io, rcov->cons_id, 0);

    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Readpair coverage histogram");
	    return;
	}
    case REG_COMPLEMENT:
    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("readpair_coverage COMPLEMENT LENGTH contig %d\n", contig);
#endif
	    if (update_obj_readpair_coverage(io, rcov, contig) == 0) {
		display_readpair_coverage(io, rcov);
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
	    printf("readpair_coverage QUIT DELETE \n");
#endif
	    readpair_coverage_shutdown(io, rcov);
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
		readpair_coverage_list(io, c, rcov, contig);
		end_message();
		break;
#endif
	    case 0: /* Remove */
		readpair_coverage_shutdown(io, rcov);
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
	  jdata->generic.data = (void *)rcov->cons_id;
  
	case TASK_DISPLAY_TICKS:
	    {
		/* int *ticks = (int *)jdata->generic.data; */
		char cmd[100];
		int win_num;
#ifdef DEBUG
		printf("Readpair_Coverage display_ticks\n");
#endif
		/* get readpair_coverage window */
		win_num = get_consistency_win_num(c, rcov->id);

		sprintf(cmd, "%s delete all", rcov->ruler->window);
		Tcl_Eval(c->interp, cmd);
		plot_readpair_coverage_ruler(c->interp, rcov, 
					     c->win_list[win_num]->canvas,
					     c->win_list[win_num]->world);
		break;
	    }
	case TASK_CANVAS_SCROLLY: 
	    {
		char *scroll = (char *)jdata->generic.data;

		win_num = get_consistency_win_num(c, rcov->id);

		/* scroll the ruler */
		readpair_coverageScrollY(c->interp, rcov, 
					c->win_list[win_num]->canvas,
					c->win_list[win_num]->world->visible,
					scroll);
		
		/* scroll the readpair_coverage plot */
		consistencyScrollY(c->interp, rcov->c_win, 
				   c->win_list[win_num]->scroll,
				   c->win_list[win_num]->world->visible, 
				   c->win_list[win_num]->canvas, scroll);
		
		break;
	    }
	case TASK_CANVAS_REDRAW:
#ifdef DEBUG
	    printf("readpair_coverage TASK_CANVAS_REDRAW\n");
#endif
	    display_readpair_coverage(io, rcov);
	    break;
	case TASK_CANVAS_RESIZE:
	  {
	      char cmd[1024];
	      int win_num;
	      /* resize readpair_coverage window */
#ifdef DEBUG
	    printf("Readpair_Coverage TASK_CANVAS_RESIZE\n");
#endif
	    win_num = get_consistency_win_num(c, rcov->id);
	    sprintf(cmd, "%s delete all", rcov->ruler->window);
	    Tcl_Eval(c->interp, cmd);
	    plot_readpair_coverage_ruler(c->interp, rcov, 
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
		printf("READPAIR_COVERAGE TASK_CANVAS_WORLD\n");
#endif
		win_num = get_consistency_win_num(c, rcov->id);
		CanvasToWorld(c->win_list[win_num]->canvas, cx, 0, &wx, &wy);
		tw->basex = wx - c->contig_offset[tw->cnum].offset;
		break;
	    }
	case TASK_CONS_WORLD:
	    {
		d_box *cw = (d_box *)jdata->generic.data;

		win_num = get_consistency_win_num(c, rcov->id);
		cw->x1 = c->win_list[win_num]->world->total->x1;
		cw->y1 = c->win_list[win_num]->world->total->y1;
		cw->x2 = c->win_list[win_num]->world->total->x2;
		cw->y2 = c->win_list[win_num]->world->total->y2;
		
		break;
	    }
	case TASK_CONS_JOIN:
	    {
#ifdef DEBUG
	      printf("readpair_coverage JOIN_TO contig %d join %d\n", 
		     contig, jdata->join.contig); 
#endif
	      readpair_coverage_join(io, c, rcov, (int)jdata->generic.data);
	      break;
	    }
	case TASK_CONS_CURSOR_DELETE:
	    {
#ifdef DEBUG
	      printf("READPAIR_COVERAGE cursor_delete\n");
#endif
	      Tcl_VarEval(c->interp, rcov->c_win, " delete cursor_y", NULL);
	      break;
	    }
	case TASK_CANVAS_CURSOR_Y:
	    {    
		char *label;
		int *cy = (int *)jdata->generic.data;
		double wx, wy;
		int win_num;

		win_num = get_consistency_win_num(c, rcov->id);
		CanvasToWorld(c->win_list[win_num]->canvas, 0, *cy, &wx, &wy);
		wy = rcov->t_max - wy + rcov->t_min;

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
int readpair_coverage_reg(GapIO *io, 
			 Tcl_Interp *interp,
			 char *frame,
			 char *rcov_win,
			 int cons_id,
			 ruler_s *ruler) 
{
    obj_reading_coverage *rcov;
    obj_consistency_disp *c;
    int id;
    int i, j;
    int start, end, length;
    char *val;

    c = result_data(io, cons_id, 0);

    /* first check that there are enough available windows */
    if (c->num_wins >= MAX_NUM_WINS)
      return -1;

    if (NULL == (rcov = (obj_reading_coverage*)xmalloc(sizeof(obj_reading_coverage))))
	return -1;

    if (NULL == (rcov->histogram1 = (int **)xmalloc(c->num_contigs*sizeof(int*))))
	return -1;
    if (NULL == (rcov->min = (int *)xmalloc(c->num_contigs*sizeof(int))))
	return -1;
    if (NULL == (rcov->max = (int *)xmalloc(c->num_contigs*sizeof(int))))
	return -1;

    id = register_id();
    
    rcov->id = id;
    rcov->cons_id = cons_id;
    strcpy(rcov->c_win, rcov_win);
    strcpy(rcov->frame, frame);
    rcov->linewidth = get_default_int(interp, gap_defs,
				      "READPAIR_COVERAGE.LINEWIDTH");
    val = get_default_string(interp, gap_defs, "READPAIR_COVERAGE.COLOUR"); 
    strcpy(rcov->colour1, val);

    rcov->t_max = INT_MIN;
    rcov->t_min = INT_MAX;
    rcov->ruler = ruler;
    
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

	if (NULL==(rcov->histogram1[i]=(int *)xmalloc((length+1)*sizeof(int)))) {
	    return -1;
	}
	/* initialise histogram and min and max values */
	for (j = 0; j <= length; j++) {
	    rcov->histogram1[i][j] = 0;
	}
	rcov->max[i] = INT_MIN;
	rcov->min[i] = INT_MAX;

        calc_readpair_coverage(io, c->contigs[i], start, end, 
			      rcov->histogram1[i], 
			      &rcov->min[i], &rcov->max[i]);
	
	/* find total min and max over all contigs */
        if (rcov->t_max < rcov->max[i]) {
	    rcov->t_max = rcov->max[i];
	}
	rcov->t_min = 0;
    }

    if (rcov->t_max == INT_MIN) {
      vmessage("No read pairs within contigs have been found\n");
      readpair_coverage_shutdown(io, rcov);
      return -2;
    }

    add_consistency_window(io, c, rcov_win, 'b', id, c->orig_total->x1,
			   rcov->t_min, c->orig_total->x2, rcov->t_max);

    display_readpair_coverage(io, rcov);

    for (i = 0; i < c->num_contigs; i++) {
        contig_register(io, c->contigs[i], readpair_coverage_callback, 
			(void *)rcov, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_CURSOR_NOTIFY | REG_GENERIC, 
			REG_TYPE_READPAIR_COVERAGE);
    }
#ifdef DEBUG
    printf("READPAIR_COVERAGE ID %d\n", id);
#endif
    return id;
}
