#include <tk.h>
#include <string.h>

#include "IO.h"
#include "io-reg.h"
#include "canvas_box.h"
#include "gap_canvas_box.h"
#include "newgap_cmds.h"
#include "template_display.h"
#include "gap_globals.h"
#include "misc.h"
#include "qual.h"
#include "quality_plot.h"
#include "tclXkeylist.h"
#include "fort.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "edUtils.h"

/*
 * Maintains quality plots with the contig registration scheme.
 * template_quality refers to the quality plot within the template display
 * quality is the stand alone, single quality plot
 */

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void template_quality_callback(GapIO *io, int contig, void *fdata,
				      reg_data *jdata);

static void quality_callback(GapIO *io, int contig, void *fdata,
				      reg_data *jdata);


/*
 *---------------------------------------------------------------------------
 * Internal functions
 *---------------------------------------------------------------------------
 */

void
glevel(char current, 
       int Y0, 
       int YP1, 
       int YP2, 
       int YM1, 
       int YM2, 
       int *y1, 
       int *y2)
{
    switch(current) {
    case R_GOOD_GOOD_EQ:
        *y1 = Y0;
        *y2 = Y0;
	break;
    case R_GOOD_BAD:
    case R_GOOD_NONE:
        *y1 = Y0;
        *y2 = YM1;
	break;
    case R_NONE_GOOD:
    case R_BAD_GOOD:
        *y1 = Y0;
        *y2 = YP1;
	break;
    case R_BAD_BAD:
    case R_BAD_NONE:
    case R_NONE_BAD:
    case R_NONE_NONE:
        *y1 = YP1;
        *y2 = YM1;
	break;
    case R_GOOD_GOOD_NE:
        *y1 = YP2;
        *y2 = YM2;
	break;
    default:
        verror(ERR_FATAL, "quality_plot","incorrect value to glevel()");
    }
} /* end glevel */


char *
quality_colour(Tcl_Interp *interp,
	       int y1, int y2, 
	       int Y0, int YM1, int YM2, int YP1, int YP2)
{
    char fieldName[100];
        
    if ((y1 == Y0) && (y2 == Y0)) {
	strcpy(fieldName, "TEMPLATE.QUALITY.BOTH_COLOUR");
    } else if ((y1 == Y0) && (y2 == YM1)) {
	strcpy(fieldName, "TEMPLATE.QUALITY.PLUS_COLOUR");
    } else if ((y1 == Y0) && (y2 == YP1)) {
	strcpy(fieldName, "TEMPLATE.QUALITY.MINUS_COLOUR");
    } else if ((y1 == YP1) && (y2 == YM1)) {
	strcpy(fieldName, "TEMPLATE.QUALITY.BAD_COLOUR");
    } else if ((y1 == YP2) && (y2 == YM2)) {
	strcpy(fieldName, "TEMPLATE.QUALITY.DISAGREE_COLOUR");
    }
 
    return get_default_astring(interp, gap_defs, fieldName);
}

/*
 * each code has its own level and we draw rectangles for segments of
 * the same code
 * The coordinate system is xmin to xmax, and ymin to ymax which here
 * are set to 1, sequence_length; -2, 2 
 */
void
plot_quality(Tcl_Interp *interp,
	     char *seq, 
	     int seq_len,
             char *q_win,
	     GapIO *io,
	     int offset)
{
    int Y0  = 0;
    int YP1 = 1;
    int YP2 = 2;
    int YM1 = -1;
    int YM2 = -2;
    int i;
    int x1, x2;
    int y1, y2;
    char current;
    char cmd[1024];
    char *colour;

    i = 0;
    current = seq[i];
    x1 = i;
    while (i < seq_len) {
        if (current != seq[i]) {
            /* printf("not equal cur %c seq(%d) %c \n", current, i, seq[i]); */

            glevel(current, Y0, YP1, YP2, YM1, YM2, &y1, &y2);
	    x2 = i;

	    colour = quality_colour(interp, y1, y2, Y0, YM1, YM2, YP1, YP2);
	    /* HACK - do something about calc of quality y position */
	    sprintf(cmd, "%s create rectangle %d %d %d %d -fill %s "
		    "-outline %s -tag {quality S}",
		    q_win, x1+offset, 24 + (y1 * 6), x2+offset, 24 + (y2 * 6),
		    colour, colour);
	    xfree(colour);
	    Tcl_Eval(interp, cmd);
            current = seq[i];
	    x1 = i;
        } 
        i++;
    }
    /* we have reached the end so finish off last rectangle */
    glevel(current, Y0, YP1, YP2, YM1, YM2, &y1, &y2);
    x2 = i;

    colour = quality_colour(interp, y1, y2, Y0, YM1, YM2, YP1, YP2);
    sprintf(cmd, "%s create rectangle %d %d %d %d -fill %s "
	    "-outline %s -tag {quality S}",
	    q_win, x1+offset, 24 + (y1 * 6), x2+offset, 24 + (y2 * 6), colour, colour);
    xfree(colour);
    Tcl_Eval(interp, cmd);
}

/*
 * Lists the quality to the output window
 */
static void quality_list(GapIO *io, char *qual, int contig, int length) {
    f_int l;
    f_int start = 1;
    f_int end;
    f_int linlen = 60;
    f_int kbout = 0;

    vfuncheader("quality listing");
    l = length;
    end = length;
    
    vmessage("Contig %s (#%d)\n", 
	     get_contig_name(io, contig),
	     io_clnbr(io, contig));
    fmtdb_(qual, &l, &start, &end, &linlen, &kbout, l);
    vmessage("\n");

}


/*
 * Provides a summary of quality values
 */
static void quality_info(GapIO *io, int contig, int len, char *qual, 
			 int header) {
    int i, both_eq = 0, plus = 0, minus = 0, bad_both = 0, both_ne = 0;

    if (header)
	vfuncheader("quality summary");

    /* was: dbscsm_(q->qual, &len, len); */

    both_eq = 0, plus = 0, minus = 0, bad_both = 0, both_ne = 0;
    for (i = 0; i <len; i++) {
	switch(qual[i]) {
	case R_GOOD_GOOD_EQ:
	    both_eq++;
	    break;
	    
	case R_GOOD_BAD:
	case R_GOOD_NONE:
	    plus++;
	    break;
	    
	case R_BAD_GOOD:
	case R_NONE_GOOD:
	    minus++;
	    break;
		
	case R_BAD_BAD:
	case R_BAD_NONE:
	case R_NONE_BAD:
	case R_NONE_NONE:
	    bad_both++;
	    break;
	    
	default: /* R_NONE_NONE */
	    both_ne++;
	}
	}
	
    vmessage("Contig %s (#%d)\n", 
	     get_contig_name(io, contig), io_clnbr(io, contig));
    vmessage("%6.2f OK on both strands and they agree(a)\n",
	     ((float)(100 * both_eq)) / len);
    vmessage("%6.2f OK on plus strand only(b,d)\n",
		 ((float)(100 * plus)) / len);
    vmessage("%6.2f OK on minus strand only(c,e)\n",
	     ((float)(100 * minus)) / len);
    vmessage("%6.2f Bad on both strands(f,g,h,j)\n",
	     ((float)(100 * bad_both)) / len);
    vmessage("%6.2f OK on both strands but they disagree(i)\n\n",
	     ((float)(100 * both_ne)) / len);
}


/*
 * template quality
 * Updates an obj_t_qual structure.
 * This includes recalculating the length, reallocing the memory used,
 * and recalculating the consensus sequence.
 *
 * Returns 0 for success, -1 for failure.
 */
static int update_obj_t_qual(GapIO *io, obj_t_qual *q, int contig) {
    int i, ret, newl;

    for (i = 0; i < q->num_contigs; i++) {
	if (q->quality[i].contig == contig) {
	    break;
	}
    }
    newl = ABS(io_clength(io, q->quality[i].contig));

    if (q->quality[i].length < newl) {
	if (q->quality[i].qual) {
	    if (NULL == (q->quality[i].qual = 
			 (char *)xrealloc(q->quality[i].qual, newl)))
		return -1;
	}
    }

    q->quality[i].length = newl;
    ret = calc_quality(q->quality[i].contig, 1, q->quality[i].length, 
		       q->quality[i].qual, 
		       q->cons_cutoff, q->qual_cutoff, database_info, 
		       (void *)io);

    return ret;
}

void
template_display_quality(GapIO *io, 
			 obj_t_qual *q, 
			 c_offset *contig_offset) 
{
    int i;
    char cmd[1024];
    obj_template_disp *t;
    
    t = result_data(io, q->template_id, 0);

    sprintf(cmd, "%s delete quality", q->window);
    Tcl_Eval(q->interp, cmd);

    for (i = 0; i < q->num_contigs; i++) {
	plot_quality(q->interp, q->quality[i].qual, q->quality[i].length, 
		     q->window, io, 
		     contig_offset[q->quality[i].contig].offset);
    }

    scaleSingleCanvas(q->interp, t->world, t->canvas, q->window, 'x', "all");
    template_update_cursors(io, t, 0);
}

/*
 * template quality plot
 * Removes the quality calculation (and unplots etc).
 */
static void template_quality_shutdown(GapIO *io, obj_t_qual *q) {
    int i;
    char cmd[1024];
    obj_template_disp *t;
  
    /*
     * if the template display has already been deleted, then the win_list
     * will have been freed anyway
     */
    if (NULL != (t = result_data(io, q->template_id, 0)))
	deleteWindow(t->win_list, &t->num_wins, q->window);

    for (i = 0; i < q->num_contigs; i++) {
	contig_deregister(io, q->quality[i].contig, template_quality_callback, 
			  (void *)q);
    }
    sprintf(cmd, "DeleteTemplateQualityPlot %s %s\n", q->frame, q->window);
    Tcl_Eval(q->interp, cmd);

    for (i = 0; i < q->num_contigs; i++) {
	if (q->quality[i].qual)
	    xfree(q->quality[i].qual);
    }
    xfree(q->quality);
    xfree(q);
}

static void
template_quality_renumber(GapIO *io,
			  obj_t_qual *q, 
			  int old_contig,
			  int new_contig)
{
    int i;
    for (i = 0; i < q->num_contigs; i++) {
	if (ABS(q->quality[i].contig) == old_contig) {
	    q->quality[i].contig = q->quality[i].contig > 0 ? new_contig : 
		-new_contig;
	    break;
	}
    }
}

static void
template_quality_join(GapIO *io,
		      obj_t_qual *q, 
		      int old_contig,
		      int new_contig)
{
    int i, length;

    for (i = 0; i < q->num_contigs; i++) {
	if (ABS(q->quality[i].contig) == old_contig) {
	    length = q->num_contigs - i - 1;
	    memmove(&(q->quality)[i], &(q->quality)[i+1], 
		    length * sizeof(c_qual));
	    q->num_contigs--;
	    break;
	}
    }
    
}


/*
 * The callback from the contig registration scheme.
 */
static void template_quality_callback(GapIO *io, int contig, void *fdata,
				      reg_data *jdata) {
    obj_t_qual *q = (obj_t_qual *)fdata;
    obj_template_disp *t;
    
    t = result_data(io, q->template_id, 0);
    
    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
#ifdef DEBUG
	    printf("quality QUERY_NAME\n");
#endif
	    sprintf(jdata->name.line, "Calculate quality");
	    return;
	}
    case REG_JOIN_TO: 

	{
#ifdef DEBUG
	    printf("quality JOIN TO contig %d join %d num_contigs %d\n", 
		   contig, jdata->join.contig, q->num_contigs);
#endif
	    /* remove window */
	    if (q->num_contigs == 1) {
		template_quality_shutdown(io, q);
		return;
	    } else {
		/* 
		 * HACK - not tested since always have 1 contig in template
		 * quality plot
		 */
		template_quality_join(io, q, contig, jdata->join.contig);
	    }
	}
	/* NOTE - flow through to update_obj_t_qual() */

    case REG_COMPLEMENT:
    case REG_LENGTH:
        {
	    reg_generic gen;
#ifdef DEBUG
	    printf("quality COMPLEMENT LENGTH contig %d\n", contig);
#endif
	    /* 
	     * need to update ALL the quality displays because the 
	     * contig_offset will have changed
	     */
	    if (update_obj_t_qual(io, q, contig) == 0) {
		gen.job = REG_GENERIC;
		gen.task = TASK_CANVAS_REDRAW;
		type_notify(io, REG_TYPE_QUALITY, (reg_data *)&gen, 1);
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
	    printf("quality QUIT DELETE \n");
#endif
	    template_quality_shutdown(io, q);
	    return;
	}

    case REG_GET_OPS:
	{
#ifdef DEBUG
	    printf("quality GET_OPS\n");
#endif
	    jdata->get_ops.ops =
		"Information\0"
		"List\0"
		    "SEPARATOR\0"
		"Remove\0";
	    return;
	}
    case REG_INVOKE_OP:
	{
#ifdef DEBUG
	    printf("quality INVOKE_OP %d\n", jdata->invoke_op.op);
#endif
	    switch (jdata->invoke_op.op) {
	    case 0: { /* Information */
		int header = 1, i;
		for (i = 0; i < q->num_contigs; i++) {
		    if (i) 
			header = 0;
		    start_message();
		    quality_info(io, q->quality[i].contig, 
				  q->quality[i].length,
				 q->quality[i].qual, header);
		    end_message(q->window);
		}
		break;
		}
	    case 1: { /* List Quality */
		int i;
		for (i = 0; i < q->num_contigs; i++) {
		    quality_list(io, q->quality[i].qual, q->quality[i].contig,
				 q->quality[i].length);
		}
		break;
	    }
	    case 2: /* Remove */
		template_quality_shutdown(io, q);
		break;
	    }
	    return;
	}
    case REG_PARAMS:
	{
	    puts("REG_PARAMS");
	    return;
	}
    case REG_NUMBER_CHANGE:
	{
#ifdef DEBUG
	    printf("quality REG_NUMBER_CHANGE contig %d number %d\n",
		   contig, jdata->number.number);
#endif
	    template_quality_renumber(io, q, contig, jdata->number.number);
	    template_display_quality(io, q, t->contig_offset);
	    return;
	}
    case REG_GENERIC:
	switch (jdata->generic.task) {
	    case TASK_CANVAS_REDRAW:
#ifdef DEBUG
	    printf("quality TASK_CANVAS_REDRAW\n");
#endif
	    template_display_quality(io, q, t->contig_offset);

	}
    }
}

/*
 *---------------------------------------------------------------------------
 * External functions
 *---------------------------------------------------------------------------
 */

/*
 * Registers and initialises a quality buffer for a particular contig.
 */
int template_quality_reg(GapIO *io, 
			 Tcl_Interp *interp,
			 int *contig_array,
			 int num_contigs,
			 float cons_cutoff, 
			 int qual_cutoff, 
			 char *frame,
			 char *win_quality,
			 int template_id) 
{
    obj_template_disp *t;
    obj_t_qual *q;
    int id;
    int i;
    reg_generic gen;
    win winfo;

    t = result_data(io, template_id, 0);

    /* first check that there are enough available windows */
    if (t->num_wins >= MAX_NUM_WINS)
	return -1;

    if (NULL == (q = (obj_t_qual *)xmalloc(sizeof(obj_t_qual))))
       return -1;

    if (NULL == (q->quality = (c_qual *)xmalloc(num_contigs * 
						sizeof(c_qual)))){
	xfree(q);
	return -1;
    }

    id = register_id();

     /* add to window list */
    winfo.window = win_quality;
    winfo.scroll = 'x';
    winfo.id = id;

    gen.job = REG_GENERIC;
    gen.task = TASK_WINDOW_ADD;
    gen.data = (void *)&winfo;
    result_notify(io, template_id, (reg_data *)&gen, 0); 

    q->cons_cutoff = cons_cutoff;
    q->qual_cutoff = qual_cutoff;
    q->template_id = template_id;
    q->num_contigs = num_contigs;
    q->interp = interp;

    strcpy(q->window, win_quality);
    strcpy(q->frame, frame);

    for (i = 0; i < num_contigs; i++) {
	int ret;

	q->quality[i].contig = contig_array[i];
	q->quality[i].length = io_clength(io, q->quality[i].contig);

	if (NULL==(q->quality[i].qual=(char *)xmalloc(q->quality[i].length))) {
	    return -1;
	}
	ret = calc_quality(q->quality[i].contig, 1, q->quality[i].length, 
			   q->quality[i].qual, q->cons_cutoff, q->qual_cutoff,
			   database_info, (void *)io);

	if (ret == -1) 
	    return -1;
	quality_info(io, q->quality[i].contig, q->quality[i].length,
		     q->quality[i].qual, 0);
    }

    template_display_quality(io, q, t->contig_offset);

    for (i = 0; i < q->num_contigs; i++) {
	contig_register(io, q->quality[i].contig, template_quality_callback, 
			(void *)q, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_GENERIC, REG_TYPE_QUALITY);
    }

    return id;
}

/*
 * stand alone quality display
 */
void
display_quality(GapIO *io, obj_qual *q) {
    char cmd[1024];

    sprintf(cmd, "%s delete all", q->window);
    Tcl_Eval(q->interp, cmd);

    plot_quality(q->interp, q->quality.qual, q->quality.length, q->window, 
		 io, q->quality.start);

    q->world->total->x1 = q->quality.start;
    q->world->total->x2 = q->quality.end;
    q->world->total->y1 = 2;
    q->world->total->y2 = -2;

    memcpy(q->world->visible, q->world->total, sizeof(d_box));

    SetCanvasCoords(q->interp, q->world->visible->x1, q->world->visible->y1, 
		    q->world->visible->x2, q->world->visible->y2, q->canvas);

    draw_single_ruler(q->interp, q->ruler, q->canvas, q->ruler->start,
		      q->ruler->end, 1);

    scaleCanvas(q->interp, q->win_list, q->num_wins, "all", 
		q->world->visible, q->canvas);
    scrollRegion(q->interp, q->win_list, q->num_wins, q->world->total, 
		 q->canvas);

    /* remove all current zooming info */
    freeZoom(&q->zoom);

    /* add first zoom */
    pushZoom(&q->zoom, q->world->visible);
}

static int 
update_obj_qual(GapIO *io, obj_qual *q, int contig) {
    int ret, newl;

    newl = ABS(io_clength(io, q->quality.contig));

    if (q->quality.length < newl) {
	if (q->quality.qual) {
	    if (NULL == (q->quality.qual = 
			 (char *)xrealloc(q->quality.qual, newl)))
		return -1;
	}
    }

    /* if need to replot, lose original lreg to rreg info */
    q->quality.start = 1;
    q->quality.end = newl;
    q->quality.length = newl;

    /* reset ruler from and to */
    q->ruler->start = 1;
    q->ruler->end = newl;
 
   ret = calc_quality(q->quality.contig, 1, q->quality.length, 
		       q->quality.qual, 
		       q->cons_cutoff, q->qual_cutoff, database_info, 
		       (void *)io);

    return ret;
}

/*
 * stand alone quality plot
 * Removes the quality calculation (and unplots etc).
 */
static void quality_shutdown(GapIO *io, obj_qual *q) {
    char cmd[1024];

    contig_deregister(io, q->quality.contig, quality_callback, (void *)q);

    delete_contig_cursor(io, q->quality.contig, q->cursor->id, 0);

    sprintf(cmd, "DeleteQualDisplay %s %s\n", q->frame, q->window);
    Tcl_Eval(q->interp, cmd);

    if (q->quality.qual)
	xfree(q->quality.qual);

    if (q->canvas)
	xfree(q->canvas);
    
    if (q->world->visible)
	xfree(q->world->visible);
    
    if (q->world->total)
	xfree (q->world->total);
    
    if (q->world)
	xfree(q->world);

    if (q->xhair.colour) free(q->xhair.colour);

    free_ruler_struct(q->ruler);

    freeZoom(&q->zoom);
    free_win_list(q->win_list, q->num_wins);
    
    xfree(q);
}

static void
quality_renumber(GapIO *io,
		 obj_qual *q, 
		 int old_contig,
		 int new_contig)
{
    if (q->quality.contig == old_contig) {
	q->quality.contig = q->quality.contig > 0 ? new_contig : -new_contig;
    }
}

/*
 * The callback from the contig registration scheme.
 */
static void quality_callback(GapIO *io, int contig, void *fdata,
			     reg_data *jdata) {
    obj_qual *q = (obj_qual *)fdata;

    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Calculate quality");
	    return;
	}
    case REG_JOIN_TO:
#ifdef DEBUG
	printf("quality JOIN TO contig %d join %d\n", contig, jdata->join.contig);
#endif
	q->quality.contig = jdata->join.contig;
	/* NOTE - flow through to update_obj_qual() */

    case REG_COMPLEMENT:
    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("quality COMPLEMENT LENGTH contig %d\n", contig);
#endif
	    if (update_obj_qual(io, q, contig) == 0) {
		display_quality(io, q);
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
	    printf("quality QUIT DELETE \n");
#endif
	    quality_shutdown(io, q);
	    return;
	}

    case REG_GET_OPS:
	{
	    jdata->get_ops.ops =
		"Information\0"
		"List\0"
		    "SEPARATOR\0"
		"Remove\0";
	    return;
	}
    case REG_INVOKE_OP:
	{
	    switch (jdata->invoke_op.op) {
	    case 0: /* Information */
		start_message();
		quality_info(io, q->quality.contig, q->quality.length, 
			     q->quality.qual, 1);
		end_message(q->window);
		break;
	    case 1: /* List Quality */
		quality_list(io, q->quality.qual, q->quality.contig,
				 q->quality.length);
		break;
	    case 2: /* Remove */
		quality_shutdown(io, q);
		break;
	    }
	    return;
	}
    case REG_PARAMS:
	{
	    puts("REG_PARAMS");
	    return;
	}
    case REG_CURSOR_NOTIFY:
	{
	    canvas_cursor_refresh(q->interp, io, q->quality.contig,
				  jdata->cursor_notify.cursor, q->cursor,
				  q->canvas, q->win_list, q->num_wins,
				  q->id, 0, &q->cursor_visible, q->world, 1);
	    return;
	}
    case REG_NUMBER_CHANGE:
	{
#ifdef DEBUG
	    printf("quality REG_NUMBER_CHANGE contig %d number %d\n",
		   contig, jdata->number.number);
#endif
	    quality_renumber(io, q, contig, jdata->number.number);
	    display_quality(io, q);
	    return;
	}
    case REG_GENERIC:
	switch (jdata->generic.task) {
	case TASK_CANVAS_SCROLLX: 
	    {
		char *scroll = (char *)jdata->generic.data;
		
		canvasScrollX(q->interp, q->window, q->win_list, q->num_wins,
			      q->world->visible, q->canvas, scroll);
		break;
	    }
	case TASK_CANVAS_ZOOMBACK:
	    {
		canvasZoomback(q->interp, q->canvas, q->window, q->world,
			   q->win_list, q->num_wins, &q->zoom);
		/* redraw ruler ticks by redrawing ruler */
		draw_single_ruler(q->interp, q->ruler, q->canvas, 
				  q->ruler->start, q->ruler->end, 1); 
		scaleSingleCanvas(q->interp, q->world, q->canvas, 
				  q->ruler->window, 'x', "all");
	     
		canvas_cursor_refresh(q->interp, io, q->quality.contig,
				      q->cursor, q->cursor,
				      q->canvas, q->win_list, q->num_wins,
				      q->id, 0, &q->cursor_visible, q->world, 0);
		
		break;
	    }
	case TASK_CANVAS_ZOOM:
	    {
		s_zoom *szoom = (s_zoom *)jdata->generic.data;
		canvasZoom(q->interp, q->canvas, q->window, q->world,
			   q->win_list, q->num_wins, &q->zoom, szoom->zoom,
			   szoom->scroll);
		/* redraw ruler ticks by redrawing ruler */
		draw_single_ruler(q->interp, q->ruler, q->canvas, 
				  q->ruler->start, q->ruler->end, 1); 
		scaleSingleCanvas(q->interp, q->world, q->canvas, 
				  q->ruler->window, 'x', "all");
		
		canvas_cursor_refresh(q->interp, io, q->quality.contig,
				      q->cursor, q->cursor,
				      q->canvas, q->win_list, q->num_wins,
				      q->id, 0, &q->cursor_visible, q->world, 0);
		
		break;
	    }
	case TASK_CANVAS_CURSOR_X:
	    {    
		char *label;
		int *cx = (int *)jdata->generic.data;
		double wx, wy;

		CanvasToWorld(q->canvas, *cx, 0, &wx, &wy);

		label = get_default_string(q->interp, gap_defs,
					   "QUALITY.CURSOR");

		canvasCursorX(q->interp, q->canvas, q->frame, label, 
			      q->xhair.colour, q->xhair.width, *cx, wx, 
			      q->win_list, q->num_wins);

		break;
	    }
	case TASK_CANVAS_RESIZE:
	    {
		char scroll_args[20];
		/* resize template display window */
		resizeCanvas(q->interp, q->window, q->win_list, q->num_wins,
			     q->world->visible, q->world->total, q->canvas);
		sprintf(scroll_args, "scroll 0 units");
		canvasScrollX(q->interp, q->window, q->win_list, q->num_wins,
			      q->world->visible, q->canvas, scroll_args);
	    break;
	    }
	case TASK_CANVAS_WORLD:
	    {
		task_world_t *tw = (task_world_t *)jdata->generic.data;
		int cx = tw->canvasx;
		double wx, wy;

		CanvasToWorld(q->canvas, cx, 0, &wx, &wy);
		tw->basex = wx;
		break;
	    }
	}
    }
}


/*
 * Registers and initialises a quality buffer for a particular contig.
 */
int quality_reg(GapIO *io, 
		Tcl_Interp *interp,
		int contig, 
		int lreg,
		int rreg,
		float cons_cutoff, 
		int qual_cutoff, 
		char *frame,
		char *win_quality,
		ruler_s *ruler,
		cursor_s xhair) 
{
    obj_qual *q;
    int id;
    int ret;

    if (NULL == (q = (obj_qual *)xmalloc(sizeof(obj_qual))))
	return -1;

    id = register_id();
    q->id = id;

    q->cons_cutoff = cons_cutoff;
    q->qual_cutoff = qual_cutoff;
    q->quality.contig = contig;
    q->quality.start = lreg;
    q->quality.end = rreg;
    q->quality.length = rreg - lreg + 1;
    q->ruler = ruler;
    q->xhair = xhair;
    q->cursor = create_contig_cursor(io, contig, 0, id);
    q->cursor_visible = 0;
    q->interp = interp;
    
    strcpy(q->window, win_quality);
    strcpy(q->frame, frame);
			   
    /* create list of windows in the quality display */
    if (NULL == (q->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;
    q->num_wins = 0;
    addWindow(q->win_list, &q->num_wins, q->window, 'x', id);
    addWindow(q->win_list, &q->num_wins, q->ruler->window, 'x', id);
    
    if (NULL == (q->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (q->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (q->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (q->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    initCanvas(q->interp, q->canvas, q->window);
    createZoom(&q->zoom);

    if (NULL == (q->quality.qual = (char *)xmalloc(q->quality.length+1))) {
	return -1;
    }
    ret = calc_quality(q->quality.contig, lreg, rreg,
		       q->quality.qual, q->cons_cutoff, q->qual_cutoff, 
		       database_info, (void *)io);

    if (ret == -1) 
	return -1;
    display_quality(io, q);

    contig_register(io, q->quality.contig, quality_callback, (void *)q, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
		    REG_CURSOR_NOTIFY | REG_NUMBER_CHANGE | REG_GENERIC, 
		    REG_TYPE_QUALITY);
    canvas_cursor_refresh(q->interp, io, q->quality.contig,
			  q->cursor, q->cursor,
			  q->canvas, q->win_list, q->num_wins, q->id,
			  0, &q->cursor_visible, q->world, 1);
    
    quality_info(io, q->quality.contig, q->quality.length, q->quality.qual, 0);

    return id;
}
