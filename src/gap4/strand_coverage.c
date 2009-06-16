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
#include "tclXkeylist.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "consen.h"
#include "qual.h" 
#include "consistency_display.h"
#include "text_output.h"
#include "ruler_display.h"
#include "template.h"
#include "strand_coverage.h"

#define FORWARD 0
#define REVERSE 1

/*
 * Maintains strand_coverage with the contig registration scheme.
 * The histogram is drawn as part of a consistency display
 */

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void strand_coverage_callback(GapIO *io, int contig, void *fdata,
				      reg_data *jdata);


/*
 *---------------------------------------------------------------------------
 * Internal functions
 *---------------------------------------------------------------------------
 */

void plot_strand_coverage(Tcl_Interp *interp,
			  int *histogram,
			  int seq_len,
			  char *c_win,
			  GapIO *io,
			  int offset,
			  int width,
			  char *colour,
			  int y_offset)
{
    int i, start;
    char cmd[1024];
    
    start = 0;

    for (i = 1; i <= seq_len; i++) {
#ifdef DEBUG
      printf("h[%d]=%d\n", i, histogram[i]);
#endif
      /* initialise start */
      if (histogram[i] && !start)
	  start = i;

      if (!histogram[i] && histogram[i-1] && start) {
	/* awful hack to deal with this special case */
	if (start == 1)
	    start = 0;

	sprintf(cmd, "%s create line %d %d %d %d -fill %s -width %d "
		"-capstyle round ", 
		c_win, start+offset, y_offset, 
		i-1+offset, y_offset, colour, width);
	
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);
	start = 0;
      }
    }
    /* draw last line */
    if (start) {
	sprintf(cmd, "%s create line %d %d %d %d -fill %s -width %d "
		"-capstyle round ", 
		c_win, start+offset-1, y_offset, 
		i-1+offset-1, y_offset, colour, width);
	
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);
    }
}

/*
 * plot a line where there isn't information ie plotting the problem areas
 */
void plot_strand_problems(Tcl_Interp *interp,
			  int *histogram,
			  int seq_len,
			  char *c_win,
			  GapIO *io,
			  int offset,
			  int width,
			  char *colour,
			  int y_offset)
{
    int i, start;
    char cmd[1024];
    
    start = 0;

    for (i = 1; i <= seq_len; i++) {
#ifdef DEBUG
      printf("h[%d]=%d\n", i, histogram[i]);
#endif
      /* initialise start */
      if (!histogram[i] && !start)
	  start = i;

      if (histogram[i] && !histogram[i-1] && start) {
	/* awful hack to deal with this special case */
	if (start == 1)
	    start = 0;
	sprintf(cmd, "%s create line %d %d %d %d -fill %s -width %d " 
		"-capstyle round ", 
		c_win, start+offset, y_offset, 
		i-1+offset, y_offset, colour, width);
	
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);
	start = 0;
      }
    }
    /* draw last line */
    if (start) {
	sprintf(cmd, "%s create line %d %d %d %d -fill %s -width %d " 
		"-capstyle round ", 
		c_win, start+offset-1, y_offset, 
		i-1+offset-1, y_offset, colour, width);
	
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);
    }
}

int calc_strand_coverage(GapIO *io,
			 int contig,
			 int start,
			 int end,
			 int *forward,
			 int *reverse)
{
    int j;
    char *qual;
    int ret;

    if (NULL==(qual=(char *)xmalloc(end-start+1))) {
      return -1;
    }

    ret = calc_quality(contig, start, end, qual, 
		       consensus_cutoff, quality_cutoff, database_info, 
		       (void *)io);
    
    for (j = start; j <= end; j++) {
      int tmp = j-start;

      switch (qual[tmp]) {
      case R_GOOD_GOOD_NE:
      case R_GOOD_GOOD_EQ:
      case R_GOOD_BAD:
      case R_BAD_GOOD:
      case R_BAD_BAD:
	/*
	 * NOTE: putting in "j-start+1" here instead of tmp was crashing
	 * MS Visual C++ 5.0. Hence the introduction of a temporary variable.
	 * I cannot see why this crashes the Microsoft compiler, but I guess
	 * it's just bugged.
	 */
	forward[tmp]++;
	reverse[tmp]++;
	break;
	    
      case R_GOOD_NONE:
      case R_BAD_NONE:
	forward[tmp]++;
	break;
	
      case R_NONE_GOOD:
      case R_NONE_BAD:
	reverse[tmp]++;
	break;
	
      default: /* R_NONE_NONE */
	;
      }
    }
    xfree(qual);
    return 0;
}

/*
 * strand_coverage histogram
 */
void display_strand_coverage(GapIO *io, obj_strand_coverage *scov) {
    char cmd[1024];
    int i;
    obj_consistency_disp *c;
    int length;
    int win_num;

    c = result_data(io, scov->cons_id, 0);

    sprintf(cmd, "%s delete all", scov->c_win);
    Tcl_Eval(c->interp, cmd);

    win_num = get_consistency_win_num(c, scov->id);
    for (i = 0; i < c->num_contigs; i++) {
        if (c->num_contigs > 1) {
	    length = ABS(io_clength(io, c->contigs[i]));
	} else {
	    length = c->end - c->start + 1;
	}
	if (scov->problems == 1) {
	  if (scov->strand == 1 || scov->strand == 3) { 
	    plot_strand_coverage(c->interp, scov->forward[i], length, 
				 scov->c_win, io, 
				 c->start+c->contig_offset[c->contigs[i]].offset, 
				 scov->linewidth, scov->colour1, 
				 scov->forward_offset);
	  }
	  if (scov->strand == 2 || scov->strand == 3) { 
	    plot_strand_coverage(c->interp, scov->reverse[i], length, 
				 scov->c_win, io, 
				 c->start+c->contig_offset[c->contigs[i]].offset, 
				 scov->linewidth, scov->colour2, 
				 scov->reverse_offset);
	  }
	} else {
	  if (scov->strand == 1 || scov->strand == 3) { 
	    plot_strand_problems(c->interp, scov->forward[i], length, 
				 scov->c_win, io, 
				 c->start+c->contig_offset[c->contigs[i]].offset, 
				 scov->linewidth, scov->colour1, 
				 scov->forward_offset);
	  }
	  if (scov->strand == 2 || scov->strand == 3) { 
	    plot_strand_problems(c->interp, scov->reverse[i], length, 
				 scov->c_win, io, 
				 c->start+c->contig_offset[c->contigs[i]].offset, 
				 scov->linewidth, scov->colour2, 
				 scov->reverse_offset);
	  }
	}
    }

    scaleCanvas(c->interp, &c->win_list[win_num], 1, "all", 
		c->win_list[win_num]->world->visible, 
		c->win_list[win_num]->canvas); 

    scrollRegion(c->interp, &c->win_list[win_num], 1, 
		 c->win_list[win_num]->world->total, 
		 c->win_list[win_num]->canvas);

    consistency_update_cursors(io, c, 0);
}

/*
 * Removes the strand_coverage (and unplots etc).
 */
static void strand_coverage_shutdown(GapIO *io, obj_strand_coverage *scov) {
    char cmd[1024];
    int i;
    obj_consistency_disp *c;
    int win_num;

    if (NULL != (c = result_data(io, scov->cons_id, 0))) {
	win_num = get_consistency_win_num(c, scov->id);
	delete_consistency_window(c, win_num);
    }
    for (i = 0; i < c->num_contigs; i++) {
      contig_deregister(io, c->contigs[i], strand_coverage_callback, 
			(void *)scov);
    }
    sprintf(cmd, "DeleteStrandCoverage %d %s %s %d\n", *handle_io(io), 
	    scov->frame, scov->c_win, scov->cons_id);
    if (TCL_ERROR == Tcl_Eval(c->interp, cmd)) 
      printf("strand_coverage_shutdown: %s\n", Tcl_GetStringResult(c->interp));

    if (scov->forward) {
        for (i = 0; i < c->num_contigs; i++) {
	  xfree(scov->forward[i]);
	}
	xfree(scov->forward);
    }
    if (scov->reverse) {
        for (i = 0; i < c->num_contigs; i++) {
	  xfree(scov->reverse[i]);
	}
	xfree(scov->reverse);
    }
    xfree(scov);
    if (c->num_wins == 0)
        consistency_shutdown(io, c);
}

/*
 * Updates an obj_strand_coverage structure.
 * This includes recalculating the length, reallocing the memory used,
 * and recalculating the consensus sequence.
 *
 * Returns 0 for success, -1 for failure.
 */
static int update_obj_strand_coverage(GapIO *io, 
					obj_strand_coverage *scov, 
					int contig) {
    int i, j, ret, newl;
    obj_consistency_disp *c;
    int start, end;

#ifdef DEBUG
    printf("UPDATE_OBJ_STRAND_COVERAGE\n");
#endif
    c = result_data(io, scov->cons_id, 0);
    
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

    if (scov->forward[i]) {
      if (NULL == (scov->forward[i] = 
		   (int *)xrealloc(scov->forward[i], (newl+1)*sizeof(int))))
	return -1;
    }
    if (scov->reverse[i]) {
      if (NULL == (scov->reverse[i] = 
		   (int *)xrealloc(scov->reverse[i], (newl+1)*sizeof(int))))
	return -1;
    }

    for (j = 0; j <= newl; j++) {
      scov->forward[i][j] = 0;
      scov->reverse[i][j] = 0;
    }

    ret = calc_strand_coverage(io, c->contigs[i], start, end, 
				 scov->forward[i], scov->reverse[i]);

    return ret;
}

#if 0
static void strand_coverage_list(GapIO *io, 
				 obj_consistency_disp *c,
				 obj_strand_coverage *scov,
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

	/* HACK - to tidy up */
	for (j = 0; j < length-1; j++) {
	  vmessage("%d %d\n", c->start+j, scov->forward[i][j]);
	  vmessage("%d %d\n", c->start+j, scov->reverse[i][j]);
	}
	vmessage("\n");
    }
}
#endif

static void strand_coverage_join(GapIO *io,
				   obj_consistency_disp *c,
				   obj_strand_coverage *scov, 
				   int index)
{
    int length;

    length = c->num_contigs - index - 1;

    if (scov->forward[index])
	xfree(scov->forward[index]);
    if (scov->reverse[index])
	xfree(scov->reverse[index]);
    memmove(&(scov->forward)[index], &(scov->forward)[index+1], 
	    length * sizeof(int*));
    memmove(&(scov->reverse)[index], &(scov->reverse)[index+1], 
	    length * sizeof(int*));

}

/*
 * The callback from the contig registration scheme.
 */
static void strand_coverage_callback(GapIO *io, int contig, void *fdata,
				       reg_data *jdata) {
    obj_strand_coverage *scov = (obj_strand_coverage *)fdata;
    obj_consistency_disp *c;
    int win_num;

    c = result_data(io, scov->cons_id, 0);

    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Strand coverage");
	    return;
	}
    case REG_COMPLEMENT:
    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("strand_coverage COMPLEMENT LENGTH contig %d\n", contig);
#endif
	    if (update_obj_strand_coverage(io, scov, contig) == 0) {
		display_strand_coverage(io, scov);
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
	    printf("strand_coverage QUIT DELETE \n");
#endif
	    strand_coverage_shutdown(io, scov);
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
		strand_coverage_list(io, c, scov, contig);
		end_message();
		break;
#endif
	    case 0: /* Remove */
		strand_coverage_shutdown(io, scov);
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
	  jdata->generic.data = (void *)scov->cons_id;
	case TASK_CANVAS_REDRAW:
#ifdef DEBUG
	    printf("strand_coverage TASK_CANVAS_REDRAW\n");
#endif
	    display_strand_coverage(io, scov);
	    break;
	case TASK_CANVAS_WORLD:
	    {
		task_world_t *tw = (task_world_t *)jdata->generic.data;
		int cx = tw->canvasx;
		double wx, wy;
#ifdef DEBUG
		printf("STRAND_COVERAGE TASK_CANVAS_WORLD\n");
#endif
		win_num = get_consistency_win_num(c, scov->id);
		CanvasToWorld(c->win_list[win_num]->canvas, cx, 0, &wx, &wy);
		tw->basex = wx - c->contig_offset[tw->cnum].offset;
		break;
	    }
	case TASK_CONS_WORLD:
	    {
		d_box *cw = (d_box *)jdata->generic.data;

		win_num = get_consistency_win_num(c, scov->id);
		cw->x1 = c->win_list[win_num]->world->total->x1;
		cw->y1 = c->win_list[win_num]->world->total->y1;
		cw->x2 = c->win_list[win_num]->world->total->x2;
		cw->y2 = c->win_list[win_num]->world->total->y2;
		
		break;
	    }
	case TASK_CONS_JOIN:
	    {
#ifdef DEBUG
	      printf("strand_coverage JOIN_TO contig %d join %d\n", 
		     contig, jdata->join.contig); 
#endif
	      strand_coverage_join(io, c, scov, (int)jdata->generic.data);
	      break;
	    }
	}
    }
}


/*
 * Registers and initialises a quality buffer for a particular contig.
 */
int strand_coverage_reg(GapIO *io, 
			Tcl_Interp *interp,
			char *frame,
			char *win,
			int cons_id,
			int strand,
			int problems) 
{
    obj_strand_coverage *scov;
    obj_consistency_disp *c;
    int id;
    int i, j;
    int start, end, length;
    char *val;

    c = result_data(io, cons_id, 0);

    /* first check that there are enough available windows */
    if (c->num_wins >= MAX_NUM_WINS)
	return -1;

    if (NULL == (scov = (obj_strand_coverage*)xmalloc(sizeof(obj_strand_coverage))))
	return -1;

    if (NULL == (scov->forward = (int **)xmalloc(c->num_contigs*sizeof(int*))))
	return -1;

    if (NULL == (scov->reverse = (int **)xmalloc(c->num_contigs*sizeof(int*))))
	return -1;

    id = register_id();
    
    scov->id = id;
    scov->cons_id = cons_id;
    strcpy(scov->c_win, win);
    strcpy(scov->frame, frame);
    scov->forward_offset = get_default_int(interp, gap_defs,
				      "STRAND_COVERAGE.FORWARD_OFFSET");
    scov->reverse_offset = get_default_int(interp, gap_defs,
				      "STRAND_COVERAGE.REVERSE_OFFSET");

    scov->linewidth = get_default_int(interp, gap_defs,
				      "STRAND_COVERAGE.LINEWIDTH");
    val = get_default_string(interp, gap_defs, "STRAND_COVERAGE.COLOUR1"); 
    strcpy(scov->colour1, val);
    val = get_default_string(interp, gap_defs, "STRAND_COVERAGE.COLOUR2"); 
    strcpy(scov->colour2, val);
    scov->strand = strand;
    scov->problems = problems;

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

	if (NULL==(scov->forward[i]=(int *)xmalloc((length+1)*sizeof(int)))) {
	    return -1;
	}
	if (NULL==(scov->reverse[i]=(int *)xmalloc((length+1)*sizeof(int)))) {
	    return -1;
	}
	/* initialise histogram and min and max values */
	for (j = 0; j <= length; j++) {
	    scov->forward[i][j] = 0;
	    scov->reverse[i][j] = 0;
	}

        calc_strand_coverage(io, c->contigs[i], start, end, scov->forward[i],
			     scov->reverse[i]);
    }

    add_consistency_window(io, c, win, 'x', id, c->orig_total->x1,
			   0, c->orig_total->x2, 0);

    display_strand_coverage(io, scov);

    for (i = 0; i < c->num_contigs; i++) {
        contig_register(io, c->contigs[i], strand_coverage_callback, 
			(void *)scov, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_CURSOR_NOTIFY | REG_GENERIC, 
			REG_TYPE_STRAND_COVERAGE);
    }

#ifdef DEBUG
    printf("STRAND_COVERAGE ID %d\n", id);
#endif
    return id;
}
