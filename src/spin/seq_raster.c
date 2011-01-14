#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "seq_results.h"
#include "seq_reg.h"
#include "seq_raster.h"
#include "tkRaster.h"
#include "misc.h"
#include "tcl_utils.h"
#include "text_output.h"
#include "spin_globals.h"

/* FIXME: should move to os.h - can use rint on unix machines */
#define ROUND(x)	((x) < 0.0 ? ceil((x) - 0.5) : floor((x) + 0.5))

#define NUM_COL 10
static char *colour[NUM_COL];

void init_raster_colour(Tcl_Interp *interp)
{
    colour[0] = get_default_string(interp, spin_defs, "RASTER.COLOUR.0");
    colour[1] = get_default_string(interp, spin_defs, "RASTER.COLOUR.1");
    colour[2] = get_default_string(interp, spin_defs, "RASTER.COLOUR.2");
    colour[3] = get_default_string(interp, spin_defs, "RASTER.COLOUR.3");
    colour[4] = get_default_string(interp, spin_defs, "RASTER.COLOUR.4");
    colour[5] = get_default_string(interp, spin_defs, "RASTER.COLOUR.5");
    colour[6] = get_default_string(interp, spin_defs, "RASTER.COLOUR.6");
    colour[7] = get_default_string(interp, spin_defs, "RASTER.COLOUR.7");
    colour[8] = get_default_string(interp, spin_defs, "RASTER.COLOUR.8");
    colour[9] = get_default_string(interp, spin_defs, "RASTER.COLOUR.9");
}

char * get_raster_colour(void) 
{
    static int col_index = 0;

    col_index = (col_index + 1) % NUM_COL;
    return(colour[col_index]);
}

/*
 * replots all visible plots in the raster widget at the current zoom level
 */
void ReplotAllCurrentZoom(Tcl_Interp *interp, 
			  char *raster_win)
{
    char cmd[1024];
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    int num_results;
    double wx0, wy0, wx1, wy1;

#ifdef START_DEBUG
    printf("START ReplotAllCurrentZoom\n");
#endif
    Tcl_GetCommandInfo(interp, raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    num_results = seq_num_results();

    /* check that there is at least one result to replot! */
    if (num_results <= 0) {
	tk_RasterClear(raster);
	return;
    }
    tk_RasterClear(raster);

    GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);

#ifdef DEBUG
    printf("&&&&&&COORDS %f %f %f %f\n", wx0, wy0, wx1, wy1);
#endif
    /* current extents of the raster */
    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
#ifdef DEBUG
    printf("&&&&&&WORLD %f %f %f %f\n", wx0, wy0, wx1, wy1);
#endif

    if ((wx1 == 0) && (wy1 == 0)) {
#ifdef DEBUG
	printf("all results are hidden !\n");
#endif
	return;
    }
    /* draw rulers */
    sprintf(cmd, "rasterHRuler %s %f %f ", raster_win, wx0, wx1);
    if (TCL_OK != Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "ReplotAllCurrentZoom", "%s\n", Tcl_GetStringResult(interp));

    sprintf(cmd, "rasterVRuler %s %f %f", raster_win, wy0, wy1);
    if (TCL_OK != Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "ReplotAllCurrentZoom", "%s \n", Tcl_GetStringResult(interp));

    /* replot results ONCE only */
    RasterCallPlotFunc(raster, RASTER_REPLOT_ZOOM,
		      (int)wx0, (int)wy0, (int)wx1, (int)wy1);

    tk_RasterRefresh(raster);
}


int GetRasterWindowSize(Tcl_Interp *interp,
			char *raster_win,
			double *x0, double *y0, double *x1, double *y1)
{
    int i;
    int win_list_argc;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    double wx0, wy0, wx1, wy1;
    int retval = -1;
    char **win_list_argv = NULL;

    *x0 = DBL_MAX;
    *y0 = DBL_MAX;
    *x1 = 0;
    *y1 = 0;

    win_list_argv = GetRasterWindowList(interp, raster_win, &win_list_argc);

    for (i = 0; i < win_list_argc; i++) {
	if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
	    goto cleanup;
	raster = (Tk_Raster*)info.clientData;
	
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

	if (wx0 < *x0) *x0 = wx0;
	if (wy0 < *y0) *y0 = wy0;
	if (wx1 > *x1) *x1 = wx1;
	if (wy1 > *y1) *y1 = wy1;
    }
#ifdef DEBUG
    printf("GetRasterWindowSize %f %f %f %f\n", *x0, *y0, *x1, *y1);
#endif
    retval = 0;


cleanup:
    if(win_list_argv) Tcl_Free((char *)win_list_argv);
    return retval;
}

/*
 * replots all visible plots in each raster in a raster window
 * and sets their new scroll size to be that of the largest scroll size
 */
int
ReplotAllRasterWindow(Tcl_Interp *interp, 
		      char *raster_win)
{
    char cmd[1024];
    char *window;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    int num_results;
    double x0, x1, y0, y1;
    double wx0, wy0, wx1, wy1;
    int win_list_argc;
    int i;
    int retval = -1;
    char **win_list_argv = NULL;

 
#ifdef START_DEBUG
    printf("START ReplotAllRasterWindow %s \n", raster_win);
#endif

    win_list_argv = GetRasterWindowList(interp, raster_win, &win_list_argc);

    /* get window scroll region for raster_win */
    GetRasterWindowSize(interp, raster_win, &x0, &y0, &x1, &y1);
    num_results = seq_num_results();

    for (i = 0; i < win_list_argc; i++) {

	window = win_list_argv[i];
	if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
	    goto cleanup;
	raster = (Tk_Raster*)info.clientData;

	/* set new x scroll region */
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
	RasterSetWorldScroll(raster, x0, wy0, x1, wy1);

	/* clear current raster */
	tk_RasterClear(raster);
	
	if ((wx1 == 0) && (wy1 == 0)) {
#ifdef DEBUG
	    printf("all results are hidden !\n");
#endif
	    retval = 0;
	    goto cleanup;
	}
	/* draw vertical ruler */
	sprintf(cmd, "rasterVRuler %s %f %f", window, wy0, wy1);
	
	if (TCL_OK != Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "ReplotAllRasterWindow", "%s \n", Tcl_GetStringResult(interp));
	    
	/* replot results ONCE only */
	RasterCallPlotFunc(raster, RASTER_REPLOT_ZOOM,
			  (int)wx0, (int)wy0, (int)wx1, (int)wy1);
	tk_RasterRefresh(raster);
    }

    /* draw horizontal ruler */
    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	goto cleanup;
    raster = (Tk_Raster*)info.clientData;
    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1); 

    sprintf(cmd, "rasterHRuler %s %f %f ", raster_win, wx0, wx1);
    if (TCL_OK != Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "ReplotAllRasterWindow", "%s\n", Tcl_GetStringResult(interp));
    retval = 0;


cleanup:
    if(win_list_argv) Tcl_Free((char *)win_list_argv);
    return retval;
}


/*
 * replot raster plot when a window resize or zooming has taken place
 */
void ReplotAllZoom(Tcl_Interp *interp,
		   char *raster_win)
{
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    int num_results;
    char cmd[1024];
    double wx0, wy0, wx1, wy1;
    double x0, y0, x1, y1;
    double ty0, ty1;

#ifdef START_DEBUG

    printf("ReplotAllZoom\n"); 
#endif
    Tcl_GetCommandInfo(interp, raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    /* current extents of the raster */
    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
    GetRasterCoords(raster, &x0, &ty0, &x1, &ty1);

    y1 = rasterY(raster, ty0);
    y0 = rasterY(raster, ty1);

    /* draw rulers */
    sprintf(cmd, "rasterHRuler %s %f %f ", raster_win, x0, x1);
    if (TCL_OK != Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "ReplotAllZoom", "%s\n", Tcl_GetStringResult(interp));

    sprintf(cmd, "rasterVRuler %s %f %f", raster_win, y0, y1);
    if (TCL_OK != Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "ReplotAllZoom", "%s \n", Tcl_GetStringResult(interp));

    /* check that there is at least one result to replot! */
    num_results = seq_num_results();
    tk_RasterClear(raster);
    if (num_results <= 0) {
	return;
    }
    /* replot results ONCE only */
    RasterCallPlotFunc(raster, RASTER_REPLOT_SLIVER,
		       (int)wx0, (int)wy0, (int)wx1, (int)wy1);

    tk_RasterRefresh(raster);
}

int seq_raster_shutdown(int r_id,
			RasterResult *result)
{
    seq_reg_generic gen;
    int i;
    int num;

#ifdef DEBUG
    printf("seq_raster_shutdown status %d\n", result->status);
#endif

    /* remove all plots in raster */
    gen.job = SEQ_GENERIC;
    gen.task = TASK_RASTER_QUIT_RESULTS;
    seq_result_notify(r_id, (seq_reg_data *)&gen, 1);

    /* remove the raster */
    /* need to deregister from all seq_num the raster is registered with */
    for (i = 0; i < result->num_seq_id; i++) {
	num = GetSeqNum(result->seq[i].seq_id);
	seq_deregister(num, seq_raster_callback, result);
	delete_cursor(num, result->cursor[i]->id, 0);
    }
    return result->status;
}

void raster_shutdown(Tcl_Interp *interp,
		     char *raster_win,
		     RasterResult *result)
{
    char *tmp;

#ifdef DEBUG
    printf("raster_shutdown %s\n", raster_win);
#endif

    tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
    if (TCL_OK != Tcl_VarEval(interp, "removeRaster ", raster_win, " ", tmp, 
			      NULL))
	verror(ERR_WARN, "raster_shutdown",  "%s\n", Tcl_GetStringResult(interp));

    xfree(result->seq);
    xfree(result->cursor);
    xfree(result);
}

RasterResult *raster_id_to_result(int raster_id)
{
    seq_reg_info info;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(raster_id, (seq_reg_data *)&info, 0);
    if (!info.result) {
#ifdef DEBUG
	printf("NO raster result info\n");
#endif
	return NULL;
    }
    return (RasterResult *)info.result;
}

RasterResult *raster_name_to_result(Tcl_Interp *interp, 
				    char * raster_win)
{
  Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
  return (raster_id_to_result(atoi(Tcl_GetStringResult(interp))));
}

/* invoked from tcl move cursor binding command */
int seq_raster_move_cursor(int raster_id,
			   Tk_Raster *raster,
			   int cursor_id,
			   int pos,
			   int direction)
{
    seq_cursor_notify cn;
    seq_reg_info info;
    RasterResult *result;
    cursor_t *cursor;
    double wx, wy;
    int seq_num = -1;
    double x0, y0, x1, y1;

    /* convert raster coords into world coords */
    RasterToWorld (raster, pos, pos, &wx, &wy);

    /* ensure cursors don't go off the raster */
    RasterGetWorldScroll(raster, &x0, &y0, &x1, &y1);

    /* have to turn world y upside down after RasterToWorld */
    wy = rasterY(raster, wy);

    if (direction == HORIZONTAL) {
	if (wx < x0) {
	    wx = x0;
	}
	if (wx > x1) {
	    wx = x1;
	}
    } else {
 	if (wy < y0) {
	    wy = y0;
	}
	if (wy > y1) {
	    wy = y1;
	}
   }
    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(raster_id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return -1;
    }
    result = (RasterResult *)info.result;
    cursor = find_cursor(&seq_num, cursor_id, -1);
#ifdef DEBUG
    printf("seq_raster_move_cursor pos %d wx %f cursorpos %d cursor_id %d %d seq_num %d\n", pos, wx, 
	   cursor->abspos, cursor_id, cursor->id, seq_num);
    printf("seq_raster_move_cursor pos %d wx %f\n", pos, wx);
#endif

    result->cursor_array[cursor->id].prev_pos = cursor->abspos;
#ifdef REMOVE
    cursor->prev_pos = cursor->abspos;
#endif
    if (direction == HORIZONTAL) {
	cursor->abspos = ROUND(wx);
    } else {
	cursor->abspos = ROUND(wy);
    }
    cursor->job = CURSOR_MOVE;

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = cursor;
    seq_notify(seq_num, (seq_reg_data *)&cn); 
    return 0;
}

/*
 * return the first found seqed cursor in direction HORIZONTAL or VERTICAL
 * (called from tcl)
 */
int seq_raster_find_edcursor(int raster_id, 
			     Tk_Raster *raster, 
			     int pos,
			     int direction, 
			     int *seq_id)
{
    int i = 0;
    RasterResult *result;
    *seq_id = -1;

    if (NULL == (result = raster_id_to_result(raster_id))) 
	return -1;

    for (i = 0; i < result->num_seq_id; i++) {

#ifdef DEBUG
	printf("seq_raster_find_edcursor: i %d seq_id %d private %d direction %d\n",
	       i, result->seq[i].seq_id, result->cursor[i]->private, result->seq[i].direction);
#endif
	if (result->cursor[i]->private && 
	    result->seq[i].direction == direction) {
	    *seq_id = result->seq[i].seq_id;
	    return result->cursor[i]->id;
	}
    }
    /* no editor cursor displayed */
    for (i = 0; i < result->num_seq_id; i++) {
	if (result->seq[i].direction == direction) {
	    *seq_id = result->seq[i].seq_id;
	    break;
	}
    }
    return -1;
}


/*
 * ---------------------------------------------------------------------------
 * Cursor manipulation
 * ---------------------------------------------------------------------------
 */

/*
 * Updates the scroll region to ensure that the cursor is visible.
 * The absolute cursor position (in bases) on the display is 'cursorx'.
 *
 * Returns whether we've scrolled the canvas.
 */
int raster_cursor_show(Tcl_Interp *interp, 
		       Tk_Raster *raster,
		       char *raster_win,
		       int cursor,
		       int direction)
{
    int x1, x2, dx;
    char cmd[1024];
    double fract;
    double wx0, wy0, wx1, wy1;
    double sx0, sy0, sx1, sy1;
    char *r_win;
    char *raster_stem;
    double w0, w1, s0, s1;
    int id;

#ifdef DEBUG
    printf("raster_cursor_show %s\n", raster_win);
#endif

    GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);
    RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);

    if (direction == HORIZONTAL) {
	w0 = wx0;
	w1 = wx1;
	s0 = sx0;
	s1 = sx1;
    } else {
	w0 = wy0;
	w1 = wy1;
	s0 = sy0;
	s1 = sy1;
	cursor = rasterY(raster, cursor);
    }
    
    /* Is it visible? */
    if (cursor >= w0 && cursor <= w1) {
#ifdef DEBUG
	printf("cursor is visible \n");
#endif
	return 0;
    }

    dx = w1 - w0;
    x1 = (double)cursor - dx / 2;
    if (x1 < s0)
	x1 = s0;
    if (x1 > s1 - dx)
	x1 = s1 - dx;
    x2 = x1 + dx;

    /* Compute xview fraction and scroll */
    fract = (double)(x1 - s0) / (s1 - s0);
#ifdef DEBUG
    printf("cursor %d sx0 %f sx1 %f wx0 %f wx1 %f fract %f x1 %f\n",
	   cursor, s0, s1, w0, w1, fract, x1);
#endif
    Tcl_VarEval(interp, "winfo parent ", raster_win, NULL);
    r_win = strdup(Tcl_GetStringResult(interp));

    Tcl_VarEval(interp, "GetRasterStem ", r_win, NULL);
    raster_stem = strdup(Tcl_GetStringResult(interp));

    /* HACK - to fix ruler_h */
    if (direction == HORIZONTAL) {
	/*
	sprintf(cmd, "after idle scrollXCmd %s %s %s.ruler_h moveto %f", 
		r_win, raster_stem, r_win, fract);
	*/
	/*
	 * 10.03.00 removed "after idle" because this causes a bug if you
	 * bring up a plot, zoom, then bring up an seqed plot (in nip4).
	 * The seqed plot is now scrolled to the left hand side rather than
	 * where the cursor is. This is because the seqed is originally
	 * brought up with the cursor at position 1 and only later repositioned
	 * at the correct cursor position. The after idle command ensured that
	 * the scroll command was executed until after the new position of the
	 * had been determined.
	 */

	sprintf(cmd, "scrollXCmd %s %s %s.ruler_h moveto %f", 
		r_win, raster_stem, r_win, fract);
    } else {
	Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
	id = atoi(Tcl_GetStringResult(interp));
	/*
	sprintf(cmd, "after idle scrollYCmd %s %s.ruler_v%d moveto %f", 
		raster_win, r_win, id, fract);
	*/
	sprintf(cmd, "scrollYCmd %s %s.ruler_v%d moveto %f", 
		raster_win, r_win, id, fract);
    }

    if (TCL_ERROR == (Tcl_Eval(interp, cmd)))
	verror(ERR_WARN, "raster_cursor_show", "%s\n", Tcl_GetStringResult(interp));

    free(r_win);
    free(raster_stem);
    return 1;
}

/*
 * Deletes an existing cursor - calls canvas_cursor_delete in tcl
 */
int raster_cursor_delete(Tcl_Interp *interp,
			 char *raster_win,
			 cursor_t *cursor,
			 int seq_id)
{
    double wx0, wy0, wx1, wy1;
    double sx0, sy0, sx1, sy1;
    int win_list_argc;
    int id_list_argc;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    int i;
    RasterResult *result;
    char cmd[1024];
    int retval           = -1;
    char **id_list_argv  = NULL;
    char **win_list_argv = NULL;



#ifdef DEBUG
    printf("raster_cursor_delete\n");
#endif

    /* HACK 
     * When deleting a window, the win_list in tcl is not updated until after
     * this code is executed, so have to check raster_id_to_result
     */
    win_list_argv = GetRasterWindowList(interp, raster_win, &win_list_argc);
    id_list_argv = GetRasterIdList(interp, raster_win, &id_list_argc);

    for (i = 0; i < win_list_argc; i++) {
	if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
	    goto cleanup;
	raster = (Tk_Raster*)info.clientData;
	GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);
	RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);

	if (NULL == (result = raster_id_to_result(atoi(id_list_argv[i]))))
	    continue;

	result->cursor_array[cursor->id].prev_pos = -1;
	if (result->cursor_array[cursor->id].visible[cursor->direction] == 0)
	    continue;

	SetDrawEnviron(interp, raster, 
		       result->cursor_array[cursor->id].env);

	if (cursor->direction == HORIZONTAL) {
	    RasterDrawLine(raster, cursor->abspos, wy0, 
			   cursor->abspos, wy1);
	    
#ifdef DEBUG
	    printf("^^^^^^^^RasterDrawLine H %d %f %f \n", cursor->abspos, wy0, wy1);
#endif
	    Tcl_VarEval(interp, "winfo parent ", result->raster_win, 
			NULL);
	    sprintf(cmd, "%s.buttons.pos1 configure -text {}", 
		    Tcl_GetStringResult(interp));
	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		printf("raster_cursor_delete: %s\n", Tcl_GetStringResult(interp));
	    }
	} else {
	    RasterDrawLine(raster, wx0, rasterY(raster,cursor->abspos),
			   wx1, rasterY(raster, cursor->abspos));
#ifdef DEBUG
	    printf("^^^^^^^^RasterDrawLine V %d %f %f \n", cursor->abspos, wx0, wx1);
#endif
	    Tcl_VarEval(interp, "winfo parent ", result->raster_win, 
			    NULL);
	    sprintf(cmd, "%s.buttons.pos2 configure -text {}", 
		    Tcl_GetStringResult(interp));
	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		printf("raster_cursor_delete: %s\n", Tcl_GetStringResult(interp));
	    }
	}

	result->cursor_array[cursor->id].visible[cursor->direction] = 0;

	tk_RasterRefresh(raster);
    }
    retval = 0;


cleanup:
    if(win_list_argv) Tcl_Free((char *)win_list_argv);
    if(id_list_argv)  Tcl_Free((char *)id_list_argv);
    return retval;
}

void raster_cursor_remove(Tcl_Interp *interp, 
			  Tk_Raster *raster,
			  cursor_t *cursor,
			  RasterResult *result,
			  int direction)
{
    double wx0, wy0, wx1, wy1;
    double sx0, sy0, sx1, sy1;

#ifdef DEBUG
    printf("raster_cursor_remove\n");
#endif

    if (!result->cursor_array[cursor->id].visible[direction]) {
	return;
    }

    GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);
    RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);
    
    if (result->cursor_array[cursor->id].env < 0) {
	result->cursor_array[cursor->id].env = 
	    raster_init_env(interp, raster, cursor);
    }
    SetDrawEnviron(interp, raster, result->cursor_array[cursor->id].env);

    if (direction == HORIZONTAL && cursor->direction == HORIZONTAL) {
	RasterDrawLine(raster, cursor->abspos, wy0, cursor->abspos, wy1);
#ifdef DEBUG
	printf("@@@@@@@@RasterDrawLine H %d %d %f %f dir %d\n", cursor->abspos, result->cursor_array[cursor->id].prev_pos, wy0, wy1, cursor->direction);
#endif

    } else if (direction == VERTICAL && cursor->direction == VERTICAL){
	RasterDrawLine(raster, wx0, rasterY(raster, cursor->abspos), 
		       wx1, rasterY(raster, cursor->abspos));
#ifdef DEBUG
	printf("@@@@@@@@RasterDrawLine V %d %d %f %f dir %d\n", cursor->abspos, result->cursor_array[cursor->id].prev_pos, wy0, wy1, cursor->direction);
#endif

    }

#ifdef DEBUG
    printf("VIS1 id %d dir %d vis 0\n", cursor->id, direction);
#endif
    result->cursor_array[cursor->id].visible[direction] = 0;

    tk_RasterRefresh(raster);
}
/*
 * Moves and existsing cursor - calls canvas_cursor_move in tcl
 *
 * Returns whether we've scrolled the canvas.
 */
int raster_cursor_move(Tcl_Interp *interp, 
		       Tk_Raster *raster,
		       cursor_t *cursor,
		       RasterResult *result,
		       int cursor_show,
		       int direction)
{
    double wx0, wy0, wx1, wy1;
    double sx0, sy0, sx1, sy1;
    int r = 0;
    char cmd[1024];

#ifdef DEBUG
    printf("raster_cursor_move cursor_id %d visible %d raster %s\n", 
	   cursor->id, result->cursor_array[cursor->id].visible[direction],
	   result->raster_win);
#endif
    GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);
    RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);

    if (result->cursor_array[cursor->id].env < 0) {
	result->cursor_array[cursor->id].env = 
	    raster_init_env(interp, raster, cursor);
    }
    SetDrawEnviron(interp, raster, result->cursor_array[cursor->id].env);

    /* draw over previous cursor */
    if (direction == HORIZONTAL && cursor->direction == HORIZONTAL) {
	if (result->cursor_array[cursor->id].visible[HORIZONTAL]) {
	    RasterDrawLine(raster, result->cursor_array[cursor->id].prev_pos, wy0, result->cursor_array[cursor->id].prev_pos, wy1);
#ifdef DEBUG
	    printf("**********RasterDrawLine H %d %f %f \n", result->cursor_array[cursor->id].prev_pos, wy0, wy1);
#endif
	}
    } else if (direction == VERTICAL && cursor->direction == VERTICAL) {
	if (result->cursor_array[cursor->id].visible[VERTICAL]) {
	    RasterDrawLine(raster, wx0, rasterY(raster, result->cursor_array[cursor->id].prev_pos), 
			   wx1, rasterY(raster, result->cursor_array[cursor->id].prev_pos));
#ifdef DEBUG
	    printf("**********RasterDrawLine V %d %f %f \n", result->cursor_array[cursor->id].prev_pos, wx0, wx1);
#endif
	}
    }
    if (direction == HORIZONTAL && cursor->direction == HORIZONTAL) {
	RasterDrawLine(raster, cursor->abspos, wy0, cursor->abspos, wy1);
	Tcl_VarEval(interp, "winfo parent ", result->raster_win, NULL);
	sprintf(cmd, "%s.buttons.pos1 configure -text %d", 
		Tcl_GetStringResult(interp), cursor->abspos);
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	    printf("raster_cursor_move: %s\n", Tcl_GetStringResult(interp));
	}

#ifdef DEBUG
    printf("!!!!!!!!!!RasterDrawLine H %d %f %f \n", cursor->abspos, wy0, wy1);
#endif
    } else if (direction == VERTICAL && cursor->direction == VERTICAL) {
	RasterDrawLine(raster, wx0, rasterY(raster,cursor->abspos), 
		       wx1, rasterY(raster, cursor->abspos));
	Tcl_VarEval(interp, "winfo parent ", result->raster_win, NULL);
	sprintf(cmd, "%s.buttons.pos2 configure -text %d", 
		Tcl_GetStringResult(interp), cursor->abspos);
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	    printf("raster_cursor_move: %s\n", Tcl_GetStringResult(interp));
	}
#ifdef DEBUG
    printf("!!!!!!!!!!RasterDrawLine V %d %f %f \n", cursor->abspos, wx0, wx1);
#endif
    }

#ifdef DEBUG
    printf("VIS2 id %d dir %d vis 1\n", cursor->id, direction);
#endif
    result->cursor_array[cursor->id].visible[direction] = 1;

    /* RasterSetCursorPos(raster, seqed_id, pos); */

    tk_RasterRefresh(raster);

    /* Make sure the cursor is still visible */
    
    if (cursor_show) 
	r = raster_cursor_show(interp, raster, result->raster_win, 
			       cursor->abspos, direction);
			       
    return r; 
}

int raster_init_env(Tcl_Interp *interp,
		    Tk_Raster *raster,
		    cursor_t *cursor)
{
    char *opts[7];
    int index;
    int fg_pixel, bg_pixel;

    /* raster need colour environment */
    if (NULL == (opts[1] = (char *)xmalloc((strlen(cursor->colour)+1) * sizeof(char))))
	return -1;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;

    opts[0] = "-fg";
    strcpy(opts[1], cursor->colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", cursor->line_width);
    opts[4] = "-function";
    opts[5] = "xor";
    opts[6] = NULL;
    index = CreateDrawEnviron(interp, raster, 6, opts);
    SetDrawEnviron(interp, raster, index);

    bg_pixel = GetBgPixel(raster);
    fg_pixel = GetFgPixel(interp, raster, index);
    SetFgPixel(interp, raster, index, fg_pixel ^ bg_pixel);

    xfree(opts[1]);
    xfree(opts[3]);
    return index;
}

void raster_update_cursor(RasterResult *result,
			  cursor_t *cursor,
			  int seq_id,
			  Tk_Raster *raster,
			  int cursor_show,
			  int direction)
{
    int index = -1, i;

    for (i = 0; i < result->num_seq_id; i++) {
	if (result->seq[i].seq_id == seq_id && 
	    result->cursor[i]->direction == direction) {
	    index = i;
	    break;
	}
    }
    if (index >= 0)
	raster_cursor_refresh(result->interp, raster, cursor, 
			      result->cursor[index], seq_id, result, 
			      cursor_show, direction);
}

void remove_all_raster_cursors(Tcl_Interp *interp,
			       Tk_Raster *raster,
			       RasterResult *result)
{
    int i;
    cursor_t *gc;

    for (i = 0; i < result->num_seq_id; i++) {
	/* need to loop through ALL cursors */
	for (gc = result->cursor[i]; gc; gc=gc->next) {
	    raster_cursor_remove(interp, raster, gc, result,
				 result->seq[i].direction);
	}    
    }
}

void raster_update_cursors(RasterResult *result,
			   Tk_Raster *raster)
{
    int i;
    cursor_t *gc;

    for (i = 0; i < result->num_seq_id; i++) {
	/* need to loop through ALL cursors */
	for (gc = result->cursor[i]; gc; gc = gc->next) {
	    raster_update_cursor(result, gc, result->seq[i].seq_id, 
				 raster, 0, gc->direction);
	    result->cursor_array[gc->id].prev_pos = gc->abspos;
#ifdef REMOVE
	    gc->prev_pos = gc->abspos;
	    raster_update_cursor(result, result->cursor[i], result->seq[i].seq_id, 
				 raster, 0, result->seq[i].direction);
	    result->cursor[i]->prev_pos = result->cursor[i]->abspos;
#endif
	}
    }
}

void update_raster_cursor(int new_id,
			  int old_id) 
{
    RasterResult *new_result;
    RasterResult *old_result;
    int i;
    int num_cursors;

    if (NULL == (new_result = raster_id_to_result(new_id))) 
	return;
    if (NULL == (old_result = raster_id_to_result(old_id)))
	return;

    /* HACK - to tidy up */
    num_cursors = get_num_cursors();
    for (i = 0; i < num_cursors; i++) {
	if (old_result->cursor_array[i].env > -1) {
	    new_result->cursor_array[i].env = -1;
	}
    }
}
/*
 * Refreshes the cursor drawing. This involves moving, creating or deleting
 * as the position and ref counts change.
 * Returns whether we've scrolled the canvas and updates the visibility
 * pointer.
 */
int raster_cursor_refresh(Tcl_Interp *interp, 
			  Tk_Raster *raster, 
			  cursor_t *changed_c, 
			  cursor_t *raster_c,
			  int seq_id,
			  RasterResult *result,
			  int cursor_show,
			  int direction)
{
    int r;

    /* raster_c = NULL; */

    /* Check for deletion */
    if (changed_c->job & CURSOR_DELETE) {

	/* check the cursor has been referenced */
	/* if (changed_c->refs > 0) { */

	if (result->cursor_array[changed_c->id].visible[direction]) {
	    raster_cursor_delete(interp, result->raster_win, changed_c, 
				 seq_id);
	}
#ifdef DEBUG
	    printf("VIS3 id %d dir %d vis 0\n", changed_c->id, direction);
#endif
#ifdef DEBUG
	    printf("DELETED \n");
#endif
	    /* } */
	return 0;
    }
    
    /* If the cursor is our own with only one reference, don't draw it */

    /* 
     * commented the following code out to quickly overcome bugs in the
     * communication of n/sip and gap to ensure that the cursors were always
     * drawn
     */
#ifdef DEBUG
    if (changed_c == raster_c && changed_c->refs <= 1) {

	/* If it's still visible, then remove it */
	if (result->cursor_array[changed_c->id].visible[direction]) {

	    raster_cursor_delete(interp, result->raster_win, changed_c, 
				 seq_id);

#ifdef DEBUG
	    printf("VIS4 id %d dir %d vis 0\n", changed_c->id, direction);
#endif
	    result->cursor_array[changed_c->id].visible[direction] = 0;
#ifdef DEBUG
	    printf("DELETED 2\n");
#endif
	}
	return 0;
    }
#endif
 
    /* Move it, creating it in the process if neccesary */
    r = raster_cursor_move(interp, raster, changed_c, result, cursor_show,
			   direction);
    return r;
}

void UpdateScaleBars(Tcl_Interp *interp,
		     double old_xscroll, 
		     double new_xscroll,
		     double old_yscroll, 
		     double new_yscroll,
		     char *raster_new,
		     int id_old,
		     int id_new,		     
		     int job) /* 0 = add, 1 = remove */
{
    double value;
    char cmd[1024];
    char *tmp;
    int current_xmag_value, current_ymag_value;
    char *r_win;
    double xscroll, yscroll;
    
    if (old_xscroll == -1 && old_yscroll == -1)
	return;
    if (new_xscroll == -1 && new_yscroll == -1)
	return;

#ifdef DEBUG
    printf("UpdateScaleBars\n");
    printf("old zoom %d new zoom %d\n", GetRasterZoom(id_old),
	   GetRasterZoom(id_new));
#endif

#ifdef REMOVE
    /* need to check if plot is zoomable */
    if (GetRasterZoom(id_old) <= 0 || GetRasterZoom(id_new) <= 0) {
#ifdef DEBUG
	printf("NOT ZOOMABLE\n");
#endif
	return;
    }
#endif

    Tcl_VarEval(interp, "winfo parent ", raster_new, NULL);
    r_win = strdup(Tcl_GetStringResult(interp));

    /* 
     * find the current scale value
     */
    tmp = get_default_string(interp, tk_utils_defs, w("RASTER.SCALEX.WIN"));
    Tcl_VarEval(interp, r_win, tmp, " get", NULL);
    current_xmag_value = atoi(Tcl_GetStringResult(interp));

#ifdef DEBUG
    printf("old %f new %f\n", old_xscroll, new_xscroll);
    printf("old %f new %f\n", old_yscroll, new_yscroll);
#endif

    if (job == 0) {
	/* adding a plot */
	xscroll = RasterSetXMag(interp, new_xscroll, current_xmag_value);
	value = RasterGetXMag(interp, old_xscroll, xscroll);
#ifdef DEBUG
	printf("current %d xscroll %f value %f\n", current_xmag_value,
	       xscroll, value);
#endif
	sprintf(cmd, "%s%s set %d", r_win, tmp, (int)ROUND(value));
	Tcl_Eval(interp, cmd);
    } else {
	/* removing a plot */
	value = RasterGetXMag(interp, old_xscroll, new_xscroll);
	sprintf(cmd, "%s%s set %d", r_win, tmp, (int)ROUND(value));
	Tcl_Eval(interp, cmd);
    }

    /* 
     * find the current scale value
     */
    tmp = get_default_string(interp, tk_utils_defs, w("RASTER.SCALEY.WIN"));
    Tcl_VarEval(interp, r_win, tmp, " get", NULL);
    current_ymag_value = atoi(Tcl_GetStringResult(interp));

    if (job == 0) {
	/* adding a plot */
	yscroll = RasterSetYMag(interp, new_yscroll, current_ymag_value);
	value = RasterGetYMag(interp, old_yscroll, yscroll);
#ifdef DEBUG
	printf("current %d yscroll %f value %f\n", current_ymag_value,
	       yscroll, value);
#endif
	sprintf(cmd, "%s%s set %d", r_win, tmp, (int)ROUND(value));
	Tcl_Eval(interp, cmd);
    } else {
	/* removing a plot */
	value = RasterGetYMag(interp, old_yscroll, new_yscroll);
	sprintf(cmd, "%s%s set %d", r_win, tmp, (int)ROUND(value));
	Tcl_Eval(interp, cmd);
    }

    free(r_win);
}

void UpdateZoomList(Tcl_Interp *interp,
		    double o_wx0, double o_wy0, double o_wx1, double o_wy1,
		    double n_wx0, double n_wy0, double n_wx1, double n_wy1,
		    char *parent,
		    int job)
{
    double x0, y0, x1, y1;
    char cmd[1024];

    x0 = (MIN(o_wx0, n_wx0));
    y0 = (MIN(o_wy0, n_wy0));
    x1 = (MAX(o_wx1, n_wx1));
    y1 = (MAX(o_wy1, n_wy1));

#ifdef DEBUG
    printf("x0 %f x1 %f y0 %f y1 %f\n n_wx0 %f n_wx1 %f n_wy0 %f n_wy1 %f\n",
	   x0, x1, y0, y1, n_wx0, n_wx1, n_wy0, n_wy1);
#endif

    /* 
     * if any of the extents have changed, update the zoom list for new
     * raster window
     */
    if (x0 != n_wx0 || y0 != n_wy0 || x1 != n_wx1 || y1 != n_wy1) {
	sprintf(cmd, "update_zoom_list %s %d {%d %d %d %d}",
		parent, job, (int)x0, (int)y0, (int)x1, (int)y1);
	
	if (TCL_ERROR == (Tcl_Eval(interp, cmd))) {
	    printf("UpdateZoomList %s\n", Tcl_GetStringResult(interp));
	}
    }
    if (job == 1) {
	sprintf(cmd, "update_zoom_list %s %d {%d %d %d %d}",
		parent, job, (int)x0, (int)y0, (int)x1, (int)y1);
	
	if (TCL_ERROR == (Tcl_Eval(interp, cmd))) {
	    printf("UpdateZoomList %s\n", Tcl_GetStringResult(interp));
	}
    }
}

/*
 * find the largest x and y values of all results in a single raster
 */
int FindRasterSize(int raster_id,
		   d_line **max_size)
{
    seq_reg_generic gen;
    d_line *max;

    gen.job = SEQ_GENERIC;

    gen.task = TASK_RASTER_MAX_SIZE;
    gen.result = NULL;
    seq_result_notify(raster_id, (seq_reg_data *)&gen, 0);

    max = (d_line *)gen.result;
    /* max will be NULL if seq_result_notify failed to find anything */
    if (max == NULL)
	return -1;

#ifdef DEBUG
    printf("FindRasterSize %d %f %f %f %f \n", raster_id, max->x0, max->y0, max->x1, max->y1);
#endif
    if ((max->x0 == DBL_MAX/2) && (max->y0 == DBL_MAX/2) && 
	(max->x1 == -DBL_MAX/2) && (max->y1 == -DBL_MAX/2)) {
#ifdef DEBUG
	printf("Raster no longer exists \n");
#endif
	*max_size = max;
	return -1;
    }
    *max_size = max;
    return 0;
}

int SetRasterWindowSize(Tcl_Interp *interp,
			char *raster_win)
{
    int i;
    int win_list_argc;
    double x0, y0, x1, y1;
    double wx0, wy0, wx1, wy1;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    int  retval = -1;
    char **win_list_argv = NULL;


#ifdef DEBUG
    printf("START SetRasterWindowSize %s\n", raster_win);
#endif

    win_list_argv = GetRasterWindowList(interp, raster_win, &win_list_argc);   

    x0 = DBL_MAX;
    y0 = DBL_MAX;
    x1 = 0;
    y1 = 0;

    for (i = 0; i < win_list_argc; i++) {
#ifdef DEBUG
	printf("i %d argv %s \n", i, win_list_argv[i]);
#endif
	if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
	    goto cleanup;
	raster = (Tk_Raster*)info.clientData;

	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
	if (wx0 < x0) x0 = wx0;
	if (wy0 < y0) y0 = wy0;
	if (wx1 > x1) x1 = wx1;
	if (wy1 > y1) y1 = wy1;
    }
#ifdef DEBUG
    printf("FINAL %f %f %f %f \n", x0, y0, x1, y1);
#endif

    for (i = 0; i < win_list_argc; i++) {
	if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
	    goto cleanup;
	raster = (Tk_Raster*)info.clientData;
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
	RasterResetWorldScroll(raster);
	RasterSetWorldScroll(raster, x0, wy0, x1, wy1);
    }
    retval = 0;


cleanup:
    if(win_list_argv) Tcl_Free((char *)win_list_argv);
    return retval;
}

int add_seq_to_raster(RasterResult *result,
		      int seq_id,
		      int seq_num,
		      int direction,
		      int line_width,
		      void (*func)(int seq_num, void *fdata, 
				   seq_reg_data *jdata))
{
    int i;
    static int max_seq = MAX_NUM_SEQ;

    /* allocate a new slot */
    if (result->num_seq_id >= max_seq) {
	max_seq *= 2;
	if (NULL == (result->seq = (seq_id_dir *)realloc(result->seq, 
							 max_seq * 
							 sizeof(seq_id_dir))))
	    return -1;
	if (NULL == (result->cursor = 
		     (cursor_t **)realloc(result->cursor, 
					  max_seq * sizeof(cursor_t *))))
	    return -1;

    }

    i = result->num_seq_id;
    /* add new seq_id */
    result->seq[i].seq_id = seq_id;
    result->seq[i].direction = direction;

    /* need to register with seq_id */
    result->cursor[i] = 
	create_cursor(seq_num, 0, NULL, line_width, 1, direction);
    result->cursor_array[result->cursor[i]->id].env = -1;
    result->cursor_array[result->cursor[i]->id].prev_pos = -1;

    seq_register(seq_num, func, (void *)result, SEQ_RASTER, result->id);
    
    /* increment num of seq_id */
    result->num_seq_id++;
    return 0;
}

/*
 * removing plots may leave plots which should not have a vertical ruler
 */
void RemoveVRuler(Tcl_Interp *interp,
		  char *raster_win,
		  int id)
{
    char *r_win;
    char cmd[1024];

#ifdef DEBUG
    printf("RemoveVRuler %s\n", raster_win);
#endif

    Tcl_VarEval(interp, "winfo parent ", raster_win, NULL);
    r_win = (Tcl_GetStringResult(interp));
    
    sprintf(cmd, "%s.ruler_v%d delete all", r_win, id); 
    Tcl_Eval(interp, cmd);
}

void UpdateVRuler(Tcl_Interp *interp,
		  char *raster_win,
		  double wy0,
		  double wy1)
{
    char cmd[1024];
    
#ifdef DEBUG
    printf("UpdateVRuler %s\n", raster_win);
#endif

    sprintf(cmd, "rasterVRuler %s %f %f", raster_win, wy0, wy1);
    if (TCL_OK != Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "UpdateVRuler", "%s \n", Tcl_GetStringResult(interp));
}

/*
 * reset the world (zoom) parameters of results if the scroll region has 
 * changed such that the world is larger than the scroll
 */
int ReSetRasterWindowWorld(Tcl_Interp *interp,
			   char *raster_old,
			   double orig_y1,
			   int raster_type)
{
    int i;
    double wx0, wy0, wx1, wy1;
    double cx0, cy0, cx1, cy1;
    double sx0, sy0, sx1, sy1;
    int win_list_argc;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    double diff_x, diff_y;
    int zoom;
    char cmd[1024];
    int retval = -1;
    char *parent_old = NULL;
    char **id_list_argv = NULL;
    char **win_list_argv = NULL;


#ifdef DEBUG
    printf("ReSetRasterWindowWorld %s\n", raster_old);
#endif

    if (Tcl_VarEval(interp, "winfo parent ", raster_old, NULL) != TCL_OK)
	/* NB - parent may not exist if we've removed the last result */
	return 0;

    parent_old = strdup(Tcl_GetStringResult(interp));

    win_list_argv = GetRasterWindowList(interp, raster_old, &win_list_argc);
    id_list_argv = GetRasterIdList(interp, raster_old, &win_list_argc);
  

    for (i = 0; i < win_list_argc; i++) {

	if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
	    goto cleanup;
	raster = (Tk_Raster*)info.clientData;

	RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);
	GetRasterCoords(raster, &cx0, &cy0, &cx1, &cy1);
	GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);

	diff_x = wx1 - wx0;
	diff_y = wy1 - wy0;
	
#ifdef DEBUG
	printf("reset world SCROLL %s %f %f %f %f \n", win_list_argv[i], sx0, sy0, sx1, sy1);
	printf("reset world COORDS %s %f %f %f %f \n", win_list_argv[i], wx0, wy0, wx1, wy1);
#endif	

	/* 
	 * if world is larger than scroll: shift world so that it falls within
	 * the current zoom but maintain at current zoom level
	 */
	/* for x:
	 *
	 *    wx0   wx1 wx0 wx1 wx0     wx1
	 *     |--A--|   |-B-|   |---C---|
	 *  |----------|------|---------------|
	 *           sx0     sx1
	 *
	 * case A: wx0 is less than scroll region therefore, move A to sx0, ie
	 * set wx0 to sx0 and wx1 to the difference 
	 *
	 *             wx0   wx1
   	 *             |--A--|
	 *  |----------|------|---------------|
	 *           sx0     sx1
	 *
	 * case B: contained with scroll region therefore OK
	 *
	 * case C: move C so that wx1 is at sx1 and wx0 to the difference
	 *    wx0   wx1 wx0 wx1 wx0     wx1
	 *            |---C---|
	 *  |----------|------|---------------|
	 *           sx0     sx1
	 *
	 * in this case wx0 is still too big therefore set it to be sx0
	 */
	    
	/* for y:
	 *              firstly, need to convert sy0 & sy1 onto world scale
	 *    -                   case A: wy0 < orig_y - sy1, so move A so
	 *    |                           wy0 is at sy0 and wy1 to the
	 *    | - wy0                     difference
	 *    | |
	 *    | A                 case B: convert wy0 & wy1 onto world scale of
	 *    | |                         new plot 
	 *    | - wy1
	 * sy0-  (orig_y - sy1)   case C: wy0 > orig_y -sy0 so move C so 
	 *    | - wy0                     wy1 is at sy1 and wy0 to the  
	 *    | |                         difference
	 *    | B
	 *    | |
	 *    | - wy1
	 * sy1-  (orig_y - sy0)
	 *    |
	 *    | - wy0
	 *    | |
	 *    | |
	 *    | |
	 *    | C
	 *    | |
	 *    | |
	 *    | |
	 *    | - wy1
	 *    - orig_y1 (on raster before removal of result)
	 * 
	 */

	/* case A */
	if (wx0 > sx1) {
	    wx0 = sx1 - diff_x;
	    wx1 = sx1;
	}
	/* case B */
	if (wx0 < sx0) {
	    wx0 = sx0;
	    wx1 = sx0 + diff_x;
	}

	/* case B */
	if ((wy0 > orig_y1 - sy1) && (wy0 < orig_y1 - sy0) &&
	    (wy1 > orig_y1 - sy1) && (wy1 < orig_y1 - sy0) ) {
	    wy0 = sy1 - (orig_y1 - wy0);
	    wy1 = sy1 - (orig_y1 - wy1);
	} else {
	    /* case A */
	    if (wy0 < orig_y1 - sy1) {
		wy0 = sy0;
		wy1 = sy0 + diff_y;
	    }
	    /* case C */
	    if (wy0 > orig_y1 - sy0) {
		wy0 = sy1 - diff_y;
		wy1 = sy1;
	    }
	}

	/* 
	 * if the current world still is larger than the scroll, then
	 * set world so it fits within the scroll region
	 */
	if (wx1 > sx1) wx1 = sx1;
	if (wx0 < sx0) wx0 = sx0;
	if (wx1 < sx0) wx1 = sx1;
	if (wx0 > sx1) wx0 = sx0;

	if (wy1 > sy1) wy1 = sy1;
	if (wy0 < sy0) wy0 = sy0;
	if (wy1 < sy0) wy1 = sy1;
	if (wy0 > sy1) wy0 = sy0;

	/* HACK - to deal with removing a 2d plot and leaving only 1d behind */
	zoom = GetRasterZoom(atoi(id_list_argv[i]));

	if (!zoom) {
	    wy0 = sy0;
	    wy1 = sy1;
	} 

#ifdef DEBUG
	printf("setraster coords %f %f %f %f\n", wx0, wy0, wx1, wy1);
#endif
	SetRasterCoords(raster, wx0, wy0, wx1, wy1);
    }

    if ((raster_type == SEQ_DOT || raster_type == SEQ_E_DOT) && 
	win_list_argc > 0) {
	/* SIP */
	sprintf(cmd, "change_zoom_list %s {%d %d %d %d} {%d %d %d %d}", 
		parent_old, (int)cx0, (int)cy0, (int)cx1, (int)cy1, 
		(int)wx0, (int)wy0, (int)wx1, (int)wy1);
	if (TCL_ERROR == (Tcl_Eval(interp, cmd))) {
	    printf("ReSetRasterWindowWorld %s\n", Tcl_GetStringResult(interp));
	}
    }
    retval = 0;


cleanup:
    if(win_list_argv) Tcl_Free((char *)win_list_argv);
    if(id_list_argv)  Tcl_Free((char *)id_list_argv);
    if(parent_old)    free(parent_old);
    return retval;
}

/*
 * return 0 if the raster contains results with configure.zoom of 0 or 1
 * return 1 if the raster contains any results with configure.zoom of 2
 */
int GetRasterZoom(int raster_id)
{
    seq_reg_generic gen;
    int result;

    gen.job = SEQ_GENERIC;
    gen.result = NULL;
    gen.task = TASK_RASTER_ZOOM;
    seq_result_notify(raster_id, (seq_reg_data *)&gen, 0);
    
    result = (int)gen.result;
    return result;
}

cursor_t *find_raster_result_cursor(RasterResult *result,
				    int seq_id,
				    int direction)
{
    int i;

    for (i = 0; i < result->num_seq_id; i++) {
	if (result->seq[i].seq_id == seq_id && 
	    result->seq[i].direction == direction) {
	    return result->cursor[i];
	}
    }
    return NULL;

}

/* return the number of results in a window  */
int GetWindowNumResults(Tcl_Interp *interp,
			char *raster_win)
{
    int i;
    int id_list_argc;
    char **id_list_argv;
    int cnt = 0;
    RasterResult *raster_result;

    id_list_argv = GetRasterIdList(interp, raster_win, &id_list_argc);

    for (i = 0; i < id_list_argc; i++) {

	if (NULL ==(raster_result = raster_id_to_result(atoi(id_list_argv[i]))))
	    continue;
    
	cnt += raster_result->num_results;
    }
    ckfree((char *)id_list_argv);  /* 7/1/99 johnt - use ckfree as allocated with ckalloc */
    return cnt;
}


void delete_seq_from_raster(int seq_id, 
			    int seq_num,
			    RasterResult *result,
			    void (*func)(int seq_num, void *fdata, 
					 seq_reg_data *jdata))
{
    int i, index = 0;
    
    seq_deregister(seq_num, func, (RasterResult *)result);
    for (i = 0; i < result->num_seq_id; i++) {
	if (seq_id == result->seq[i].seq_id) {

	    result->cursor_array[result->cursor[i]->id].env = -2;

#ifdef DEBUG
	    printf("VIS5 id %d dir 0 vis 0\n", result->cursor[i]->id);
	    printf("VIS6 id %d dir 1 vis 0\n", result->cursor[i]->id);
#endif
	    result->cursor_array[result->cursor[i]->id].visible[HORIZONTAL] = 0;
	    result->cursor_array[result->cursor[i]->id].visible[VERTICAL] = 0;
	    result->cursor_array[result->cursor[i]->id].prev_pos = -1;

	    delete_cursor(seq_num, result->cursor[i]->id, 0);
	    index = i;
	    break;
	}
    }
    /* only move if not the last sequence in array */
    if (index < result->num_seq_id - 1) {
	memmove(&result->seq[index], 
		&result->seq[index + 1], 
		(result->num_seq_id - index - 1) * sizeof(seq_id_dir));
	memmove(&result->cursor[index], 
		&result->cursor[index + 1], 
		(result->num_seq_id - index - 1) * sizeof(cursor_t*));
    }
    result->num_seq_id--;
}

/*
 * Deletes an element of the raster seq[i] array. This is called when we have
 * the same sequence in both directions (otherwise we call delete_seq_from_
 * raster).
 */
static void delete_seq_from_direction(int r_seq_id,
				      RasterResult *result)
{
    result->cursor_array[result->cursor[r_seq_id]->id].env = -2;

    result->cursor_array[result->cursor[r_seq_id]->id].visible[HORIZONTAL] = 0;
    result->cursor_array[result->cursor[r_seq_id]->id].visible[VERTICAL] = 0;
    result->cursor_array[result->cursor[r_seq_id]->id].prev_pos = -1;

    delete_cursor(GetSeqNum(result->seq[r_seq_id].seq_id),
		  result->cursor[r_seq_id]->id, 0);

    /* only move if not the last sequence in array */
    if (r_seq_id < result->num_seq_id - 1) {
	memmove(&result->seq[r_seq_id], 
		&result->seq[r_seq_id + 1], 
		(result->num_seq_id - r_seq_id - 1) * sizeof(seq_id_dir));
	memmove(&result->cursor[r_seq_id], 
		&result->cursor[r_seq_id + 1], 
		(result->num_seq_id - r_seq_id - 1) * sizeof(cursor_t*));
    }
    result->num_seq_id--;
}

static void remove_redundant_sequences(RasterResult *raster, 
				       seq_result *result)
{
    int horiz_seq_id = result->seq_id[0];
    int vert_seq_id = result->seq_id[1];
    int horiz_r_id = -1, vert_r_id = -1;
    int horiz_r_count, vert_r_count;
    int i, num_funcs, num_elements;
    seq_result **data, *tmp_result;
    out_raster *output;
    int del1, del2;

    num_elements = seq_num_results();
    if (num_elements == 0)
	return;

    /*
     * Find elements in raster->seq[] that correspond to the sequences used
     * by this seq_result.
     */
    for (i = 0; i < raster->num_seq_id; i++) {
	if (raster->seq[i].seq_id == horiz_seq_id)
	    horiz_r_id = i;
	if (raster->seq[i].seq_id == vert_seq_id)
	    vert_r_id = i;
    }

    /*
     * Find how many cases horiz_r_id and vert_r_id are used within this
     * raster, as well as horiz_r_seq and vert_r_seq.
     */
    data = (seq_result **)xmalloc(num_elements * sizeof(seq_result *));
    search_reg_data(comparison2, (void **)data, &num_funcs);
    horiz_r_count = vert_r_count = 0;

    for (i = 0; i < num_funcs; i++) {
	tmp_result = (seq_result *)data[i];
	output = tmp_result->output;
	if (strcmp(output->raster_win, raster->raster_win) == 0) {
	    if (tmp_result->seq_id[HORIZONTAL] == horiz_seq_id)
		horiz_r_count++;
	    if (tmp_result->seq_id[VERTICAL] == vert_seq_id)
		vert_r_count++;
	}
    }

    /*
     * If horiz_r_count or vert_r_count are 0 then we can now remove
     * this 'raster sequence' from the raster. (The sip result will have
     * already been removed - hence 0 instead of 1.)
     */
    del1 = del2 = -1;
    if (horiz_r_count == 0)
	del1 = horiz_r_id;
    if (vert_r_count == 0) {
	if (del1 == -1)
	    del1 = vert_r_id;
	else {
	    if (vert_r_id > del1) {
		del2 = del1;
		del1 = vert_r_id;
	    } else {
		del2 = vert_r_id;
	    }
	}
    }

    if (del1 != -1)
	delete_seq_from_direction(del1, raster);
    if (del2 != -1)
	delete_seq_from_direction(del2, raster);

    xfree(data);
}


/*
 * return an array of raster names contained in a window
 */
char **GetRasterWindowList(Tcl_Interp *interp,
			 char *raster_win,
			 int *num_windows)
{
    char *parent;
    int win_list_argc;
    char **win_list_argv;
    
    Tcl_VarEval(interp, "GetRasterParent ", raster_win, NULL);
    parent = strdup(Tcl_GetStringResult(interp));

    if (TCL_ERROR == (Tcl_VarEval(interp, "GetRasterWinList ", parent, NULL))) {
	printf("GetRasterWindowList: %s\n", Tcl_GetStringResult(interp));
	free(parent);
	return NULL;
    }
    if (Tcl_SplitList(interp, Tcl_GetStringResult(interp),
		      &win_list_argc, &win_list_argv) != TCL_OK) {
	free(parent);
	return NULL;
    }

    *num_windows = win_list_argc;
    free(parent);
    return win_list_argv;
}

/*
 * return an array of raster ids contained in a window
 */
char **GetRasterIdList(Tcl_Interp *interp,
		     char *raster_win,
		     int *num_windows)
{
    char *parent;
    int id_list_argc;
    char **id_list_argv;
    
    Tcl_VarEval(interp, "GetRasterParent ", raster_win, NULL);
    parent = strdup(Tcl_GetStringResult(interp));

    if (TCL_ERROR == (Tcl_VarEval(interp, "GetRasterIdList ", parent, NULL))) {
	printf("GetRasterIdList: %s\n", Tcl_GetStringResult(interp));
	free(parent);
	return NULL;
    }
    if (Tcl_SplitList(interp, Tcl_GetStringResult(interp),
		      &id_list_argc, &id_list_argv) != TCL_OK) {
	free(parent);
	return NULL;
    }

    *num_windows = id_list_argc;
    free(parent);
    return id_list_argv;
}


void AddResultToRaster(RasterResult *result)
{
    result->num_results++;
}

int DeleteResultFromRaster(RasterResult *result)
{
    seq_reg_delete del;
   
    result->num_results--;
    if (result->num_results == 0) {

	del.job = SEQ_QUIT;
	seq_result_notify(result->id, (seq_reg_data *)&del, 0);
	return 0;
    }
    return result->num_results;
}


/*
 * convert world y into raster y
 */
double rasterY(Tk_Raster *raster, 
	       double y)
{
    double x0, y0, x1, y1;

    RasterGetWorldScroll(raster, &x0, &y0, &x1, &y1);
    return (y1 - y + y0);
}

void FindRasterResultY0(Tk_Raster *raster,
			int raster_id,
			config *configure,
			int num_results,
			double *y0,
			double *tick_height)
{
    double sx0, sy0, sx1, sy1;
    double wx0, wy0, wx1, wy1;
    double tick_ht;
    double tx, ty;
    double pos;

    /* world scroll limits */
    RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);

    /* world limits */
    GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);
#ifdef DEBUG
    printf("RASTER WORLD %f %f %f %f\n", wx0, wy0, wx1, wy1);
#endif
    /* height of plot */
    if (configure->height <= 1.0) {
	tick_ht = (wy1 - wy0) * configure->height;
    } else {
	/* find tick height in world coords */
	RasterToWorld(raster, 0, 0, &tx, &ty);
	RasterToWorld(raster, 0, configure->height, &tx, &tick_ht);
	tick_ht = tick_ht - ty;
    }

    if (configure->scroll == 0) {
	/* always to remain in the same place */
	if (configure->y_direction == '-') {
	    pos = ((wy1 - wy0) * configure->position) + wy0;
	} else {
	    pos = wy1 - ((wy1 - wy0) * configure->position);
	}
    } else {
	/* moves as move y scroll */
	if (configure->y_direction == '-') {
	    pos = ((sy1 - sy0) * configure->position) + sy0;
	} else {
	    pos = sy1 - ((sy1 - sy0) * configure->position);
	}
    }
#ifdef DEBUG
    printf("sy0 %f sy1 %f position %f \n", sy0, sy1, configure->position);
    printf("pos %f \n", pos);
#endif
    /* position of plot */
    if (configure->zoom == 1 && num_results == 1) {
	/* 
	 * HACK to deal with tick plot on its own to ensure it starts at the
	 * base line and the crosshair coords are correct
	 */
	if (configure->y_direction == '-') {
	    pos = wy0;
	} else {
	    pos = wy1;
	}
    } else {
	/* y0 = sy1 - pos + sy0; */
    }
    
    *y0 = sy1 - pos + sy0;
    *tick_height = tick_ht;

#ifdef DEBUG
    printf("pos %f tick_ht %f \n", pos, tick_ht);
    printf("y0 %f\n", *y0);
#endif
}

int SeqReSetRasterWindowSize(Tcl_Interp *interp,
			     char *raster_win,
			     int raster_type)
{
    int i;
    int win_list_argc;
    double x0, y0, x1, y1;
    d_line *max_size;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    int retval = TCL_ERROR;
    char **win_list_argv = NULL;
    char **id_list_argv = NULL;



#ifdef DEBUG
    printf("START  SeqReSetRasterWindowSize %s\n", raster_win);
#endif

    win_list_argv = GetRasterWindowList(interp, raster_win, &win_list_argc);
    id_list_argv = GetRasterIdList(interp, raster_win, &win_list_argc);

    x0 = DBL_MAX/2;
    y0 = DBL_MAX/2;
    x1 = -DBL_MAX/2;
    y1 = -DBL_MAX/2;

    for (i = 0; i < win_list_argc; i++) {
	FindRasterSize(atoi(id_list_argv[i]), &max_size);

	if (max_size->x0 < x0) x0 = max_size->x0;
	if (max_size->y0 < y0) y0 = max_size->y0;
	if (max_size->x1 > x1) x1 = max_size->x1;
	if (max_size->y1 > y1) y1 = max_size->y1;

	xfree(max_size);
    }
#ifdef DEBUG
    printf("FINAL %f %f %f %f \n", x0, y0, x1, y1);
#endif

    for (i = 0; i < win_list_argc; i++) {
	FindRasterSize(atoi(id_list_argv[i]), &max_size);
#ifdef DEBUG
	printf("max size %s %f %f %f %f\n", win_list_argv[i], max_size->x0, 
	       max_size->y0, max_size->x1,max_size->y1);
#endif
	if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
	    goto cleanup;
	raster = (Tk_Raster*)info.clientData;
	RasterResetWorldScroll(raster);

    if (raster_type == SEQ_DOT || raster_type == SEQ_E_DOT) {
	RasterSetWorldScroll(raster, x0, y0, x1, y1);
    } else {
	RasterSetWorldScroll(raster, x0, max_size->y0, x1, max_size->y1);
	/* 
	 * added 2/3/98 to reset world values to be the same as the scroll ie
	 * lose all previous zoom info.
	 */
	SetRasterCoords(raster, x0, max_size->y0, x1, max_size->y1);
    }

#ifdef DEBUG
	{
	    double tx0, ty0, tx1, ty1;
	printf("RESET SCROLL %s %f %f %f %f\n", win_list_argv[i], x0, max_size->y0, x1, max_size->y1);
	GetRasterCoords(raster, &tx0, &ty0, &tx1, &ty1);
	printf("RESET COORDS %s %f %f %f %f\n", win_list_argv[i], tx0, ty0, tx1, ty1);
	}
#endif
	/* HACK goes wrong when add base comp to gene search */
	/* but need it when remove gene search from base comp */
	/* SetRasterCoords(raster, x0, y0, x1, y1); */
	xfree(max_size);
    }
    retval = TCL_OK;


cleanup:
    if(win_list_argv) Tcl_Free((char *)win_list_argv);
    if(id_list_argv)  Tcl_Free((char *)id_list_argv);
    return retval;
}

int comparison3(void *v_result,
		int type)
{
    if (type == SEQ_PLOT_PERM) {
	return 1;
    }
    return 0;
}

int comparison_raster(void *v_result,
		      int type)
{
    if (type == SEQ_RASTER) {
	return 1;
    }
    return 0;
}

/*
 * move a result from one raster to another
 * scale the plot being added so that it is at the same x and y scale as
 * the host raster
 */
void SeqUpdateResultWindow(Tcl_Interp *interp,
			   char *raster_old,
			   char *raster_new,
			   int old_id,
			   int new_id,
			   int result_id,
			   int job)
{
    int num_elements;
    seq_result *result;
    out_raster *output;
    char *opts[5];
    Tcl_CmdInfo cmd_info;
    Tk_Raster *rasternew, *rasterold;
    double wx0, wy0, wx1, wy1;
    double o_wx0, o_wy0, o_wx1, o_wy1;
    d_line *max_size;
    double new_xscroll, old_xscroll;
    double new_yscroll, old_yscroll;
    char *parent_old, *parent_new;
    RasterResult *new_result;
    RasterResult *old_result;
    int i, j, index = 0;
    char cmd[1024];
    double orig_y1;
    int line_width;
    d_line *dim;
    seq_reg_info info;
    int old_zoom, new_zoom;

#ifdef DEBUG
    printf("SeqUpdateResultWindow\n");
#endif

    num_elements = seq_num_results();
    if (num_elements == 0)
	return;

    /* get Raster structures */
    if (Tcl_GetCommandInfo(interp, raster_old, &cmd_info) == 0) {
	return;
    }
    rasterold = (Tk_Raster*)cmd_info.clientData;

    if (Tcl_GetCommandInfo(interp, raster_new, &cmd_info) == 0) 
	return;
    rasternew = (Tk_Raster*)cmd_info.clientData;
    
    /* get raster result structures */
    old_result = raster_id_to_result(old_id);
    new_result = raster_id_to_result(new_id);

    /* get seq result structure */
    result = seq_id_to_result(result_id);

    /* increment num_results in new_result */
    AddResultToRaster(new_result);

    /* need to do this before start changing the raster_win of old plot */
    old_zoom = GetRasterZoom(old_id);

    output = result->output;
    strcpy(output->raster_win, raster_new);

    if (NULL == (opts[1] = (char *)xmalloc((strlen(GetRasterColour(interp, 
								   rasterold, 
								   output->env_index))+1) 
					   * sizeof(char))))
	return;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return;
    
    opts[0] = "-fg";
    strcpy(opts[1], GetRasterColour(interp, rasterold, output->env_index));
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", GetRasterLineWidth(interp, rasterold, 
					      output->env_index));
    opts[4] = NULL;

    output->env_index = CreateDrawEnviron(interp, rasternew, 4, opts);
    xfree(opts[1]);
    xfree(opts[3]);
    RasterInitPlotFunc(rasternew, SeqRasterPlotFunc);  

    /* 
     * get the current scroll region for the new and old rasters
     */
    RasterGetWorldScroll(rasterold, &o_wx0, &o_wy0, &o_wx1, &o_wy1);
    /* need to correct the scale */
    o_wy0 = (o_wy0 - output->sf_c) / output->sf_m;
    o_wy1 = (o_wy1 - output->sf_c) / output->sf_m;
    RasterGetWorldScroll(rasternew, &wx0, &wy0, &wx1, &wy1);
    
    if(wx0 >= (DBL_MAX / 2)) wx0 = o_wx0;
    if(wy0 >= (DBL_MAX / 2)) wy0 = o_wy0;
    if(wx1 <= (-DBL_MAX / 2)) wx1 = o_wx1;
    if(wy1 <= (-DBL_MAX / 2)) wy1 = o_wy1;
    RasterSetWorldScroll(rasternew, wx0, wy0, wx1, wy1);

#ifdef DEBUG
    printf("OLD %s %f %f %f %f\n", raster_old, o_wx0, o_wy0, o_wx1, o_wy1);
    printf("NEW %s %f %f %f %f\n", raster_new, wx0, wy0, wx1, wy1);
#endif

    new_xscroll = wx1 - wx0;
    old_xscroll = o_wx1 - o_wx0;
    new_yscroll = wy1 - wy0;
    old_yscroll = o_wy1 - o_wy0;
    orig_y1 = o_wy1;

    if (job == ADD) {
	/* add new raster to a window already containing rasters */

	/* reset scale factors */
	output->sf_m = 1.0;
	output->sf_c = 0.0;
	FindRasterSize(new_id, &max_size);
	RasterSetWorldScroll(rasternew, max_size->x0, max_size->y0, 
			     max_size->x1, max_size->y1);
	
	SeqAddRasterToWindow(interp, raster_new, result->graph);
        UpdateVRuler(interp, raster_new, max_size->y0, max_size->y1);
	xfree(max_size);
    } else if (job == NEW) {
	/* add new raster to a new window */

	/* reset scale factors */
	output->sf_m = 1.0;
	output->sf_c = 0.0;
	
	/* find size of result */
	FindRasterSize(new_id, &max_size);

	/* set world coords */
	SetRasterCoords(rasternew, max_size->x0, max_size->y0, max_size->x1, 
			max_size->y1);

	 RasterSetWorldScroll(rasternew, max_size->x0, max_size->y0, max_size->x1, max_size->y1); 
	xfree(max_size);
    } else {
	double p2, q2;
	double m, c;

	/* find dimensions of result */
	info.job = SEQ_RESULT_INFO;
	info.op = DIMENSIONS;
	info.result = NULL;
	seq_result_notify(result_id, (seq_reg_data *)&info, 0);
	if (!info.result) {
	    return;
	}
	dim = (d_line *)info.result;

	/* superimpose, update an exising raster */
	p2 = (wy1 - wy0) * (dim->y0 - o_wy0) / (o_wy1 - o_wy0) + wy0;
	q2 = (wy1 - wy0) * (dim->y1 - o_wy0) / (o_wy1 - o_wy0) + wy0;

	m = (p2 - q2) / (dim->y0 - dim->y1);
	c = p2 - (m * dim->y0);
	
	/* 
	 * kfs 19.2.98, I think I need to reset sf_m and sf_c and things seem
	 * to work. This was to fix the case of superimposing frame 1 of
	 * base pref onto base comp and then moving frame 1 of base pref
	 * onto frame 2 of base pref. Previously sf_m and sf_c would be their
	 * previous value whereas they should end up as 1.0 and 0.0
	 */
#ifdef DEBUG
	printf("result_id %d\n", result_id);
	printf("before sf_c %f sf_m %f\n", result->sf_c, result->sf_m);
#endif

	output->sf_c = 0.0;
	output->sf_m = 1.0;

	output->sf_c = (m * output->sf_c) + (c);
	output->sf_m *= m;


	/* find size of result */
	FindRasterSize(new_id, &max_size);

	/* enlarge x scroll region if necessary */
	RasterSetWorldScroll(rasternew, max_size->x0, max_size->y0, 
			     max_size->x1, max_size->y1);

	xfree(max_size);

	/* add old cursor to new raster */
	update_raster_cursor(new_id, old_id);

	line_width = get_default_int(interp, spin_defs, 
				     w("SEQ.CURSOR.LINE_WIDTH"));

	for (i = 0; i < old_result->num_seq_id; i++) {
	    for (j = 0; j < new_result->num_seq_id; j++) {
		if ((old_result->seq[i].seq_id != new_result->seq[j].seq_id) ||
		    (old_result->seq[i].direction != 
		    new_result->seq[j].direction)) {
		    index = i;
		} else {
		    index = -1;
		    break;
		}
	    }
	    if (index > -1) {
		add_seq_to_raster(new_result, old_result->seq[index].seq_id,
				  GetSeqNum(old_result->seq[index].seq_id),
				  old_result->seq[index].direction, 
				  line_width, seq_raster_callback);
	    }
	}
    }  

    /* 
     * need to replot all the rasters of the new raster window in case the
     * scroll region has changed 
     */
    Tcl_VarEval(interp, "winfo parent ", raster_new, NULL);
    parent_new = strdup(Tcl_GetStringResult(interp));
    Tcl_VarEval(interp, "winfo parent ", raster_old, NULL);
    parent_old = strdup(Tcl_GetStringResult(interp));

    /* only need to do this if rasters are in different windows */
    if (strcmp(parent_new, parent_old) != 0) { 
	ReplotAllRasterWindow(interp, raster_new);
    }

    /* update the old raster plot and remove the old result */

    /* want to remove result from raster */
    old_result->num_results--;

    /* deregister raster_old for seq_num of result if necessary */
    SeqDeregisterRasterWindow(result->seq_id[HORIZONTAL], old_result, raster_old);

    if (result->seq_id[VERTICAL] > -1) {
	SeqDeregisterRasterWindow(result->seq_id[VERTICAL], old_result, raster_old);
    }

    /* Remove 'raster sequences' if appropriate */
    if (old_result->num_results > 0)
	remove_redundant_sequences(old_result, result);

    /* 
     * if removed a graph plot from a tick plot, need to reset the pixel 
     * size of the raster plot
     */
    /* 
     * removed (30/9/97) because moving gene search A onto B and then back
     * again, goes wrong if the rasters have been shrunk
     */
    /* ReSetRasterSize(interp, old_id, raster_old); */

    /*
     * need this if removing a larger scroll region so that the remaining
     * scroll regions in the raster are rescaled 
     */
    SeqReSetRasterWindowSize(interp, raster_old, result->graph);
    /* 
     * also need to reset the world for each of the rasters in the old 
     * raster window if necessary ie if the the scroll region has changed such
     * that the world is outside this
     */
    ReSetRasterWindowWorld(interp, raster_old, orig_y1, result->graph);
    ReplotAllRasterWindow(interp, raster_old);

    RasterGetWorldScroll(rasterold, &o_wx0, &o_wy0, &o_wx1, &o_wy1);

    if (GetRasterZoom(old_id) <= 0) {
	RemoveVRuler(interp, raster_old, old_id);
    }

    if (result->graph == SEQ_DOT || result->graph == SEQ_E_DOT) {
	if (job == SUPER) {
	    UpdateZoomList(interp, o_wx0, o_wy0, o_wx1, o_wy1, 
			   wx0, wy0, wx1, wy1, parent_new, 0);
	}

	RasterGetWorldScroll(rasterold, &o_wx0, &o_wy0, &o_wx1, &o_wy1);
	new_xscroll = o_wx1 - o_wx0;
	new_yscroll = o_wy1 - o_wy0;
	
#ifdef DEBUG
	printf("new_xscroll %f old_xscroll %f new_yscroll %f old_yscroll %f\n",
	       new_xscroll, old_xscroll, new_yscroll, old_yscroll);
#endif

	if (new_xscroll < old_xscroll || new_yscroll < old_yscroll) {
	    /* 
	     * if window has still has results, then remove current extents if the
	     * zoom range has decreased, ie a plot of larger extents has been 
	     * removed 
	     */
	    if (GetWindowNumResults(interp, raster_old) > 0) {
		sprintf(cmd, "update_zoom_list %s %d {%d %d %d %d}", 
			parent_old, job, (int)o_wx0, (int)o_wy0, (int)o_wx1, (int)o_wy1);
		
		if (TCL_ERROR == (Tcl_Eval(interp, cmd))) {
		    printf("UpdateZoomList %s\n", Tcl_GetStringResult(interp));
		}
		
	    }
	}
	
    } else {
	new_zoom = GetRasterZoom(new_id);
	
	if (job == SUPER || job == ADD) {
	    /* update scale bars of new raster - only do x */
	    if (old_zoom >= 0 && new_zoom >= 0) {
		UpdateScaleBars(interp, old_xscroll, new_xscroll,
				new_yscroll, new_yscroll, raster_new, old_id,
				new_id, 0); 
	    }
	}
	/* reset scale bars of old raster */
	if (old_zoom >= 0 && new_zoom >= 0) {
	    UpdateScaleBars(interp, 1.0, 1.0, 1.0, 1.0, raster_old, old_id, 
			    new_id, 1);
	}
	

    }

    /* 
     * I know this is a bit of a hack, but I have to delete the raster last
     * even though I decremented result->num_results earlier
     */
    if (old_result->num_results == 0) {
	seq_reg_delete del;
	del.job = SEQ_QUIT;
	seq_result_notify(old_result->id, (seq_reg_data *)&del, 0);
    }

    free(parent_new);
    free(parent_old);
}

/*
 * move a raster to another
 */
void SeqUpdateRasterWindow(Tcl_Interp *interp,
			   char *raster_old,
			   char *raster_new,
			   int new_id,
			   int old_id,
			   int job)
{
    int num_elements;
    seq_result **data;
    seq_result *result = NULL;
    int num_funcs;
    int i;
    out_raster *output;
    char *opts[5];
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    Tk_Raster *rasterold;
    d_line *max_size = NULL;
    double new_xscroll, old_xscroll;
    double new_yscroll, old_yscroll;
    double o_wy0, o_wy1;
    double n_wy0, n_wy1;
    int line_width;
    int num_results = 0;
    RasterResult *new_result;
    RasterResult *old_result;
    seq_reg_info r_info;
    char *parent_old, *parent_new;

#ifdef DEBUG
    printf("SeqUpdateRasterWindow\n");
#endif

    num_elements = seq_num_results();
    if (num_elements == 0)
	return;

    data = (seq_result **)xmalloc(num_elements * sizeof(seq_result *));
 	
    search_reg_data(comparison2, (void **)data, &num_funcs);

    opts[0] = "-fg";
    opts[2] = "-linewidth";
    opts[4] = NULL;

    if (Tcl_GetCommandInfo(interp, raster_new, &info) == 0) 
	return;
    raster = (Tk_Raster*)info.clientData;
    
    if (Tcl_GetCommandInfo(interp, raster_old, &info) == 0) {
	return;
    }
    rasterold = (Tk_Raster*)info.clientData;

    FindRasterSize(old_id, &max_size);
    old_xscroll = max_size->x1 - max_size->x0;
    old_yscroll = max_size->y1 - max_size->y0;
    o_wy0 = max_size->y0;
    o_wy1 = max_size->y1;
    xfree(max_size);

    FindRasterSize(new_id, &max_size);
    new_xscroll = max_size->x1 - max_size->x0;
    new_yscroll = max_size->y1 - max_size->y0;
    n_wy0 = max_size->y0;
    n_wy1 = max_size->y1;
    xfree(max_size);

#ifdef DEBUG
    printf("OLD %f %f NEW %f %f\n", o_wy0, o_wy1, n_wy0, n_wy1);
    printf("OLD %s NEW %s \n", raster_old, raster_new);
#endif

    old_result = raster_id_to_result(old_id);
    new_result = raster_id_to_result(new_id);

    /* update the old raster info with new raster info */
    max_size = NULL;
    for (i = 0; i < num_funcs; i++) {
	result = data[i];
	output = result->output;

	if (strcmp(output->raster_win, raster_old) == 0) {
	    num_results++;

	    strcpy(output->raster_win, raster_new);

	    if (NULL == (opts[1] = (char *)xmalloc((strlen(GetRasterColour(interp, rasterold, output->env_index))+1) 
						   * sizeof(char))))
		return;
	    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
		return;

	    strcpy(opts[1], GetRasterColour(interp, rasterold, output->env_index));
	    sprintf(opts[3], "%d", GetRasterLineWidth(interp, rasterold, output->env_index));

	    RasterInitPlotFunc(raster, SeqRasterPlotFunc);

	    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);

	    if (job == SUPER) {
		double p2, q2;
		double m, c;		
		d_line *dim;

		/* update an exising raster */

		/* find dimensions of result */
		r_info.job = SEQ_RESULT_INFO;
		r_info.op = DIMENSIONS;
		r_info.result = NULL;
		seq_result_notify(result->id, 
				  (seq_reg_data *)&r_info, 0);
		if (!r_info.result) {
		    return;
		}
		dim = (d_line *)r_info.result;
		/* update an exising raster */

		p2 = (n_wy1 - n_wy0) * (dim->y0 - o_wy0) / 
		    (o_wy1 - o_wy0) + n_wy0;
		q2 = (n_wy1 - n_wy0) * (dim->y1 - o_wy0) / 
		    (o_wy1 - o_wy0) + n_wy0;
		
		if (dim->y0 - dim->y1 != 0) {
		    m = (p2 - q2) / (dim->y0 - dim->y1);
		} else {
		    m = 0;
		}
		c = p2 - (m * dim->y0);

		output->sf_c = (m * output->sf_c) + c;
		output->sf_m *= m;
 	    }
	    /* find the scroll size for all results in this raster */
	    if (max_size)
		xfree(max_size);	    
	    FindRasterSize(new_id, &max_size);
	    RasterSetWorldScroll(raster, max_size->x0, max_size->y0, 
				 max_size->x1, max_size->y1);

	    xfree(opts[1]);
	    xfree(opts[3]);
	}
    }

   /*
    * add required number of results to raster result
    */
    for (i = 0; i < num_results; i++) {
	AddResultToRaster(new_result);
    }
    
    if (job == ADD) {
	/* add raster to a new raster in window */
	SeqAddRasterToWindow(interp, raster_new, result->graph);
	ReplotAllRasterWindow(interp, raster_new);
	
    } else if (job == NEW) {
	/* add raster to new window */
	SetRasterCoords(raster, max_size->x0, max_size->y0, 
			max_size->x1, max_size->y1);
	ReplotAllRasterWindow(interp, raster_new);

    } else {
	/* superimpose */
	int index = -1;
	int j;

	/* add old cursor to new raster */
	update_raster_cursor(new_id, old_id);

	line_width = get_default_int(interp, spin_defs, 
				     w("SEQ.CURSOR.LINE_WIDTH"));
	for (i = 0; i < old_result->num_seq_id; i++) {
	    for (j = 0; j < new_result->num_seq_id; j++) {
		if (old_result->seq[i].seq_id != new_result->seq[j].seq_id ||
		    old_result->seq[i].direction !=
		    new_result->seq[j].direction) {
		    index = i;
		} else {
		    index = -1;
		    break;
		}
	    }
	    if (index > -1) {
		add_seq_to_raster(new_result, old_result->seq[index].seq_id,
				  GetSeqNum(old_result->seq[index].seq_id),
				  old_result->seq[index].direction, 
				  line_width, seq_raster_callback);
	    }
	}
    }

    if (!(result->graph == SEQ_DOT || result->graph == SEQ_E_DOT)) {
	/* 
	 * need to replot all the rasters of the new raster window in case the
	 * scroll region has changed 
	 */
	Tcl_VarEval(interp, "winfo parent ", raster_new, NULL);
	parent_new = strdup(Tcl_GetStringResult(interp));
	Tcl_VarEval(interp, "winfo parent ", raster_old, NULL);
	parent_old = strdup(Tcl_GetStringResult(interp));

	/* only need to do this if rasters are in different windows */
	if (strcmp(parent_new, parent_old) != 0) { 
	    ReplotAllRasterWindow(interp, raster_new);
	}
	free(parent_new);
	free(parent_old);
    }

   /*
    * delete required number of results to raster result
    */
    for (i = 0; i < num_results; i++) {
	DeleteResultFromRaster(old_result);
    }

    /*
     * need this if removing a larger scroll region so that the remaining
     * scroll regions in the raster are rescaled 
     */
    SeqReSetRasterWindowSize(interp, raster_old, result->graph);
    /* 
     * also need to reset the world for each of the rasters in the old 
     * raster window if necessary ie if the the scroll region has changed such
     * that the world is outside this
     */
    ReSetRasterWindowWorld(interp, raster_old, o_wy1, result->graph);
    ReplotAllRasterWindow(interp, raster_old);
    
    if (max_size)
	xfree(max_size);
    xfree(data);
}


/*
 * add a new raster (+ results) to an existing window
 * need to:
 * find the current world value for existing window
 * set the world of raster_new to this
 * set all scroll region of all the rasters in the window to this value 
 */
int SeqAddRasterToWindow(Tcl_Interp *interp,
			 char *raster_win,
			 int raster_type)
{
    int i;
    int win_list_argc;
    Tcl_CmdInfo info;
    Tk_Raster *rasterorig = NULL;
    Tk_Raster *rasternew;
    double n_wx0, n_wy0, n_wx1, n_wy1;
    double o_wx0, o_wy0, o_wx1, o_wy1;
    double o_x0, o_y0, o_x1, o_y1;
    double y0_offset = 0.0, y1_offset = 1.0;
    double new_xscroll = -1, orig_xscroll = -1;
    double new_yscroll = -1, orig_yscroll = -1;
    int raster_id_new, raster_id_orig = -1;
    char *parent;
    int n_rasters = 0;
    char *raster_orig = NULL;
    int orig_zoom = -1, new_zoom;
    int retval = TCL_ERROR;
    char **win_list_argv = NULL;


    /* find all rasters in a window */
    win_list_argv = GetRasterWindowList(interp, raster_win, &win_list_argc);

    /* find the world coords of any of the rasters already in a window */
    for (i = 0; i < win_list_argc; i++) {
	if (strcmp(raster_win, win_list_argv[i]) != 0) {
	    
	    if (Tcl_GetCommandInfo(interp, win_list_argv[i], &info) == 0) 
		goto cleanup;
	    rasterorig = (Tk_Raster*)info.clientData;

	    raster_orig = win_list_argv[i];
	    GetRasterCoords(rasterorig, &o_x0, &o_y0, &o_x1, &o_y1);
	    RasterGetWorldScroll(rasterorig, &o_wx0, &o_wy0, &o_wx1, &o_wy1);
	    orig_xscroll = o_wx1 - o_wx0;
	    orig_yscroll = o_wy1 - o_wy0;
	    y0_offset = (o_y0 - o_wy0) / (o_wy1 - o_wy0);
	    y1_offset = (o_y1 - o_wy0) / (o_wy1 - o_wy0);
	    Tcl_VarEval(interp, "GetRasterId ", win_list_argv[i], NULL);
	    raster_id_orig = atoi(Tcl_GetStringResult(interp));

#ifdef DEBUG
	    printf("o_x0 %f o_x1 %f o_y0 %f o_y1 %f\n", 
		   o_x0, o_x1, o_y0, o_y1);
	    printf("y0_offset %f y1_offset %f y0 %f y1 %f wy0 %f wy1 %f\n", 
		   y0_offset, y1_offset, o_y0, o_y1, o_wy0, o_wy1);
#endif
	    n_rasters++;
	    break;
	}	
    }

    /* set the new raster's world coords ie zoom level */
    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	goto cleanup;
    rasternew = (Tk_Raster*)info.clientData;
    RasterGetWorldScroll(rasternew, &n_wx0, &n_wy0, &n_wx1, &n_wy1);

#ifdef DEBUG
    printf("new %f %f %f %f\n", n_wx0, n_wx1, n_wy0, n_wy1);
#endif

    new_xscroll = n_wx1 - n_wx0;
    new_yscroll = n_wy1 - n_wy0;

    if (win_list_argc == 1) {
	/* only 1 raster in window */
	o_x0 = n_wx0;
	o_x1 = n_wx1;
	o_y0 = n_wy0;
	o_y1 = n_wy1;
    }

    /* HACK: other windows in raster haven't been plotted yet */
    if (o_x0 == 0.0 && o_x1 == 0.0) {
	o_x0 = n_wx0;
	o_x1 = n_wx1;
	o_y0 = n_wy0;
	o_y1 = n_wy1;
	y0_offset = 0.0;
	y1_offset = 1.0;
	new_xscroll = n_wx1 - n_wx0;
	new_yscroll = n_wy1 - n_wy0;
    }
    
    if (raster_type == SEQ_DOT || raster_type == SEQ_E_DOT) {
	SetRasterCoords(rasternew, o_x0, o_y0, o_x1, o_y1);

	if (n_rasters > 0) {
	    /* SIP */
	    /* scale all plots in x & y to the original and update zoom list */
	    Tcl_VarEval(interp, "winfo parent ", raster_orig, NULL);
	    parent = strdup(Tcl_GetStringResult(interp));

	    RasterGetWorldScroll(rasterorig, &o_wx0, &o_wy0, &o_wx1, &o_wy1);
	    SeqReSetRasterWindowSize(interp, raster_orig, raster_type);
	    ReplotAllRasterWindow(interp, raster_orig);
    
	    UpdateZoomList(interp, n_wx0, n_wy0, n_wx1, n_wy1, 
			   o_wx0, o_wy0, o_wx1, o_wy1, parent, 0);
	    free(parent);
	}
    } else {
	o_y0 = (n_wy1 - n_wy0) * y0_offset + n_wy0;
	o_y1 = (n_wy1 - n_wy0) * y1_offset + n_wy0;
	
	SetRasterCoords(rasternew, o_x0, o_y0, o_x1, o_y1);
	
	/* update the world scroll region for all rasters in this window */
	SetRasterWindowSize(interp, raster_win);
	
#ifdef DEBUG
	printf("new_xscroll %f orig_xscroll %f\n", new_xscroll, orig_xscroll);
	printf("new_yscroll %f orig_yscroll %f\n", new_yscroll, orig_yscroll);
#endif
	Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
	raster_id_new = atoi(Tcl_GetStringResult(interp));

	if (raster_id_orig > -1) {
	    orig_zoom = GetRasterZoom(raster_id_orig);
	}
	new_zoom = GetRasterZoom(raster_id_new);
	
	if (n_rasters > 0 && orig_zoom >= 0 && new_zoom >= 0) {
	    UpdateScaleBars(interp, new_xscroll, orig_xscroll,
			    new_yscroll, new_yscroll, raster_win, 
			    raster_id_orig, raster_id_new, 0);
	}
    }
    retval = TCL_OK;


cleanup:
    if(win_list_argv) Tcl_Free((char *)win_list_argv);
    return retval;
}

void seq_raster_callback(int seq_num, void *obj, seq_reg_data *jdata) 
{
    RasterResult *result = (RasterResult *) obj;
    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "raster plot");
	    return;
	}
    case SEQ_DELETE:	
	{
	    int seq_id;
	    int *seq_nums;
	    int num_seqs;

	    seq_nums = result_to_seq_nums(result->id, &num_seqs);
	    xfree(seq_nums);	    

	    /* if only 1 sequence is present then QUIT the display */
	    if (num_seqs > 1) {
		/* delete sequence */
		seq_id = GetSeqId(seq_num);
		delete_seq_from_raster(seq_id, seq_num, result, 
				   seq_raster_callback);
		break;
	    } /* else flow through to SEQ_QUIT */
	}
    case SEQ_QUIT:
	{
	    /* quit raster display */
	    /* 
	     * HACK to prevent re-entry into this routine when deleting a 
	     * result in the raster which may invoke a QUIT event if the 
	     * result is the last one in the raster
	     */
#ifdef DEBUG
	    printf("raster SEQ_QUIT %d\n", result->status);
#endif
	    if (result->status) {
		return;
	    } else {
		result->status = 1;
	    }
	    seq_raster_shutdown(result->id, result);
	    raster_shutdown(result->interp, result->raster_win, result); 
	    return;
	}
    case SEQ_GET_OPS:
	jdata->get_ops.ops = "Remove\0"; 
	return;

    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* remove */
	    /* 
	     * HACK to prevent re-entry into this routine when deleting a 
	     * result in the raster which may invoke a QUIT event if the 
	     * result is the last one in the raster
	     */
	    if (result->status) {
		return;
	    } else {
		result->status = 1;
	    }
	    seq_raster_shutdown(result->id, result);
	    raster_shutdown(result->interp, result->raster_win, result); 
	    return;
	}
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case RESULT:
	    jdata->info.result = (void *)result;
	    break;
	case WIN_NAME:
	    {
		char *r_win = result->raster_win;
		jdata->info.result = (void *)r_win;
	    }
	break;
	}
	break;
    case SEQ_GENERIC: {
	int num_elements;
	seq_result **data;
	seq_result *sq_result;
	int num_funcs;
	int i;

	num_elements = seq_num_results();
	if (num_elements == 0)
	    return;
	
 	if (NULL == (data = (seq_result **)xmalloc(num_elements * 
						  sizeof(seq_result *))))
	    return;
	
	/* find all plotted results */
	search_reg_data(comparison3, (void **)data, &num_funcs);
	
	switch (jdata->generic.task) {
	case TASK_RASTER_QUIT_RESULTS: 
	    {
		out_raster *output;
		seq_reg_quit info;
		info.job = SEQ_QUIT;    
		
		/* find number of results in raster */
		for (i = 0; i < num_funcs; i++) {
		    sq_result = data[i];
		    output = sq_result->output;
		    
		    if (strcmp(output->raster_win, result->raster_win) == 0){
			seq_result_notify(sq_result->id, 
					  (seq_reg_data *)&info, 0);
		    }
		}
		break;
	    }
	case TASK_RASTER_MAX_SIZE: 
	    {
		out_raster *output;
		d_line *max;
		d_line *dim;
		seq_reg_info info;
		double x0, y0, x1, y1;
		double t_x0, t_y0, t_x1, t_y1;
		double max_x0, max_y0, max_x1, max_y1;
		double sf_m, sf_c;
                double max_max_x0, max_max_y0, max_max_x1, max_max_y1;

		if (NULL == (max = (d_line *)xmalloc(sizeof(d_line)))) {
		    xfree(data);
		    return;
		}
    
		max->x0 = DBL_MAX/2;
		max->y0 = DBL_MAX/2;
		max->x1 = -DBL_MAX/2;
		max->y1 = -DBL_MAX/2;

		x0 = DBL_MAX/2;
		y0 = DBL_MAX/2;
		x1 = -DBL_MAX/2;
		y1 = -DBL_MAX/2;

		max_x0 = DBL_MAX/2;
		max_y0 = DBL_MAX/2;
		max_x1 = -DBL_MAX/2;
		max_y1 = -DBL_MAX/2;

		t_x0 = DBL_MAX/2;
		t_y0 = DBL_MAX/2;
		t_x1 = -DBL_MAX/2;
		t_y1 = -DBL_MAX/2;

                max_max_x0 = 0;
                max_max_x1 = 0;
                max_max_y0 = 0;
                max_max_y1 = 0;

		/* find number of results in raster */
		for (i = 0; i < num_funcs; i++) {
		    sq_result = data[i];
		    if (sq_result->output == NULL) {
			jdata->info.result = (void *)max;
			xfree(data);
			return;
		    }

		    output = sq_result->output;
		    sf_m = output->sf_m;
		    sf_c = output->sf_c;

		    /* find dimensions of result */
		    info.job = SEQ_RESULT_INFO;
		    info.op = DIMENSIONS;
		    info.result = NULL;
		    seq_result_notify(sq_result->id, 
				      (seq_reg_data *)&info, 0);
		    if (!info.result) {
			verror(ERR_WARN, "seq_raster_callback", 
			       "No result found\n");
			xfree(data);
			return;
		    }
		    dim = (d_line *)info.result;
		    
#ifdef DEBUG
		    printf("output %s raster %s\n", output->raster_win, 
			   result->raster_win);
#endif
                    /*
                     * If the raster is going in a new window, then output->raster_win
                     * and result->raster_win will never match. Therefore need to find
                     * the maximums from any that do exist, so that maximums of the new
                     * window are not left as DBL_MAX.
                     */
                    if(dim->x0 > max_max_x0) max_max_x0 = dim->x0;
                    if(dim->x1 > max_max_x1) max_max_x1 = dim->x1;
                    if(dim->y0 > max_max_y0) max_max_y0 = dim->y0;
                    if(dim->y1 > max_max_y1) max_max_y1 = dim->y1;

		    if (strcmp(output->raster_win, result->raster_win) == 0){

			if (output->configure && (output->configure[0]->zoom 
						  == 0 ||
						  output->configure[0]->zoom 
						  == 1)) {
			    x0 = dim->x0;
			    x1 = dim->x1;
			    
			    t_y0 = sf_m * dim->y0 + sf_c;
			    t_y1 = sf_m * dim->y1 + sf_c;
			} else {
			
			    x0 = dim->x0;
			    x1 = dim->x1;
			
			    y0 = sf_m * dim->y0 + sf_c;
			    y1 = sf_m * dim->y1 + sf_c;
			}
			if (x0 < max->x0) max->x0 = x0;
			if (y0 < max->y0) max->y0 = y0;
			if (x1 > max->x1) max->x1 = x1;
			if (y1 > max->y1) max->y1 = y1;

			if (t_x0 < max_x0) max_x0 = t_x0;
			if (t_x1 > max_x1) max_x1 = t_x1;
			if (t_y0 < max_y0) max_y0 = t_y0;
			if (t_y1 > max_y1) max_y1 = t_y1;
		    }
		}

		if (max->x0 == DBL_MAX/2) max->x0 = max_x0;
		if (max->x1 == -DBL_MAX/2) max->x1 = max_x1;
		if (max->y0 == DBL_MAX/2) max->y0 = max_y0;
		if (max->y1 == -DBL_MAX/2) max->y1 = max_y1;

		/*
                 * If output->raster_win and result->raster_win
                 * never matched, max->x0 etc. will still equal
                 * DBL_MAX. In which case we need to set them to
                 * the maximums of all the windows that were found.
                 */
                if(max->x0 > max_max_x0) max->x0 = max_max_x0;
                if(max->x1 > max_max_x1) max->x1 = max_max_x1;
                if(max->y0 > max_max_y0) max->y0 = max_max_y0;
                if(max->y1 > max_max_y1) max->y1 = max_max_y1;

#ifdef DEBUG
		printf("MAX %f %f %f %f\n", max->x0, max->y0, max->x1, max->y1);
#endif
		jdata->info.result = (void *)max;
		break;
	    }
	case TASK_RASTER_WINDOW_SIZE:
	    {
		seq_reg_info info;
		d_point *max;
		d_point *pt;
		out_raster *output;

		info.job = SEQ_RESULT_INFO;
		info.op = WIN_SIZE;
		info.result = NULL;

		if (NULL == (max = (d_point *)xmalloc(sizeof(d_point)))) {
		    xfree(data);
		    return;
		}

		max->x = 0;
		max->y = 0;

		for (i = 0; i < num_funcs; i++) {
		    sq_result = data[i];
		    output = sq_result->output;
		    
		    if (strcmp(output->raster_win, result->raster_win) == 0){
			seq_result_notify(sq_result->id, 
					  (seq_reg_data *)&info, 0);
			pt = (d_point *)info.result;
#ifdef DEBUG			
			printf("id %d x %d y %f\n", i, pt->x, pt->y);
#endif
			if (pt->x > max->x) max->x = pt->x;
			if (pt->y > max->y) max->y = pt->y;
		    }
		}
		jdata->info.result = (void *)max;
		break;
	    }
	case TASK_RASTER_ZOOM: 
	    {
		/* 
		 * return 0 if all the results have configure.zoom of 0 or 1
		 * return 1 if any of the results have configure.zoom of 2
		 */
		out_raster *output;
		int zoom = -1;

		for (i = 0; i < num_funcs; i++) {
		    sq_result = data[i];
		    output = sq_result->output;

		    if (output == NULL) {
		      zoom = -1;
		      jdata->info.result = (void *)zoom;
		      break;
		    }
		    if (!output->raster_win) {
			zoom = 0;
			break;
		    }
		    
		    if (!output->configure) {
			zoom = -1;
			break;
		    }
		    
		    if (strcmp(output->raster_win, result->raster_win) == 0){
			if (output->configure[0]->zoom == 2) {
			    zoom = 1;
			    break;
			} else {
			    zoom = 0;
			}
		    }
		}
		jdata->info.result = (void *)zoom;
		break; /* TASK_RASTER_ZOOM */
	    }
	case TASK_RASTER_UPDATE_ALL: {
	    /* contains new raster */
	    update_struct *update = (update_struct *)jdata->generic.data; 
	    
	    /*
	     * result->raster_win is the old raster
	     * update->raster is the new raster
	     */
	    SeqUpdateRasterWindow(result->interp, update->raster,
				  result->raster_win, 
				  result->id, update->old_id,
				  update->job);
	    
	    break; /* TASK_RASTER_UPDATE_ALL */
	} 
	case TASK_RASTER_UPDATE_SINGLE: 
	    {
		/* contains old raster & result id */
		update_struct *update = (update_struct *)jdata->generic.data; 
		
		SeqUpdateResultWindow(result->interp, update->raster, 
				      result->raster_win, update->old_id, 
				      result->id, update->id, update->job);
		break; /* TASK_RASTER_UPDATE_SINGLE */
	    }
	} /* switch */
	xfree(data);
	break;
    } /* case GENERIC */
    case SEQ_CURSOR_NOTIFY:
	{
	    cursor_t *cursor = (cursor_t *)jdata->cursor_notify.cursor;
	    cursor_t *gc;
	    Tk_Raster *raster;
	    Tcl_CmdInfo info;
	    int i;
	    int direction = -1;

	    Tcl_GetCommandInfo(result->interp, result->raster_win, &info);
	    raster = (Tk_Raster*)info.clientData;
	    
	    for (i = 0; i < result->num_seq_id; i++) {
		for (gc = result->cursor[i]; gc; gc = gc->next) {
		/* gc = result->cursor[i]; */

		    if (GetSeqId(seq_num) == result->seq[i].seq_id &&
			gc == cursor && 
			gc->direction == result->seq[i].direction) {

			direction = cursor->direction;
			
			if (direction != -1) {
			    raster_update_cursor(result, cursor, GetSeqId(seq_num), 
						 raster, 1, direction);
			} else {
			    printf("ERROR: SEQ_CURSOR_NOTIFY no direction found\n");
			}
		    } 
		}
	    }
	    result->cursor_array[cursor->id].prev_pos = cursor->abspos;
#ifdef REMOVE
	    cursor->prev_pos = cursor->abspos;
#endif
	    break;
	}
    }
}

/*
 * return an array of nip results in a particular raster 
 */
seq_result ** seq_get_raster_results(char *raster_win, 
				     int *num_results)
{
    int num_elements;
    seq_result **data;
    seq_result **result, *first;
    int num_funcs;
    int i;
    int cnt= 0;
    out_raster *output;

    num_elements = seq_num_results();
    if (num_elements == 0)
	return NULL;
    
    if (NULL == (data = (seq_result **)xmalloc(num_elements * 
					      sizeof(seq_result *))))
	return NULL;

    /* find all plotted results */
    search_reg_data(comparison3, (void **)data, &num_funcs);

    if (NULL == (result = (seq_result **)xmalloc(num_funcs * 
						 sizeof(seq_result *) +
						 num_funcs *
						 sizeof(seq_result))))
	return NULL;
    first = (seq_result *)(((char *)result) +
			   num_funcs * sizeof(seq_result *));
    for (i = 0; i < num_funcs; i++) {
	result[i] = &first[i];
    }	
    
    for (i = 0; i < num_funcs; i++) {

	result[cnt] = data[i];
	output = result[cnt]->output;
		    
	if (strcmp(output->raster_win, raster_win) == 0){
	    cnt++;
	}
    }	    

    xfree(data);
    *num_results = cnt;
    return result;
}

int seq_raster_reg(Tcl_Interp *interp,
		   char *raster_win, 
		   seq_id_dir *seq_array,
		   int num_seq_id)
{
    RasterResult *raster_result;
    int id;
    int i;
    int line_width;
    int *seq_counts_h, *seq_counts_v, nseq;

    if (NULL == (raster_result = (RasterResult *)xmalloc(sizeof(RasterResult))))
	return -1;

    if (NULL == (raster_result->cursor = (cursor_t **)xmalloc(MAX_NUM_SEQ * 
						      sizeof(cursor_t*))))
	return -1;

    id = get_reg_id();
  
    raster_result->op_func = seq_raster_callback;

    sprintf(raster_result->raster_win, "%s%d", raster_win, id);
    raster_result->interp = interp;
    raster_result->ed_cursor = -1;
    raster_result->id = id;
    raster_result->num_seq_id = num_seq_id;
    raster_result->seq = seq_array;
    raster_result->num_results = 0;
    raster_result->status = 0;

    /* initialised raster cursor array */
    for (i = 0; i < NUM_CURSORS; i++) {
	raster_result->cursor_array[i].env = -2;
	raster_result->cursor_array[i].visible[HORIZONTAL] = 0;
	raster_result->cursor_array[i].visible[VERTICAL] = 0;
	raster_result->cursor_array[i].prev_pos = -1;
    }

    line_width = get_default_int(interp, spin_defs, 
				 w("SEQ.CURSOR.LINE_WIDTH"));

    /* need to register with each seq_id in a raster */
    nseq = NumSequences();
    if (NULL == (seq_counts_h = (int *)xmalloc(nseq * sizeof(int))))
	return -1;
    if (NULL == (seq_counts_v = (int *)xmalloc(nseq * sizeof(int))))
	return -1;
    for (i = 0; i < nseq; i++) {
	seq_counts_h[i] = 0;
	seq_counts_v[i] = 0;
    }

    for (i = 0; i < num_seq_id; i++) {
	int snum;

	snum = GetSeqNum(seq_array[i].seq_id);
	if (seq_array[i].direction == HORIZONTAL) {
	    raster_result->cursor[i] = 
		create_cursor(snum, 0, NULL, line_width, ++seq_counts_h[snum],
			      seq_array[i].direction);
	} else {
	    raster_result->cursor[i] = 
		create_cursor(snum, 0, NULL, line_width, ++seq_counts_v[snum],
			      seq_array[i].direction);
	}
	raster_result->cursor_array[raster_result->cursor[i]->id].env = -1;
    }

    for (i = 0; i < num_seq_id; i++) {
	seq_register(GetSeqNum(seq_array[i].seq_id), seq_raster_callback, 
		     (void *)raster_result, SEQ_RASTER, id);
    }

    xfree(seq_counts_h);
    xfree(seq_counts_v);
    return id;
}

/*
 * return raster window which conforms to seq_id, type and frame
 */
static char *seq_get_seq_type(int seq_id,
			      int type,
			      int frame)
{
    int num_elements;
    seq_result **data;
    int num_funcs;
    int i;
    out_raster *output;

    num_elements = seq_num_results();
    if (num_elements == 0) {
	return NULL;
    }
    if (NULL == (data = (seq_result **)xmalloc(num_elements * 
					      sizeof(seq_result *))))
	return NULL;
    
    /* find all plotted results */
    search_reg_data(comparison3, (void **)data, &num_funcs);
    for (i = num_funcs-1; i >= 0; i--) {
#ifdef DEBUG
	printf("seqid %d %d type %d %d frame %d %d\n",
	       data[i]->seq_id, seq_id, data[i]->type, type,
	       data[i]->frame, frame);
#endif
	if (data[i]->seq_id[0] == seq_id && data[i]->type & type) {
	    if (data[i]->frame == 0 || frame == data[i]->frame) {
		output = data[i]->output;
		xfree(data);
		if (output == NULL)
		    return NULL;
		return output->raster_win; 
	    }
	}
    }	    
    xfree(data);
    return NULL;
}

/*
 * returns the appropriate raster name for either a new raster or the 
 * default one associated with a sequence (NIP)
 */
char *get_raster_frame_graph(Tcl_Interp *interp,
			     int seq_id,
			     int type,
			     int frame)
{
    char *val;
    char *win, *raster_win, *raster_frame;
    int seq_num;
    seq_id_dir *seq_array;
    int raster_id;

    /* 
     * if type is set, find a raster window containing a plot of that type
     */
    if (type != -1) {
	win = seq_get_seq_type(seq_id, type, frame);
	if (win != NULL) {
	    raster_win = strdup(win);
	    return raster_win;
	}
    }
    seq_num = GetSeqNum(seq_id);
	
    /* find default raster_win for sequence */
    if (NULL == (raster_frame = GetRaster(seq_num))) {
	/* if no default defined, then create one */
	
	if (NULL == (raster_frame = (char *)xmalloc(1024 * sizeof(char))))
	    return NULL;
	
	Tcl_VarEval(interp, "GetRasterWindow", NULL);
	strcpy(raster_frame, Tcl_GetStringResult(interp));
	
	/* define default raster for sequence */
	SetRaster(seq_num, raster_frame);
    }
    if (NULL == (raster_win = (char *)xmalloc(1024 * sizeof(char))))
	return NULL;
    strcpy(raster_win, raster_frame);
    val = get_default_string(interp, tk_utils_defs, w("RASTER.R.WIN"));
    if (NULL == (seq_array = (seq_id_dir *)xmalloc(MAX_NUM_SEQ*sizeof(seq_id_dir))))
	return NULL;
    
    seq_array[0].direction = HORIZONTAL;
    seq_array[0].seq_id = seq_id;
    
    sprintf(raster_win, "%s%s", raster_win, val);
    raster_id = seq_raster_reg(interp, raster_win, seq_array, 1);
    sprintf(raster_win, "%s%d", raster_win, raster_id);
    return raster_win;

}

/*
 * returns the appropriate raster name for either a new raster or the 
 * default one associated with a sequence pair (SIP)
 */
int get_raster_frame_dot(Tcl_Interp *interp,
			 int seq_id_h,
			 int seq_id_v,
			 char *raster_win)
{
    int i, j;
    int num_elements;
    RasterResult **data = NULL;
    RasterResult *raster_result = NULL;
    int num_funcs;
    int h_found = 0;
    int v_found = 0;
    seq_id_dir *seq_array;
    char *val;
    int raster_id;
    int p_seq_id_h, p_seq_id_v, p_result_id;
    int line_width;

    raster_win[0] = 0;
    num_elements = seq_num_results();
    if (num_elements > 0) {
	
	if (NULL == (data = (RasterResult **)xmalloc(num_elements * 
						     sizeof(RasterResult *))))
	    return -1;
	
	/* find all rasters */
	search_reg_data(comparison_raster, (void **)data, &num_funcs);
    
	/* 
	 * loop through all results in raster_win to see if any are registered with
	 * seq_id
	 */
	p_seq_id_h = GetParentalSeqId(GetSeqNum(seq_id_h));
	p_seq_id_v = GetParentalSeqId(GetSeqNum(seq_id_v));
	for (i = 0; i < num_funcs; i++) {
	    h_found = v_found = 0;
	    raster_result = data[i];
	    for (j = 0; j < raster_result->num_seq_id; j++) {
		
		/* look at parental seq_id */
		p_result_id = GetParentalSeqId(GetSeqNum(raster_result->seq[j].seq_id));

		if (raster_result->seq[j].direction == HORIZONTAL && 
		    p_result_id == p_seq_id_h) {
		    h_found = 1;
		}
		if (raster_result->seq[j].direction == VERTICAL &&
		    p_result_id == p_seq_id_v) {
		    v_found = 1;
		}
	    }
	    if (h_found && v_found) {
		strcpy(raster_win, raster_result->raster_win);
		break;
	    }
	}	    
	/* 
	 * found raster_win but also need to check if the raster is regisered
	 * with the current seq_id's ie if using ranges of seqs
	 */
	if (raster_win[0]) {
	    h_found = v_found = 0;
	    for (i = 0; i < num_funcs; i++) {
		raster_result = data[i];
		for (j = 0; j < raster_result->num_seq_id; j++) {
		    if (raster_result->seq[j].direction == HORIZONTAL && 
			raster_result->seq[j].seq_id == seq_id_h) {
			h_found = 1;
		    }
		    if (raster_result->seq[j].direction == VERTICAL &&
			raster_result->seq[j].seq_id == seq_id_v) {
			v_found = 1;
		    }
		}
	    }	    
	    /* if failed to find current seq_id's registered, then register */
	    line_width = get_default_int(interp, spin_defs, 
					 w("SEQ.CURSOR.LINE_WIDTH"));
	    if (!h_found) {
		add_seq_to_raster(raster_result, seq_id_h, 
				  GetSeqNum(seq_id_h), HORIZONTAL, 
				  line_width, seq_raster_callback);
	    }

	    if (!v_found) {
		add_seq_to_raster(raster_result, seq_id_v, 
				  GetSeqNum(seq_id_v), VERTICAL, 
				  line_width, seq_raster_callback);
	    }	
	    xfree(data);
	    return 0;
	}
    }

    /* failed to find existing raster, therefore create new window & raster */

    Tcl_VarEval(interp, "GetRasterWindow", NULL);
    strcpy(raster_win, Tcl_GetStringResult(interp));

    val = get_default_string(interp, tk_utils_defs, w("RASTER.R.WIN"));
    if (NULL == (seq_array = (seq_id_dir *)xmalloc(MAX_NUM_SEQ*sizeof(seq_id_dir)))) {
	xfree(data);
	return -1;
    }

    seq_array[0].direction = HORIZONTAL;
    seq_array[0].seq_id = seq_id_h;
    seq_array[1].direction = VERTICAL;
    seq_array[1].seq_id = seq_id_v;

    sprintf(raster_win, "%s%s", raster_win, val);
    
    raster_id = seq_raster_reg(interp, raster_win, seq_array, 2);

    sprintf(raster_win, "%s%d", raster_win, raster_id);
    xfree(data);
    return 0;
}

int raster_select_cursor_graph(int raster_id,
			       Tk_Raster *raster,
			       char *raster_win,
			       int pos,
			       int max_dist)
{
    int cursor_id = -1;
    int i;
    int closest = INT_MAX;
    int id, cursor_pos;
    int diff;
    RasterResult *raster_result;

    if (NULL == (raster_result = raster_id_to_result(raster_id)))
	return -1;

    for (i = 0; i < raster_result->num_seq_id; i++) {

	id = find_nearest_cursor(raster, 
				 GetSeqNum(raster_result->seq[i].seq_id), 
				 pos, max_dist,raster_result->seq[i].direction,
				 &cursor_pos);
	if (id != -1) {
	    diff = abs(cursor_pos - pos);
		if (diff < closest) {
		    closest = diff;
		    cursor_id = id;
		}
	}
    }
    return cursor_id;
}

int raster_select_cursor_dot(int raster_id,
			     Tk_Raster *raster,
			     char *raster_win,
			     int rx,
			     int ry,
			     int max_dist,
			     int *cursor_id_h,
			     int *cursor_id_v)
{
    int cursor_id = -1;
    int i, j, k, l;
    int closest = INT_MAX;
    int cursor_pos, dir;
    int diff;
    RasterResult *raster_result;
    cursor_res *id;
    int cnt = 0;
    seq_result **result;
    int *found_h, *found_v;
    int num_results, direction = -1;
    int *id_h, *id_v;

    /* initialise cursor_id's */
    *cursor_id_h = -1;
    *cursor_id_v = -1;

    raster_result = raster_id_to_result(raster_id);
    if (NULL == (id = (cursor_res *)xmalloc(raster_result->num_seq_id * 
					    sizeof(cursor_res))))
	return -1;

    for (i = 0; i < raster_result->num_seq_id; i++) {
	if (raster_result->seq[i].direction == HORIZONTAL) {
	    id[cnt].cursor_id = find_nearest_cursor(raster, 
						    GetSeqNum(raster_result->seq[i].seq_id), 
						    rx, max_dist, 
						    raster_result->seq[i].direction, 
						    &cursor_pos);
	    dir = HORIZONTAL;
	    id[cnt].seq_id = raster_result->seq[i].seq_id;
	    id[cnt].direction = HORIZONTAL; 
	} else {

	    id[cnt].cursor_id = find_nearest_cursor(raster, 
						    GetSeqNum(raster_result->seq[i].seq_id), 
						    ry, max_dist, 
						    raster_result->seq[i].direction, 
						    &cursor_pos);
	    dir = VERTICAL;
	    id[cnt].seq_id = raster_result->seq[i].seq_id;
	    id[cnt].direction = VERTICAL;
	}
	if (id[cnt].cursor_id != -1) {
	    if (raster_result->seq[i].direction == HORIZONTAL) {
		diff = abs(cursor_pos - rx);
	    } else {
		diff = abs(cursor_pos - ry);
	    }
	    if (diff < closest) {
		closest = diff;
		cursor_id = id[cnt].cursor_id;
		direction = dir;
	    }
	    cnt++;
	}
    }
    
    /* only no or one result */
    if (cnt < 2) {
	if (direction == HORIZONTAL) {
	    *cursor_id_h = cursor_id;
	} else if (direction == VERTICAL) {
	    *cursor_id_v = cursor_id;
	}
	xfree(id);
	return 0;
    }

    result = seq_get_raster_results(raster_win, &num_results);
    k = l = 0;
    
    if (NULL == (id_h = (int *)xmalloc(cnt * num_results * sizeof(int))))
	return -1;

    if (NULL == (id_v = (int *)xmalloc(cnt * num_results * sizeof(int))))
	return -1;

    if (NULL == (found_h = (int *)xmalloc(cnt * num_results * sizeof(int))))
	return -1;

    if (NULL == (found_v = (int *)xmalloc(cnt * num_results * sizeof(int))))
	return -1;

    {
	int len = cnt * num_results;
	
	for (i = 0; i < len; i++) {
	    found_h[i] = -1;
	    found_v[i] = -1;
	}
    }

    for (i = 0; i < cnt; i++) {
	/* need to find cursor pairs */
	for (j = 0; j < num_results; j++) {
#ifdef DEBUG
	    printf("result[%d] h=%d v=%d id[%d] %d\n", j, 
		   result[j]->seq_id[HORIZONTAL], result[j]->seq_id[VERTICAL],
		   i, id[i].seq_id);
#endif		   
	    if (id[i].direction == HORIZONTAL && 
		 result[j]->seq_id[HORIZONTAL] == id[i].seq_id) {
		found_h[k] = j;
		id_h[k++] = id[i].cursor_id;
	    }
	    if (id[i].direction == VERTICAL && 
		result[j]->seq_id[VERTICAL] == id[i].seq_id) {
		found_v[l] = j;
		id_v[l++] = id[i].cursor_id;
	    }
	}
    }
    for (i = 0; i < k; i++) {
	for (j = 0; j < l; j++) {
	    if (found_h[i] == found_v[j] && found_h[i] != -1) {
		*cursor_id_h = id_h[i];
		*cursor_id_v = id_v[j];
		break;
	    }
	}
    }

    xfree(id);
    xfree(result);
    xfree(id_h);
    xfree(id_v);
    xfree(found_h);
    xfree(found_v);
    return 0;
}

int SeqDeregisterRasterWindow(int seq_id,
			      RasterResult *old_result,
			      char *raster_win)
{
    int i;
    int num_elements;
    seq_result **data;
    seq_result *result;
    int num_funcs;
    int found = 0;
    out_raster *output;
    int result_cnt = 0;

#ifdef DEBUG
    printf("SeqDeregisterRasterWindow\n");
#endif
    num_elements = seq_num_results();
    if (num_elements == 0)
	return -1;
	
    if (NULL == (data = (seq_result **)xmalloc(num_elements * 
					      sizeof(seq_result *))))
	return -1;
    
    /* find all plotted results */
    search_reg_data(comparison3, (void **)data, &num_funcs);
	
    /* 
     * loop through all results in raster_win to see if any are registered with
     * seq_id
     */
    for (i = 0; i < num_funcs; i++) {
	result = data[i];
	output = result->output;
		  
	if (strcmp(output->raster_win, raster_win) == 0){
	    result_cnt++;
	    if ((result->seq_id[HORIZONTAL] == seq_id) || 
		(result->seq_id[VERTICAL] == seq_id)) {
		found = 1;
	    }
	}
    }	    
    /* 
     * if no results are registered with seq_id, then deregister raster_win
     * from this seq_id
     */
    if (!found && result_cnt > 0) {
	delete_seq_from_raster(seq_id, GetSeqNum(seq_id), old_result,
			       seq_raster_callback);
    }
    xfree(data);
    return 0;
}

int init_graph_raster(Tcl_Interp *interp, 
		      int seq_id,  
		      int result_id,
		      char *raster_win, 
		      int raster_id,
		      config *configure,
		      char *colour,
		      int line_width)
{
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    out_raster *output; 
    char *opts[5];
    int seq_num;
    seq_cursor_notify cn;
    RasterResult *raster_result;
    cursor_t *cursor;
    int superimpose = 1;
    seq_result *result;
    r_graph *data;

    seq_num = GetSeqNum(seq_id);
    if (NULL == (result = result_data(result_id, seq_num))) {
	return -1;
    }
    output = result->output;
    data = result->data;
    
   if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	return -1;

    raster = (Tk_Raster*)info.clientData;

    RasterInitPlotFunc(raster, SeqRasterPlotFunc);
    
    /* need to check if superimposing result on another plot */
    raster_result = raster_id_to_result(raster_id);
    if (raster_result->num_results == 0) {
	superimpose = 0;
    }

    if (NULL == (opts[1] = (char *)xmalloc((strlen(colour)+1) * 
					   sizeof(char))))
	return -1;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;

    strcpy(output->raster_win, raster_win);
    output->raster_id = raster_id;
    output->interp = interp;
    output->hidden = 0;

    opts[0] = "-fg";
    strcpy(opts[1], colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = NULL;

    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);

    if (NULL == (output->configure = (config **)xmalloc(sizeof(config*))))
	return -1;

    output->configure[0] = configure;
    output->n_configure = 1;
    output->scroll = 'b';
    output->sf_m = 1.0;
    output->sf_c = 0.0;

    if (superimpose) {
	SeqSuperimposeResult(interp, output->raster_win, result_id, 
			     data->dim.x0, data->dim.y0, 
			     data->dim.x1, data->dim.y1);
    } else {
      /* 
       * update world scroll region for this raster, if necessary to be the 
       * largest range 
       */
	RasterSetWorldScroll(raster, data->dim.x0, data->dim.y0, 
			     data->dim.x1, data->dim.y1);
	SeqAddRasterToWindow(interp, raster_win, result->graph);
    }

    raster_result = raster_id_to_result(raster_id);
    cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);

    /* move cursor to start position if no cursor is yet present */
    if (raster_result->cursor_array[cursor->id].prev_pos == -1 
	&& data->dim.x0 > raster_result->cursor_array[cursor->id].prev_pos) {
	cursor->abspos = data->dim.x0;
    }

    AddResultToRaster(raster_result);
    ReplotAllCurrentZoom(interp, output->raster_win);
    xfree(opts[1]);
    xfree(opts[3]);

    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(seq_num, (seq_reg_data *)&cn); 

  return 0;
}

int seq_gene_search_plot(Tcl_Interp *interp,
			 int result_id,
			 int seq_num,
			 char *raster_win,
			 char *colour,
			 int line_width)
{
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    char *opts[5];
    double min, max;
    int superimpose;
    RasterResult *raster_result;
    seq_result *result;
    config *configure;
    out_raster *output;
    gene_search *data;

    if (NULL == (output = (out_raster *)xmalloc(sizeof(out_raster))))
	return -1;

    if (NULL == (result = result_data(result_id, seq_num)))
	return -1;

    result->output = output;
    data = result->data;

    output->scroll = 'b';

    output->sf_m = 1.0;
    output->sf_c = 0.0;

    if (NULL == (configure = (config *)xmalloc(sizeof(config))))
	return -1;

    if (NULL == (output->configure = (config **)xmalloc(sizeof(config*))))
	return -1;

    configure->position = 0.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;    

    output->configure[0] = configure;
    output->n_configure = 1;
 
    superimpose = 1;

    if (NULL == (opts[1] = (char *)xmalloc(100 * sizeof(char))))
	return -1;
    
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;
    
    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	return -1;
    raster = (Tk_Raster*)info.clientData;

    RasterInitPlotFunc(raster, SeqRasterPlotFunc);

    strcpy(output->raster_win, raster_win);

    output->interp = interp;
    output->hidden = 0;

    /* update the world scroll region for all rasters in this window */
    /* need to check if superimposing result on another plot */
    raster_result = raster_name_to_result(interp, raster_win);
    if (raster_result->num_results == 0) {
	superimpose = 0;
    }

    if (!superimpose) {
	RasterSetWorldScroll(raster, data->dim.x0, data->dim.y0, 
			     data->dim.x1, data->dim.y1);
    } 

    opts[0] = "-fg";
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = NULL;

    strcpy(opts[1], colour);
    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);

    if (superimpose) {
	min = DBL_MAX;
	max = -DBL_MAX;
	if (min > data->dim.y0)
	    min = data->dim.y0;
	if (max < data->dim.y1)
	    max = data->dim.y1;
    }
    if (superimpose) {
	SeqSuperimposeResult(interp, output->raster_win, result_id, 
			     data->dim.x0, min, 
			     data->dim.x1, max);
    } else {
	SeqAddRasterToWindow(interp, raster_win, result->graph);
    }
    ReplotAllCurrentZoom(interp, output->raster_win);
    xfree(opts[1]);
    xfree(opts[3]);

    return 0;
}

int init_gene_search_raster(Tcl_Interp *interp, 
			    int num_items,
			    char **win_list, 
			    char **id_list,
			    int seq_id, 
			    char **result_id,
			    char **colour_list,
			    int line_width)
{
    char *seq;
    int seq_len;
    int sequence_type;
    int seq_num;
    seq_cursor_notify cn;
    RasterResult *raster_result;
    cursor_t *cursor;
    int i;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    raster_result = raster_id_to_result(atoi(id_list[0]));
    cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);
    
    for (i = 0; i < num_items; i++) {
	seq_gene_search_plot(interp, atoi(result_id[i]), seq_num, win_list[i], 
			     colour_list[i], line_width);
    }
    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);
    cn.job = SEQ_CURSOR_NOTIFY;

    for (i = 0; i < num_items; i++) {
	raster_result = raster_id_to_result(atoi(id_list[i]));
	cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);
	cn.cursor = cursor;
	cn.cursor->job = CURSOR_MOVE;
	seq_notify(seq_num, (seq_reg_data *)&cn); 
	AddResultToRaster(raster_result);
    }
    return 0;

}

int init_stick_raster(Tcl_Interp *interp, 
		      int seq_id,  
		      int result_id,
		      char *raster_win, 
		      int raster_id,
		      config *configure,
		      char *colour,
		      int line_width,
		      int tick_ht)
{
    out_raster *output; 
    int seq_num;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    char *opts[5];
    seq_cursor_notify cn;
    RasterResult *raster_result;
    cursor_t *cursor;
    int superimpose = 1;
    seq_result *result;
    stick *data;

    if (NULL == (output = (out_raster *)xmalloc(sizeof(out_raster))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    if (NULL == (result = result_data(result_id, seq_num)))
	return -1;

    data = result->data;
    result->output = (void *)output;
    
    /* if given a raster name, then draw a tick at start each string */    
    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	return -1;
    raster = (Tk_Raster*)info.clientData;

    /* HACK - don't know where to put this yet! */
    RasterInitPlotFunc(raster, SeqRasterPlotFunc);

    /* need to check if superimposing result on another plot */
    raster_result = raster_id_to_result(raster_id);
    
    if (raster_result->num_results == 0) {
	superimpose = 0;
    }

    if (NULL == (opts[1] = (char *)xmalloc((strlen(colour)+1) * 
					   sizeof(char))))
	return -1;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;

    opts[0] = "-fg";
    strcpy(opts[1], colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = NULL;

    strcpy(output->raster_win, raster_win);
    output->interp = interp;
    output->hidden = 0;
    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);
    
   if (NULL == (output->configure = (config **)xmalloc(sizeof(config*))))
	return -1;

    output->configure[0] = configure;
    output->scroll = 'x';
    output->sf_m = 1.0;
    output->sf_c = 0.0;

    if (superimpose) {
	SeqSuperimposeResult(interp, output->raster_win, result_id, 
			     data->ap_array[0].dim.x0, 
			     data->ap_array[0].dim.y0, 
			     data->ap_array[0].dim.x1, 
			     data->ap_array[0].dim.y1);

    } else {
      /* 
       * update world scroll region for this raster, if necessary to be the 
       * largest range 
       */
	RasterSetWorldScroll(raster, data->ap_array[0].dim.x0, 
			     data->ap_array[0].dim.y0, 
			     data->ap_array[0].dim.x1, 
			     data->ap_array[0].dim.y1);
	SeqAddRasterToWindow(interp, raster_win, result->graph);
    }

    raster_result = raster_id_to_result(raster_id);
    cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);

    /* move cursor to start position if no cursor is yet present */
    if (raster_result->cursor_array[cursor->id].prev_pos == -1 
	&& data->ap_array[0].dim.x0 > raster_result->cursor_array[cursor->id].prev_pos) {
	cursor->abspos = data->ap_array[0].dim.x0;
    }

    AddResultToRaster(raster_result);
    ReplotAllCurrentZoom(interp, output->raster_win);

    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);
    cn.job = SEQ_CURSOR_NOTIFY;

    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(seq_num, (seq_reg_data *)&cn); 

    xfree(opts[1]);
    xfree(opts[3]);

    return 0;
}

int init_dot_plot(Tcl_Interp *interp, 
		  int seq_id_h,  
		  int seq_id_v,
		  int result_id,
		  char *name,
		  char *raster_win, 
		  int raster_id,
		  char **opts,
		  int n_opts,
		  int plot_style,
		  d_line dim)
{
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    out_raster *output; 
    seq_cursor_notify cn;
    RasterResult *raster_result;
    cursor_t *cursor;
    seq_result *result;
    int seq1_num, seq2_num;
    config *configure;

    if (NULL == (output = (out_raster *)xmalloc(sizeof(out_raster))))
	return -1;

    seq1_num = GetSeqNum(seq_id_h);
    seq2_num = GetSeqNum(seq_id_v);

    result = result_data(result_id, seq1_num);
    result->output = output;

    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	return -1;

    raster = (Tk_Raster*)info.clientData;

    RasterInitPlotFunc(raster, SeqRasterPlotFunc);
    raster_result = raster_id_to_result(raster_id);

    output->name = strdup(name);
    strcpy(output->raster_win, raster_win);
    output->interp = interp;
    output->hidden = 0;
    output->env_index = CreateDrawEnviron(interp, raster, n_opts, opts);
    output->plot_style = plot_style;

    if (NULL == (configure = (config *)xmalloc(sizeof(config))))
	return -1;

    configure->position = 0.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    if (NULL == (output->configure = (config **)xmalloc(sizeof(config*))))
	return -1;

    output->configure[0] = configure;
    output->n_configure = 1;

    output->sf_m = 1.0;
    output->sf_c = 0.0;
    
    /* 
     * if the raster already has a result plotted in it, then scale new result
     * to fit on scale of previous
     */
    if (raster_result->num_results > 0) {
	double o_wx0, o_wy0, o_wx1, o_wy1;
	char *parent;

	/* scale all plots in x & y to the original and update zoom list */
	Tcl_VarEval(interp, "winfo parent ", raster_win, NULL);
	parent = strdup(Tcl_GetStringResult(interp));

	RasterGetWorldScroll(raster, &o_wx0, &o_wy0, &o_wx1, &o_wy1);
	SeqReSetRasterWindowSize(interp, raster_win, result->graph);
	ReplotAllRasterWindow(interp, raster_win);
    
	UpdateZoomList(interp, dim.x0,  dim.y0, dim.x1, 
		       dim.y1, o_wx0, o_wy0, o_wx1, o_wy1, 
		       parent, 0);
	free(parent);
    } else {

	/* set world scroll region for this raster */
	RasterSetWorldScroll(raster, dim.x0,  dim.y0,
			     dim.x1, dim.y1);

	SeqAddRasterToWindow(interp, output->raster_win, result->graph);
	
	ReplotAllCurrentZoom(interp, output->raster_win);
    }
    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */

    Tcl_VarEval(interp, "update idletasks ", NULL);
    cursor = find_raster_result_cursor(raster_result, GetSeqId(seq1_num), 
				       HORIZONTAL);
    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(seq1_num, (seq_reg_data *)&cn); 

    cursor = find_raster_result_cursor(raster_result, GetSeqId(seq2_num), 
				       VERTICAL);
    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(seq2_num, (seq_reg_data *)&cn); 

    AddResultToRaster(raster_result);

    return 0;
}

