#include <string.h>

#include "tkRaster.h"
#include "cli_arg.h"
#include "raster_structs.h"
#include "seq_reg.h"
#include "xalloc.h"
#include "seq_raster.h"
#include "tcl_utils.h"
#include "text_output.h"
#include "raster_globals.h"
#include "ruler_tick.h"

/*
 * replot all the matches in the raster
 */
int 
RasterReplotZoom(ClientData clientData, 
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[]) 
{    
    replot_arg args;

    cli_args a[] = {
	{"-raster",   ARG_STR, 1, NULL, offsetof(replot_arg, raster)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ReplotAllZoom(interp, args.raster);
    return TCL_OK;
}

/*
 * replot all the matches in the raster
 */
int 
RasterReplotAll(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{    
    replot_arg args;
    cli_args a[] = {
	{"-raster", ARG_STR, 1, NULL, offsetof(replot_arg, raster)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ReplotAllCurrentZoom(interp, args.raster);
    return TCL_OK;
}

/*
 * set the colour and width of the matches in the raster of result 'id' 
 */
int 
RasterConfig(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{    
    result_config_arg args;
    out_raster *output;
    seq_reg_info info;
    char *opts[5];
    Tk_Raster *raster;
    Tcl_CmdInfo cmd_info;

    cli_args a[] = {
	{"-index", ARG_INT, 1, NULL, offsetof(result_config_arg, id)},
	{"-fill",  ARG_STR, 1, NULL, offsetof(result_config_arg, colour)},
	{"-width", ARG_INT, 1, NULL, offsetof(result_config_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    if (NULL == (opts[1] = (char *)xmalloc((strlen(args.colour)+1) * 
					   sizeof(char))))
	return TCL_OK;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return TCL_OK;

    /* 
     * update output parameters with new colour and line width
     */
    info.job = SEQ_RESULT_INFO;
    info.op = OUTPUT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    output = (out_raster *)info.result;

    opts[0] = "-fg";
    strcpy(opts[1], args.colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", args.line_width);
    opts[4] = NULL;

    Tcl_GetCommandInfo(interp, output->raster_win, &cmd_info);
    raster = (Tk_Raster*)cmd_info.clientData;

    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);

    /* 
     * need to replot everything for line width to take effect - otherwise
     * there is no way of removing the previous points
     */
    ReplotAllCurrentZoom(interp, output->raster_win);
    xfree(opts[1]);
    xfree(opts[3]);
    return TCL_OK;
}

/*
 * get the colour of the matches in the raster from result id
 */
int 
RasterGetConfig(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{    
    result_info_arg args;
    out_raster *output;
    Tk_Raster *raster;
    Tcl_CmdInfo cmd_info;
    seq_reg_info info;
    DrawEnvironment* drawenvptr;
    int result;
    
    cli_args a[] = {
	{"-index", ARG_INT, 1, NULL, offsetof(result_info_arg, id)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = OUTPUT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    output = (out_raster *)info.result;

    Tcl_GetCommandInfo(interp, output->raster_win, &cmd_info);
    raster = (Tk_Raster*)cmd_info.clientData;

    Tcl_ResetResult(interp);

    DrawEnvIndex (interp, raster, output->env_index, (ClientData*)&drawenvptr);
    result = SetDrawEnv (interp, raster, drawenvptr);

    vTcl_SetResult(interp, "{-fill %s} {-width %d}", 
		   GetRasterColour(interp, raster, output->env_index),
		   GetRasterLineWidth(interp, raster, output->env_index));
    return TCL_OK;
}

/*
 * move a result id from raster_old to raster_new and replot
 */
int
UpdateRasterWindow(ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, 
		   char *argv[])
{
    update_win_arg args;
    seq_reg_generic gen;
    update_struct update;

    cli_args a[] = {
	{"-old_id", ARG_INT, 1,NULL, offsetof(update_win_arg, old_id)},
	{"-new_id", ARG_INT, 1,NULL, offsetof(update_win_arg, new_id)},
	{"-old", ARG_STR, 1, NULL, offsetof(update_win_arg, raster_old)},
	{"-new", ARG_STR, 1, NULL, offsetof(update_win_arg, raster_new)},
	{"-result_id", ARG_INT, 1, "-1", offsetof(update_win_arg, result_id)},
	{"-job", ARG_STR, 1, NULL, offsetof(update_win_arg, job)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    gen.job = SEQ_GENERIC;

    if (strcmp(args.job, "SUPER") == 0) {
	update.job = SUPER;
    } else if (strcmp(args.job, "ADD") == 0)  {
	update.job = ADD;
    } else if (strcmp(args.job, "NEW") == 0) {
	update.job = NEW;
    } else {
	verror(ERR_WARN, "UpdateRasterWindow", "No such job \n");
	return TCL_OK;
    }

    if (args.result_id == -1) {
	update.raster = args.raster_old;
	update.id = args.result_id;
	update.old_id = args.old_id;
	gen.data = (void *)&update;

	gen.task = TASK_RASTER_UPDATE_ALL;
	seq_result_notify(args.new_id, (seq_reg_data *)&gen, 0);
    } else {
	update.raster = args.raster_old;
	update.id = args.result_id;
	update.old_id = args.old_id;
	gen.data = (void *)&update;

	gen.task = TASK_RASTER_UPDATE_SINGLE;
	seq_result_notify(args.new_id, (seq_reg_data *)&gen, 0);
    }
    return TCL_OK;
}

int
RasterResults(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[])
{
    result_arg args;
    seq_reg_generic gen;

    cli_args a[] = {
	{"-id",     ARG_INT, 1, NULL, offsetof(result_arg, id)},
	{"-option", ARG_STR, 1, NULL, offsetof(result_arg, option)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = SEQ_GENERIC;
    gen.result = NULL;

    if (strcmp(args.option, "zoom") == 0) {
	int zoom;

	gen.task = TASK_RASTER_ZOOM;
	seq_result_notify(args.id, (seq_reg_data *)&gen, 0);

	/* result = (char *)gen.result; */
	zoom = (int)gen.result;
	vTcl_SetResult(interp, "%d", zoom);

    } else if (strcmp(args.option, "number") == 0) {
	RasterResult *raster_result;

	raster_result = raster_id_to_result(args.id);
	if (raster_result) {
	    vTcl_SetResult(interp, "%d", raster_result->num_results);
	} else {
	    vTcl_SetResult(interp, "%d", 0);
	}
    } else {
	verror(ERR_WARN, "RasterResults", "option unknown \n");
    }

    return TCL_OK;

}

int
RulerTicks(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    ticks_arg args;
    double firstTick;
    double step;
    int numTicks;

    cli_args a[] = {
	{"-min",       ARG_DOUBLE, 1, "0.0", offsetof(ticks_arg, min)},
	{"-max",       ARG_DOUBLE, 1, "0.0", offsetof(ticks_arg, max)},
	{"-num_ticks", ARG_INT,    1, "0", offsetof(ticks_arg, num_ticks)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ruler_ticks(args.min, args.max, args.num_ticks, &firstTick, &step, &numTicks);

    if (step >= 1) {
	vTcl_SetResult(interp, "%g %d %d", firstTick, (int)step, numTicks);
    } else {
	vTcl_SetResult(interp, "%g %f %d", firstTick, step, numTicks);
    }
    return TCL_OK;
}

int
RasterMoveCursor(ClientData clientData, 
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[]) 
{
    raster_mcursor_arg args;
    Tcl_CmdInfo info;
    Tk_Raster *raster;

    cli_args a[]= {
	{"-id",    ARG_INT, 1, NULL, offsetof(raster_mcursor_arg, raster_id)},
	{"-raster",ARG_STR, 1, NULL, offsetof(raster_mcursor_arg, raster)},
	{"-pos",   ARG_INT, 1, NULL, offsetof(raster_mcursor_arg, pos)},
	{"-cursor_id", ARG_INT, 1, NULL, offsetof(raster_mcursor_arg, cursor_id)},
	{"-direction",ARG_INT, 1, "-1", offsetof(raster_mcursor_arg, direction)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    if (Tcl_GetCommandInfo(interp, args.raster, &info) == 0) 
	return TCL_ERROR;
    raster = (Tk_Raster*)info.clientData;
    if (args.direction == -1) 
	args.direction = HORIZONTAL;

    seq_raster_move_cursor(args.raster_id, raster, args.cursor_id, args.pos,
			   args.direction);
    return TCL_OK;
}

int
RasterFindEdCursor(ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, 
		   char *argv[]) 
{
    raster_fcursor_arg args;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    int cursor_id;
    int seq_id;

    cli_args a[]= {
	{"-id",    ARG_INT, 1, NULL, offsetof(raster_fcursor_arg, raster_id)},
	{"-raster",ARG_STR, 1, NULL, offsetof(raster_fcursor_arg, raster)},
	{"-pos",   ARG_INT, 1, NULL, offsetof(raster_fcursor_arg, pos)},
	{"-direction",ARG_INT, 1, "-1", offsetof(raster_fcursor_arg, direction)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    if (Tcl_GetCommandInfo(interp, args.raster, &info) == 0) 
	return TCL_ERROR;
    raster = (Tk_Raster*)info.clientData;

    if (args.direction == -1) 
	args.direction = HORIZONTAL;

    cursor_id = seq_raster_find_edcursor(args.raster_id, raster, args.pos,
					 args.direction, &seq_id);

#ifdef DEBUG
    printf("raster_find_edcursor cursor %d seq %d\n", cursor_id, seq_id);
#endif
    vTcl_SetResult(interp, "%d %d", cursor_id, seq_id);
    return TCL_OK;
}

/* set up tcl commands which call C procedures */
/*****************************************************************************/
/*                              RasterUtils_Init                             */
/*****************************************************************************/
int
RasterUtils_Init(Tcl_Interp *interp) {

    Tcl_CreateCommand(interp, "raster_replot_zoom", RasterReplotZoom, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_replot_all", RasterReplotAll, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_config", RasterConfig, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_getconfig", RasterGetConfig, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "update_raster_window", UpdateRasterWindow, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_results", RasterResults, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ruler_ticks", RulerTicks, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_move_cursor", RasterMoveCursor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_find_edcursor", 
		      RasterFindEdCursor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
}
