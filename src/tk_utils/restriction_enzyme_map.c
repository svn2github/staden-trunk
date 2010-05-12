#include <stdio.h>
#include <tcl.h>
#include <string.h>

#include "renz_utils.h"
#include "canvas_box.h"
#include "text_output.h"
#include "xalloc.h"
#include "tcl_utils.h"
#include "misc.h"

void renz_shutdown(R_Enz *r_enzyme, 
		   int num_enzymes,
		   R_Match *match,
		   CanvasPtr *canvas,
		   WorldPtr *world,
		   StackPtr *zoom)
{
    int i, j;

    /* Free memory */
    if (r_enzyme) {
	for (i = 0; i < num_enzymes; i++) {
	    xfree(r_enzyme[i].name);
	    for (j = 0; j < r_enzyme[i].num_seq; j++) {
		xfree(r_enzyme[i].seq[j]);
	    }
	    xfree(r_enzyme[i].seq);
	    xfree(r_enzyme[i].cut_site);
	}
	xfree(r_enzyme);
    }

    xfree(match);

    if (canvas)
	xfree(canvas);
    
    if (world->visible)
	xfree(world->visible);

    if (world->total)
	xfree (world->total);
    
    if (world)
	xfree(world);

    freeZoom(&zoom);
}

/*
 * tcl interface to plotting lines on restriction enzyme display
 * stand alone and template displays
 */
void
PlotStickMap(Tcl_Interp *interp,
	     char *win_name,
	     int unpadded_cut_site,
	     int cut_site,
	     int xoffset,
	     int yoffset,
	     int tick_ht,
	     int tick_wd,
	     char *colour,
	     int index,
	     int from,
	     int to)
{
    char cmd[1024];
    int xcut;

    /* Force cut sites off the end of the contig to be displayed. */
    xcut = MAX(cut_site, from);
    xcut = MIN(to, xcut);
    
    sprintf(cmd,"%s create line %d %d %d %d -fill %s -width %d "
	    "-tag {S renz re_%d pos_%d}",
	    win_name, xcut+xoffset, yoffset, xcut+xoffset, 
	    yoffset + tick_ht, colour, tick_wd, index, unpadded_cut_site);
    Tcl_Eval(interp, cmd);
}

void
plot_renz_matches(Tcl_Interp *interp,
		  char *re_win,
		  char *names_win,
		  int text_offset,
		  char *t_colour,
		  int yoffset,
		  int num_enzymes,
		  R_Enz *r_enzyme,
		  ruler_s *ruler,
		  int sequence_len,
		  int num_matches,
		  R_Match *match,
		  tick_s *tick,
		  char *frame,
		  WorldPtr *world,
		  CanvasPtr *canvas,
		  win **win_list,
		  int num_wins,
		  StackPtr **zoom)
{
    char cmd[1024];
    int item;
    int offset;
    int y_offset;
    int t_offset;
    int i;
    int cut_site;

    sprintf(cmd, "%s delete all", re_win);
    Tcl_Eval(interp, cmd);

    sprintf(cmd, "%s delete all", names_win);
    Tcl_Eval(interp, cmd);

    t_offset = text_offset;
    y_offset = yoffset;

    /* for each selected enzyme, display its name and cut sites */
    for (item = 0; item < num_enzymes; item++) {

        offset = yoffset;
	sprintf(cmd,"%s create text 10 %d -text %s -anchor w -fill %s "
		"-font enzyme_font -tag {S re_%d}",
		names_win, t_offset, r_enzyme[item].name, t_colour, 
		item);
	Tcl_Eval(interp, cmd);
	sprintf(cmd, "%s create line %d %d %d %d -tag contig -fill %s",
		re_win, ruler->start, y_offset, ruler->end, y_offset, ruler->colour);
	Tcl_Eval(interp, cmd);

	for (i = 0; i < num_matches; i++) {
	    if (item == match[i].enz_name) {
	
		offset = (item * tick->ht) + yoffset;
		/* cut_site = match[i].cut_pos; */
		cut_site = ruler->start - 1 + match[i].padded_cut_pos;

		PlotStickMap(interp, re_win,
			     ruler->start - 1 + match[i].cut_pos, cut_site,
			     0, offset, 
			     tick->ht, tick->line_width, tick->colour, item,
			     ruler->start, ruler->end);
	    }
	}
	y_offset += tick->ht;
	t_offset += tick->ht;
    }
    /* draw final line */
    sprintf(cmd, "%s create line %d %d %d %d -tag contig -fill %s",
	    re_win, ruler->start, y_offset, ruler->end, y_offset, ruler->colour);
    Tcl_Eval(interp, cmd);

    /* add back selection rectangles */
    if (TCL_ERROR == Tcl_VarEval(interp, "ReSelectRect ", frame, " ",
				 names_win, NULL))
	verror(ERR_WARN, "plot_renz_matches", "%s\n",
	       Tcl_GetStringResult(interp));

    world->total->x1 = ruler->start;
    world->total->x2 = ruler->end;
    world->total->y1 = 1;
    world->total->y2 = y_offset;

    memcpy(world->visible, world->total, sizeof(d_box));
    world->visible->y2 = canvas->height;

    SetCanvasCoords(interp, world->visible->x1, world->visible->y1, 
		    world->visible->x2, world->visible->y2, canvas);

    draw_single_ruler(interp, ruler, canvas, ruler->start, ruler->end, 1);
    scaleCanvas(interp, win_list, num_wins, "all", world->visible, 
		canvas);
    scrollRegion(interp, win_list, num_wins, world->total, canvas);

    /* remove all current zooming info */
    freeZoom(zoom);

    /* add first zoom */
    pushZoom(zoom, world->visible);
}

