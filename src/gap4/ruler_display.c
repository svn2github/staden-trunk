#include <tk.h>
#include <math.h>

#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "misc.h" 
#include "IO.h"
#include "gap_globals.h"
#include "canvas_box.h"
#include "ruler_display.h"
#include "template_display.h"
#include "ruler_tick.h"

int display_ruler(Tcl_Interp *interp, 
		  GapIO *io,
		  CanvasPtr *canvas,    
		  c_offset *contig_offset,
		  int *contig_array,
		  int num_contigs,
		  int disp_ruler,
		  int disp_ticks,
		  ruler_s *ruler,
		  char *frame,
		  PlotRec **ruler_coord)
{
    char cmd[1024];
    int i;
    PlotRec *CArray;
    int depth;

    if (!disp_ruler) {
	return 0;
    }

    if (NULL == (CArray = (PlotRec *)xmalloc(num_contigs * sizeof(PlotRec)))) {
	return -1;
    }

    /* remove current ruler before drawing the new ruler */
    Tcl_VarEval(interp, ruler->window, " delete contig", NULL);
    Tcl_VarEval(interp, ruler->window, " delete tag", NULL);
    Tcl_VarEval(interp, ruler->window, " delete tick", NULL);
    /*
     * set up CArray 
     */
    for (i = 0; i < num_contigs; i++) {
	CArray[i].num = contig_array[i];

	CArray[i].l.x1 = 1 + contig_offset[contig_array[i]].offset;
	CArray[i].l.x2 = io_clength(io, contig_array[i]) 
	    + contig_offset[contig_array[i]].offset;

	CArray[i].colour = ruler->colour;

	if (NULL == (CArray[i].type = (char *)xmalloc(40))) {
	    verror(ERR_FATAL, "display_ruler", "out of memory");
	    return -1;
	}		

	sprintf(CArray[i].type, "{contig c_%d num_%d hl_%d S}", 
		i+1, contig_array[i], contig_array[i]);
	strcpy(CArray[i].arrow, "none");
    }

    CalcYDepth(num_contigs, CArray, num_contigs, &depth);

    for (i = 0; i < num_contigs; i++) {
	CArray[i].l.y1 = CArray[i].l.y1 * ruler->offset;
	CArray[i].l.y2 = CArray[i].l.y2 * ruler->offset;
    }
#ifdef DEBUG
    for (i = 0; i < num_contigs; i++) {
	printf("%d %d %d %d %s\n", CArray[i].l.x1, CArray[i].l.y1,
	       CArray[i].l.x2, CArray[i].l.y2, CArray[i].colour);
    }
#endif

    plot_lines(interp, CArray, num_contigs, ruler->window, ruler->line_width);

    *ruler_coord = CArray;

    if (!disp_ticks) {
	sprintf(cmd, "RulerWindowSize %d %s %s ", disp_ticks, frame,
		ruler->window);
	Tcl_Eval(interp, cmd);

	return(0);
    }

    for (i = 0; i < num_contigs; i++) {
      display_ruler_ticks(interp, canvas, 
			    contig_offset[contig_array[i]].offset, 
			    CArray[i].l.y1, ruler, 1, CArray[i].l.x2 -
			    CArray[i].l.x1 + 1);
    }

    sprintf(cmd, "RulerWindowSize %d %s %s ", 1, frame, ruler->window);
    Tcl_Eval(interp, cmd);

    return(0);
}

/*
 * draw a vertical ruler
 */
void display_ruler_v(Tcl_Interp *interp,
		     CanvasPtr *canvas,
		     ruler_s *ruler,
		     double wy1,
		     double wy2)
{
    char cmd[1024];

    Tcl_VarEval(interp, ruler->window, " delete all", NULL);
    sprintf(cmd, "%s create line %d %f %d %f\n", ruler->window, ruler->offset,
	    wy1, ruler->offset, wy2);
    Tcl_Eval(interp, cmd);

    display_ruler_ticks_v(interp, canvas, ruler, wy1, wy2);
}

