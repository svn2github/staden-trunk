#include <tk.h>
#include "IO.h"
#include "newgap_cmds.h"
#include "template_display.h"
#include "gap_globals.h"
#include "misc.h"
#include "qual.h"

/* frame name for tk */
#define FRAME_LEN 30

void
glevel(char current,
       float Y0,
       float YP1,
       float YP2,
       float YM1,
       float YM2,
       float *yfrom,
       float *yto)
{
    switch(current) {
    case R_GOOD_GOOD_EQ:
        *yfrom = Y0;
        *yto = Y0;
	break;
    case R_GOOD_NONE:
        *yfrom = Y0;
        *yto = YM1;
	break;
    case R_NONE_GOOD:
        *yfrom = Y0;
        *yto = YP1;
	break;
    case R_BAD_BAD:
        *yfrom = YP1;
        *yto = YM1;
	break;
    case R_GOOD_GOOD_NE:
        *yfrom = YP2;
        *yto = YM2;
	break;
    case R_GOOD_BAD:
        *yfrom = Y0;
        *yto = YM1;
	break;
    case R_BAD_GOOD:
        *yfrom = Y0;
        *yto = YP1;
	break;
    case R_BAD_NONE:
        *yfrom = YP1;
        *yto = YM1;
	break;
    case R_NONE_BAD:
        *yfrom = YP1;
        *yto = YM1;
	break;
    case R_NONE_NONE:
	break;
        *yfrom = YP1;
        *yto = YM1;
    default:
        verror(ERR_FATAL, "quality_plot","incorrect value to glevel()");
    }
} /* end glevel */

/*
 * Removes a quality plot from the item display
 */
void plot_quality_shutdown(char *win_name) {
    Tcl_VarEval(GetInterp(), "DeleteQualPlot ", win_name, NULL);
}

/*
 * each code has its own level and we draw rectangles for segments of
 * the same code
 * The coordinate system is xmin to xmax, and ymin to ymax which here
 * are set to 1, sequence_length; -2, 2
 */
void
plot_quality(char *seq,
	     int seq_len,
             char *win_name,
	     GapIO *io,
	     int init)
{
    float Y0  = 0.;
    float YP1 = 1.;
    float YP2 = 2.;
    float YM1 = -1.;
    float YM2 = -2.;
    int xmin = 0;
    int xmax = seq_len;
    float ymin = YM2;
    float ymax = YP2;
    int i;
    int xfrom, xto;
    float yfrom, yto;
    char current;
    char cmd[1024];

    char q_name[FRAME_LEN];
    char *scroll = "x";

    sprintf(q_name, "%s.quality", win_name);

    if (init) {
	/* add a new canvas to current view */
	/* view = cur_view; */
	sprintf(cmd, "CreateQualPlot %d %s %d %d", *handle_io(io), win_name,
		xmin, xmax);
	Tcl_Eval(GetInterp(), cmd);
	sprintf(cmd, "%s,disp_quality", win_name);
    } else {
	/* Delete the old plot */
 	Tcl_VarEval(GetInterp(), q_name, " delete all", NULL);
    }

    i = 1;
    current = seq[i];
    xfrom = i;
    while (i < seq_len) {
        if (current != seq[i]) {
            /* printf("not equal cur %c seq(%d) %c \n", current, i, seq[i]); */
            glevel(current, Y0, YP1, YP2, YM1, YM2, &yfrom, &yto);
	    xto = i /* - 1 */;

	    sprintf(cmd,"plot_rec %s %d %.20f %d %.20f "
		    "%.20f %.20f %.20f %.20f %.20f %.20f %.20f",
		    q_name, xfrom, yfrom, xto, yto,
		    Y0, YP1, YP2, YM1, YM2,
		    ymin, ymax);
#ifdef DEBUG
	    printf("%s %s \n", GetInterp()->result, cmd);
#endif
	    Tcl_Eval(GetInterp(), cmd);
            current = seq[i];
	    xfrom = i;
        }
        i++;
    }
    /* we have reached the end so finish off last rectangle */
    glevel(current, Y0, YP1, YP2, YM1, YM2, &yfrom, &yto);
    xto = i;
    sprintf(cmd,"plot_rec %s %d %.20f %d %.20f "
            "%.20f %.20f %.20f %.20f %.20f %.20f %.20f",
            q_name, xfrom, yfrom, xto, yto,
            Y0, YP1, YP2, YM1, YM2,
            ymin, ymax);
    Tcl_Eval(GetInterp(), cmd);

    /* dreadful hack necessary because C doesn't know about the Consts that
     * tcl knows about
     */
    sprintf(cmd, "SetQualityCanvas %s %s %s %d %d", q_name, win_name, scroll, xmin, xmax);
    Tcl_Eval(GetInterp(), cmd);

}
