#if !defined(POSTSCRIPT_H)
#define POSTSCRIPT_H

#include <stdio.h>
#include <stdlib.h>
#include <io_lib/Read.h>

#include "tk.h"

typedef struct {
    /* Dimensions in pixels. */
    int		page_height, page_width;
    char	*orientation; /* portrait or landscape. */
    int		top_margin, bottom_margin, left_margin, right_margin;
    char	*font_type;
    int		font_size;
    int		panel_height, panel_width; /* Height and width of each printed segment. */
    int	      	panel_sep; /* Vertical separation between printed segments. */
    int	       	n_panel; /* The number of panels that will fit on the page. */
} PS_OPTIONS;

typedef struct {
    int		lw; /* Line width. */
    char	*colour_str;
    XColor	*XCol;
    float	rgb[3];
    char	*dash_str;
    int		*dash; /* Last element contains the offset. */
    int		n_dash;
} LINE_OPTIONS;

typedef struct {
    int	x, y;
} PS_POINT;

typedef struct {
    char	*text;
    int		x, y;
} PS_TEXT;

int	ps_configure(PS_OPTIONS *ps_options, int argc, char **argv);
int	free_ps_options(PS_OPTIONS *ps_options);
int	ps_configure_line(Tcl_Interp *interp, Tk_Window tkwin,
			  LINE_OPTIONS *options, int argc, char **argv);
int	free_ps_line(LINE_OPTIONS *options);
FILE   	*ps_fopen(char *filename, PS_OPTIONS ps_options);
void   	ps_newpage(FILE *ps_file, char *label, int page_number, PS_OPTIONS ps_options);
void   	ps_finishpage(FILE *ps_file);
void   	xfree_ps_text(PS_TEXT *ps_text, int n);
void   	ps_draw_lines(FILE *ps_file, LINE_OPTIONS options, PS_POINT *p, int NPoints);
void   	ps_draw_text(FILE *ps_file, PS_TEXT *t, int nt, float rgb[3], char align);
int	char_to_ps_text(PS_TEXT *ps_text, char c, int x_pos, int y_pos);
int 	int_to_ps_text(PS_TEXT *ps_text, int num, int x_pos, int y_pos);
int    	inch_to_pixel(double inch);
double 	pixel_to_inch(int pixel);

#endif

