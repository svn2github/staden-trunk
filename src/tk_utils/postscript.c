#include <staden_config.h>

#include <string.h>
#include <ctype.h>
#include "misc.h"
#include "postscript.h"
#include "cli_arg.h"
#include "split.h"

int ps_configure(PS_OPTIONS *ps_options, int argc, char **argv) {
    cli_args	a[] = {
	{"-page_height",      ARG_INT, 1, "842",       offsetof(PS_OPTIONS, page_height)},
	{"-page_width",       ARG_INT, 1, "595",       offsetof(PS_OPTIONS, page_width)},
	{"-orientation",      ARG_STR, 1, "Portrait",  offsetof(PS_OPTIONS, orientation)},
	{"-top_margin",       ARG_INT, 1, "72",        offsetof(PS_OPTIONS, top_margin)},
	{"-bottom_margin",    ARG_INT, 1, "72",        offsetof(PS_OPTIONS, bottom_margin)},
	{"-left_margin",      ARG_INT, 1, "72",        offsetof(PS_OPTIONS, left_margin)},
	{"-right_margin",     ARG_INT, 1, "72",        offsetof(PS_OPTIONS, right_margin)},
	{"-panel_height",     ARG_INT, 1, "103",       offsetof(PS_OPTIONS, panel_height)},
	{"-panel_width",      ARG_INT, 1, "453",       offsetof(PS_OPTIONS, panel_width)},
	{"-panel_separation", ARG_INT, 1, "13",        offsetof(PS_OPTIONS, panel_sep)},
	{"-n_panel",          ARG_INT, 1, "6",         offsetof(PS_OPTIONS, n_panel)},
	{"-font",             ARG_STR, 1, "Helvetica", offsetof(PS_OPTIONS, font_type)},
	{"-font_size",        ARG_INT, 1, "12",        offsetof(PS_OPTIONS, font_size)},
	{NULL,                0,       0, NULL,        0}
    };

    if(-1 == parse_args(a, ps_options, argc, argv)) {
	return TCL_ERROR;
    }

    ps_options->orientation = strdup(ps_options->orientation);
    ps_options->font_type   = strdup(ps_options->font_type);

    return TCL_OK;
}

int free_ps_options(PS_OPTIONS *ps_options) {
    if(ps_options->orientation != NULL) {
	free(ps_options->orientation);
	ps_options->orientation = NULL;
    }
    if(ps_options->font_type != NULL) {
	free(ps_options->font_type);
	ps_options->font_type = NULL;
    }
    return 0;
}

int ps_configure_line(Tcl_Interp *interp, Tk_Window tkwin,
		      LINE_OPTIONS *options, int argc, char **argv) {
    char		**dash_list;
    int			i;

    cli_args	a[] = {
	{"-line_width", ARG_INT, 1, "0",     offsetof(LINE_OPTIONS, lw)},
	{"-colour",     ARG_STR, 1, "black", offsetof(LINE_OPTIONS, colour_str)},
	{"-dash",       ARG_STR, 1, "0",     offsetof(LINE_OPTIONS, dash_str)},
	{NULL,          0,       0, NULL,    0}
    };

    if(-1 == parse_args(a, options, argc, argv)) {
	return TCL_ERROR;
    }

    /* Convert *colour_str to rgb[3] */
    options->XCol = Tk_GetColor(interp, tkwin, options->colour_str);
    options->rgb[0] = ((float) options->XCol->red)   / 65535.0;
    options->rgb[1] = ((float) options->XCol->green) / 65535.0;
    options->rgb[2] = ((float) options->XCol->blue)  / 65535.0;

    /* Convert dash_str into an integer array. */
    dash_list = split(options->dash_str, " \t\n");
    if(NULL == (options->dash = (int *) xmalloc(strlen(options->dash_str) * sizeof(int)))) {
	return TCL_ERROR;
    }
    i = 0;
    while(dash_list[i]) {
	options->dash[i] = atoi(dash_list[i]);
	i++;
    }
    options->n_dash = i;
    if(NULL == (options->dash = (int *) xrealloc(options->dash, options->n_dash * sizeof(int) + 1))) {
	return TCL_ERROR;
    }
    split_xfree(dash_list);

    return TCL_OK;
}

int free_ps_line(LINE_OPTIONS *options) {
    if(options->XCol != NULL) {
	Tk_FreeColor(options->XCol);
	options->XCol = NULL;
    }
    if(options->dash != NULL) {
	xfree(options->dash);
	options->dash = NULL;
    }
    return 0;
}

FILE *ps_fopen(char *filename, PS_OPTIONS ps_options) {
    FILE	*ps_file;

    ps_file = fopen(filename, "w");
    if(ps_file == NULL) {
	return NULL;
    }

    fprintf(ps_file, "%%!PS-Adobe-3.0\n");
    fprintf(ps_file, "%%%%Creator:\ttrace_print\n");
    switch(tolower(*(ps_options.orientation))) {
    case 'l' :
	fprintf(ps_file, "%%%%Orientation:\tLandscape\n");
	break;
    case 'p' : default:
	fprintf(ps_file, "%%%%Orientation:\tPortrait\n");
	break;
    }
    fprintf(ps_file, "%%%%BeginProlog\n");
    fprintf(ps_file, "/t {translate} def\n");
    fprintf(ps_file, "/r {rotate} def\n");
    fprintf(ps_file, "/m {moveto} def\n");
    fprintf(ps_file, "/rm {rmoveto} def\n");
    fprintf(ps_file, "/l {lineto} def\n");
    fprintf(ps_file, "/rl {rlineto} def\n");
    fprintf(ps_file, "/s {show} def\n");
    fprintf(ps_file, "/rgb {setrgbcolor} def\n");
    fprintf(ps_file, "/lw {setlinewidth} def\n");
    fprintf(ps_file, "/st {stroke} def\n");
    fprintf(ps_file, "/n {newpath} def\n");
    fprintf(ps_file, "/rep {repeat} def\n");
    fprintf(ps_file, "/dash {setdash} def\n");
    /* The next four work on the last string on the postscript stack. */
    fprintf(ps_file, "/ln {stringwidth} def\n");
    fprintf(ps_file, "/l_h {ln exch -0.5 mul exch rm} def\n"); /* Move left by half the stringwidth. */
    fprintf(ps_file, "/l_f {ln exch -1 mul exch rm} def\n"); /* Move left by a full stringwidth. */
    fprintf(ps_file, "/r_h {ln exch 0.5 mul exch rm} def\n"); /* Move right by half the stringwidth. */
    fprintf(ps_file, "%%%%EndProlog\n");
    fprintf(ps_file, "%%%%BeginSetup\n");
    /*
    fprintf(ps_file, "<</PageSize [%d %d]>> setpagedevice\n",
	    ps_options.page_width, ps_options.page_height);
    */
    fprintf(ps_file, "/%s findfont %d scalefont setfont\n", ps_options.font_type, ps_options.font_size);
    fprintf(ps_file, "%%%%EndSetup\n");

    return ps_file;
}

void ps_newpage(FILE *ps_file, char *label, int page_number, PS_OPTIONS ps_options) {
    fprintf(ps_file, "%%%%Page: %s %d\n", label, page_number);
    fprintf(ps_file, "%%%%BeginPageSetup\n");
    switch(tolower(*(ps_options.orientation))) {
    case 'l' :
	fprintf(ps_file, "90 r 0 %d t\n", (-1 * ps_options.page_height));
	break;
    case 'p' : default:
	break;
    }
    /* Translate to top left-hand corner. */
    fprintf(ps_file, "%d %d t\n", ps_options.left_margin, ps_options.page_height - ps_options.top_margin);
    fprintf(ps_file, "%%%%EndPageSetup\n");

    fprintf(ps_file, "0 0 m\n");
    fprintf(ps_file, "(%s) s\n", label);
}

void ps_finishpage(FILE *ps_file) {
    fprintf(ps_file, "showpage\n");
}

void xfree_ps_text(PS_TEXT *ps_text, int n) {
    int	i;

    for(i = 0; i < n; i++) {
	xfree(ps_text[i].text);
    }
    xfree(ps_text);
}

void ps_draw_lines(FILE *ps_file, LINE_OPTIONS options, PS_POINT *p, int NPoints) {
    int	i;

    fprintf(ps_file, "n\n");

    fprintf(ps_file, "%d %d m\n", p[0].x, p[0].y);
    for(i = NPoints - 1; i >= 1; i--) {
	fprintf(ps_file, "%d %d\n", p[i].x - p[i - 1].x, p[i].y - p[i - 1].y);
    }
    fprintf(ps_file, "%d {rl} rep\n", NPoints - 1);

    fprintf(ps_file, "%d lw\n", options.lw);
    fprintf(ps_file, "%04.2f %04.2f %04.2f rgb\n", options.rgb[0], options.rgb[1], options.rgb[2]);

    fprintf(ps_file, "[");
    for(i = 0; i < options.n_dash - 1; i++) {
	fprintf(ps_file, "%d ", options.dash[i]);
	i++;
    }
    fprintf(ps_file, "] %d dash\n", options.dash[i]);

    fprintf(ps_file, "st\n");
}

void ps_draw_text(FILE *ps_file, PS_TEXT *t, int nt, float rgb[3], char align) {
    int	i, len;

    if(rgb[0] != -1.0) {
	fprintf(ps_file, "%04.2f %04.2f %04.2f rgb\n", rgb[0], rgb[1], rgb[2]);
    }
    for(i = 0; i < nt; i++) {
	fprintf(ps_file, "%d %d m\n", t[i].x, t[i].y);
	switch(align) {
	case 'c':		/* Center */
	    fprintf(ps_file, "(%s) l_h\n", t[i].text);
	    break;
	case 'r':		/* Right */
	    fprintf(ps_file, "(%s) l_f", t[i].text);
	    break;
	case 'f':		/* Center on first character of string. */
	    len = strlen(t[i].text);
	    fprintf(ps_file, "(%c) l_h\n", t[i].text[len - 1]);
	    break;
	case 'e':		/* Center on last character of string. */
	    fprintf(ps_file, "(%s) l_f\n", t[i].text);
	    len = strlen(t[i].text);
	    fprintf(ps_file, "(%c) r_h\n", t[i].text[len - 1]);
	    break;
	case 'l': default:	/* Left. This is the default in postscript. */
	    break;
	}

	fprintf(ps_file, "(%s) s\n", t[i].text);
    }
}

int char_to_ps_text(PS_TEXT *ps_text, char c, int x_pos, int y_pos) {
    ps_text->text = (char *) xmalloc(2 * sizeof(char));
    if(NULL == ps_text->text) {
	return -1;
    }

    sprintf(ps_text->text, "%c", c);
    ps_text->x = x_pos;
    ps_text->y = y_pos;

    return 0;
}

int int_to_ps_text(PS_TEXT *ps_text, int num, int x_pos, int y_pos) {
    ps_text->text = (char *) xmalloc(30 * sizeof(char));
    if(NULL == ps_text->text) {
	return -1;
    }

    sprintf(ps_text->text, "%d", num);
    ps_text->x = x_pos;
    ps_text->y = y_pos;

    return 0;
}

int inch_to_pixel(double inch) {
    return ((int) (72.0 * inch));
}

double pixel_to_inch(int pixel) {
    return (((double) pixel) / 72.0);
}
