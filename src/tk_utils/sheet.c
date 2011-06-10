/*
 * The Sheet widget, used internally by the Editor widget and by the sheet tk
 * widget.
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <string.h>

/* Function prototypes missing on windows implementation of Xlib.h, tk.h pulls them in */
#ifdef _WIN32
#include <tk.h>
#endif
#include <X11/Xlib.h>

#include "xalloc.h"
#include "sheet.h"
#include "misc.h"

/* ---- Local prototypes ---- */

static void redisplayRegion(Sheet *sw, XRectangle *expose);
static void destroy_array(sheet_array a);
static sheet_array create_array (int r, int c, size_t s);
static void move_array (sheet_array a, sheet_array b);
static void extend_array (sheet_array *a, int r, int c);
/* static char *get_array_element(sheet_array a, int r, int c); */

/* ---- Local defines ---- */
#define FONT_WIDTH(W) ((W)->font_width)
#define FONT_HEIGHT(W) ((W)->fm.linespace)

#define GET_ARRAY_CELL(A,R,C)\
    ( &A->base[(R * A->cols + C)*A->size] )

/* ---- External functions ---- */

/*
 * Clears a sheet to spaces and default inks. Used when creating or
 * resizing a sheet.
 */
void sheet_clear(Sheet *sw) {
    int c, r;
    sheet_ink ink_base;
    sheet_paper paper_base;

    for (r = 0; r < sw->rows; r++) {
	ink_base = (sheet_ink)GET_ARRAY_CELL(sw->ink, r, 0);
	paper_base = (sheet_paper)GET_ARRAY_CELL(sw->paper, r, 0);

	memset(paper_base, ' ', sw->columns);

	for (c = 0; c < sw->columns; c++, ink_base++)
	    ink_base->sh = sh_default;
    }
}

/*
 * Destroys a sheet (deallocates memory allocated by it, but does not
 * deallocate the sheet itself).
 */
void sheet_destroy(Sheet *sw) {
    if (sw->paper)
	destroy_array(sw->paper);
    if (sw->ink)
	destroy_array(sw->ink);
    if (sw->dbl_buffer)
	Tk_FreePixmap(sw->display, sw->dbl_buffer);
}

/*
 * Updates the local arrays used in a sheet widget.
 * NOTE: The new size is already containing in the Sheet pointer, and
 * the old size is passed as arguments.
 */
void sheet_resize(Sheet *sw, int rows, int columns) {
    if (!sw->rows || !sw->columns) {
	return;
    }

    if (rows != sw->rows || columns != sw->columns) {
	if (NULL == sw->paper)
	    sw->paper = create_array(sw->rows, sw->columns,
				     sizeof(sheet_paper_struct));
	else
	    extend_array (&sw->paper, sw->rows, sw->columns);

	if (NULL == sw->ink)
	    sw->ink = create_array(sw->rows, sw->columns,
				   sizeof(sheet_ink_struct));
	else
	    extend_array (&sw->ink,   sw->rows, sw->columns);

	sheet_clear(sw);

	if (sw->dbl_buffer)
	    Tk_FreePixmap(sw->display, sw->dbl_buffer);
	if (Tk_IsMapped(sw->tkwin))
	    sw->dbl_buffer = Tk_GetPixmap(sw->display, Tk_WindowId(sw->tkwin),
					  sw->width_in_pixels,
					  sw->height_in_pixels,
					  Tk_Depth(sw->tkwin));
	else
	    sw->dbl_buffer = (Pixmap)NULL;
    }
}

/*
 * (re)configures a sheet. This contains such things are allocating foreground
 * and background colours.
 */
void sheet_config(Sheet *sw, Pixel light, Pixel fg, Pixel bg, Pixel ifg) {
    unsigned long valuemask;
    XGCValues values;

    sw->light = light;
    sw->foreground = fg;
    sw->background = bg;
    sw->indel_foreground = ifg;

    /*
     * Change GC colours
     */
    valuemask = GCFont | GCGraphicsExposures | GCForeground | GCBackground;
    values.graphics_exposures = False;

    values.foreground = fg;
    values.background = bg;
    values.font = Tk_FontId(sw->bold_font);
    sw->normgcB = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.font = Tk_FontId(sw->font);
    sw->normgc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = bg;
    values.background = bg;
    sw->sparegc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = light;
    values.background = bg;
    sw->greygc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = bg;
    values.background = fg;
    sw->whitegc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = ifg;
    values.background = bg;
    sw->indelgc = Tk_GetGC(sw->tkwin, valuemask, &values);
}


/*
 * Initialises a sheet, returns -1 for error
 */
int sheet_create(Sheet *sw, Pixel light, Pixel fg, Pixel bg, Pixel ifg) {
    unsigned long valuemask;
    XGCValues values;


    sw->paper = NULL;
    sw->ink = NULL;
    sw->light = light;
    sw->foreground = fg;
    sw->background = bg;
    sw->indel_foreground = ifg;
    sw->cursor_row = sw->cursor_column = -1;
    sw->display_cursor = 1;
    sw->window = 0;
    sw->yflip = 0;
    sw->dbl_buffer = 0;
    sw->hollow_cursor = 0;

    sheet_resize(sw, 0, 0);

    /*
     * GCs and things
     */
    valuemask = GCFont | GCGraphicsExposures | GCForeground | GCBackground ;
    values.graphics_exposures = False;

    values.foreground = sw->foreground;
    values.background = sw->background;
    values.font = Tk_FontId(sw->bold_font);
    sw->normgcB = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.font = Tk_FontId(sw->font);
    sw->normgc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = sw->indel_foreground;
    values.background = sw->background;
    sw->indelgc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = sw->background;
    values.background = sw->background;
    sw->sparegc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = sw->foreground;
    values.background = sw->background;
    sw->greygc = Tk_GetGC(sw->tkwin, valuemask, &values);

    values.foreground = sw->background;
    values.background = sw->foreground;
    sw->whitegc = Tk_GetGC(sw->tkwin, valuemask, &values);

#define grey_width 2
#define grey_height 2
    if (DisplayPlanes(sw->display, DefaultScreen(sw->display)) == 1) {
	static char grey_bits[] = { 0x01 , 0x02 };

	sw->grey_stipple =
	    XCreateBitmapFromData(sw->display,
				  RootWindowOfScreen(Tk_Screen(sw->tkwin)),
				  grey_bits,
				  grey_width,
				  grey_height);
	XSetFillStyle(sw->display, sw->greygc, FillOpaqueStippled);
	XSetStipple(sw->display, sw->greygc, sw->grey_stipple);
#undef grey_width
#undef grey_height
    } /* else already allocated by X default mechanism */

    return 0;
}

void sheet_display(Sheet *sw) {
    XRectangle expose;

    if (!Tk_IsMapped(sw->tkwin))
	return;

    sw->window = Tk_WindowId(sw->tkwin);

    expose.x = 0;
    expose.y = 0;
    expose.width  = sw->width_in_pixels;
    expose.height = sw->height_in_pixels;

    redisplayRegion(sw, &expose);
}

/* ---- Local functions ---- */
int binary_op(int src, int dst, int op)
{
    switch (op & HOP_MASK) {
	case 0: return 0;
	case 1: return ~ (src | dst);
	case 2: return ~src & dst;
	case 3: return ~src;
	case 4: return src & ~dst;
	case 5: return ~dst;
	case 6: return src ^ dst;
	case 7: return ~(src & dst);
	case 8: return src & dst;
	case 9: return ~src ^ ~dst;
	case 10: return dst;
	case 11: return ~src|dst;
	case 12: return src;
	case 13: return src|~dst;
	case 14: return src|dst;
	case 15: return 1;
    }

    /* Default case should never occur */
    return 0;
}

static void destroy_array(sheet_array a)
{
    xfree (a->base);
    xfree ((char *)a);
}

static sheet_array create_array (int r, int c, size_t s)
{
    sheet_array new = (sheet_array) xcalloc (1,sizeof(sheet_array_struct));
    if (new != NULL) {
	new->base = (char *) xcalloc (r*c,s);
	if (new->base == NULL) {
	    xfree ((char *)new);
	    new = NULL;
	} else {
	    new->rows = r;
	    new->cols = c;
	    new->size = s;
	}
    }
    return new;
}

static void move_array (sheet_array a, sheet_array b)
{
    size_t r,c;
    unsigned int i;

    c = min (a->cols*a->size, b->cols*b->size);
    r = min (a->rows, b->rows);
    for (i=0; i<r; i++)
	memcpy(
	    (char *) GET_ARRAY_CELL(b,i,0),
	    (char *) GET_ARRAY_CELL(a,i,0),
	    c
	);

}

static void extend_array (sheet_array *a, int r, int c)
/*
** Extending stategy: For rows.
** Do rows really need extending?
** Yes:
**    Will twice old rows do?
**    Yes:
**       Use twice old rows.
**    No:
**	 Use new rows plus EXTEND_ROWS_GUESS (A wild quess that's enough)
** No:
**    Don't extend
**
** Extending stategy: For columns.
** As for rows.
*/
{
    int newr;
    int newc;
    sheet_array b;

#define EXTEND_ROWS_GUESS 5
#define EXTEND_COLS_GUESS 5
    newr = (r<=(*a)->rows)
	?(*a)->rows
	:(r<=(*a)->rows*2)
	    ?(*a)->rows*2
	    :r+EXTEND_ROWS_GUESS;

    newc = (c<=(*a)->cols)
	?(*a)->cols
	:(c<=(*a)->cols*2)
	    ?(*a)->cols*2
	    :c+EXTEND_COLS_GUESS;

    if (newr!=(*a)->rows || newc!=(*a)->cols) {
	b = create_array(newr,newc,(*a)->size);
	move_array(*a,b);
	destroy_array(*a);
	*a = b;
    }
}

/*
static char * get_array_element(sheet_array a, int r, int c)
{
    if (r < 0 || c < 0)
	return NULL;
    if (a->rows > r || a->cols > c)
	return NULL;

    return GET_ARRAY_CELL(a,r,c);
}
*/

static GC setGC(Sheet *sw, sheet_ink ink_base)
{
    unsigned long valuemask;
    XGCValues values;

    valuemask = GCFont | GCGraphicsExposures | GCForeground | GCBackground ;
    values.font = ink_base->sh & sh_bold
	? Tk_FontId(sw->bold_font)
	: Tk_FontId(sw->font);
    values.graphics_exposures = False;

    if (ink_base->sh & sh_inverse) {
	if (ink_base->sh&sh_bg)
	    values.foreground = ink_base->bg;
	else
	    values.foreground = sw->background;
	if (ink_base->sh&sh_fg)
	    values.background = ink_base->fg;
	else
	    values.background = sw->foreground;
    } else {
	if (ink_base->sh&sh_fg)
	    values.foreground = ink_base->fg;
	else
	    values.foreground = sw->foreground;
	if (ink_base->sh&sh_bg)
	    values.background = ink_base->bg;
	else
	    values.background = sw->background;
    }

    /* FIXME:
     * We never Tk_FreeGC so the reference count increases indefinitely
     */
    return Tk_GetGC(sw->tkwin, valuemask, &values);
}


static void _repaint_colour(Sheet *sw, int c, int r, int l, sheet_ink ink, char *s)
{
    sheet_ink_struct my_ink;
    GC bg_gc;
    XGCValues values;
    unsigned long mask;
    unsigned long bg;

    my_ink = *ink;
    if (ink->sh & sh_light) {
	my_ink.sh = (my_ink.sh | sh_fg) /* & ~sh_bg */ ;
	my_ink.fg = sw->light;
    }

    if (ink->sh & sh_inverse)
	bg = (ink->sh & sh_fg) ? ink->fg : sw->foreground;
    else
	bg = (ink->sh & sh_bg) ? ink->bg : sw->background;

    /* Create or change the GC for filling the background */
    values.foreground = bg;
    values.function = GXcopy;
    values.fill_style = FillSolid;
    mask = GCForeground  | GCFunction | GCFillStyle;
    bg_gc = Tk_GetGC(sw->tkwin, mask, &values);

#if 0
    /* fill the background */
    XFillRectangle(sw->display,
		   sw->window,
		   bg_gc,
		   (int) COL_TO_PIXEL(sw,c),
		   (int) ROW_TO_PIXEL(sw,r-1),
		   FONT_WIDTH(sw) * l,
		   FONT_HEIGHT(sw)
		   );
    Tk_FreeGC(sw->display, bg_gc);

    /* finally draw the string */
    sw->sparegc = setGC(sw, &my_ink);
    Tk_DrawChars(sw->display,sw->window,sw->sparegc,sw->font,s,l,
		 (int) COL_TO_PIXEL(sw,c),(int) ROW_TO_BASELINE_PIXEL(sw,r));

    /* and underline if required */
    if (ink->sh & sh_select || ink->sh & sh_underline) {
	XDrawLine(sw->display,
		  sw->window,
		  sw->sparegc,
		  (int) COL_TO_PIXEL(sw,c),
		  (int) ROW_TO_BASELINE_PIXEL(sw,r),
		  (int) COL_TO_PIXEL(sw,c+l)-1,
		  (int) ROW_TO_BASELINE_PIXEL(sw,r)
		  );
    }

    /* add ruler tick if required */
    if (ink->sh & sh_tick) {
	int i;
	for (i = 0; i < l; i++) {
	    XDrawRectangle(sw->display,
			   sw->window,
			   sw->sparegc,
			   (int) COL_TO_PIXEL(sw,c+i) + COL_TO_PIXEL(sw,1)/2-1,
			   (int) ROW_TO_BASELINE_PIXEL(sw,r)+1,
			   1, 2);
	}
    }
#else
 {
     /* Double buffered */
     GC copygc;

     Pixmap p = sw->dbl_buffer;
     mask = GCFunction|GCGraphicsExposures|GCForeground;
     values.function = GXcopy;
     values.graphics_exposures = False;
     values.foreground = sw->foreground;
     copygc = Tk_GetGC(sw->tkwin, mask, &values);

     /* fill the background */
     XFillRectangle(sw->display,
		    p,
		    bg_gc,
		    (int) COL_TO_PIXEL(sw,c),
		    (int) ROW_TO_PIXEL(sw,r-1),
		    FONT_WIDTH(sw) * l,
		    FONT_HEIGHT(sw)
		    );
     Tk_FreeGC(sw->display, bg_gc);

     /* finally draw the string */
     sw->sparegc = setGC(sw, &my_ink);
     Tk_DrawChars(sw->display,p,sw->sparegc,sw->font,s,l,
		  (int) COL_TO_PIXEL(sw,c),(int) ROW_TO_BASELINE_PIXEL(sw,r));

     /* and underline if required */
     if (ink->sh & sh_select || ink->sh & sh_underline) {
	 int row = ROW_TO_BASELINE_PIXEL(sw,r);
	 if (sw->fm.descent > 1)
	     row++;
	 XDrawLine(sw->display,
		   p,
		   (ink->sh & sh_indel) ? sw->indelgc : sw->sparegc,
		   (int) COL_TO_PIXEL(sw,c),
		   row,
		   (int) COL_TO_PIXEL(sw,c+l)-1,
		   row);
     }

     if (ink->sh & (sh_box | sh_box_alt)) {
	 XDrawRectangle(sw->display,
			p,
			sw->sparegc,
			(int) COL_TO_PIXEL(sw,c),
			(int) ROW_TO_PIXEL(sw,r-1),
			sw->font_width * l - 1,
			sw->fm.linespace - 1);
     }

     /* add ruler tick if required */
     if (ink->sh & sh_tick) {
	 int i;
	 for (i = 0; i < l; i++) {
	     XDrawRectangle(sw->display,
			    p,
			    sw->sparegc,
			    (int)COL_TO_PIXEL(sw,c+i) + COL_TO_PIXEL(sw,1)/2-1,
			    (int)ROW_TO_BASELINE_PIXEL(sw,r)+1,
			    1, 2);
	 }
     }

     /* Add carets */
     if (ink->sh & sh_caret_l) {
	 int i;
	 GC gc = (ink->sh & sh_indel) ? sw->indelgc : sw->sparegc;
	 
	 for (i = 0; i < l; i++) {
	     int x = (int)COL_TO_PIXEL(sw, c+i);
	     int y = (int)ROW_TO_BASELINE_PIXEL(sw,r)+1;
	     XDrawLine(sw->display, p, gc, x,   y-1, x+2, y+1);
	 }
     }

     /* add insertion caret, split in two halves */
     if (ink->sh & sh_caret_r) {
	 int i;
	 GC gc = (ink->sh & sh_indel) ? sw->indelgc : sw->sparegc;

	 for (i = 0; i < l; i++) {
	     int x = (int)COL_TO_PIXEL(sw, c+i+1);
	     int y = (int)ROW_TO_BASELINE_PIXEL(sw,r)+1;
	     XDrawLine(sw->display, p, gc, x-1, y-1, x-3, y+1);
	 }
     }
     

     XCopyArea(sw->display, p, sw->window, copygc,
	       (int) COL_TO_PIXEL(sw,c),
	       (int) ROW_TO_PIXEL(sw,r-1),
	       FONT_WIDTH(sw) * l,
	       FONT_HEIGHT(sw),
	       (int) COL_TO_PIXEL(sw,c),
	       (int) ROW_TO_PIXEL(sw,r-1));
     Tk_FreeGC(sw->display, copygc);
 }
#endif

    Tk_FreeGC(sw->display, sw->sparegc);
}


static void _repaint_monochrome(Sheet *sw, int c, int r, int l, sheet_ink ink, char *s)
{

    GC fg_gc;
    GC bg_gc;

#define L  ( ink->sh & sh_light )
#define I  ( ink->sh & sh_inverse )
#define BG ( ink->sh & (sh_bg | sh_fg) )
    /*
    ** bg_determination
    */
    bg_gc = ( I && !L ) ? sw->normgc :
	(I || (!L && BG)) ? sw->greygc :
	    sw->whitegc;
    /*
    ** fg_determination
    */
    fg_gc = ( !I && !L ) ? sw->normgc :
	(!I || (!L && BG)) ? sw->greygc :
	    sw->whitegc;
#undef L
#undef I
#undef BG

    XFillRectangle(
		   sw->display,
		   sw->window,
		   bg_gc,
		   (int) COL_TO_PIXEL(sw,c),
		   (int) ROW_TO_PIXEL(sw,r-1),
		   FONT_WIDTH(sw) * l,
		   FONT_HEIGHT(sw)
		   );

   /* 7/1/99 johnt - changed to use TK_DrawChars for Windows support */
    Tk_DrawChars(sw->display,sw->window,fg_gc,sw->font,s,l,
		(int) COL_TO_PIXEL(sw,c),(int) ROW_TO_BASELINE_PIXEL(sw,r));

    if (ink->sh & sh_select || ink->sh & sh_underline) {
	sw->sparegc = setGC(sw, ink);
	XDrawLine(
		  sw->display,
		  sw->window,
		  fg_gc,
		  (int) COL_TO_PIXEL(sw,c),
		  (int) ROW_TO_BASELINE_PIXEL(sw,r),
		  (int) COL_TO_PIXEL(sw,c+l)-1,
		  (int) ROW_TO_BASELINE_PIXEL(sw,r)
		  );
    }

     if (ink->sh & (sh_box | sh_box_alt)) {
	 XDrawRectangle(sw->display,
			sw->window,
			fg_gc,
			(int) COL_TO_PIXEL(sw,c),
			(int) ROW_TO_PIXEL(sw,r-1),
			sw->font_width * l - 1,
			sw->fm.linespace - 1);
     }

     /* add ruler tick if required */
     if (ink->sh & sh_tick) {
	 int i;
	 for (i = 0; i < l; i++) {
	     XDrawRectangle(sw->display,
			    sw->window,
			    sw->sparegc,
			    (int) COL_TO_PIXEL(sw,c+i) + COL_TO_PIXEL(sw,1)/2-1,
			    (int) ROW_TO_BASELINE_PIXEL(sw,r)+1,
			    1, 2);
	 }
     }
}



static void _repaint(Sheet *sw, int c, int r, int l, sheet_ink ink, char *s)
{
    if (!Tk_IsMapped(sw->tkwin))
	return;

    /* Incase it hasn't been created yet due to an unmapped window */
    if (!sw->dbl_buffer)
	sw->dbl_buffer = Tk_GetPixmap(sw->display, Tk_WindowId(sw->tkwin),
				      sw->width_in_pixels,
				      sw->height_in_pixels,
				      Tk_Depth(sw->tkwin));

    if (!sw->window)
	sw->window = Tk_WindowId(sw->tkwin);

    /* Belt and braces mode, just incase! */
    if (!sw->window)
	return;

    /*
    printf("Repaint %p => (%d,%d) len %d, ink->sh=%d, text=%.*s\n",
	   sw, c, r, l, ink->sh, l, s);
    */

    if (ink->sh==sh_default) {
	XGCValues values;
	unsigned long mask;
	GC bg_gc;

	/* Create or change the GC for filling the background */
	values.foreground = sw->background;
	values.function = GXcopy;
	values.fill_style = FillSolid;
	mask = GCForeground  | GCFunction | GCFillStyle;
	bg_gc = Tk_GetGC(sw->tkwin, mask, &values);

#if 0
	/* fill the background */
	XFillRectangle(sw->display,
		       sw->window,
		       bg_gc,
		       (int) COL_TO_PIXEL(sw,c),
		       (int) ROW_TO_PIXEL(sw,r-1),
		       FONT_WIDTH(sw) * l,
		       FONT_HEIGHT(sw)
		       );
	Tk_FreeGC(sw->display, bg_gc);

	/* finally draw the string */
	Tk_DrawChars(sw->display,sw->window,sw->normgc,sw->font,s,l,
	(int) COL_TO_PIXEL(sw,c),(int) ROW_TO_BASELINE_PIXEL(sw,r));
#else
 {
     /* Double buffered */
     GC copygc;

     Pixmap p = sw->dbl_buffer;
     mask = GCFunction|GCGraphicsExposures|GCForeground;
     values.function = GXcopy;
     values.graphics_exposures = False;
     values.foreground = sw->foreground;
     copygc = Tk_GetGC(sw->tkwin, mask, &values);

     /* fill the background */
     XFillRectangle(sw->display,
		    p,
		    bg_gc,
		    (int) COL_TO_PIXEL(sw,c),
		    (int) ROW_TO_PIXEL(sw,r-1),
		    FONT_WIDTH(sw) * l,
		    FONT_HEIGHT(sw)
		    );
     Tk_FreeGC(sw->display, bg_gc);

     /* finally draw the string */
     Tk_DrawChars(sw->display,p,sw->normgc,sw->font,s,l,
		  (int) COL_TO_PIXEL(sw,c),(int) ROW_TO_BASELINE_PIXEL(sw,r));

     XCopyArea(sw->display, p, sw->window, copygc,
	       (int) COL_TO_PIXEL(sw,c),
	       (int) ROW_TO_PIXEL(sw,r-1),
	       FONT_WIDTH(sw) * l,
	       FONT_HEIGHT(sw),
	       (int) COL_TO_PIXEL(sw,c),
	       (int) ROW_TO_PIXEL(sw,r-1));
     Tk_FreeGC(sw->display, copygc);
 }
#endif

    } else {
	if (DisplayPlanes(sw->display,DefaultScreen(sw->display))==1)
	    _repaint_monochrome(sw,c,r,l,ink,s);
	else
	    _repaint_colour(sw,c,r,l,ink,s);
    }

}


static void redrawCursor(Sheet *sw, Boolean draw)
{
    SheetRow    r = sw->cursor_row;
    SheetColumn c = sw->cursor_column;
    sheet_ink ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c);
    sheet_paper paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
    sheet_ink_struct ink;

    /* check cursor is on screen */
    if (r < 0 || r > sw->rows-1 ||
	c < 0 || c > sw->columns-1 ) return;

    ink.fg = ink_base->fg;
    ink.bg = ink_base->bg;
    if (draw) {
	if (sw->hollow_cursor)
	    ink.sh = ink_base->sh | sh_box;
	else
	    ink.sh = ink_base->sh | sh_inverse;
    } else {
	ink.sh = ink_base->sh;
    }

    /*
    cursor_was = sw->display_cursor;
    sw->display_cursor = 0;
    */
    _repaint(sw, c, r, 1, &ink, paper_base);
    /*
    sw->display_cursor = cursor_was;
    */
}

static void repaintText(Sheet *sw, int c, int r, int l)
{
    sheet_ink ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c);
    sheet_ink ink_peek;
    sheet_paper paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
    sheet_paper paper_peek;
    SheetColumn c_peek;
    int i;

    while (l > 0) {
	/* find stretch where all hilight the same */
	ink_peek = ink_base;
	ink_peek++;
	paper_peek = paper_base;
	paper_peek++;
	c_peek = c;
	c_peek++;
	i = 1;
	l--;
/*#define implies(A,B) ((B)|!(A)) */
#define implies(A,B) (!(A)||(B))
	while ( (l > 0) &&
	    (ink_peek->sh == ink_base->sh) &&
	    implies(ink_base->sh&sh_fg,ink_peek->fg==ink_base->fg) &&
	    implies(ink_base->sh&sh_bg,ink_peek->bg==ink_base->bg) ) {
	    ink_peek++;
	    paper_peek++;
	    c_peek++;
	    i++;
	    l--;
	}

	_repaint(sw, c, r, i, ink_base, paper_base);

	paper_base = paper_peek;
	ink_base = ink_peek;
	c = c_peek;
    }

}

static void redisplayRegion(Sheet *sw, XRectangle *expose)
{
    int tlc,brc;
    int tlr,brr,r;

    if (sw->columns <= 0 || sw->rows <= 0)
	return;

    tlc = PIXEL_TO_COL(sw,expose->x);
    tlr = PIXEL_TO_ROW(sw,expose->y);
    brc = PIXEL_TO_COL(sw,expose->x+expose->width-1);
    brr = PIXEL_TO_ROW(sw,expose->y+expose->height-1);

    if (tlr > brr) {
	int tmp = tlr;
	tlr = brr;
	brr = tmp;
    }
    /* Adjust for potentially partial lines */
    tlr--;
    brr++;

    if (tlc < 0) tlc = 0;
    if (tlr < 0) tlr = 0;
    if (brc < 0) brc = 0;
    if (brr < 0) brr = 0;
    if (tlc >= sw->columns) tlc = sw->columns-1;
    if (tlr >= sw->rows)    tlr = sw->rows-1;
    if (brc >= sw->columns) brc = sw->columns-1;
    if (brr >= sw->rows)    brr = sw->rows-1;

    for (r=tlr;r<=brr;r++) {
	repaintText(sw, tlc, r, brc-tlc+1);
    }

    if (sw->display_cursor &&
	sw->cursor_row >= tlr &&
	sw->cursor_row <= brr &&
	sw->cursor_column >= tlc &&
	sw->cursor_column <= brc)
    {
	/* better redraw cursor */
	redrawCursor(sw,True);
    }
}



/* The Xaw interface */




void XawSheetPutText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l,
		     char *s)
/*
** Put plain text
*/
{
    int i;
    sheet_ink ink_base;
    sheet_paper paper_base;
    char *sp;

/*    printf("Printing @ %d,%d %s\n", c, r, s); */

    if (r>=0 && r<sw->rows &&
	c+l>0 && c<sw->columns &&
	l > 0) {
	if (c<0) { l += c; s -= c; c = 0; }
	if (c+l>sw->columns) l = sw->columns - c;
	for (
	    i = 0, sp = s,
	    ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c),
	    paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
	    i < l;
	    i++, ink_base++, paper_base++, sp++) {
	    ink_base->sh = sh_default;
	    *paper_base = *sp;
	}
	if (Tk_IsMapped(sw->tkwin)) {
	    _repaint(sw, c, r, l, (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c), s);

	    if (sw->display_cursor &&
		sw->cursor_row == r &&
		sw->cursor_column >= c &&
		sw->cursor_column < c+l)
	    {
		/* better redraw cursor */
		redrawCursor(sw,True);
	    }
	}
    }
}


void XawSheetPutJazzyText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l, char *s, sheet_ink ink_list)
/*
** Put multi-coloured text
*/
{
    int i;
    sheet_ink ink_base;
    sheet_paper paper_base;
    char *sp;

    if (r>=0 && r<sw->rows &&
	c+l>0 && c<sw->columns &&
	l > 0) {
	if (c<0) { l += c; s -= c; c = 0; }
	if (c+l>sw->columns) l = sw->columns - c;
	for (
	    i = 0, sp = s,
	    ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c),
	    paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
	    i < l;
	    i++, ink_base++, ink_list++, paper_base++, sp++) {
	    ink_base->fg = ink_list->fg;
	    ink_base->bg = ink_list->bg;
	    ink_base->sh = ink_list->sh;
	    *paper_base = *sp;
	}
	if (Tk_IsMapped(sw->tkwin)) {
	    repaintText(sw, c, r, l);

	    if (sw->display_cursor &&
		sw->cursor_row == r &&
		sw->cursor_column >= c &&
		sw->cursor_column < c+l)
	    {
		/* better redraw cursor */
		redrawCursor(sw,True);
	    }
	}
    }
}


void XawSheetPutHilightText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l, char *s)
/*
** Put text using default hilights
*/
{
    int i;
    sheet_ink ink_base;
    sheet_paper paper_base;
    char *sp;

    if (r>=0 && r<sw->rows &&
	c+l>0 && c<sw->columns &&
	l > 0) {
	if (c<0) { l += c; s -= c; c = 0; }
	if (c+l>sw->columns) l = sw->columns - c;
	for (
	    i = 0, sp = s,
	    ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c),
	    paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
	    i < l;
	    i++, ink_base++, paper_base++, sp++) {
	    ink_base->sh = sw->default_sh;
	    ink_base->fg = sw->default_fg;
	    ink_base->bg = sw->default_bg;
	    *paper_base = *sp;
	}
	if (Tk_IsMapped(sw->tkwin)) {
	    _repaint(sw, c, r, l, (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c), s);
	    if (sw->display_cursor &&
		sw->cursor_row == r &&
		sw->cursor_column >= c &&
		sw->cursor_column < c+l)
	    {
		/* better redraw cursor */
		redrawCursor(sw,True);
	    }
	}
    }
}

void XawSheetHilightText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l, Pixel fg, Pixel bg, SheetHilight h)
/*
** Hilight already draw text
*/
{
sheet_ink ink_base;
sheet_paper paper_base;

/*
Hilights currently supported:

    sh_default		yes
    sh_fg		yes
    sh_bg		yes
    sh_underline	yes
    sh_inverse		yes
    sh_light		yes
    sh_tick		yes
    sh_bold		yes
    sh_italic		no
*/

    if (r>=0 && r<sw->rows &&
	c+l>0 && c<sw->columns &&
	l > 0) {
	int i;

	if (c<0) { l += c; c = 0; }
	if (c+l>sw->columns) l = sw->columns - c;
	for (
	    i = 0,
	    ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c),
	    paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
	    i < l;
	    i++, ink_base++, paper_base++)
	{
	    if (h==sh_default) {
	        ink_base->sh = sh_default;
	    } else {
		if (h & sh_fg) ink_base->fg  = fg;
		if (h & sh_bg) ink_base->bg  = bg;
		ink_base->sh |= h;
	    }
	}
	repaintText(sw, (int)c, (int)r, (int)l);
    }
}

void XawSheetUnhilightText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l, Pixel fg, Pixel bg, SheetHilight h)
/*
** Remove hilighting from text
*/
{
sheet_ink ink_base;
sheet_paper paper_base;


    if (r>=0 && r<sw->rows &&
	c+l>0 && c<sw->columns &&
	l > 0) {
	int i;

	if (c<0) { l += c; c = 0; }
	if (c+l>sw->columns) l = sw->columns - c;
	for (
	    i = 0,
	    ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c),
	    paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
	    i < l;
	    i++, ink_base++, paper_base++)
	{
	    if (h==sh_default) {
	    } else {
		if (h & sh_fg) ink_base->fg  = fg;
		if (h & sh_bg) ink_base->bg  = bg;
		ink_base->sh &= !h&sh_mask;
	    }
	}
	repaintText(sw, (int)c, (int)r, (int)l);
    }
}

void XawSheetOpHilightText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l, SheetHilight h, int op)
/*
** Perform boolean operations on text
*/
{
sheet_ink ink_base;
sheet_paper paper_base;


    if (r>=0 && r<sw->rows &&
	c+l>0 && c<sw->columns &&
	l > 0) {
	int i;

	if (c<0) { l += c; c = 0; }
	if (c+l>sw->columns) l = sw->columns - c;
	for (
	    i = 0,
	    ink_base = (sheet_ink) GET_ARRAY_CELL(sw->ink,r,c),
	    paper_base = (sheet_paper) GET_ARRAY_CELL(sw->paper,r,c);
	    i < l;
	    i++, ink_base++, paper_base++)
	{
	    ink_base->sh = binary_op(h,ink_base->sh,op)&sh_mask;
	}
	repaintText(sw, (int)c, (int)r, (int)l);
	if (sw->display_cursor &&
	    sw->cursor_row == r &&
	    sw->cursor_column >= c &&
	    sw->cursor_column < c+l)
	{
	    /* better redraw cursor */
	    redrawCursor(sw,True);
	}

    }
}

void XawSheetPositionCursor(Sheet *sw, SheetColumn c, SheetRow r)
{
    if (Tk_IsMapped(sw->tkwin) && sw->display_cursor) {
	sw->window = Tk_WindowId(sw->tkwin);
	redrawCursor(sw,False);
    }
    sw->cursor_column = c;
    sw->cursor_row = r;
    if (Tk_IsMapped(sw->tkwin) && sw->display_cursor)
	redrawCursor(sw,True);
}

void XawSheetDisplayCursor(Sheet *sw, Boolean b)
{
    if (sw->display_cursor^b) {/*state change*/
	sw->display_cursor = b;
	if (Tk_IsMapped(sw->tkwin)) redrawCursor(sw, b);
    }
}

/* 7/1/99 johnt - removed XawSheetColourNameToPixel - no longer used, and doesn't compile under windows */
#ifdef NOTDEF
Pixel XawSheetColourNameToPixel(Sheet *sw, char *c)
{
    XColor rgb_db_def, hardware_def;
    Colormap cmap;
    Status s;

    cmap = DefaultColormap(sw->display, DefaultScreen(sw->display));
    s = XAllocNamedColor(sw->display, cmap, c, &rgb_db_def, &hardware_def);

    return hardware_def.pixel;
}
#endif


void XawSheetSetHilight(Sheet *sw, Pixel fg, Pixel bg, SheetHilight h)
{
    if (h & sh_fg) sw->default_fg  = fg;
    if (h & sh_bg) sw->default_bg  = bg;
    sw->default_sh = h;

}


void XawSheetDrawLine(Sheet *sw, int x1, int y1, int x2, int y2) {
    XDrawLine(sw->display,
	      sw->window,
	      sw->greygc,
	      (int) COL_TO_PIXEL(sw, x1),
	      (int) ROW_TO_PIXEL(sw, y1),
	      (int) COL_TO_PIXEL(sw, x2),
	      (int) ROW_TO_PIXEL(sw, y2)
	      );
    XDrawLine(sw->display,
	      sw->window,
	      sw->normgc,
	      (int) COL_TO_PIXEL(sw, x1),
	      (int) ROW_TO_PIXEL(sw, y1)+1,
	      (int) COL_TO_PIXEL(sw, x2),
	      (int) ROW_TO_PIXEL(sw, y2)+1
	      );
    XDrawLine(sw->display,
	      sw->window,
	      sw->greygc,
	      (int) COL_TO_PIXEL(sw, x1),
	      (int) ROW_TO_PIXEL(sw, y1)+2,
	      (int) COL_TO_PIXEL(sw, x2),
	      (int) ROW_TO_PIXEL(sw, y2)+2
	      );
}
