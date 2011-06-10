#ifndef _TK_SheetP_h
#define _TK_SheetP_h

/* Shouldn't need either of these, except for a X11/Xlib.h */
#include <tk.h>
#include "intrinsic_type.h"

typedef int SheetRow;
typedef int SheetColumn;
typedef int SheetHilight;

/* hilights available */
#define sh_default	(0L)
#define sh_fg		(1L<<0)
#define sh_bg		(1L<<1)
#define sh_underline	(1L<<2)
#define sh_inverse	(1L<<3)
#define sh_light	(1L<<4)
#define sh_tick  	(1L<<5)
#define sh_bold 	(1L<<6)
#define sh_italic	(1L<<7)
#define sh_select       (1L<<8)
#define sh_box          (1L<<9)
#define sh_box_alt      (1L<<10) /* alternative to allow neighbouring boxes */
#define sh_caret_l  	(1L<<11)
#define sh_caret_r  	(1L<<12)
#define sh_indel  	(1L<<13)
#define sh_mask         ((1L<<14) - 1)

/* hilight operations */
#define HOP_MASK 0xF
#define HOP_SRC 0xC
#define HOP_DST 0xA
#define HOP_SET 0xD
#define HOP_CLR 0x2
#define HOP_TOG 0x6
#define HOP_AND(S,D) ((S & D) & HOP_MASK)
#define HOP_OR(S,D) ((S | D) & HOP_MASK)
#define HOP_NOT(S)  ((!S) & HOP_MASK)

/* define unique representation types not found in <X11/StringDefs.h> */

typedef struct {
    SheetRow rows;
    SheetColumn cols;
    char *base;
    size_t size;
} *sheet_array, sheet_array_struct;

typedef struct {
    Pixel fg;
    Pixel bg;
    SheetHilight sh;
} *sheet_ink, sheet_ink_struct, XawSheetInk;

typedef char *sheet_paper, sheet_paper_struct;

typedef struct {
    SheetRow row;
    SheetColumn column;
    Boolean visible;
} sheet_cursor;

typedef struct {
    /* Some tk specifics - should remove these somehow */
    Display *display;
    Tk_Window tkwin;
    Drawable window;
    Tk_Font        font;
    Tk_Font        bold_font;
    Tk_FontMetrics fm;
    int            font_width;

    /* resources */
    Pixel foreground;
    Pixel background;
    Pixel indel_foreground;
    sheet_cursor   cursor;
    SheetRow       rows;
    SheetColumn    columns;
    Boolean        display_cursor;
    SheetRow       cursor_row;
    SheetColumn    cursor_column;
    int            yflip;

    /* private state */
    sheet_array    paper;
    sheet_array    ink;
    int		   border_width;
    int            width_in_pixels;
    int            height_in_pixels;
    int		   hollow_cursor;
    GC             normgc;
    GC             normgcB;
    GC             greygc;
    GC             whitegc;
    GC		   indelgc;
    GC             sparegc;
    Pixel	   default_fg;
    Pixel	   default_bg;
    Pixel          light;
    Pixmap         grey_stipple;
    SheetHilight   default_sh;
    Pixmap	   dbl_buffer;
} Sheet;

extern void sheet_destroy(Sheet *sw);
extern int sheet_create(Sheet *sw, Pixel light, Pixel fg, Pixel bg, Pixel ifg);
extern void sheet_config(Sheet *sw, Pixel light, Pixel fg, Pixel bg, Pixel ifg);
extern void sheet_resize(Sheet *sw, int old_rows, int old_columns);
extern void sheet_display(Sheet *sw);
extern void sheet_clear(Sheet *sw);

/*
 * A couple of handy utility functions for obtaining base coordinates
 */
#define PIXEL_TO_COL(W,P) \
    (((long)(P) - (long)(W)->border_width) / \
     (long)(W)->font_width)
#define PIXEL_TO_ROW(W,P) \
    ((W)->yflip \
        ? ((W)->rows - (((long)(P) - (long)(W)->border_width) / (long)(W)->fm.linespace) -1) \
        : (((long)(P) - (long)(W)->border_width) / (long)(W)->fm.linespace))

#define COL_TO_PIXEL(W,C) \
    ((W)->font_width * (C) + (W)->border_width)
#define ROW_TO_PIXEL(W,R) \
    ((W)->yflip \
        ? ((W)->fm.linespace * ((W)->rows - ((R)+1)-1) + (W)->border_width) \
        : ((W)->fm.linespace * ((R)+1) + (W)->border_width))

#define ROW_TO_BASELINE_PIXEL(W,R) \
    ((W)->yflip \
        ? ((W)->fm.linespace * ((W)->rows - (R)-1) + (W)->border_width + \
           (W)->fm.ascent) \
        : ((W)->fm.linespace * (R) + (W)->border_width + \
           (W)->fm.ascent))
#define BASELINE_PIXEL_TO_ROW(W,P) \
    ((W)->yflip \
        ? (W)->rows - (((long)(P) - (long)(W)->border_width - \
            (long)(W)->fm.ascent) / (long)(W)->fm.linespace) -1\
        : (((long)(P) - (long)(W)->border_width - \
            (long)(W)->fm.ascent) / (long)(W)->fm.linespace))

void XawSheetPutText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l,
		     char *s);
void XawSheetPutJazzyText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l,
			  char *s, sheet_ink ink_list);
void XawSheetDisplayCursor(Sheet *sw, Boolean b);
void XawSheetPositionCursor(Sheet *sw, SheetColumn c, SheetRow r);
void XawSheetOpHilightText(Sheet *sw, SheetColumn c, SheetRow r, Dimension l, SheetHilight h, int op);

#endif /* _TK_SheetP_h */
