#ifndef _TK_SEQED_H
#define _TK_SEQED_H

#include <tk.h>
#include "sheet.h"
#include "tkSeqedNames.h"
#include "seq_reg.h"
#include "renz_utils.h"

#define TKSHEET(se)   ((tkSheet *)(se))

#define MAX_DISPLAY_WIDTH 200
#define MAX_DISPLAY_HEIGHT 20
#define MAX_TRANS_LINES   7
#define MAX_LINES         12

#define AUTO    0
#define FRAME1p 1
#define FRAME2p 2
#define FRAME3p 3
#define FRAME1c 4
#define FRAME2c 5
#define FRAME3c 6
#define SEQ 7
#define RULER 8
#define COMP 9
#define RENZ 10
#define AUTO_C 11

typedef struct _tkSeqed {
#   include "tkSheet_struct.h"
    int rid; /* registration id */
    char *xScrollCmd;
    char *yScrollCmd;
    int heightmin;        /* min display height */
    int heightmax;       /* max display height */
    int displayWidth;                /* current screen width */
    int displayHeight;                /* current screen height */
    int cursorPos;
    int cursorSeq;
    char cursorCol[10];
    int displayPos;
    int displayYPos;
    int extent_left;    /* always 0 */
    int extent_right;   /* length of sequence */
    int num_lines;     /* number of display lines */
    char *sequence;                  /* sequence */
    int seq_len;                     /* sequence length */
    char *seq_name;                  /* sequence name */
    int seq_id;            /* registered sequence id */
    int sequence_type;               /* 0 = linear, 1 = circular */
    int rulerDisplayed;              /* ruler display switch */
    int complementDisplayed;         /* complemented seq display switch */
    int translationDisplayed;         /* translation display switch */
    int autoDisplayed;         /* automatic translation display switch */
    int renzDisplayed;         /* restriction enzymes display switch */

    int trans[MAX_TRANS_LINES]; /* which translation lines to display */
    int trans_mode;    /* 3 or 1 letter code */
    int trans_lines;   /* number of translation lines */
    int trans_lines_p;   /* number of translation lines */
    int trans_lines_c;   /* number of translation lines */
    int renz_lines;   /* number of restriction enzyme lines */
    int auto_trans_lines; /* number of automatic translation lines */
    int auto_c_trans_lines; /* number of automatic complemented translation lines */
    int lines[MAX_LINES]; 
    int anchor_pos;     /* pos of SEQ in seqed which remains const as scroll */
    cursor_t *cursor;
    int prev_pos;
    int cursor_visible;
    R_Enz *r_enzyme;
    int num_enzymes;
} tkSeqed;

int Seqed_Init(Tcl_Interp *interp);
void bell(void);

#endif
