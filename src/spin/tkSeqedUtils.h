#ifndef _TKSEQEDUTILS_H_
#define _TKSEQEDUTILS_H_

/*
 * Useful distances
 * (treat as symbolic rather than actual distances)
 */
#define D_screen     80
#define D_halfScreen 40
#define D_character   1

/*
typedef struct _codon {
    char seq[3];
    int pos[3];
} codon;
*/

typedef struct _region {
    int start;        /* start of coding region */
    int end;          /* end of coding region */
    int num_char;     
    int line_num;     /* position to write translation: if overlaps */
    int join;         /* BOOLEAN: 1 = join to previous exon */
    int complement;   /* BOOLEAN: 0 = complement */
    Pixel colour;
} region;

void initSeqed(tkSeqed *se);
void setDimensions(tkSeqed *se);
void seqed_set_h_sb_pos(tkSeqed *se, int pos);
void seqed_set_v_sb_pos(tkSeqed *se, int pos);
int seqed_add_sequence(tkSeqed *se, int length, char *sequence, 
		       char *seq_name, int sequence_type, int seq_id, int x, 
		       int y);
void seqedTranslateAdd(Tcl_Interp *interp, tkSeqed *se,  int mode);

void seqed_positionCursor(tkSeqed *se, int seq, int pos);
void seqed_showCursor(tkSeqed *se, int seq, int pos);

void seqed_setDisplayPos2(tkSeqed *se, int pos, int keep_x_pos);
void seqed_redisplay_seq(tkSeqed *se, int pos, int keep_x_pos);
void seqed_add_renzyme(tkSeqed *se, char *filename, char *list, int num_items);
void seqedSeqInfo(tkSeqed *se, int x, int y);
void seqedTranslateDelete(tkSeqed *se, int mode);
void seqedTransMode(tkSeqed *se, int mode);
void seqed_setCursorPos(tkSeqed *se, int pos);
void set_lines(tkSeqed *se, int set_display_pos, int keep_x_pos);
void seqed_save(tkSeqed *se, char *filename, int from, int to, 
		int line_length);
int seqed_search(tkSeqed *se, char *string, int direction, int strand,
		 double per_match, int new_search, int use_iub_code);
void reset_anchor(tkSeqed *se);
#endif

