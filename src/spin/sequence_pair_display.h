#ifndef _SEQUENCE_DISPLAY_H_
#define _SEQUENCE_DISPLAY_H_

typedef struct seqdispresult {
    void (*op_func)(int seq_num,
		     void *obj,
		     seq_reg_data *data);
    int seq_id[2];
    int result_id;
    Tcl_Interp *interp;
    char seq_disp_win[1024];
    char *colour;
    int index;
    cursor_t *cursor[2];
    int cursor_visible[2];
    int prev_pos[2];
} seq_pair_disp_result;

int seq_disp_reg(Tcl_Interp *interp, char *seq_disp_win, int seq_id_h,
		 int seq_id_v, int cursor_id_h, int cursor_id_v,
		 int result_id, int wx, int wy);

int update_seqs(Tcl_Interp *interp, char *window1, char *window2,
		char *win_diff, char *seq1, char *seq2, int seq1_len,
		int seq2_len, int seq1_left, int seq2_left, int width, int type);

int spin_list_alignment ( char *seq1, char *seq2, char *name1, char *name2,
		     int pos1, int pos2, char *title, int type);

void SequencePairDisplay(Tcl_Interp *interp, char *raster_win, 
			 int result_index, int seq_id_h, int seq_id_v);
void DestroySequencePairDisplay(Tcl_Interp *interp, int result_index);

#endif
