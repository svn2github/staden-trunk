#ifndef _SEQED_H_
#define _SEQED_H_

#include <tcl.h>

#include "tkSeqed.h"

typedef struct seqedresult {
    void (*op_func)(int seq_num,
		     void *obj,
		     seq_reg_data *data);
    int seq_id;
    Tcl_Interp *interp;
    char seqed_win[1024];
    char *colour;
    int index;
} SeqedResult;

tkSeqed * seqed_id_to_se(Tcl_Interp *interp, int seqed_id);

int seqed_seq_length(Tcl_Interp *interp, int seqed_id);


/*
 * add sequence to seqed widget
 */
int add_seq_seqed(Tcl_Interp *interp, char *sequence, char *seqed_win,
		  int seq_num);

/*
 * update sequence position
 */
void update_seqed(Tcl_Interp *interp, char *seqed_win, int pos);
void seqed_set_cursor_pos(Tcl_Interp *interp, int seqed_id, int pos);
#endif
