#if !defined(RENZYME_SEARCH_H)
#define RENZYME_SEARCH_H

#define MAXMATCHES 10000

#include "renzyme_box.h"
#include "read_sequence.h"
#include "editor.h"

/* for each cut */
typedef struct rmatch {
    unsigned short enz_name; /* position of enzyme in RENZYMES array */
    int cut_pos1;  /* 3'cut position in bases */
    int cut_pos2;  /* 5'cut position in bases */
} R_MATCH;

int RenzSearch(ClientData cData, Tcl_Interp *interp, int argc, char **argv);
int RenzSearchSelect(ClientData cData, Tcl_Interp *interp, int argc, char **argv);
int REnzyme_search_Init(Tcl_Interp *interp);
RENZYMES *get_selected_renzyme (int num_sel, char **sel);
void free_sites ( SITES *ss );
int insert_fragment_in_editor (EDITOR_RECORD *editor_record, EDIT *edit, int add_to_undo);
int editor_complement_shutdown (EDITOR_RECORD *er, EDIT *edit, int add_to_undo);
SEQUENCE *create_fragment_from_sequence (int seq_num, int start, int end, char *renz_start, char *renz_end);

int sequences_notify_graphic (int seq_num, int pos);
SITES *init_sites (void);
SITE *init_site (void);
SITES *realloc_sites ( SITES *ss, int num_site);
int max_rec_seq_length (RENZYMES *sel_rs);

extern  RENZYMES *renzymes;
extern  EDITOR_RECORDS *editor_records;


#endif
