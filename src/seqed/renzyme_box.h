#if !defined(RENZYME_BOX_H)
#define RENZYME_BOX_H

typedef struct R_ENZYME_ {
    char  *name;
    int   ID;
    char  *rec_seq;
    char  *rec_seq_text;  
    int   cut_pos_1;
    int   cut_pos_2;
    float av_frag_size;
    char  *prototype;
    char  *supplier_codes;
    int   stock_level;
} RENZYME;

typedef struct R_ENZYMES_ {
    RENZYME **renzyme;
    int capacity;
    int used;
    int max_rec_seq; /* maximum recoginition sequence length for checking edits*/
} RENZYMES;

int REnzyme_box_Init(Tcl_Interp *interp);
int GetRenzInfo(ClientData cData, Tcl_Interp *interp, int argc, char **argv);
int SaveRenzInfo(ClientData cData, Tcl_Interp *interp, int argc, char **argv);
RENZYME *init_renzyme (void);
RENZYMES *init_renzymes (void);
void free_renzyme (RENZYME *renzyme);
void free_renzymes ( RENZYMES *rs );
RENZYMES *realloc_renzymes ( RENZYMES *rs, int num_ren);

extern  RENZYMES *renzymes;
#endif
