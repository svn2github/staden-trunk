
typedef struct _r_line{
    int num;
    int same;
} r_line;


void seqed_delete_renzyme(tkSeqed *se);
int seqed_redisplay_renzyme(tkSeqed *se, int pos, int keep_x_pos);
int seqedREnzyme(tkSeqed *se, char *filename, char *list, int num_items, int width);
int seqed_write_renzyme(char *sequence,
			int sequence_type,
			R_Enz *r_enzyme,
			int num_enzymes,
			int pos, 
			int line_length,
			int name_overlap,
			char ***alines,
			int *max_renz_lines,
			int *num_lines);

int seqed_get_max_lines(void);
int seqed_get_max_name_length(void);
void free_lines(void);
void free_r_enzyme(R_Enz *r_enzyme, int num_enzymes);
