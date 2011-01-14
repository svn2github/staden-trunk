#ifndef _SEQUENCE_FORMATS_H_
#define _SEQUENCE_FORMATS_H_

#define DNA 1
#define PROTEIN 2

#define loc_len 1024
#define number_keys 63
#define number_quas 70

typedef struct base_pos {
  int start_pos;
  int end_pos;
  char type_range[2];
  struct base_pos *next;
}BasePos;

typedef struct {
  BasePos *loca;
  char type_loca[3];
  int id;
  char *cdsexpr;
  char *qualifier[number_quas];
} Featcds;



#ifdef _MSC_VER
# ifdef BUILDING_SEQ_UTILS_DLL
#  define SEQ_UTILS_EXPORT  extern __declspec(dllexport)
# else
#  define SEQ_UTILS_EXPORT  extern __declspec(dllimport)
# endif
#else
# define SEQ_UTILS_EXPORT extern
#endif


SEQ_UTILS_EXPORT char feat_quas[number_quas][20];
SEQ_UTILS_EXPORT char feat_key[number_keys][16];
SEQ_UTILS_EXPORT char genetic_code_ft[16][10];


int parse_feat(char *locexpr, Featcds **key_index, int i);
int read_cds_pos(char *locexpr, int *start_pos, int *end_pos);
int read_cds_pos_join(BasePos **head, char *locexpr);
BasePos *add_list_item(BasePos **head, BasePos *entry, int start_pos, int end_pos, char *type_range);

int get_seq_ft ( Featcds **key_index, char **seq, int max_len, 
                 int *seq_len, char *file_name, char *entry_name_in, 
                 char **identifier, int *err);

int get_seq ( char **seq, int max_len, int *seq_len, char *file_name,
  char *entry_name_in);

/* check sequence content to see if it looks ok: 85% a,c,g,t is DNA
   98% protein codes is protein, else is crap 
   returns 1 for dna, 2 for protein, 0 for anything else */
int get_seq_type ( char *seq, int seq_len );

int get_identifiers (char *file_name, 
		     char ***list,
		     int *num_identifiers);

void free_key_index(Featcds **key_index);
int purify_range (char *range);

#endif

