typedef struct {
  int seq_num;
} seq_ops_arg;

typedef struct {
    char *seqed_win;
    int pos;
} seq_update_arg;
 
typedef struct {
    int id;
    char *position;
} result_pos_arg;

typedef struct {
    int win_len;
    int a;
    int c;
    int g;
    int t;
    int start;
    int end;
    int seq_id;
    int strand;
} pbc_arg;

typedef struct {
    char *codon_table;
    int win_len;
    int start;
    int end;
    int option;
    int seq_id;
    int strand;
} gene_arg;

typedef struct {
    char *codon_table;
    double error;
    int start;
    int end;
    int seq_id;
    int strand;
} author_arg;

typedef struct {
    int start;
    int end;
    int seq_id;
    int strand;
} trna_arg;

typedef struct {
    int seq_id;
    char *result_id;
} nip_list_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    int strand;
  int result_id;
} s_codon_arg;

typedef struct {
  double min;
  double max;
  int num_ticks;
} ticks_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    char *donor;
    char *acceptor;
    int strand;
} splice_arg;

typedef struct {
    int strand;
    float match;
    char *string;
    int use_iub_code;
    int start;
    int end;
    int seq_id;
} string_arg;

typedef struct {
    int start;
    int end;
    int seq_id;
    char *wt_matrix;
    int strand;
} wtmatrix_arg;

typedef struct {
    char *filename;
    char *frame;
    char *win_name;
    char *plot;
    char *win_ruler;
    char *inlist;
    int num_items;
    int text_offset;
    char *text_fill;
    int tick_ht;
    int tick_wd;
    char *tick_fill;
    int cursor_wd;
    char *cursor_fill;
    int yoffset;
    int seq_id;
    int start;
    int end;
} nip_renz_arg;

typedef struct {
    int id;
    int enzyme;
} nip_enz_name_arg;

typedef struct {
    int seqed_id;
    int pos;
} nip_seqed_cursor_arg;

typedef struct {
    int id;
    int cx;
} nip_world_arg;

typedef struct {
    int seq_id;
    char *filename;
    char *con_file;
    int concat;
    int format;
    int strand;
    char *table;
    char *range;
} codon_usage_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    int line_length;
    int size;
    int feat;
    char *selcds;
} nip_trans_arg;

typedef struct {
    int seq_id;
} length_arg;

typedef struct {
    char *codon_table;
} codon_arg;

typedef struct {
    int seq_id;
} create_cursor_arg;

typedef struct {
    int result_id;
    int print_opt;
} nip_enz_info_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    int min_orf;
    int strand;
} trans_ft_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    int min_orf;
    int strand;
    char *filename;
} trans_fasta_arg;

typedef struct {
    int code_type;
    char *c_table;
} genetic_code_arg;

typedef struct {
    char *filename;
} read_enz_arg;


