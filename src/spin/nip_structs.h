typedef struct {
  int seq_num;
} seq_ops_arg;

typedef struct {
    int id;
    char *option;
} result_info_arg;

typedef struct {
    int result_id;
} result_id_arg;

typedef struct {
    int seq_id;
} seq_id_arg;

typedef struct {
    int seq_num;
    int option;
    int data;
} invoke_arg;

typedef struct {
    int library;
    int entry_mode;
    char *file;
    char *entry;
    char *sequence;
} set_seq_arg;

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
    char *raster;
    int raster_id;
    int seq_id;
    int result_id;
    char *colour;
    int line_width;
} pbc_arg;

typedef struct {
    char *codon_table;
    int win_len;
    int start;
    int end;
    int option;
    char *raster1;
    char *raster2;
    char *raster3;
    char *colour1;
    char *colour2;
    char *colour3;
    int raster_id1;
    int raster_id2;
    int raster_id3;
    int seq_id;
    char *result_id;
    int line_width;
} gene_arg;

typedef struct {
    char *codon_table;
    double error;
    int start;
    int end;
    int option;
    char *raster1;
    char *raster2;
    char *raster3;
    char *colour1;
    char *colour2;
    char *colour3;
    int raster_id1;
    int raster_id2;
    int raster_id3;
    int seq_id;
    int line_width;
} author_arg;

typedef struct {
    int start;
    int end;
    char *raster;
    char *colour;
    int line_width;
    int tick_ht;
    int raster_id;
    int seq_id;
    int result_id;
} trna_arg;

typedef struct {
    char *rasters;
    char *raster_ids;
    int seq_id;
    char *result_id;
    int start;
    int end;
    char *strand;
    char *colours;
    int line_width;
    float tick_ht;
} s_codon_arg;

typedef struct {
  double min;
  double max;
  int num_ticks;
} ticks_arg;

typedef struct {
    char *raster;
    int raster_id;
    char *result_id;
    int seq_id;
    int start;
    int end;
    char *colour1;
    char *colour2;
    char *colour3;
    int line_width;
    float tick_ht;
    char *donor;
    char *acceptor;
} splice_arg;

typedef struct {
    char *strand;
    float match;
    char *string;
    int use_iub_code;
    int start;
    int end;
    char *raster;
    int raster_id;
    char *result_id;
    int seq_id;
    char *colour;
    int line_width;
    int tick_ht;
} string_arg;

typedef struct {
    int start;
    int end;
    char *raster;
    int raster_id;
    char *result_id;
    int seq_id;
    char *colour;
    int line_width;
    int tick_ht;
    char *wt_matrix;
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
    char *xscroll;
    char *yscroll;
} nip_scroll_arg;

typedef struct {
    int seq_id;
    int id;
    int x1;
    int y1;
    int x2;
    int y2;
    char *scroll;
} nip_zoom_arg;

typedef struct {
    int id;
    int cx;
} nip_cursor_arg;

typedef struct {
    int id;
} nip_resize_arg;

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
    char *strand;
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
    int raster_id;
    char *raster;
    int pos;
} raster_scursor_arg;

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
    int seq_id;
    int id;
} delete_cursor_arg;

typedef struct {
    int seq_id;
    int id;
    int pos;
} cursor_notify_arg;

typedef struct {
    int cursorid;
    int seq_id;
} qc_arg;

typedef struct {
    char *raster;
    char *seq_id;
    char *seq_num;
} reg_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
} set_range_arg;

typedef struct {
    int result_id;
    int print_opt;
} nip_enz_info_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    int min_orf;
    char *strand;
} trans_ft_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    int min_orf;
    char *strand;
    char *filename;
} trans_fasta_arg;

typedef struct {
    int code_type;
    char *c_table;
} genetic_code_arg;

typedef struct {
    char *filename;
} read_enz_arg;

typedef struct {
    char *raster;
    char *raster_id;
    int seq_id;
    char *result_id;
    char *colour;
    int line_width;
    float tick_ht;
} plot_arg;

