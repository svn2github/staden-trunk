typedef struct {
    char *raster;
} replot_all_arg;

typedef struct {
    int seq_num;
    int direction;
} active_seq_arg;

typedef struct {
    int win_len;
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int type_h;
    int type_v;
    int use_av_comp;
} find_prob_arg;

typedef struct {
    int win_len;
    int num_matches;
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int type_h;
    int type_v;
    int use_av_comp;
} find_score_arg;

typedef struct {
    int pt_x;
    int pt_y;
    int id;
    int match;
} nearest_match_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int match;
    int mismatch;
    int start_gap;
    int cont_gap;
} align_seqs_arg;

typedef struct {
    char *window;
    int min_x;
    int min_y;
    int max_x;
    int max_y;
} raster_zoom_arg;

typedef struct {
    char *type;
} result_type_arg;

typedef struct {
    int seq_num;
} result_info_arg;

typedef struct {
    int id;
    char *colour;
    int line_width;
} result_config_arg;

typedef struct {
    int seq_num;
    int option;
    int rf;
    int start;
    int end;
} invoke_arg;

typedef struct {
    int seq_id;
    int rf;
    int start;
    int end;
} translate_arg;

typedef struct {
    int library;
    int entry_mode;
    char *file;
    char *entry;
    char *sequence;
    int direction;
} set_seq_arg;

typedef struct {
    char *direction;
    int result_index;
} seq_len_arg; 

typedef struct {
    char *window1;
    char *window2;
    char *label1;
    char *label2;
    int win_len;
    int result_index;
} insert_seq_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int result_id;
    char *raster;
    int raster_id;
    char *colour;
    int line_width;
} plot_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int win_len;
    int min_match;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int result_id;
    char *raster;
    int raster_id;
    char *colour;
    int line_width;
} similar_spans_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int word_len;
} identity_arg;

typedef struct {
    char *window;
    int raster_id;
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int win_len;
    int min_match;
    int word_len;
    float min_sd;
    char *colour;
    int line_width;
} quick_scan_arg;

typedef struct {
    char *file;
} archive_list_arg;

typedef struct {
    char *file;
    int type;
} set_score_matrix_arg;

typedef struct {
    int type;
} get_score_matrix_arg;

typedef struct {
    int seq_id;
} length_arg;

typedef struct {
    char *type;
} raster_win_arg;

typedef struct {
    int result_id;
} result_id_arg;

typedef struct {
    int seq_id;
} seq_id_arg;

typedef struct {
    int seq_num;
} seq_num_arg;

typedef struct {
    char *name;
} seq_name_arg;

typedef struct {
    char *type;
} sip_seq_names_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
} set_range_arg;

typedef struct {
    char *window;
    int raster_id;
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    char *colour;
    int line_width;
    int num_align;
    float score_align;
    float match;
    float transition;
    float transversion;
    float start_gap;
    float cont_gap;
} sim_arg;

typedef struct {
    int seq_id;
    int f1;
    int f2;
    int f3;
    int all;
    int start;
    int end;
} translate_seq_arg;

typedef struct {
    int seq_id;
    int origin;
} rotate_arg;

typedef struct {
    int seq_id;
    char *result_id;
} list_arg;
