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
    int strand;
} align_seqs_arg;

typedef struct {
    char *type;
} result_type_arg;

typedef struct {
    int id;
    char *colour;
    int line_width;
} result_config_arg;

typedef struct {
    int seq_id;
    int rf;
    int start;
    int end;
} translate_arg;

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
    char *element;
    char *container;
    int element_id;
    int container_id;
    int seq_id_h;
    int seq_id_v;
    int result_id;
    char *element_type;
    char *colour;
    int line_width;
    Tcl_Obj *results;
} splot_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int win_len;
    int min_match;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int strand;
    int score;
} similar_spans_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int word_len;
    int strand;
} identity_arg;

typedef struct {
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
    int strand;
} quick_scan_arg;

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
    int seq_num;
} seq_num_arg;

typedef struct {
    char *name;
} seq_name_arg;

typedef struct {
    char *type;
} sip_seq_names_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int start_h;
    int end_h;
    int start_v;
    int end_v;
    int num_align;
    float score_align;
    float match;
    float transition;
    float transversion;
    float start_gap;
    float cont_gap;
    int strand;
} sim_arg;
