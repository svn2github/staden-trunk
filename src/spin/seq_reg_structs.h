typedef struct {
    int raster_id;
    int result_id;
} seq_result_names_arg;

typedef struct {
    char *win_diff;
    char *win_1;
    char *win_2;
    int left1;
    int left2;
    int win_len;
    int result_index;
} update_seqs_arg;

typedef struct {
    char *window;
    int seq_id_h;
    int seq_id_v;
    int cursor_id_h;
    int cursor_id_v;
    int result_id;
    int wx;
    int wy;
} display_arg;

typedef struct {
    int seqdisp_id;
    int direction;
    int pos;
} move_cursor_arg;

typedef struct {
    int raster_id;
    char *raster;
    int rx;
    int ry;
} raster_scursor_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
} raster_frame_arg;

typedef struct {
    char *raster;
    char *seq_id;
    char *seq_num;
} reg_arg;

typedef struct {
    int seq_id;
    char *type;
    int frame;
} raster_win_arg;

typedef struct {
    char *seqed_win;
    int seq_id;
} seq_disp_arg;

typedef struct {
    char *raster;
    int seq_id_h;
    int seq_id_v;
} find_result_arg;

typedef struct {
    int seq_num;
} seq_get_arg;

typedef struct {
    int seq_num;
    int option;
    int data;
} seq_invoke_arg;

typedef struct {
    int id;
    int option;
    int rf;
} invoke_arg;

typedef struct {
    int id;
    char *option;
} update_arg;

typedef struct {
    int seq_num;
    int line_width;
    int direction;
    int private;
} create_cursor_arg;

typedef struct {
    int seq_num;
    int id;
    int private;
} delete_cursor_arg;

typedef struct {
    int seq_num;
    int id;
    int pos;
    int direction;
} cursor_notify_arg;

typedef struct {
    int seq_num;
    int id;
    int ref;
    int direction;
} cursor_ref_arg;

typedef struct {
    int cursorid;
    int seq_num;
} qc_arg;

typedef struct {
    int id;
} get_ops_arg;

typedef struct {
    int id;
    char *option;
    int direction;
} result_info_arg;

typedef struct {
    char *file;
} archive_list_arg;

typedef struct {
    int seq_id;
} seq_id_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
    int format;
    char *file;
} file_save_arg;

typedef struct {
    int seq_id;
    int start;
    int end;
} set_range_arg;

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
    char *rid;
    int seq_id;
} sender_arg;

typedef struct {
    int seq_id;
    int direction;
} active_seq_arg;

typedef struct {
    int library;
    int entry_mode;
    char *file;
    char *entry;
    char *sequence;
    int direction;
    char *name;
} set_seq_arg;

typedef struct {
    int raster_id;
    int seq_id;
    int direction;
    int line_width;
} seqtoraster_arg;

typedef struct {
    int result_id;
} result_id_arg;

typedef struct {
    int pt_x;
    int pt_y;
    int id;
    int match;
} nearest_match_arg;

typedef struct {
    int seq_id;
    char *result_id;
} list_arg;

typedef struct {
    int seq_id_h;
    int start_h;
    int end_h;
    int seq_id_v;
    int start_v;
    int end_v;
    int graph;
    char *data;
} emboss_arg;

typedef struct {
    int seq_id_h;
    int seq_id_v;
    int result_id;
    int graph;
    char *name;
    char *raster;
    int raster_id;
    char *colour;
    int line_width;
} plot_arg;








