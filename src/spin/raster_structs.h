typedef struct {
    char *raster;
} replot_arg;

typedef struct {
    char *window;
    int min_x;
    int min_y;
    int max_x;
    int max_y;
} raster_zoom_arg;

typedef struct {
    int id;
} result_info_arg;

typedef struct {
    int id;
    char *colour;
    int line_width;
} result_config_arg;

typedef struct {
    char *raster;
    int pos;
    int raster_id;
    int seqed_id;
    char *colour;
} move_arg;

typedef struct {
    int seqed_id;
} pos_arg;

typedef struct {
    char *window;
} win_arg;

typedef struct {
    char *raster;
    int pos;
} cursor_arg;

typedef struct {
    int old_id;
    int new_id;
    char *raster_old;
    char *raster_new;
    int result_id;
    char *job;
} update_win_arg;

typedef struct {
    int id;
} id_arg;

typedef struct {
    int id;
    char *option;
} result_arg;

typedef struct {
  double min;
  double max;
  int num_ticks;
} ticks_arg;

typedef struct {
    int raster_id;
    char *raster;
    int pos;
    int cursor_id;
    int direction;
} raster_mcursor_arg;

typedef struct {
    int raster_id;
    char *raster;
    int pos;
    int direction;
} raster_fcursor_arg;

