#ifndef _NIP_RESULTS_H_
#define _NIP_RESULTS_H_
#include <tcl.h>

#define TASK_NIP_RENZ_INFO    0

typedef struct out_canvas_ {
  Tcl_Interp *interp;
  cursor_t *cursor;
  int cursor_visible;
} out_canvas;

#endif
