#ifndef _TMAN_INTERFACE_H
#define _TMAN_INTERFACE_H

#include "tman_display.h"
#include "dstring.h"

#define TRACE_TYPE_SEQ	0
#define TRACE_TYPE_CON	1
#define TRACE_TYPE_DIFF	2
#define TRACE_TYPE_MINI 3
#define TRACE_TYPE_POS_CONTROL 4
#define TRACE_TYPE_NEG_CONTROL 5

typedef struct tman_dc_ {
    DisplayContext *dc;
    int type;
    tg_rec seq;
    int pos;
    int derivative_seq; /* zero if none */
    int derivative_offset; /* diff in base position from der.seq */
    edview *xx;
} tman_dc;

extern tman_dc *find_free_edc(void);

extern void tman_reposition_traces(edview *xx, int pos, int mini_trace);

extern DisplayContext *tman_manage_trace(
					 char *format,
					 char *rawDataFile,
					 int baseNum,
					 int leftCutOff,
					 int cutLength,
					 int complemented,
					 int baseSpacing,
					 char *traceTitle,
					 edview *xx,
					 tg_rec seq,
					 int allow_dup,
					 int mini_trace
					 );

extern void tman_unhighlight(tman_dc *edc);

extern void tman_highlight(tman_dc *edc);

extern tman_dc *find_edc(DisplayContext *dc);

extern int tman_get_lock(edview *xx);

extern void tman_set_lock(edview *xx, int val);

extern void tman_shutdown_traces(edview *xx, int limit_to);

extern void tman_init_problem_traces(char *spec);

extern void tman_problem_traces(edview *xx, int pos);

void cons_edc_trace(edview *xx, int start, int end, int strand, int match,
		    int exception);

Read *cons_trace(edview *xx, int start, int end, int strand,
		 int match, int exception);

DisplayContext *diff_edc_traces(edview *xx, tman_dc *ed1, tman_dc *ed2);

void tman_set_lock(edview *xx, int val);

int tman_get_lock(edview *xx);

void deleteTrace(edview *xx, char *path);

void diffTrace(edview *xx, char *path1, char *path2);

DisplayContext *diff_traces(edview *xx, int seq1, int seq2, int pos);

/*
 * From a single scroll of trace named 'path', this scrolls all other
 * traces in the same trace display. Ie a real trace->trace lock mode.
 */
void edScrollTraces(edview *xx, char *path, char *command);

/*
 * Called when the user has selected to automatically perform trace
 * differencing on their reference traces (both the wildtype/negative control
 * and optionally a positive control).
 */
int auto_diff(edview *xx, int seq, int pos);

/*
 * Converts editor consensus position to trace coordinate (bases)
 */
int tman_get_trace_position(edview *xx, tman_dc *dc, int pos, int *end);

#endif /* _TMAN_INTERFACE_H */
