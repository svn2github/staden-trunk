#ifndef TK_TRACEIO_H
#define TK_TRACEIO_H

/*
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int trace_load(DNATrace *t, char *file, char *format);

/*
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int trace_save(DNATrace *t, char *file, char *format);

/*
 * Initialises the tracePos and tracePosE arrays.
 */
void trace_init_pos(DNATrace *t);

/*
 * Insert a base leftwards of 'pos'. Position 0 is left of first base.
 */
void trace_insert(DNATrace *t, int pos, char base);

/*
 * Delete a base leftwards of pos
 */
void trace_delete(DNATrace *t, int pos);

/*
 * Creates a Trace structure from a memory copy of a Read structure
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int trace_memory_load(DNATrace *t, Read *r);

/*
 * Deallocates and destroys a trace allocated by the trace_load function.
 */
void trace_unload(DNATrace *t);

#endif /* TK_TRACEIO_H */
