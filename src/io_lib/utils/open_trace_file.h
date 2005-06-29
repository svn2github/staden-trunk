#ifndef _OPEN_TRACE_FILE_H_
#define _OPEN_TRACE_FILE_H_

#include "mFILE.h"

/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in RAWDATA (which is a colon separated list).
 *
 * Returns a mFILE pointer when found.
 *           NULL otherwise.
 */
mFILE *open_trace_mfile(char *file, char *relative_to);
FILE *open_trace_file(char *file, char *relative_to);


#endif /* _OPEN_TRACE_FILE_H_ */
