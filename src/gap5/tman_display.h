#ifndef _tman_display_h
#define _tman_display_h

#include "os.h"
#include <io_lib/Read.h>
#include "tkTrace.h"

#define MAXCONTEXTS 1000
#define FILE_NAME_LENGTH 4096

typedef struct {
    int used;
    char file[FILE_NAME_LENGTH];
    char path[1024];
    DNATrace *tracePtr;
    int complemented;
    int mini_trace; /* 1 == in editor, 0 == full trace display */
} DisplayContext;

extern DisplayContext *manageTrace(edview *xx,
				   char *format,
				   char *rawDataFile,
				   int baseNum,
				   int leftCutOff,
				   int cutLength,
				   int complemented,
				   int baseSpacing,
				   char *traceTitle,
				   int allow_dup,
				   int seq
				   );

extern void repositionSeq(edview *xx, DisplayContext *dc, int baseNum);

extern DisplayContext *getTDisplay(edview *xx, char *file,
				   int allow_dup, int mini_trace, int *exists);

extern void deleteTrace(edview *xx, char *path);

extern void diffTrace(edview *xx, char *path1, char *path2);

extern DisplayContext *trace_path_to_dc(char *path);

extern void deleteTraceDisplay(edview *xx, DisplayContext *dc);

#endif /* _tman_display_h */
