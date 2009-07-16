#ifndef _tman_display_h
#define _tman_display_h

#include <io_lib/Read.h>
#include "os.h"
#include "IO.h"         /* IMPORT - FILE_NAME_LENGTH */
#include "tkTrace.h"

#define MAXCONTEXTS 1000

typedef struct {
    int used;
    char file[FILE_NAME_LENGTH];
    char path[1024];
    DNATrace *tracePtr;
    int complemented;
    int mini_trace; /* 1 == in editor, 0 == full trace display */
} DisplayContext;

extern DisplayContext *manageTrace(EdStruct *xx,
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

extern void repositionSeq(EdStruct *xx, DisplayContext *dc, int baseNum);

extern DisplayContext *getTDisplay(EdStruct *xx, char *file,
				   int allow_dup, int mini_trace, int *exists);

extern void deleteTrace(EdStruct *xx, char *path);

extern void diffTrace(EdStruct *xx, char *path1, char *path2);

extern DisplayContext *trace_path_to_dc(char *path);

extern void deleteTraceDisplay(EdStruct *xx, DisplayContext *dc);

#endif /* _tman_display_h */
