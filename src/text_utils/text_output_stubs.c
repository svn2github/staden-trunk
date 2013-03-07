/*
 * Stubs for the text_output system to be used when in a non windowing
 * environment.
 */

#include <stdio.h>
#include <stdarg.h>

#include "text_output.h"

static int header_outputted = 0;

void start_message(void) {}
void end_message(const char *parent) {}

/*
 * Usage: verror(priority, name, format, args...);
 * NB: don't pass more than 8K per call
 */
/* ARGSUSED */
__PRINTF_FORMAT__(3,4)
void verror(int priority, const char *name, const char *fmt, ...) { 
    va_list args;

    fprintf(stderr, "%s: ", name);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
}


/*
 * Usage: vmessage(format, args...);
 */
__PRINTF_FORMAT__(1,2)
void vmessage(const char *fmt, ...) {
    va_list args;

    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
}

/*
 * Adds a new header to the text output window.
 */
__PRINTF_FORMAT__(1,2)
void vfuncheader(const char *fmt, ...) {
    va_list args;

    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    fprintf(stdout, "\n");
    header_outputted = 1;
}

__PRINTF_FORMAT__(1,2)
void vfuncparams(const char *fmt, ...) {
}

/*
 * Used for grouping outputs together (such as the 2D plot results).
 * Basically we don't output a header if the last output was from this
 * group and there haven't been function headers outputted since.
 *
 * group numbers:
 * 1	2D plot matches
 * 2	Information from template display
 */
__PRINTF_FORMAT__(2,3)
void vfuncgroup(int group, const char *fmt, ...) {
    static int group_num = 0;
    va_list args;


    if (header_outputted || group != group_num) {
	va_start(args, fmt);
	vfprintf(stdout, fmt, args);
	fprintf(stdout, "\n");

	header_outputted = 0;
	group_num = group;
    }
}


/* Dummy function for text_output lib, but used in tk_utils */
/* ARGSUSED */
void log_file(const char *fn, const char *message) {}
int log_vmessage(int log) {
    return 0;
}
