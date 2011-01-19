#include <stdio.h>
#include <stdarg.h>

#include "error.h"

/*
 * Usage:
 *
 * errout(format, args...);
 */
__PRINTF_FORMAT__(1,2)
void errout(char *fmt, ...) {
    va_list args;

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);

    va_end(args);
}

/*
 * Usage:
 *
 * messout(format, args...);
 */
__PRINTF_FORMAT__(1,2)
void messout(char *fmt, ...) {
    va_list args;

    va_start(args, fmt);
    vfprintf(stdout, fmt, args);

    va_end(args);
}
