#include <stdio.h>
#include <stdarg.h>

/*
 * Usage:
 *
 * errout(format, args...);
 */
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
void messout(char *fmt, ...) {
    va_list args;

    va_start(args, fmt);
    vfprintf(stdout, fmt, args);

    va_end(args);
}
