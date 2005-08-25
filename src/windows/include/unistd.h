#ifndef _UNISTD_H_
#define _UNISTD_H_

#include <io.h>

/*
 * Mock unistd.h defined for windows
 */

/* access */
/* #define access _access */
#ifndef F_OK
#  define F_OK 0
#  define X_OK 1 /* not available */
#  define R_OK 2
#  define W_OK 4
#endif

#endif
