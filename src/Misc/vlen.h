#ifndef _VLEN_H_
#define _VLEN_H_

#include "misc.h"

extern int vflen(const char *fmt, va_list ap);
extern int flen(const char *fmt, ...) __PRINTF_FORMAT__(1,2);

#endif /* _VLEN_H_ */
