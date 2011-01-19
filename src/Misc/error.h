#ifndef _ERROR_H_
#define _ERROR_H_

#include "misc.h"

extern void messout(char *fmt, ...) __PRINTF_FORMAT__(1,2);
extern void errout(char *fmt, ...) __PRINTF_FORMAT__(1,2);

#endif /*_GAP_ERROR_H_*/
