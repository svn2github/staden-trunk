/*
 * File: error.h
 */

#ifndef _XERROR_H_
#define _XERROR_H_

#define xerr_set(e,s) ( (e) ? xerr_set_globals((e),(s),__LINE__,__FILE__) : 0 )

extern void xperror(char *s, void (*out_func)(char *name, char *str));
extern int xerr_set_globals(int errnum, char *string, int line, char *file); /* 7/1/99 johnt - changed errno to errnum as errno is a macro in Visual C++ */
extern int get_xerrnum(void);

#endif /*_XERROR_H_*/
