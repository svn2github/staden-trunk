#ifndef REG_EXP_H
#define REG_EXP_H

/*
 * History:
 * 20/07/93 jkb creation
 * 7/1/99 johnt - changed to use TCL regular expression code when possible (NO_TCL_REGEXP not defined)
 * 31/3/99 johnt - added interp arg, so that TCL REGEXPS don't crash !
 */
#ifdef NO_TCL_REGEXP
/*
 * Solaris makes it really awkward to use BSD style regular expression code
 * and we wish to steer clear of /usr/ucblib. So we will be able to handle
 * both System V and BSD styles.
 */
#  ifdef SYSV_REGEX
#     include <libgen.h>
#     define REGCMP(i,x) regcmp(x, (char *)0)
#     define REGEX(i,str, exp) (regex(exp, str) ? 1 : 0)
#     define REGFREE(i,exp) free(exp)
#  else
      extern char *re_comp(char *s);
      extern int re_exec(char *s);
#     define REGCMP(i,x) (re_comp(x) ? 0 : (char *)1)
#     define REGEX(i,str, exp) re_exec(str)
#     define REGFREE(i,exp)
#  endif /* SYSV_REGEX */

#else

#  include <tcl.h>
#  define REGCMP(i,x)      (char*)Tcl_RegExpCompile(i,x)
#  define REGEX(i,str,exp) Tcl_RegExpExec(i,(struct Tcl_RegExp_ *)exp,str,str)
#  define REGFREE(i,exp)   

#endif

#endif /* REG_EXP_H */
