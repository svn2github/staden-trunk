#ifndef _F2C_H_
#define _F2C_H_

#define VOID void
typedef int integer;
typedef int ftnlen;
typedef unsigned uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;

#ifndef min
#    define min(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef max
#    define max(a,b) ((a)>(b)?(a):(b))
#endif

int swritf_(char *sbuf, char *fmt, ...);
int jfromc_(char *data, integer *length_p, ftnlen);
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb);
int s_copy(char *a, char *b, ftnlen la, ftnlen lb);
integer i_len(char *s, ftnlen n);
integer i_sign(integer *a, integer *b);

#endif
