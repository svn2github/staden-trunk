#ifndef _ALIGN_LIB_H
#define _ALIGN_LIB_H

#include <stdio.h>
#include <sys/types.h>

extern char consen_6(FastInt B[6]);
extern void expand(char *A, char *B, int M, int N,
		   char *A2, char *B2, int *AL, int *BL, FastInt *S,
		   int keep_pads);
extern void expand_6(char *A, FastInt (*B)[6], int M, int N,
		     char *A2, FastInt (*B2)[6], int *AL, int *BL, FastInt *S,
		     int keep_pads);

#endif
