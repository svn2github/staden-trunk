#include <stdio.h>
#include <stdlib.h>

#include "align.h"

/* Turn vector into single (most likely) nucleotide */
char consen_6(FastInt B[6]) {
    FastInt a = B[0], b = 0;

    if (B[1] > a)
        a = B[1], b = 1;
    if (B[2] > a)
        a = B[2], b = 2;
    if (B[3] > a)
        a = B[3], b = 3;
    if (B[4] > a)
        a = B[4], b = 4;
    if (B[5] > a)
        a = B[5], b = 5;

    return a == 0 ? '-' : "ACGT*-"[b];
}

void expand(char *A, char *B, int M, int N,
	    char *A2, char *B2, int *AL, int *BL, FastInt *S, int keep_pads)
{
    register int i, j, op;
    char *a = A2, *b = B2;
    

    A--; B--;
    i = j = op = 0;
    /*       fpad= abs(S[0]);
       for(i=0; i<M; i++){ 
        if (S[i] !=0 )
	  lpad = i;
       }
     
       i=0;*/
  
    while (i < M || j < N)
    {
        
      if (op == 0 && *S == 0)
	{
	    op = *S++;
	    *a++ = A[++i];
	    *b++ = B[++j];
	} else {
	    if (op == 0)
		op = *S++;
	    if (op > 0)
	    {
		*a++ = '.';
		*b++ = B[++j];
		op--;
	    } else {
		*a++ = A[++i];
		*b++ = '.';
		op++;
	    }
	}
    }

    /* Strip back expanded sequences to first non pad character */
    if (!keep_pads) {
	while (*--a == '.');
	while (*--b == '.');
    } else {
	--a;
	--b;
    }

    /* & NULL terminate */
    *++a = 0;
    *++b = 0;
    *AL = (int)(a - A2);
    *BL = (int)(b - B2);


}

void expand_6(char *A, FastInt (*B)[6], int M, int N,
	      char *A2, FastInt (*B2)[6], int *AL, int *BL, FastInt *S,
	      int keep_pads)
{
    register int i,  j, k, op;
    char *a = A2;
    FastInt (*b)[6] = B2;
    extern char base_val[128];

    A--; B--;
    i = j = op = 0;
    while (i < M || j < N)
    {
	if (op == 0 && *S == 0)
	{
	    op = *S++;
	    *a++ = A[++i];
	    for (k=0; k<6; k++) {
		(*b)[k] = B[j][k];
	    }
	    b++;j++;
	} else {
	    if (op == 0)
		op = *S++;
	    if (op > 0)
	    {
		*a++ = ' ';
		for (k=0; k<6; k++) {
		    (*b)[k] = B[j][k];
		}
		b++;j++;
		op--;
	    } else {
		*a++ = A[++i];
		for (k=0; k<6; k++) {
		    (*b)[k] = base_val['*'];
		}
		b++;j++;
		op++;
	    }
	}
    }

    /* Strip back expanded sequences to first non pad character */
    if (!keep_pads) {
	while (*--a == '*');
	while ((*b)[0] == base_val['*']) b--;
    } else {
	--a;
	--b;
    }

    *AL = (int)(a - A2);
    *BL = (int)(b - B2);
}

