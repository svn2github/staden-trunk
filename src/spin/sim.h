#include "align.h"

long SIM(char *A,
	 char *B,
	 long M,
	 long N,
	 long K,
	 long V[][128],
	 long Q,
	 long R,
	 long nseq,
	 float score_align,
	 align_int **S,
	 long *start1,
	 long *start2,
	 long *end1,
	 long *end2);
