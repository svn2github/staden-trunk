#ifndef _FINISH_HASH_H_
#define _FINISH_HASH_H_

#include "align_lib.h"
#include "hash_lib.h"

#define FIN_MAXPRIMERLEN 50

double hash_compare_primer(Hash *h, char *prim, int lprim,
			   double minmat, int skip_self, int skip_strand);

double compare_primer(char *seq1, int len1, char *prim, int lprim,
		      double minmat, int skip_self, int skip_strand);

#endif /* _FINISH_HASH_H_ */
