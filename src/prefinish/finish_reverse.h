#ifndef _FINISH_REVERSE_H_
#define _FINISH_REVERSE_H_

#include "finish.h"

experiments_t *experiment_reverse(finish_t *fin, int pos,
				  int chem, int dir,
				  int prob_start, int prob_end, int *nexp_p);

#endif /* _FINISH_REVERSE_H_ */
