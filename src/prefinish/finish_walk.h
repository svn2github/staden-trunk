#ifndef _FINISH_WALK_H_
#define _FINISH_WALK_H_

#include "finish.h"
#include "template.h"

experiments_t *experiment_walk(finish_t *fin, int pos, int chem, int dir,
			       int prob_start, int prob_end, int *nexp_p,
			       int etype);

#endif /* _FINISH_WALK_H_ */
