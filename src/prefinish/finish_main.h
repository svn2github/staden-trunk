#ifndef _FINISH_MAIN_H
#define _FINISH_MAIN_H

unsigned int *classify_bases(finish_t *fin,
			     int start,
			     int end,
			     int **virtual,
			     int (*info_func)(int        job,
					      void       *mydata,
					      info_arg_t *theirdata),
			     void *info_data);

void implement_solutions(Tcl_Interp *interp, finish_t *fin);

#endif /* _FINISH_MAIN_H */
