#ifndef _SAM_H_
#define _SAM_H_

#include <tg_index.h>

int parse_bam(GapIO *io, const char *fn, tg_args *a);
int parse_sam(GapIO *io, const char *fn, tg_args *a);
char *sam_aux_stringify(char *s, int len);


#endif /* _SAM_H_ */
