#ifndef _CLONES_H_
#define _CLONES_H_

#include "IO.h"
#include "seqInfo.h"

extern int add_vector(GapIO *io, char *V, int level);
extern int add_clone(GapIO *io, char *CN, char *CV);
extern int add_template(GapIO *io, char *TN, char *SV, char *ST, char *SI,
			int clone);
extern int add_seq_details(GapIO *io, int N, SeqInfo *si);

#endif
