#ifndef _ACE_H_
#define _ACE_H_

#include <tg_index.h>

/*
 * Parses a new ACE format file passed in.
 *
 * Returns 0 on success
 *	  -1 on error
 */
int parse_ace(GapIO *io, char *ace_fn, tg_args *a);


#endif /* _ACE_H_ */
