#ifndef _ACE_H_
#define _ACE_H_

/*
 * Parses a new ACE format file passed in.
 *
 * Returns 0 on success
 *	  -1 on error
 */
int parse_ace(GapIO *io, int max_size, char *ace_fn, int no_tree,
	      int pair_reads, int merge_contigs);


#endif /* _ACE_H_ */
