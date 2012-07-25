#ifndef _TG_CHECK_H_
#define _TG_CHECK_H_

/*
 * Checks a single contig. lib_hash may be NULL, but if not then it is a
 * HacheTable keyed on library record numbers. We then validate all sequences
 * against this to check that no unknown libraries are present.
 *
 * Returns the number of errors found
 *         0 on success.
 */
int check_contig(GapIO *io, tg_rec crec, int fix, int level,
		 HacheTable *lib_hash, HacheTable *scaf_hash,
		 int *fixed, int *removed);


/*
 * Performs a thorough internal consistency check of all on disk data
 * structures. It's therefore quite slow, but can highlight algorithm
 * problems during development.
 *
 * Level 1 checks the basics: bins, contigs, etc.
 * Level 2 also checks the contents of the bins: mainly seqs and tags.
 *
 * Returns the number of errors found
 *         or 0 on success.
 */
int check_database(GapIO *io, int fix, int level);

/*
 * Ensures that the parent bin is large enough to cover this bin. Grow it
 * if necessary.
 */
void grow_bin(GapIO *io, bin_index_t *bin);

#endif /* _TG_CHECK_H_ */
