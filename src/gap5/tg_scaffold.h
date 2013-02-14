#ifndef _TG_SCAFFOLD_H_
#define _TG_SCAFFOLD_H_


/*
 * Adds a contig to the end of a scaffold.
 * "prev_gap" is the distance from the previous contig.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_add(GapIO *io, tg_rec scaffold, tg_rec contig,
		 int gap_size, int gap_type, int evidence);

/*
 * Adds a contig named ctg_name to a scaffold named scaf_name. The names are
 * looked up in the B+Tree index.
 */
int scaffold_add_by_name(GapIO *io, char *scaf_name, char *ctg_name,
			 int gap_size, int gap_type, int evidence);

/*
 * Removes a contig from a scaffold.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_remove(GapIO *io, tg_rec scaffold, tg_rec contig);

/*
 * Given a contig order and a set of current scaffolds, this updates the
 * order of entries within each scaffold to match the contig order.
 *
 * For example if we have contigs in order 1 3 5 2 6 8 4 7 9 and
 * scaffolds {1 2 3 4} {5 6 7 8 9} we would shuffle the scaffold members
 * to        {1 3 2 4} {5 6 8 7 9}
 *
 * The purpose is for integration with contig shuffling in the Contig List
 * or Contig Selector. The master contig order array is what gets shuffled
 * manually by the user and it is also the definitive order to use when
 * outputting data (so it is completely under users control whether they
 * sort by name, size or scaffold).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int update_scaffold_order(GapIO *io);

/*
 * Loads a new scaffold from an AGP file.
 * Contigs that are not listed in this will keep their old scaffold data.
 * Contigs which are listed and have an existing scaffold will be amended.
 *
 * AGP format is http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_from_agp(GapIO *io, char *fn);

/*
 * Exports Scaffold information to an AGP file
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_to_agp(GapIO *io, char *fn);

/*
 * Complements a scaffold; both complementing each contig within it and
 * reversing the order of contigs in the scaffold.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int complement_scaffold(GapIO *io, tg_rec srec);

/*
 * Search for a scaffold record by name.
 * 
 * Returns record number on success
 *         -1 on failure
 */

tg_rec scaffold_index_query(GapIO *io, char *name);

#endif /* _TG_SCAFFOLD_H_ */
