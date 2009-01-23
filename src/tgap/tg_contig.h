#ifndef _TG_CONTIG_H_
#define _TG_CONTIG_H_

#include <limits.h>

/*
 * 'get' functions - simply returns the structure member.
 *
 * <type> contig_get_XXX(contig_t **c)
 */
#define contig_get_start(c)  ((*(c))->start)
#define contig_get_end(c)    ((*(c))->end)
#define contig_get_bin(c)    ((*(c))->bin)
#define contig_get_name(c)   ((*(c))->name)
#define contig_get_length(c) (contig_get_end(c) - contig_get_start(c) + 1)


/*
 * 'set' functions all have prototype:
 *
 * int contig_set_XXX(GapIO *io, contig_t **c, <type> new_value)
 *
 * Returns 0 for success, possibly also modifying *c pointer
 *        -1 for failure
 */
int contig_set_start(GapIO *io, contig_t **c, int value);
int contig_set_end  (GapIO *io, contig_t **c, int value);
int contig_set_bin  (GapIO *io, contig_t **c, int value);
int contig_set_name (GapIO *io, contig_t **c, char *name);

/*
 * Returns the offset to apply to the root bin of the contig.
 * This takes into account whether it have been complemented or not.
 */
int contig_offset(GapIO *io, contig_t **c);

int contig_insert_base(GapIO *io, contig_t **c, int pos, char base, int conf);

int contig_delete_base(GapIO *io, contig_t **c, int pos);


contig_t *find_contig_by_name(GapIO *io, char *name);
GRec contig_index_query(GapIO *io, char *name);
int contig_index_update(GapIO *io, char *name, int name_len, GRec rec);

contig_t *contig_new(GapIO *io, char *name);

/* FIXME: move this elsewhere */
#define get_bin(io, bnum) ((bin_index_t *)cache_search((io), GT_Bin, (bnum)))

rangec_t *contig_seqs_in_range(GapIO *io, contig_t **c, int start, int end,
			       int *count);
rangec_t *contig_bins_in_range(GapIO *io, contig_t **c, int start, int end,
			       int *count);



/* ---------------------------------------------------------------------- */

/*
 * The iterator is basically a way to easily loop through all sequences
 * in a given range. It can optionally fetches the next range of data too
 * if you iterate off the edge.
 */
typedef struct {
    rangec_t *r;     /* r[] array itself and size of it */
    int nitems;
    int index;       /* current index into r[] array */
    int cnum;        /* contig number */
    int start;       /* current region fetch coords */
    int end;
    int cstart;      /* requested limits of region to fetch */
    int cend;
    int auto_extend; /* whether to extend past cstart..cend */
} contig_iterator;

#define CITER_FIRST  0
#define CITER_LAST   1
#define CITER_CSTART INT_MIN
#define CITER_CEND   INT_MAX

/*
 * Allocates and initialises a new contig_iterator struct.
 * 'whence' may be either CITER_FIRST or CITER_LAST and it controls whether
 * we start the iteration point from the beginning or end of the list.
 * The start and end parameters dictate the initial region to query. We
 * may specify them as either coordinates or use CITER_CSTART and CITER_CEND
 * as synonyms for the first and last coordinate in the contig.
 * Finally auto_extend controls whether the start..end range is just a
 * location to start iterating from (auto_extend == 1) or a hard limit
 * with no desire to iterate past that range (auto_extend == 0).
 *
 * Returns contig_iterator pointer on success
 *         NULL on failure
 */
contig_iterator *contig_iter_new(GapIO *io, int cnum, int auto_extend,
				 int whence, int start, int end);

/*
 * Dealloctes memory used by a contig_iterator including the cached ranges.
 */
void contig_iter_del(contig_iterator *ci);

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_next(GapIO *io, contig_iterator *ci);

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_prev(GapIO *io, contig_iterator *ci);

/*
 * Queries and/or creates a track of a given type from a region of the
 * contig at a require resolution (bpv = bases per value).
 *
 * Returns track_t to free with free_track on success
 *         NULL on failure
 */
track_t *contig_get_track(GapIO *io, contig_t **c, int start, int end,
			  int type, double bpv);

/*
 * Produces a postscript file containing a plot of the contig bin structure.
 */
void contig_dump_ps(GapIO *io, contig_t **c, char *fn);

#endif /* _TG_CONTIG_H_ */
