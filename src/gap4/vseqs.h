#ifndef _VSEQS_H_
#define _VSEQS_H_

#include <tcl.h>
#include "IO.h"
#include "qual.h"

/*
 * This structure defines a virtual sequence.
 * A virtual sequence is a sequence which does not exist in the database, but
 * does link in with proper database structures via the vrseq (below).
 */
typedef struct vseq_struct {
    char *seq;
    int1 *conf;
    GReadings r;
} vseq_t;

/*
 * This structure is a place holder for virtual or real sequences (vrseqs).
 * We hold our own doubly linked list of vrseqs in much the same way that
 * the left/right neighbours in real sequences form a doubly linked list.
 * Each vrseq has either a virtual sequence pointer or a real sequence
 * number, but never both.
 */
typedef struct vrseq_struct {
    struct vrseq_struct *left;		/* NULL for end of contig */
    struct vrseq_struct *right;		/* NULL for end of contig */
    vseq_t *vseq;			/* NULL if it's a real sequence */
    int rnum;				/* Reading number */
    int position;			/* Position of sequence in contig */
} vrseq_t;

/*
 * This structure holds the virtual contig. It shadows a real contig, but
 * the contents of the contig will consist of real sequences interwoven with
 * virtual sequences.
 * Instead of having left and right sequence numbers (real sequences) we
 * have left and right pointers to vrseq_t (virtual or real sequences).
 * If the 'cons' array is NULL then GET_SEQ callback request to
 * virtual_info_func will not be available for virtual sequences, otherwise 
 * virtual sequences will just inherit their sequence from the appropriate
 * section of the consensus.
 */
typedef struct {
    GapIO *io;		/* Gap IO handle */
    int contig;		/* Real contig number */
    vrseq_t *left;	/* First virtual/real sequence in contig */
    vrseq_t *right;	/* Last virtual/real sequence in contig */
    int next_rnum;	/* Next available reading number */
    Tcl_HashTable num_hash;	/* Hash table linking pointers to rnums */
    char *cons;		/* Cached consensus. Optional - may be NULL. */
} vcontig_t;

/*
 * This defines a quality distribution along the sequence.
 * It consists of a series of ranges, specified as a %age of the sequence
 * length, and start/end quality values within those ranges. The quality
 * will then be linearly interpolated between the start/end quality and range
 * positions.
 */
typedef struct {
    int q_start[10];
    int q_end[10];
    int p_start[10];
    int p_end[10];
    int n_items; /* max 10 */
} vseq_qdist_t;

/*
 * new_vcontig
 *
 * Creates a virtual contig. This creates the linked list of vrseqs,
 * which initially just contain references to the real sequences in this
 * contig.
 * The vcontig_t should destroyed by calling del_vcontig.
 *
 * Arguments:
 *	io		Gap4 IO handle
 *	db		A vcontig_t database structure
 *	contig		The contig number to copy into db
 *
 * Returns:
 *	The virtual contig (vcontig_t *) on success.
 *	NULL on failure.
 */
vcontig_t *new_vcontig(GapIO *io, int contig);

/*
 * del_vcontig
 *
 * This frees the memory allocated within (and including) the vcontig, as
 * allocated by new_vcontig and subsequent functions.
 *
 * Arguments:
 *	vc		The virtual contig pointer to free
 *
 * Returns:
 * 	None
 */
void del_vcontig(vcontig_t *vc);

/*
 * new_vrseq
 *
 * Allocates and blanks a new virtual sequence within a virtual contig.
 * Although we return a vrseq_t (virtual or real) this contains a reference
 * to the new virtual sequence.
 *
 * Arguments:
 *	vc		The virtual contig to place this virtual sequence in
 *
 * Returns:
 *	On success; A vrseq_t struct referencing the new virtual sequence.
 *	On failure; NULL
 */
vrseq_t *new_vrseq(vcontig_t *vc);

/*
 * del_vrseq
 *
 * Deletes a virtual sequence from the database. Although we pass in a vrseq,
 * we cannot delete this if it references a real sequence.
 * This also relinks any neighbours of this vrseq as appropriate.
 *
 * Arguments:
 *	vc		The virtual contig containing the sequence.
 *	vrseq		The real/virtual sequence referencing the virtual seq.
 *
 * Returns:
 *	None
 */
void del_vrseq(vcontig_t *vc, vrseq_t *vrseq);

/*
 * link_vrseq
 *
 * This links a vrseq into a virtual contig. In order to do this it needs
 * to find the correct offset within the contig and insert vrseq into the
 * doubly linked list.
 * (FIXME: For speed, we may wish to keep an array of positions so that we
 * can binary search this to identify the insertion point quickly)
 *
 * Arguments:
 *	vc		The virtual contig
 *	vrseq		The real/virtual sequence
 *	position	Where to place the real/virtual sequence
 *
 * Returns:
 *	None
 */
void link_vrseq(vcontig_t *vc, vrseq_t *vrseq, int position);

/*
 * virtual_info_func
 *
 * This is a query function suitable to pass as an argument to the consensus
 * calculate (calc_consensus) and fragment finder (find_fragments) functions.
 *
 * This is based on the database_info function found in qualIO.c, but it
 * extends this by allowing the addition of sequences that are not in the
 * database; "fabricated sequences" if you wish. The purpose of this is to
 * allow the consensus calculation and find_fragments functions to operate
 * on our expected experimental results without actually doing them.
 *
 * Arguments:
 *	job		The type of information we are requesting.
 *	mydata		Client data used by this function (through parent) 
 *	theirdata	Both input and output: where to put results of job
 *
 * Returns:
 *	No return value, but the clientdata (depth_t) is modified.
 */
int virtual_info_func(int job, void *mydata, info_arg_t *theirdata);

vrseq_t *vrseq_index2ptr(vcontig_t *vc, int num);

#endif /* _VSEQS_H_ */
