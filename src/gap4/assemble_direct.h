#ifndef _ASSEMBLE_DIRECT_H_
#define _ASSEMBLE_DIRECT_H_

#include "IO.h"
#include "align.h"
#include "seqInfo.h"

typedef struct {
    char *con_all;
    char **con_item;
    int con_len;
    int num_contigs;
} consen_info;

typedef struct {
    align_int *res;
    int start1;
    int len1;
    int start2;
    int len2;
} align_info;

/*
 * Produces the consensus sequence for all contigs in the database.
 *
 * Returns consen_info struct for success, NULL for failure.
 */
extern consen_info *all_consensus(GapIO *io, float percd);
extern void free_all_consensus(consen_info *ci);

extern char *assemble_direct(GapIO *io, int display, double max_mism,
			     char *inlist, int do_alignments, int enter_all,
			     int ignore_vec);

extern align_info *assemble_align(GapIO *io, SeqInfo *si, consen_info *ci,
				  int contig, int *pos, int dir, int tol,
				  int display, int maxpads, double max_mism,
				  int *ierr);

extern int enter_reading(GapIO *io, SeqInfo *si, int comp, align_info *ai,
			 int contig, int position);

extern int link_reading(GapIO *io, int from_rnum, int rnum,
			int contig, int position);

extern int recalc_consensus(GapIO *io, consen_info *ci, int contig, int pos,
			    int read_len, int old_clen, int new_clen);

#endif
