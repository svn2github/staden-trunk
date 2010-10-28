#ifndef _PILEUP_H_
#define _PILEUP_H_

#include "bam.h"
#include "sam.h"
#include "sam_header.h"
#include "faidx.h"
#include "bam_maqcns.h"

typedef struct pileup {
    struct pileup *next;  // A link list, for active seqs
    void *cd;		  // General purpose per-seq client-data

    bam1_t b;		  // Bam entry associated with struct
    unsigned char *b_qual;// cached bam1_qual
    unsigned char *b_seq; // cached bam1_seq
    int  b_strand;        // 0 => fwd, 1 => rev

    int  pos;             // Current unpadded position in seq
    int  nth;		  // nth base at unpadded position 'pos'
    int  seq_offset;      // Current base position in s->seq[] array.

    int  cigar_ind;       // Current location in s->alignment cigar str
    int  cigar_op;        // Current cigar operation
    int  cigar_len;       // Remaining length of this cigar op

    int  eof;		  // True if this sequence has finished
    int  qual;            // Current qual (for active seq only)
    char base;		  // Current base (for active seq only)
    char start;		  // True if this is a new sequence
} pileup_t;

int pileup_loop(samfile_t *fp,
		int (*seq_init)(void *client_data,
				samfile_t *fp,
				pileup_t *p),
		int (*func)(void *client_data,
			    samfile_t *fp,
			    pileup_t *p,
			    int depth,
			    int pos,
			    int nth),
		void *client_data);

#endif
