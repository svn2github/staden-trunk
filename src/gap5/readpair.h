#ifndef _READPAIR_H_
#define _READPAIR_H_

#include "tg_gio.h"
#include "io_lib/hash_table.h"

enum readpair_mode {
    all_all = 0,
    end_all = 1,
    end_end = 2
};

typedef struct {
    tg_rec template; /* maybe 0 if only a direct pair */
    // tg_rec library;  /* ? want this */
    tg_rec rec[2];
    int start[2], end[2];
    tg_rec contig[2];
    int mqual[2];
} read_pair_t;

int 
find_read_pairs(GapIO *io, 
		int num_contigs, 
		contig_list_t *contig_array,
		enum readpair_mode mode,
		int end_size, int min_mq, int min_freq,
		tg_rec *library, int nlibrary);

HashTable *create_lib_hash(tg_rec *library, int nlibrary);

read_pair_t *spanning_pairs(GapIO *io, int num_contigs,
			    contig_list_t *contig_array,
			    enum readpair_mode mode,
			    int end_size, int min_mq, int min_freq,
			    HashTable *lib_hash);

#endif /* _READPAIR_H_ */
