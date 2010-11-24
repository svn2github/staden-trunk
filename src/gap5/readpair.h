#ifndef _READPAIR_H_
#define _READPAIR_H_

enum readpair_mode {
    all_all = 0,
    end_all = 1,
    end_end = 2
};

int 
find_read_pairs(GapIO *io, 
		int num_contigs, 
		contig_list_t *contig_array,
		enum readpair_mode mode,
		int end_size);


#endif /* _READPAIR_H_ */
