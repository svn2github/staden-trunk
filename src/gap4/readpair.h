int check_all(GapIO *io);

int check_single(GapIO *io, int contig, f_int *relpg, f_int *lngthg,
		 f_int *idbsiz, int lreg, int rreg,
		 f_int *margl, f_int *margr, f_int *margb, f_int *margt,
		 f_int *isxmax, f_int *isymax);

int 
find_read_pairs(GapIO *io, 
		int num_contigs, 
		int *contig_array);

