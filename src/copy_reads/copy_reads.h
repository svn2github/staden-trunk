int init_copy_reads(Tcl_Interp *interp, GapIO *io_from, GapIO *io_to, 
		    int compare_mode, int mask,
		    int min_overlap, double max_percent_mismatch, int word_len,
		    int min_match, int band,
		    double align_max_mism, double min_average_qual,
		    int display_cons, int display_seq,
		    int min_contig_len, int num_contigs_from,
		    contig_list_t *contig_array_from, int num_contigs_to,
		    contig_list_t *contig_array_to, Tcl_DString *copied_reads);
