char *seq_left_end ( char seq[], int seq_length, int start, int window_length,
		    int inc);

char *seq_right_end ( char seq[], int seq_length, int end, int window_length,
		     int inc );

void get_base_comp ( char seq[], int seq_length, double base_comp[5] );
