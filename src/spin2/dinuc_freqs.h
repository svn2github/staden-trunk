#ifndef _DINUC_FREQS_H
#define _DINUC_FREQS_H

void calc_dinuc_freqs ( char seq[], int user_start, int user_end,
			double freqs[5][5]);

void calc_expected_dinuc_freqs ( char seq[], int user_start, int user_end,
				 double freqs[5][5]);


#endif
