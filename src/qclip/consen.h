#ifndef _CONSEN_H_
#define _CONSEN_H_

int end_of_good ( char *seq, 
		 int start,
		 int window_len1, 
		 int max_unknown1, 
		 int window_len2,
		 int max_unknown2);

int start_of_good ( char *seq, 
		   int start,
		   int window_len1, 
		   int max_unknown1, 
		   int window_len2,
		   int max_unknown2);
#endif
