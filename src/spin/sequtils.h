#ifndef _SEQUTILS_H_
#define _SEQUTILS_H_

void sequence_info(char *seq_name, char *sequence, int start, int end, 
		   int seq_structure, int seq_type);
double get_seq_mass (int seq_num);
#endif
