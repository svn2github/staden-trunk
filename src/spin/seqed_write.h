#ifndef _SEQED_WRITE_H_
#define _SEQED_WRITE_H_

void seqed_write_sequence(char *sequence, int pos,  int line_length,
			  char *line);

void seqed_write_ruler(int pos, int line_length, char *line);

void seqed_write_complement(char *sequence, int pos, int line_length,
			    char *line);

int seqed_write(tkSeqed *se,FILE *fp, int from, int to, int line_length);

#endif

