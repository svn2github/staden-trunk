#ifndef _READPAM_H_
#define _READPAM_H_
#include "align.h"

void set_dna_lookup(void);
void set_protein_lookup(void);

extern int **score_matrix;

void identity_dna_matrix(int ***matrix);
void identity_prot_matrix (int ***matrix);

void set_score_matrix(int **matrix);

/*
 * return current character order for protein or dna score matrix
 */
char *get_order(int type);
int 
create_pam_matrix(char *file_name_ptr,
		  int ***matrix);


#endif
