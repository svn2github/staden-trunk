#ifndef _GENETIC_CODE_H
#define _GENETIC_CODE_H

extern char  codon_to_acid1(char *codon);
extern char *codon_to_acid3(char *codon);
extern char  codon_to_cacid1(char *codon);
extern char *codon_to_cacid3(char *codon);

int write_screen_cod_table (double codon_table[4][4][4]);
void codon_table_percent(double codon_table[4][4][4]);
void init_genetic_code (void);
int read_cod_table ( FILE *in_file, double codon_table[4][4][4] );
int write_cod_table ( FILE *out_file, double codon_table[4][4][4] );
void codon_table_64(double codon_table[4][4][4], double array64[64], int job);
int legal_codon ( char *codon );
void gen_cods_from_ac ( double codon_table[4][4][4] );
void gen_cods_from_bc ( double codon_table[4][4][4], double base_comp[5] );
void scale_codon_table ( double codon_table[4][4][4], double total );
int read_genetic_code ( FILE *in_file, char code_table[5][5][5] );
int read_global_genetic_code ( FILE *in_file);
void reset_genetic_code ( char new_genetic_code[5][5][5] );
int write_screen_genetic_code (char code_table[5][5][5] );
void average_acid_comp(double codon_table[4][4][4]);
void even_cods_per_acid(double codon_table[4][4][4]);
void third_pos_comp ( double codon_table[4][4][4] );
char (*get_global_genetic_code(void))[5][5];
int *get_genetic_code_idx(int complemented);
int load_genetic_code_number(int index);
int embl_aa_three2one(char *aa);
int edit_genetic_code(char *change);
#define TO_64 1
#define FROM_64 2
#endif
