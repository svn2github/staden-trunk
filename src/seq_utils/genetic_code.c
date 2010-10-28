/*	routines for handling the genetic code

	We need to be able to have a standard code and allow changes
	to make a current code.
	We need to translate dna to amino acids.
	We need to count codons.
	We need to use indexes into the codon tables for several methods.

	The genetic code is written tcag order, our lookups are acgt.
*/
#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#include "dna_utils.h"
#include "array_arith.h"
#include "genetic_code.h"
#include "text_output.h"
#include "misc.h"

static char one_letter[] = {"ACDEFGHIKLMNPQRSTVWY*-"};

static char *three_letter[] = {
    "Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys",
    "Leu","Met","Asn","Pro","Gln","Arg","Ser","Thr","Val",
    "Trp","Tyr","***","---"};

static char genetic_code [5][5][5];

static char std_genetic_code [5][5][5] = {
    
   {{'F','F','L','L','-'},
    {'S','S','S','S','S'},
    {'Y','Y','*','*','-'},
    {'C','C','*','W','-'},
    {'-','-','-','-','-'}},

   {{'L','L','L','L','L'},
    {'P','P','P','P','P'},
    {'H','H','Q','Q','-'},
    {'R','R','R','R','R'},
    {'-','-','-','-','-'}},

   {{'I','I','I','M','-'},
    {'T','T','T','T','T'},
    {'N','N','K','K','-'},
    {'S','S','R','R','-'},
    {'-','-','-','-','-'}},

   {{'V','V','V','V','V'},
    {'A','A','A','A','A'},
    {'D','D','E','E','-'},
    {'G','G','G','G','G'},
    {'-','-','-','-','-'}},

   {{'-','-','-','-','-'},
    {'-','-','-','-','-'},
    {'-','-','-','-','-'},
    {'-','-','-','-','-'},
    {'-','-','-','-','-'}} 
};

/*	change genetic code */
void reset_genetic_code ( char new_genetic_code[5][5][5] ) {
    int i,j,k;

    for ( i=0;i<5;i++ ) {
	for ( j=0;j<5;j++ ) {
	    for ( k=0;k<5;k++ ) {
		genetic_code[i][j][k] = new_genetic_code[i][j][k];
	    }
	}
    }
}

/*	set genetic code to the standard table */
void init_genetic_code (void) {

    int i,j,k;

    for ( i=0;i<5;i++ ) {
	for ( j=0;j<5;j++ ) {
	    for ( k=0;k<5;k++ ) {
		genetic_code[i][j][k] = std_genetic_code[i][j][k];
	    }
	}
    }
}



int legal_codon ( char *codon ) {

    if ( dna_lookup [ (unsigned)codon [ 0 ] ] == 4 ) return 0;
    if ( dna_lookup [ (unsigned)codon [ 1 ] ] == 4 ) return 0;
    if ( dna_lookup [ (unsigned)codon [ 2 ] ] == 4 ) return 0;

    return 1;
}


int genetic_code_idx[5] = {2,1,3,0,4};
int cgenetic_code_idx[5] = {0,3,1,2,4};

char *three_letter_code ( char one_letter_code ) {

/*	Input character representing an amino acid in 1 letter code
	Output its three letter code equivalent
	Note we assume incoming character is upper case
*/

    unsigned int i;

    one_letter_code = toupper(one_letter_code);

    for (i=0;i<sizeof(one_letter)-1;i++){
	if ( one_letter_code == one_letter[i] ) return three_letter[i];
    }
    return "   ";
}

/************************************************************/
char codon_to_acid1 ( char *codon ) {

/*	given a codon of three bases, return the one letter code for the acid */

    return genetic_code 
	[ genetic_code_idx [ dna_lookup [ (unsigned)codon [ 0 ] ] ] ]
	    [ genetic_code_idx [ dna_lookup [ (unsigned)codon [ 1 ] ] ] ]
		[ genetic_code_idx [ dna_lookup [ (unsigned)codon [ 2 ] ] ] ];
}

/************************************************************/
char codon_to_cacid1 ( char *codon ) {

/*	given a codon of three bases, return the one letter code for the acid 
        on the complementary strand */

    return genetic_code 
	[ cgenetic_code_idx [ dna_lookup [ (unsigned)codon [ 2 ] ] ] ]
	    [ cgenetic_code_idx [ dna_lookup [ (unsigned)codon [ 1 ] ] ] ]
		[ cgenetic_code_idx [ dna_lookup [ (unsigned)codon [ 0 ] ] ] ];
}

/************************************************************/
char *codon_to_acid3 ( char *codon ) {

/*	given a codon of three bases, return the three letter code for the acid */
    char c;

    c = genetic_code    
	[ genetic_code_idx [ dna_lookup [ (unsigned)codon [ 0 ] ] ] ]
	    [ genetic_code_idx [ dna_lookup [ (unsigned)codon [ 1 ] ] ] ]
		[ genetic_code_idx [ dna_lookup [ (unsigned)codon [ 2 ] ] ] ];
    return three_letter_code(c);
}

/************************************************************/
char *codon_to_cacid3 ( char *codon ) {

/*	given a codon of three bases, return the three letter code for the acid */
    char c;

    c = genetic_code   
	[ cgenetic_code_idx [ dna_lookup [ (unsigned)codon [ 2 ] ] ] ]
	    [ cgenetic_code_idx [ dna_lookup [ (unsigned)codon [ 1 ] ] ] ]
		[ cgenetic_code_idx [ dna_lookup [ (unsigned)codon [ 0 ] ] ] ];
    return three_letter_code(c);
}


/************************************************************/
int write_cod_table ( FILE *out_file, double codon_table[4][4][4] ) {

/*	write out a codon table */

    int i,k;
    char bases[] = {"tcag"};
/*
12345678901234567890123456789012345678901234567890
      ===========================================
      F TTT  17. S TCT  16. Y TAT  10. C TGT   2.
      F TTC   3. S TCC   6. Y TAC   0. C TGC   0.
      L TTA  11. S TCA   6. * TAA   1. * TGA   2.
      L TTG   8. S TCG   2. * TAG   0. W TGG   6.
      ===========================================
      L CTT   7. P CCT   5. H CAT   4. R CGT   8.
      L CTC   5. P CCC   2. H CAC   3. R CGC   4.
      L CTA   2. P CCA   3. Q CAA   7. R CGA   3.
      L CTG   7. P CCG   1. Q CAG   7. R CGG   1.
      ===========================================
      I ATT  12. T ACT  10. N AAT  18. S AGT   3.
      I ATC   4. T ACC   4. N AAC   4. S AGC   4.
      I ATA   1. T ACA   1. K AAA  13. R AGA   3.
      M ATG   7. T ACG   2. K AAG   7. R AGG   1.
      ===========================================
      V GTT  14. A GCT   8. D GAT   7. G GGT   7.
      V GTC   2. A GCC   3. D GAC   7. G GGC   1.
      V GTA   7. A GCA   9. E GAA   6. G GGA   0.
      V GTG   1. A GCG   0. E GAG   6. G GGG   2.
      ===========================================
01234567890123456789012345678901234567890123456789012
      ===============================================
      F ttt    17 S tct    16 Y tat    10 C tgt     2
      L ttc     3 P tcc     6 H tac     0 R tgc     0
      I tta    11 T tca     6 N taa     1 S tga     2
      V ttg     8 A tcg     2 D tag     0 G tgg     6
      ===============================================
      F ctt     7 S cct     5 Y cat     4 C cgt     8
      L ctc     5 P ccc     2 H cac     3 R cgc     4
      I cta     2 T cca     3 N caa     7 S cga     3
      V ctg     7 A ccg     1 D cag     7 G cgg     1
      ===============================================
      L att    12 S act    10 * aat    18 * agt     3
      L atc     4 P acc     4 Q aac     4 R agc     4
      I ata     1 T aca     1 K aaa    13 R aga     3
      V atg     7 A acg     2 E aag     7 G agg     1
      ===============================================
      L gtt    14 S gct     8 * gat     7 W ggt     7
      L gtc     2 P gcc     3 Q gac     7 R ggc     1
      M gta     7 T gca     9 K gaa     6 R gga     0
      V gtg     1 A gcg     0 E gag     6 G ggg     2
      ===============================================
12345678901234567890123456789012345678901234567890123

*/


    for(i=0;i<4;i++) {
	fprintf(out_file,"      ===============================================\n");
	for(k=0;k<4;k++) {
	    fprintf(out_file,"      %c %c%c%c%6.0f %c %c%c%c%6.0f %c %c%c%c%6.0f %c %c%c%c%6.0f\n",
		    genetic_code[i][0][k],bases[i],bases[0],bases[k],
		    codon_table[i][0][k],
		    genetic_code[i][1][k],bases[i],bases[1],bases[k],
		    codon_table[i][1][k],
		    genetic_code[i][2][k],bases[i],bases[2],bases[k],
		    codon_table[i][2][k],
		    genetic_code[i][3][k],bases[i],bases[3],bases[k],
		    codon_table[i][3][k]);
	}
    }
	fprintf(out_file,"      ===============================================\n");
    return 1;
}

/************************************************************/
int write_screen_cod_table (double codon_table[4][4][4] ) {

/*	write out a codon table */

    int i,k;
    char bases[] = {"tcag"};
    
    for(i = 0; i < 4; i++) {
	vmessage("      ===============================================\n");
	for(k=0;k<4;k++) {
	    vmessage("      %c %c%c%c%6.0f %c %c%c%c%6.0f %c %c%c%c%6.0f %c %c%c%c%6.0f\n",
		    genetic_code[i][0][k],bases[i],bases[0],bases[k],
		    codon_table[i][0][k],
		    genetic_code[i][1][k],bases[i],bases[1],bases[k],
		    codon_table[i][1][k],
		    genetic_code[i][2][k],bases[i],bases[2],bases[k],
		    codon_table[i][2][k],
		    genetic_code[i][3][k],bases[i],bases[3],bases[k],
		    codon_table[i][3][k]);
	}
    }
	vmessage("      ===============================================\n");
    return 1;
}

/************************************************************/
int read_cod_table ( FILE *in_file, double codon_table[4][4][4] ) {

/*	read in a codon table */

    int i,k;
    char line[60];
/*
12345678901234567890123456789012345678901234567890
      ===========================================
      F TTT  17. S TCT  16. Y TAT  10. C TGT   2.
      F TTC   3. S TCC   6. Y TAC   0. C TGC   0.
      L TTA  11. S TCA   6. * TAA   1. * TGA   2.
      L TTG   8. S TCG   2. * TAG   0. W TGG   6.
      ===========================================
      L CTT   7. P CCT   5. H CAT   4. R CGT   8.
      L CTC   5. P CCC   2. H CAC   3. R CGC   4.
      L CTA   2. P CCA   3. Q CAA   7. R CGA   3.
      L CTG   7. P CCG   1. Q CAG   7. R CGG   1.
      ===========================================
      I ATT  12. T ACT  10. N AAT  18. S AGT   3.
      I ATC   4. T ACC   4. N AAC   4. S AGC   4.
      I ATA   1. T ACA   1. K AAA  13. R AGA   3.
      M ATG   7. T ACG   2. K AAG   7. R AGG   1.
      ===========================================
      V GTT  14. A GCT   8. D GAT   7. G GGT   7.
      V GTC   2. A GCC   3. D GAC   7. G GGC   1.
      V GTA   7. A GCA   9. E GAA   6. G GGA   0.
      V GTG   1. A GCG   0. E GAG   6. G GGG   2.
      ===========================================
012345678901234567890123456789012345678901234567890
*/
    if (fgets(line,55,in_file) == NULL ) return 0;

    /* read past all lines preceeding the line of ======== */
    while (1) {
	if (strncmp(&line[6], "==========", 10) == 0) {
	    break;
	}
	if (fgets(line, 55, in_file) == NULL) return 0;
    }

    if ( line[50] != '=' ) {
	for(i=0;i<4;i++) {
	    for(k=0;k<4;k++) {
		if (fgets(line,55,in_file) == NULL ) return 0;
		line[17] = line[28] = line[39] = line[50] = '\0';
		codon_table[i][0][k] = atof(&line[11]);
		codon_table[i][1][k] = atof(&line[22]);
		codon_table[i][2][k] = atof(&line[33]);
		codon_table[i][3][k] = atof(&line[44]);
	    }
	    if (fgets(line,55,in_file) == NULL ) return 0;
	}
    }
    else {
/*
      V gtg     1 A gcg     0 E gag     6 G ggg     2
      ===============================================
012345678901234567890123456789012345678901234567890123
         1         2         3         4         5
*/

	for(i=0;i<4;i++) {
	    for(k=0;k<4;k++) {
		if (fgets(line,55,in_file) == NULL ) return 0;
		line[18] = line[30] = line[42] = line[54] = '\0';
		codon_table[i][0][k] = atof(&line[11]);
		codon_table[i][1][k] = atof(&line[23]);
		codon_table[i][2][k] = atof(&line[35]);
		codon_table[i][3][k] = atof(&line[47]);
	    }
	    if (fgets(line,55,in_file) == NULL ) return 0;
	}
    }
    return 1;
}


/* 	argos average amino acid composition */
double av_protein_comp[] = { 
    8.3,/* A */
    1.7,/* C */
    5.3,/* D */
    6.2,/* E */
    3.9,/* F */
    7.2,/* G */
    2.2,/* H */
    5.2,/* I */
    5.7,/* K */
    9.0,/* L */
    2.4,/* M */
    4.4,/* N */
    5.1,/* P */
    4.0,/* Q */
    5.7,/* R */
    6.9,/* S */
    5.8,/* T */
    6.6,/* V */
    1.3,/* W */
    3.2,/* Y */
    0.0,/* * */
    0.0 /* - */
};


/*	generate codon table that gives average amino acid composition */
void gen_cods_from_ac ( double codon_table[4][4][4] ) {

    int i,j,k;
    size_t l;
    double sum_codons;

    /* loop for each amino acid type to give each of its codons
       the average amino acid frequency
    */

    for ( l=0;l<sizeof(one_letter);l++ ) {
	sum_codons = 0.0;
	for ( i=0;i<4;i++) {
	    for ( j=0;j<4;j++) {
		for ( k=0;k<4;k++) {
		    if ( genetic_code[i][j][k] == one_letter[l] ) {
			sum_codons += 1.0;
		    }
		}
	    }
	}
	if ( sum_codons > 0.0 ) {
	    sum_codons = av_protein_comp[l] / sum_codons;
	}
	else {
	    sum_codons = 0.0;
	}
	for ( i=0;i<4;i++) {
	    for ( j=0;j<4;j++) {
		for ( k=0;k<4;k++) {
		    if ( genetic_code[i][j][k] == one_letter[l] ) {
			codon_table[i][j][k] = sum_codons;
		    }
		}
	    }
	}
    }
}

/*	generate codon table assuming random base order */
void gen_cods_from_bc ( double codon_table[4][4][4], double base_comp[5] ) {

    int i,j,k,ig,jg,kg;

    for ( i=0;i<4;i++) {
	ig = genetic_code_idx [ i ];
	for ( j=0;j<4;j++) {
	    jg = genetic_code_idx [ j ];
	    for ( k=0;k<4;k++) {
		kg = genetic_code_idx [ k ];
		codon_table[ig][jg][kg] =
		    base_comp[i] * base_comp[j] * base_comp[k];
	    }
	}
    }
}

/*	codon_table to array64 */
void codon_table_64 ( double codon_table[4][4][4], double array64[64],
		     int job ) {

    int i,j,k,l;
    if ( job == TO_64 ) {

	for ( i=0,l=0;i<4;i++) {
	    for ( j=0;j<4;j++) {
		for ( k=0;k<4;k++,l++) {
		    array64[l]= codon_table[i][j][k];
		}
	    }
	}
    }
    else if ( job == FROM_64 ) {

	for ( i=0,l=0;i<4;i++) {
	    for ( j=0;j<4;j++) {
		for ( k=0;k<4;k++,l++) {
		    codon_table[i][j][k] = array64[l];
		}
	    }
	}
    }
}

/*	set sum of codon counts to total */
void scale_codon_table ( double codon_table[4][4][4], double total ) {

    double array64[64];
    codon_table_64 ( codon_table, array64, TO_64 );
    scale_double_array ( array64, 64, total );
    codon_table_64 ( codon_table, array64, FROM_64 );
}

/* convert codon table counts into a percentage */
void codon_table_percent(double codon_table[4][4][4])
{
    int i, j, k;
    size_t l;
    int sum_codons;

    for (l = 0; l < sizeof(one_letter); l++) {
	sum_codons = 0;
	for (i = 0; i < 4; i++) {
	    for (j = 0; j < 4; j++) {
		for (k = 0; k < 4; k++) {
		    if (genetic_code[i][j][k] == one_letter[l] ) {
			sum_codons += codon_table[i][j][k];
			
		    }
		}
	    }
	}
#ifdef DEBUG
	printf("sum_codons %d %c \n", sum_codons, one_letter[l]);
#endif
	for (i = 0; i < 4; i++) {
	    for (j = 0; j < 4; j++) {
		for (k = 0; k < 4; k++) {
		    if (genetic_code[i][j][k] == one_letter[l] ) {
			if (sum_codons > 0) {
#ifdef DEBUG
				printf("%f %d", codon_table[i][j][k], sum_codons);
#endif
			    codon_table[i][j][k] = 
				codon_table[i][j][k] / sum_codons * 100.0;

#ifdef DEBUG
				printf("after %f \n", codon_table[i][j][k]);
#endif

			} else {
			    codon_table[i][j][k] = 0.0;
			}
		    }
		}
	    }
	}
	
    }
}

void even_cods_per_acid(double codon_table[4][4][4])
{
  int i, j, k, num_codons;
  double sum_codons;
  size_t l;

  for (l = 0; l < sizeof(one_letter); l++) {
    num_codons = 0;
    sum_codons = 0.0;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
	for (k = 0; k < 4; k++) {
	  if (genetic_code[i][j][k] == one_letter[l] ) {
	    num_codons++;
	    sum_codons += codon_table[i][j][k];
	  }
	}
      }
    }
    if ( num_codons > 0 ) {
      sum_codons /= num_codons;
      for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	  for (k = 0; k < 4; k++) {
	    if (genetic_code[i][j][k] == one_letter[l] ) {
	      codon_table[i][j][k] = sum_codons;
	    }
	  }
	}
      }
    }
  }
}

void average_acid_comp(double codon_table[4][4][4])
{
  int i, j, k, l;
  double sum_codons;

  for (l = 0; l < 20; l++) {
    sum_codons = 0.0;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
	for (k = 0; k < 4; k++) {
	  if (genetic_code[i][j][k] == one_letter[l] ) {
	    sum_codons += codon_table[i][j][k];
	  }
	}
      }
    }
    if ( sum_codons > 0.0 ) {
      for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	  for (k = 0; k < 4; k++) {
	    if (genetic_code[i][j][k] == one_letter[l] ) {
	      codon_table[i][j][k] *= av_protein_comp[l]/sum_codons;
	    }
	  }
	}
      }
    }
  }
}

void third_pos_comp ( double codon_table[4][4][4] ) {

  /* make the base composition from the third posns of a codon table
   * and then make a codon table from this composition.
   * FIXME something fishy here so it is not used at present
   */

  int i, j, k, num_codons;
  size_t l;
  double base_comp[5], total;

  j = 5;
  for (i=0;i<j;i++) {
    base_comp[i] =  0.0;
  }
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      for (k=0;k<4;k++) {
	base_comp[k] += codon_table[i][j][k];
      }
    }
  }
  j = 5;
  for (i=0,total=0.0;i<j;i++) {
    total += base_comp[i];
  }
  if ( total > DBL_EPSILON ) {
    for (i=0;i<j;i++) {
      base_comp[i] /= total;
    }
  }
  for (l = 0; l < sizeof(one_letter); l++) {
    num_codons = 0;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        for (k = 0; k < 4; k++) {
          if (genetic_code[i][j][k] == one_letter[l] ) {
            num_codons++;
          }
        }
      }
    }
    if ( num_codons > 0 ) {
      total = 100.0/num_codons;
      for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
          for (k = 0; k < 4; k++) {
            if (genetic_code[i][j][k] == one_letter[l] ) {
              codon_table[i][j][k] = total * base_comp[k];
            }
          }
        }
      }
    }
  }
}

int read_genetic_code ( FILE *in_file, char code_table[5][5][5] ) {

    /*	read in a codon table to get genetic code*/

    int i,j,k;
    char line[60];
/*
12345678901234567890123456789012345678901234567890
      ===========================================
      F TTT  17. S TCT  16. Y TAT  10. C TGT   2.
      F TTC   3. S TCC   6. Y TAC   0. C TGC   0.
      L TTA  11. S TCA   6. * TAA   1. * TGA   2.
      L TTG   8. S TCG   2. * TAG   0. W TGG   6.
      ===========================================
      L CTT   7. P CCT   5. H CAT   4. R CGT   8.
      L CTC   5. P CCC   2. H CAC   3. R CGC   4.
      L CTA   2. P CCA   3. Q CAA   7. R CGA   3.
      L CTG   7. P CCG   1. Q CAG   7. R CGG   1.
      ===========================================
      I ATT  12. T ACT  10. N AAT  18. S AGT   3.
      I ATC   4. T ACC   4. N AAC   4. S AGC   4.
      I ATA   1. T ACA   1. K AAA  13. R AGA   3.
      M ATG   7. T ACG   2. K AAG   7. R AGG   1.
      ===========================================
      V GTT  14. A GCT   8. D GAT   7. G GGT   7.
      V GTC   2. A GCC   3. D GAC   7. G GGC   1.
      V GTA   7. A GCA   9. E GAA   6. G GGA   0.
      V GTG   1. A GCG   0. E GAG   6. G GGG   2.
      ===========================================
012345678901234567890123456789012345678901234567890
*/

    for(i=0;i<5;i++) {
	for(j=0;j<5;j++) {
	    for(k=0;k<5;k++) {
		code_table[i][j][k] = '-';
	    }
	}
    }
    if (fgets(line,55,in_file) == NULL ) return 0;

    /* read past all lines preceeding the line of ======== */
    while (1) {
	if (strncmp(&line[6], "==========", 10) == 0) {
	    break;
	}
	if (fgets(line, 55, in_file) == NULL) return 0;
    }

    if ( line[50] != '=' ) {
	for(i=0;i<4;i++) {
	    for(k=0;k<4;k++) {
		if (fgets(line,55,in_file) == NULL ) return 0;
		code_table[i][0][k] = line[6];
		code_table[i][1][k] = line[17];
		code_table[i][2][k] = line[28];
		code_table[i][3][k] = line[39];
	    }
	    if (fgets(line,55,in_file) == NULL ) return 0;
	}
    }
    else {
/*
      V gtg     1 A gcg     0 E gag     6 G ggg     2
      ===============================================
012345678901234567890123456789012345678901234567890123
         1         2         3         4         5
*/

	for(i=0;i<4;i++) {
	    for(k=0;k<4;k++) {
		if (fgets(line,55,in_file) == NULL ) return 0;
		code_table[i][0][k] = line[6];
		code_table[i][1][k] = line[18];
		code_table[i][2][k] = line[30];
		code_table[i][3][k] = line[42];
	    }
	    if (fgets(line,55,in_file) == NULL ) return 0;
	}
    }

    /* now find which codons we can still decode when they end with '-' */

    for(i=0;i<4;i++) {
	for(j=0;j<4;j++) {
	    if ( 
		(code_table[i][j][0] == code_table[i][j][1]) &&
		(code_table[i][j][0] == code_table[i][j][2]) &&
		(code_table[i][j][0] == code_table[i][j][3]) )
		 code_table[i][j][4] =  code_table[i][j][0];
	}
    }
    return 1;
}


/************************************************************/
int write_screen_genetic_code (char code_table[5][5][5] ) {

/*	write out a genetic code table */

    int i,k;
    char bases[] = {"tcag-"};
    
    for(i = 0; i < 4; i++) {
	vmessage("      ===============================================\n");
	for(k = 0; k < 4; k++) {
	    vmessage("      %c %c%c%-7c %c %c%c%-7c %c %c%c%-7c %c %c%c%-7c\n",
		     genetic_code[i][0][k],bases[i],bases[0],bases[k],
		     genetic_code[i][1][k],bases[i],bases[1],bases[k],
		     genetic_code[i][2][k],bases[i],bases[2],bases[k],
		     genetic_code[i][3][k],bases[i],bases[3],bases[k]);
	}
    }
    vmessage("      ===============================================\n");
    return 1;
}

int read_global_genetic_code (FILE *in_file) {
  return read_genetic_code(in_file, genetic_code);
}

char (*get_global_genetic_code())[5][5] {
    return genetic_code;
}

int *get_genetic_code_idx(int complemented) {
    return complemented ? cgenetic_code_idx : genetic_code_idx;
}

/*
 * Loads a genetic code with a specific index as defined in the EMBL
 * feature table /transl_table qualifier.
 *
 * Returns 0 for success,
 *	  -1 for failure.
 */
int load_genetic_code_number(int index) {
    char buf[1024];
    char *env = getenv("STADTABL");
    FILE *fp;
    int err;

    if (!env)
	return -1;
    sprintf(buf, "%s/gcodes/code_%d", env, index);

    if (NULL == (fp = fopen(buf, "r")))
	return -1;

    err = read_global_genetic_code(fp);
    fclose(fp);

    return err;
}

/*
 * Converts an amino-acid 3-letter code into a 1-letter code following the EMBL
 * naming conventions. These are pretty standard except TERM indicates the
 * stop codon which we refer to as '*'.
 *
 * Any unknown codes are returned as '-'.
 */
int embl_aa_three2one(char *aa) {
    int i;

    if (strncmp(aa, "TERM", 4) == 0) {
	return '*';
    } else {
	for (i = 0; i < sizeof(one_letter); i++) {
	    if (strncmp(three_letter[i], aa, 3) == 0) {
		return one_letter[i];
	    }
	}
    }

    return '-';
}

/*
 * Adjusts the genetic code according to an edit specified in EMBL /codon
 * qualifier format.
 * Eg (seq:"tga",aa:Trp)
 *
 * Returns 0 for success,
 *        -1 for failure.
 */
int edit_genetic_code(char *change) {
    char *na, *aa, *cp;
    int aa_char;

    /* Parse */
    if (cp = strchr(change, ':')) {
	na = cp+1;
	if (*na == '"')
	    na++;
    } else {
	return -1;
    }

    if (cp = strchr(na, ':')) {
	aa = cp+1;
	if (*aa == '"')
	    aa++;
    } else {
	return -1;
    }

    /* Check for valid nucleic acid codes */
    if (!legal_codon(na))
	return -1;

    /* Find 1 letter AA code */
    aa_char = embl_aa_three2one(aa);

    /* Edit table */
    genetic_code
	[genetic_code_idx[dna_lookup[na[0]]]]
	[genetic_code_idx[dna_lookup[na[1]]]]
	[genetic_code_idx[dna_lookup[na[2]]]] = aa_char;

    return 0;
}
