/*	routines for handling the genetic code

	We need to be able to have a standard code and allow changes
	to make a current code.
	We need to translate dna to amino acids.
	We need to count codons.
	We need to use indexes into the codon tables for several methods.

	The genetic code is written tcag order, our lookups are acgt.
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "dna_utils.h"
/* #include "readpam.h" */


/************************************************************/
char genetic_code [5][5][5] = {

    'F','F','L','L','-',
    'S','S','S','S','S',
    'Y','Y','*','*','-',
    'C','C','*','W','-',
    '-','-','-','-','-',

    'L','L','L','L','L',
    'P','P','P','P','P',
    'H','H','Q','Q','-',
    'R','R','R','R','R',
    '-','-','-','-','-',

    'I','I','I','M','-',
    'T','T','T','T','T',
    'N','N','K','K','-',
    'S','S','R','R','-',
    '-','-','-','-','-',

    'V','V','V','V','V',
    'A','A','A','A','A',
    'D','D','E','E','-',
    'G','G','G','G','G',
    '-','-','-','-','-',

    '-','-','-','-','-',
    '-','-','-','-','-',
    '-','-','-','-','-',
    '-','-','-','-','-',
    '-','-','-','-','-' 
};

int genetic_code_index[5] = {2,1,3,0,4};
int cgenetic_code_index[5] = {0,3,1,2,4};

char *three_letter_code ( char one_letter_code ) {

/*	Input character representing an amino acid in 1 letter code
	Output its three letter code equivalent
	Note we assume incoming character is upper case
*/
    char one_letter[] = {"ACDEFGHIKLMNPQRSTVWY*-"};
    char *three_letter[] = {
     "Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys",
     "Leu","Met","Asn","Pro","Gln","Arg","Ser","Thr","Val",
     "Trp","Tyr","***","---"};
    int i;

    one_letter_code = toupper(one_letter_code);

    for (i=0;i<sizeof(one_letter)-1;i++){
	if ( one_letter_code == one_letter[i] ) return three_letter[i];
    }
    return "   ";
}

/************************************************************/
int codon_table_index ( char *codon ) {

/*	given a codon of three bases, return an index into 
	a 64 element codon table */

/* check for illegal symbol. If found return 64 Calling routines please note*/

    int i;
    for(i=0;i<3;i++) {
	if ( dna_lookup [ codon [ i ] ] == 4 ) return 64;
    }

    return   genetic_code_index [ dna_lookup [ codon [ 0 ] ] ] |
            (genetic_code_index [ dna_lookup [ codon [ 1 ] ] ]) <<2 |
            (genetic_code_index [ dna_lookup [ codon [ 2 ] ] ]) <<4;
}

/************************************************************/
char codon_to_acid1 ( char *codon ) {

/*	given a codon of three bases, return the one letter code for the acid */

    return genetic_code [ genetic_code_index [ dna_lookup [ codon [ 0 ] ] ] ]
			[ genetic_code_index [ dna_lookup [ codon [ 1 ] ] ] ]
		        [ genetic_code_index [ dna_lookup [ codon [ 2 ] ] ] ];
}

/************************************************************/
char codon_to_cacid1 ( char *codon ) {

/*	given a codon of three bases, return the one letter code for the acid 
        on the complementary strand */

    return genetic_code [ cgenetic_code_index [ dna_lookup [ codon [ 2 ] ] ] ]
			[ cgenetic_code_index [ dna_lookup [ codon [ 1 ] ] ] ]
		        [ cgenetic_code_index [ dna_lookup [ codon [ 0 ] ] ] ];
}

/************************************************************/
char *codon_to_acid3 ( char *codon ) {

/*	given a codon of three bases, return the three letter code for the acid */
    char c;

    c = genetic_code    [ genetic_code_index [ dna_lookup [ codon [ 0 ] ] ] ]
			[ genetic_code_index [ dna_lookup [ codon [ 1 ] ] ] ]
		        [ genetic_code_index [ dna_lookup [ codon [ 2 ] ] ] ];
    return three_letter_code(c);
}

/************************************************************/
char *codon_to_cacid3 ( char *codon ) {

/*	given a codon of three bases, return the three letter code for the acid */
    char c;

    c = genetic_code   [ cgenetic_code_index [ dna_lookup [ codon [ 2 ] ] ] ]
		       [ cgenetic_code_index [ dna_lookup [ codon [ 1 ] ] ] ]
		       [ cgenetic_code_index [ dna_lookup [ codon [ 0 ] ] ] ];
    return three_letter_code(c);
}

/************************************************************/
int read_codon_table ( FILE *in_file, double codon_table[4][4][4] ) {

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


    for(i=0;i<4;i++) {
	if (fgets(line,55,in_file) == NULL ) return 0;
	for(k=0;k<4;k++) {
	    if (fgets(line,55,in_file) == NULL ) return 0;
	    line[17] = line[28] = line[39] = line[50] = '\0';
	    codon_table[i][0][k] = atof(&line[11]);
	    codon_table[i][1][k] = atof(&line[22]);
	    codon_table[i][2][k] = atof(&line[33]);
	    codon_table[i][3][k] = atof(&line[44]);
	}
    }
    return 1;
}

/************************************************************/
int write_codon_table ( FILE *out_file, double codon_table[4][4][4] ) {

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
		    genetic_code[k][0][i],bases[i],bases[0],bases[k],
		    codon_table[i][0][k],
		    genetic_code[k][1][i],bases[i],bases[1],bases[k],
		    codon_table[i][1][k],
		    genetic_code[k][2][i],bases[i],bases[2],bases[k],
		    codon_table[i][2][k],
		    genetic_code[k][3][i],bases[i],bases[3],bases[k],
		    codon_table[i][3][k]);
	}
    }
	fprintf(out_file,"      ===============================================\n");
    return 1;
}

/************************************************************/
void zero_codon_table ( double codon_table[4][4][4] ) {

/*	zero a codon table */

    int i,j,k;
    for(i=0;i<4;i++) {
	for(j=0;j<4;j++) {
	    for(k=0;k<4;k++) {
		codon_table[i][j][k] = 0.0;
	    }
	}
    }
}

/************************************************************/
#if 0
main() {
    char seq[] = "aatgctcgataggctcgtatatgcgcgctattatatatgcgcg";
    char codon[] = {'t','c','a'};
    int i;
    double codon_table[4][4][4];
    FILE *in_file;
    char *file_name = "codon_table";
    in_file = fopen(file_name,"r");
    i = read_codon_table(in_file,codon_table);
    i = write_codon_table(stdout,codon_table);
    (void) zero_codon_table(codon_table);
    i = write_codon_table(stdout,codon_table);
    sbl();
/*    codon[1] = 't';
    codon[2] = 't';*/
    codon[0] = 'g';
    codon[1] = 'g';
    codon[2] = 'g';
    printf(" index %d\n",codon_table_index(codon));
    printf("tca %c\n",acid_from_codon(codon));
    printf("tca %c\n",cacid_from_codon(codon));
    codon[0] = 'a';
    printf("tca %c\n",cacid_from_codon(codon));
    codon[1] = 'a';
    printf("tca %c\n",cacid_from_codon(codon));
    printf("%c\n",genetic_code[1][2][3]);
    printf("%c\n",genetic_code[0][1][4]);
    printf("%c\n",genetic_code
	   [genetic_code_index[1]]
	   [genetic_code_index[2]]
	   [genetic_code_index[3]]);
    printf("%s\n",three_letter_code('W'));
    printf("%d %d %d\n",strlen(seq),3*strlen(seq),3*(strlen(seq)/3));
    for(i=0;i<3*(strlen(seq)/3);i+=3) {
	printf("%s",three_letter_codon(&seq[i]));
    }
	printf("\n");
    for(i=0;i<3*(strlen(seq)/3);i+=3) {
	printf("%c  ",acid_from_codon(&seq[i]));
    }
	printf("\n");
}
#endif
