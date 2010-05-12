/* DNA (and protein lookup!) utility routines */
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "os.h"
#include "FtoC.h"
#include "dna_utils.h"
#include "misc.h"


/* start of new stuff */

int char_set_size = 0, *char_lookup = 0, *char_match = 0;

#define MAX_CHAR_SET_SIZE 24
#define DNA 1
static char complementary_base[256];

int dna_lookup[256] = {0};
int dna_match[256] = {0};
int iubc_lookup[256] = {0};
int hash4_lookup[256] = {0}; /* initialised along with iubc_lookup */

/* 7/1/99 johnt - must initialise globals to force export with WINNT */
int unknown_char=0;

/* hard code the unknown_char value for iubc table */
/* a = sequence, b = string (containing iubc symbols) */
#define IUBC_MATCH(a,b) (((iubc_lookup[(int)a] < 16) && iubc_match[iubc_lookup[(int)b]][iubc_lookup[(int)a]]) ? 1 : 0)
#define IUBC_MISMATCH(a,b) (((iubc_lookup[(int)a] < 16) && iubc_match[iubc_lookup[(int)b]][iubc_lookup[(int)a]]) ? 0 : 1)

/************************************************************/

void set_dna_lookup() {

/* 	set up table of values for permitted dna characters */

    int i;

    for (i=0;i<256;i++) dna_lookup[i] = 4;
    for (i=0;i<256;i++) dna_match[i] = i+256;
    for (i=0;i<256;i++) complementary_base[i] = i;

    dna_lookup['a'] = 0;
    dna_lookup['c'] = 1;
    dna_lookup['g'] = 2;
    dna_lookup['t'] = 3;
    dna_lookup['A'] = 0;
    dna_lookup['C'] = 1;
    dna_lookup['G'] = 2;
    dna_lookup['T'] = 3;
    dna_lookup['U'] = 3;
    dna_lookup['u'] = 3;

    /* the following does not know if rna or dna so we complement
       u to a, but a to t */

    complementary_base['a'] = 't';
    complementary_base['c'] = 'g';
    complementary_base['g'] = 'c';
    complementary_base['t'] = 'a';
    complementary_base['u'] = 'a';
    complementary_base['A'] = 'T';
    complementary_base['C'] = 'G';
    complementary_base['G'] = 'C';
    complementary_base['T'] = 'A';
    complementary_base['U'] = 'A';

    complementary_base['n'] = 'n';
    complementary_base['-'] = '-';
    complementary_base['b'] = 'v';
    complementary_base['d'] = 'h';
    complementary_base['h'] = 'd';
    complementary_base['k'] = 'm';
    complementary_base['m'] = 'k';
    complementary_base['r'] = 'y';
    complementary_base['s'] = 's';
    complementary_base['v'] = 'b';
    complementary_base['w'] = 'w';
    complementary_base['y'] = 'r';

    complementary_base['B'] = 'V';
    complementary_base['D'] = 'H';
    complementary_base['H'] = 'D';
    complementary_base['K'] = 'M';
    complementary_base['M'] = 'K';
    complementary_base['R'] = 'Y';
    complementary_base['S'] = 'S';
    complementary_base['V'] = 'B';
    complementary_base['W'] = 'W';
    complementary_base['Y'] = 'R';


    /* dna_match matches the valid characters in dna_lookup 
       such as T = t = U = u, otherwise any character only
       matches itself */

    for (i=0;i<256;i++) {
	if ( dna_lookup[i] != 4 ) dna_match[i] = dna_lookup[i];
    }
    /* note dna_match is not used anywhere! */
}

/************************************************************/

void set_iubc_lookup() {

/* 	set up table of index values for iubc characters
	and hash length 4 values: hash4_lookup[]
*/

    int i;

    for (i=0;i<256;i++) iubc_lookup[i] = 16;

    iubc_lookup['a'] = 0;
    iubc_lookup['c'] = 1;
    iubc_lookup['g'] = 2;
    iubc_lookup['t'] = 3;
    iubc_lookup['u'] = 3;
    iubc_lookup['A'] = 0;
    iubc_lookup['C'] = 1;
    iubc_lookup['G'] = 2;
    iubc_lookup['T'] = 3;
    iubc_lookup['U'] = 3;
    iubc_lookup['r'] = 4;
    iubc_lookup['y'] = 5;
    iubc_lookup['m'] = 6;
    iubc_lookup['k'] = 7;
    iubc_lookup['s'] = 8;
    iubc_lookup['w'] = 9;
    iubc_lookup['b'] = 10;
    iubc_lookup['d'] = 11;
    iubc_lookup['h'] = 12;
    iubc_lookup['v'] = 13;
    iubc_lookup['n'] = 14;
    iubc_lookup['-'] = 15;
    iubc_lookup['R'] = 4;
    iubc_lookup['Y'] = 5;
    iubc_lookup['M'] = 6;
    iubc_lookup['K'] = 7;
    iubc_lookup['S'] = 8;
    iubc_lookup['W'] = 9;
    iubc_lookup['B'] = 10;
    iubc_lookup['D'] = 11;
    iubc_lookup['H'] = 12;
    iubc_lookup['V'] = 13;
    iubc_lookup['N'] = 14;

    for (i=0;i<256;i++) hash4_lookup[i] = 0; /* NB: unknown = a !!! */

    hash4_lookup['a'] = 0;
    hash4_lookup['c'] = 1;
    hash4_lookup['g'] = 2;
    hash4_lookup['t'] = 3;
    hash4_lookup['u'] = 3;
    hash4_lookup['A'] = 0;
    hash4_lookup['C'] = 1;
    hash4_lookup['G'] = 2;
    hash4_lookup['T'] = 3;
    hash4_lookup['U'] = 3;

}


static int iubc_match[17][17] = {

/* table of definite matches between symbols in the left column
   and those in the top row. Ie the table is not symmetrical so
   that string symbols (like restriction site data) should be in
   the first dimension, and the symbols from the sequence being
   searched, in the second, as shown below:
   iubc_match [ iubc_lookup [ word [ j ] ] ] [ iubc_lookup [ seq [ i ] ] ];
*/

/*      a c g t r y m k s w b d h v n - ? */

/* a */ 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
/* c */ 0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
/* g */ 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
/* t */ 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
/* r */ 1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
/* y */ 0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,
/* m */ 1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
/* k */ 0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,
/* s */ 0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
/* w */ 1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,
/* b */ 0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,0,
/* d */ 1,0,1,1,1,0,0,1,0,1,0,1,0,0,0,0,0,
/* h */ 1,1,0,1,0,1,1,0,0,1,0,0,1,0,0,0,0,
/* v */ 1,1,1,0,1,0,1,0,1,0,0,0,0,1,0,0,0,
/* n */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
/* - */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
/* ? */ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

    };

/************************************************************/

int protein_lookup[256] = {0};
static int protein_match[256];


void set_protein_lookup() {

/* 	set up table of values for permitted protein characters */

    int i;

    for (i=0;i<256;i++) protein_lookup[i] = 22;
    for (i=0;i<256;i++) protein_match[i] = i+256;

    protein_lookup['a'] = 0;
    protein_lookup['b'] = 1;
    protein_lookup['c'] = 2;
    protein_lookup['d'] = 3;
    protein_lookup['e'] = 4;
    protein_lookup['f'] = 5;
    protein_lookup['g'] = 6;
    protein_lookup['h'] = 7;
    protein_lookup['i'] = 8;
    protein_lookup['k'] = 9;
    protein_lookup['l'] = 10;
    protein_lookup['m'] = 11;
    protein_lookup['n'] = 12;
    protein_lookup['p'] = 13;
    protein_lookup['q'] = 14;
    protein_lookup['r'] = 15;
    protein_lookup['s'] = 16;
    protein_lookup['t'] = 17;
    protein_lookup['v'] = 18;
    protein_lookup['w'] = 19;
    protein_lookup['x'] = 22;
    protein_lookup['y'] = 20;
    protein_lookup['z'] = 21;
    protein_lookup['A'] = 0;
    protein_lookup['B'] = 1;
    protein_lookup['C'] = 2;
    protein_lookup['D'] = 3;
    protein_lookup['E'] = 4;
    protein_lookup['F'] = 5;
    protein_lookup['G'] = 6;
    protein_lookup['H'] = 7;
    protein_lookup['I'] = 8;
    protein_lookup['K'] = 9;
    protein_lookup['L'] = 10;
    protein_lookup['M'] = 11;
    protein_lookup['N'] = 12;
    protein_lookup['P'] = 13;
    protein_lookup['Q'] = 14;
    protein_lookup['R'] = 15;
    protein_lookup['S'] = 16;
    protein_lookup['T'] = 17;
    protein_lookup['V'] = 18;
    protein_lookup['W'] = 19;
    protein_lookup['X'] = 22;
    protein_lookup['Y'] = 20;
    protein_lookup['Z'] = 21;
    protein_lookup['*'] = 23;
    protein_lookup['-'] = 24;   /* is X = - sensible? */
    for (i=0;i<256;i++) {
	if ( protein_lookup[i] < 22 ) protein_match[i] = protein_lookup[i]; 
    }
    /* note that we allow b=b and z=z ! */
    /* note also that protein_match is not used anywhere! */
}

int *get_protein_lookup(void) {
    return protein_lookup;
}

/************************************************************/

/*
 *
 * seq_type == 1 for DNA
 * seq_type != 1 for protein
 *
 */

void set_char_set(int seq_type) {
  /* note that the arrays dna_match and protein_match are not used anywhere
   * now that char_match is set to the lookup arrays
  */
    if ( DNA == seq_type ) {
	char_set_size = 5;
	char_lookup = dna_lookup;
	char_match = dna_lookup;
        unknown_char = 4;
    }
    else {
	char_set_size = 25;
	char_lookup = protein_lookup;
	char_match = protein_lookup;
        unknown_char = 22;
    }
}


/************************************************************/
int word_match( char *seq, int seq_pos, int seq_len, char *word, int word_len) {

/* Does word match seq at seq_pos */

    register int i,j;

    for ( i = seq_pos, j = 0;
	 i < seq_len && j < word_len &&
	 SEQ_MATCH( seq[i], word[j] ); 
	 i++,j++);

    return ( j == word_len ) ? 1 : 0;
}
/************************************************************/
int iubc_word_match( char *seq, int seq_pos, int seq_len, char *word, int word_len) {

/* Does word in iubc codes match seq STARTING at seq_pos */

    register int i,j;

    for ( i = seq_pos, j = 0;
	 i < seq_len && j < word_len &&
	 iubc_match [ iubc_lookup [ (unsigned) word[j] ] ]
	    [ iubc_lookup [ (unsigned) seq[i] ] ];
	 i++,j++);

    return ( j == word_len ) ? 1 : 0;
}


int iubc_word_match_padded( char *seq, int seq_pos, int seq_len, char *word, int word_len) {

/* Does word in iubc codes match seq STARTING at seq_pos */

    register int i,j;

    /* Allow for pads in seq, but not the word we are searching */
    for ( i = seq_pos, j = 0; i < seq_len && j < word_len; i++) {
	if (seq[i] == '*')
	    continue;

	if (!iubc_match [ iubc_lookup [ (unsigned) word[j] ] ]
	    [ iubc_lookup [ (unsigned) seq[i] ] ])
	    break;

	j++;
    }

    return ( j == word_len ) ? 1 : 0;
}


/************************************************************/


int match_len ( char *seq1, int seq1_start, int seq1_len, char *seq2,
	       int seq2_start, int seq2_len ) {

/*	find length of match between seq1 and seq2 
	starting at seq1_start and seq2_start */

    register int i,j;

    for ( i = seq1_start, j = seq2_start;
	 i < seq1_len && j < seq2_len &&
	 SEQ_MATCH( seq1[i], seq2[j] ); 
	 i++,j++);

    return i - seq1_start;
}


/************************************************************/
int literal_mismatch(char a, char b)
{
  
  if (a == b || toupper(a) == b || (a) == toupper(b))
    return 0;
  else
    return 1;

}

/*
 * Find all positions in seq at which string has >= min_match matching
 * characters.
 */
int iubc_inexact_match (char *seq, 
			int seq_len, 
			char *string, 
			int string_len,
			int min_match,
			int use_iub_code,
			int *match, 
			int *score, 
			int max_matches) 
{
    int i, j, k, l, mismatch, max_mismatch, n_matches;
    int *compar, slen;

    compar = (int *)xmalloc(256 * string_len * sizeof(int));
    if (NULL == compar)
	return 0;
 
    /* use iub lookup table */
    if (use_iub_code) {
      for (i = 0; i < 256; i++) {
	for (j = 0; j < string_len; j++) {
	  compar[i + j*256] = IUBC_MISMATCH(i, string[j]);
	} 
      }
    } else { 
      /* do a literal search eg "n" == "n" */
      for (i = 0; i < 256; i++) {
	for (j = 0; j < string_len; j++) {
	  compar[i + j*256] = literal_mismatch((char)i, string[j]);
	}
      }
    }

    max_mismatch = string_len - min_match + 1;
    n_matches = 0;
    l = seq_len - string_len + 1;
    slen = 256 * string_len;

    for (i=0; i<l; i++ ) {
	for (j=0, k=i, mismatch=max_mismatch; j < slen; j+=256,k++) {
	    if (compar[seq[k] + j]) {
		if (--mismatch <= 0)
		    break;
	    }
	}
	
	if (mismatch > 0) {
	    if (n_matches < max_matches) {
		match[n_matches] = i;
		score[n_matches] = string_len - (max_mismatch - mismatch);
		n_matches++;

	    } else {
		/* make positions start at 1 */
		for (i=0; i < max_matches; i++) {
		    match[i]++;
		}
		xfree(compar);
		return -1; /* out of match storage */
	    }
	}
    }

    /* make positions start at 1 */
    for (i=0; i < n_matches; i++) {
	match[i]++;
    }
    xfree(compar);

    return n_matches;
}
/************************************************************/

/*
 * Find all positions in seq at which string has >= min_match matching
 * characters.
 */
int inexact_match (char *seq, int seq_len, char *string, int string_len,
		   int min_match, int *match, int *score, int max_matches) {

    int i, j, k, l, mismatch, max_mismatch, n_matches;
    int *compar, slen;

    compar = (int *)xmalloc(256 * string_len * sizeof(int));
    if (NULL == compar)
	return 0;

    for (i = 0; i < 256; i++) {
	for (j = 0; j < string_len; j++) {
	    compar[i + j*256] = SEQ_MISMATCH(i, string[j]);
	}
    }

    max_mismatch = string_len - min_match + 1;
    n_matches = 0;
    l = seq_len - string_len + 1;
    slen = 256 * string_len;

    for (i=0; i<l; i++ ) {
	for (j=0, k=i, mismatch=max_mismatch; j < slen; j+=256,k++) {
	    if (compar[seq[k] + j]) {
		if (--mismatch <= 0)
		    break;
	    }
	}
	
	if (mismatch > 0) {
	    if (n_matches < max_matches) {
		match[n_matches] = i;
		score[n_matches] = string_len - (max_mismatch - mismatch);
		n_matches++;

	    } else {
		/* make positions start at 1 */
		for (i=0; i < max_matches; i++) {
		    match[i]++;
		}
		xfree(compar);
		return -1; /* out of match storage */
	    }
	}
    }

    /* make positions start at 1 */
    for (i=0; i < n_matches; i++) {
	match[i]++;
    }
    xfree(compar);

    return n_matches;
}

/*
 * Find the best inexact match of 'string' in seq. See inexact_match for
 * more details.
 *
 * Returns the best match value.
 */
int best_inexact_match(char *seq, int seq_len, char *string, int string_len,
		       int *match) {
    int i, j, k, l, mismatch, max_mismatch = string_len;
    int *compar, slen;

    compar = (int *)xmalloc(256 * string_len * sizeof(int));
    if (NULL == compar)
	return 0;

    for (i = 0; i < 256; i++) {
	for (j = 0; j < string_len; j++) {
	    compar[i + j*256] = SEQ_MISMATCH(i, string[j]);
	}
    }

    l = seq_len - string_len + 1;
    slen = 256 * string_len;

    for (i=0; i<l; i++ ) {
	for (j=0, k=i, mismatch=max_mismatch; j < slen; j += 256, k++) {
	    if (compar[seq[k] + j]) {
		if (--mismatch <= 0)
		    break;
	    }
	}
	
	if (mismatch > 0) {
	    max_mismatch -= mismatch;
	    if (match)
		*match = i+1;

	    if (max_mismatch == 0)
		break;
	}
    }

    xfree(compar);
    return string_len - max_mismatch;
}

/* end of new stuff */



/************************************************************/


char complement_base (char base) {
    return complementary_base[(unsigned char)base];
}

void complement_seq ( char *seq, int seq_len ) {

    int i, middle, j;
    char temp;

    middle = seq_len/2;
    for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
	temp = complementary_base [ (unsigned char) seq[i] ];
	seq[i] = complementary_base [ (unsigned char) seq[j] ];
	seq[j] = temp;
    }

    if ( seq_len % 2 )
      seq[middle] = complementary_base [ (unsigned char) seq[middle] ];
}

/************************************************************/


/*
 * These two combine to equal the complement_seq above.
 */
void reverse_dna( char *seq, int seq_len ) {
    register int i, middle, j;
    char temp;

    middle = seq_len/2;
    for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
	temp = seq[i];
	seq[i] = seq[j];
	seq[j] = temp;
    }
}

void complement_dna(char *seq, int seq_len) {
    register int i;

    for ( i=0; i<seq_len; i++ ) {
	seq[i] = complementary_base[ (unsigned char) seq[i] ];
    }
}

/************************************************************/


void copy_seq ( char *copy, char *original, int seq_len ) {

/*	copy a sequence  */

    int i;

    for ( i=0; i < seq_len; i++) *copy++ = *original++;
}

/************************************************************/
/* count identities (ignoring case) between two aligned strings */

int same_char(char c1, char c2) {
    if (toupper(c1) == toupper(c2))
	return 1;

    if ((c1 == '*' || c1 == ',' || c1 == '.') &&
	(c2 == '*' || c2 == ',' || c2 == '.'))
	return 1;

    return 0;
}

int identities ( char *seq1, char *seq2 ) {
    int j,k,seq_len;
    seq_len = strlen(seq1);
    for (k=0,j=0;k<seq_len;k++) {
	j += same_char(seq1[k], seq2[k]);
    }
    return j;
}

static int iubc_identities ( char *seq1, char *seq2 ) {
    int j,k,seq_len;
    seq_len = strlen(seq1);
    for (k = 0, j = 0; k < seq_len; k++) {
	j += IUBC_MATCH(seq2[k], seq1[k]);
    }
    return j;
}


/* routine to display a message */

void info_ ( char *fstring, int_fl s_len ) {
    char word[1024];

    Fstr2Cstr(fstring, s_len, word, 1024);
    vmessage("%s\n", word);
}


/* routine to receive an alignment from fortran and turn it into C
   strings and variables */

int forta_ ( char *fseq1, char *fseq2, int *slen, char *fname1, char *fname2,
	      int *nlen, int *fpos1, int *fpos2, char *ftitle, int *tlen,
	    int_fl seq1_l, int_fl seq2_l, int_fl name1_l, int_fl name2_l,
	    int_fl title_l)

{
    char *seq1, *seq2, *name1, *name2, *title;
    int pos1, pos2;

    if ( ! (seq1 = (char *) xmalloc ( sizeof(char) * (*slen) + 1))) {
		return -1;
	    }
    if ( ! (seq2 = (char *) xmalloc ( sizeof(char) * (*slen) + 1))) {
		return -1;
	    }
    if ( ! (name1 = (char *) xmalloc ( sizeof(char) * (*nlen) + 1))) {
		return -1;
	    }
    if ( ! (name2 = (char *) xmalloc ( sizeof(char) * (*nlen) + 1))) {
		return -1;
	    }

    if ( ! (title = (char *) xmalloc ( sizeof(char) * (*tlen) + 1))) {
		return -1;
	    }

    copy_seq ( seq1, fseq1, *slen);
    copy_seq ( seq2, fseq2, *slen);
    seq1[*slen] = '\0';
    seq2[*slen] = '\0';
    copy_seq ( name1, fname1, *nlen);
    copy_seq ( name2, fname2, *nlen);
    name1[*nlen] = '\0';
    name2[*nlen] = '\0';
    copy_seq ( title, ftitle, *tlen);
    title[*tlen] = '\0';
    pos1 = *fpos1;
    pos2 = *fpos2;
    list_alignment( seq1, seq2, name1, name2, pos1, pos2, title);
    free ( seq1 );
    free ( seq2 );
    free ( name1 );
    free ( name2 );
    free ( title );
    return 0;
}

/* routine to receive an alignment. It consists of two null terminated
   strings containing the aligned sequences, their names and the positions
   of their left ends in the original sequences. */

int list_alignment ( char *seq1, char *seq2, char *name1, char *name2,
		     int pos1, int pos2, char *title )

{
    int i,j,k,seq_len,p1,p2,line_length=60;
    char match_syms[] = " :";
    int spads1, spads2;
    int l,p11,p22;
    p11=pos1;
    p22=pos2;
    

    seq_len = strlen(seq1);

    vmessage("%s\n", title);
   
    i = identities ( seq1, seq2 );
    vmessage(" Percentage mismatch %5.1f\n",
	     100*(seq_len ? (float)(seq_len-i)/seq_len : (float)1));

    for ( i=0,p1=pos1,p2=pos2;i<seq_len;i+=line_length) {
	vmessage("        ");

	for (j=0;j<6 && p1<pos1+seq_len;j++,p1+=10) {
            spads1=0;
	    for(l=0; l<10 && j*10+l+i < seq_len; l++){
	     	if (seq1[j*10+l+i]=='.')
		    spads1++;
	    }

	    if ( seq1[p1-pos1]=='.'){
		vmessage("%10c", '-');
		p11=p11-spads1+10;
	    } else{
		vmessage("%10d",p11);
		p11=p11-spads1+10;
	    }
	} 
	
	vmessage("\n%16.16s %.*s\n                 ",
		 name1,
		 i+line_length < seq_len ? 60 : seq_len - i,
		 &seq1[i]);

       	for (k=i;k<seq_len && k<i+line_length;k++) {
	    vmessage("%c",match_syms[same_char(seq1[k], seq2[k])]);
	}

	vmessage("\n%16.16s %.*s\n        ",
		 name2,
		 i+line_length < seq_len ? 60 : seq_len - i,
		 &seq2[i]);

	for (j=0;j<6 && p2<pos2+seq_len;j++,p2+=10) {
	    spads2=0;
            for(l=0; l<10 && j*10+l+i < seq_len; l++){
		if(seq2[j*10+l+i]=='.')
		    spads2++;
	    } 
	 
	    if (seq2[p2-pos2]=='.'){
		vmessage("%10c", '-');
		p22=p22-spads2+10;
	    } else{
		vmessage("%10d",p22);
		p22=p22-spads2+10;
	    }
	}
	vmessage("\n\n");
    }

    return 0;
}
/* end changes by kfs 27/1/95 */

/* routine to receive an alignment. It consists of two null terminated
   strings containing the aligned sequences, their names and the positions
   of their left ends in the original sequences. */
int iubc_list_alignment(char *seq1, char *seq2, char *name1, char *name2,
			int pos1, int pos2, char *title )

{
    int i,j,k,seq_len,p1,p2,line_length=60;
    char sym;

    seq_len = strlen(seq1);

    vmessage("%s\n", title);

    i = iubc_identities ( seq1, seq2 );
    vmessage(" Percentage mismatch %5.1f\n",
	     100*(seq_len ? (float)(seq_len-i)/seq_len : (float)1));

    for (i = 0, p1 = pos1, p2 = pos2; i < seq_len; i+=line_length) {
	vmessage("        ");

	for (j = 0; j < 6 && p1 < pos1+seq_len; j+=1, p1+=10) {
	    vmessage("%10d", p1);
	}
	vmessage("\n%16.16s %.*s\n                 ",
		 name1,
		 i+line_length < seq_len ? 60 : seq_len - i,
		 &seq1[i]);

	for (k = i; k < seq_len && k < i+line_length; k++) {
	    if (same_char(seq1[k], seq2[k])) {
		sym = ':';
	    } else if (IUBC_MATCH(seq2[k], seq1[k])) {
		sym = '.';
	    } else {
		sym = ' ';
	    }
	    vmessage("%c", sym);
	}
	vmessage("\n%16.16s %.*s\n        ",
		 name2,
		 i+line_length < seq_len ? 60 : seq_len - i,
		 &seq2[i]);

	for (j=0;j<6 && p2<pos2+seq_len;j+=1,p2+=10) {
	    vmessage("%10d",p2);
       }
	vmessage("\n\n");
    }
    return 0;
}

int rotate_seq ( char *seq, int seq_len, int origin ) {

    /* rotate seq seq so it starts at base origin.
       note numbering: base 1 is stored in seq[0]
    */

    char *buf;
    int i, j;

    if ( origin > seq_len+1 ) return -2;
    origin = (origin-1)%seq_len + 1;
    if ( origin < 1 ) return -3;

    if (origin == 1)
	return 0;

    if ( ( NULL == ( buf = ( char* ) xmalloc ( sizeof (char) * (origin-1))))) {
	return -1;
    }

    /* save up to origin to temp buffer */

    for ( i = 0; i < origin-1; i++ ) {
	buf[i] = seq[i];
    }

    /* move origin onwards to start of input array */

    for ( j = 0; i < seq_len; j++, i++ ) {
	seq[j]= seq[i];
    }

    /* put back original left end */

    for ( i = 0; i < origin-1; i++, j++ ) {
	seq[j] = buf[i];
    }
    xfree ( buf );
    return 0;
}

/*
 * Depad the sequence (length *len) in string str.  Array depad_to_pad is
 * filled with the padded location of each depadded base. The edits to str
 * are made in-situ. depad_to_pad array may be NULL.
 *
 * Returns: Modified len and str. Fills out depad_to_pad array.
 */
void depad_seq(char *str, int *len, int *depad_to_pad)
{
    int i;
    int curr_pos = 0;
    int old_len = *len;
    int x = old_len;
    char *a = str;
    char *b = str;

    /*str[old_len] = 0;*/

    for (i = 0; i < old_len; i++) {
	if (*b != '*') {
	    *a++ = *b++;
	    if (depad_to_pad)
		depad_to_pad[curr_pos++] = i;
	} else {
	    (*len)--;
	    b++;
	}
    }

    if (depad_to_pad) {
	for (i = curr_pos; i < old_len; i++) {
	    depad_to_pad[curr_pos++] = x++;
	}
    }

    if (*len < old_len) {
	*a = 0;
    }
}

/*
 * Given a combination of A, C, G or T, all of which are 0 for not present
 * and 1 for present, this returns an ambiguity code.
 */
char bases2ambiguity(int A, int C, int G, int T) {
    return "nTGKCYSBAWRDMHVN"[((A&1)<<3)+((C&1)<<2)+((G&1)<<1)+((T&1)<<0)];
}

/*
 * Given an ambiguity code, this stores in the A, C, G and T pointers either
 * 0 or 1 indicating if this code contains that element. Unknown codes
 * are treated as N.
 */
void ambiguity2bases(char ambig, int *A, int *C, int *G, int *T) {
    char *codes = "nTGKCYSBAWRDMHVN", *cp;
    int ind = (cp = strchr(codes, ambig)) ? cp - codes : 15;

    *A = (ind>>3) & 1;
    *C = (ind>>2) & 1;
    *G = (ind>>1) & 1;
    *T = (ind>>0) & 1;
}

/*
 * As base2ambiguity, but this time 'bits' encodes A, C, G or T (as bit 3, 2,
 * 1 and 0 respectively)
 */
char basebit2ambiguity(int bits) {
    return "nTGKCYSBAWRDMHVN"[bits];
}

/*
 * As ambiguity2bases, except we return a bit-pattern instead of 4 values.
 */
int ambiguity2basebit(char ambig) {
    char *codes = "nTGKCYSBAWRDMHVN", *cp;
    return  (cp = strchr(codes, ambig)) ? cp - codes : 15;
}

/*
 * Given nucleotides (possibly ambiguity codes themselves) we return
 * the IUB ambiguity codes.
 *
 * Logically speaking, this is equivalent to
 *    return basebit2ambiguity(ambiguity2basebit(b1) | ambiguity2basebit(b2));
 */
char ambiguity_code(char b1, char b2) {
    char *codes = "nTGKCYSBAWRDMHVN", *cp;
    int i1 = (cp = strchr(codes, b1)) ? cp - codes : 15;
    int i2 = (cp = strchr(codes, b2)) ? cp - codes : 15;
    return codes[i1 | i2];
}

#if 0 /* Not any faster than the above and short code */
/*
 * Given nucleotides (possibly ambiguity codes themselves) we return
 * the IUB ambiguity codes.
 *
 * NB: This table can be generated by (and hence also replaced by):
 *
 *  for (i = 0; i < 14; i++) {
 *	ambig1 = "ACGTRYMWSKDHVBN"[i];
 *	bit1 = ambiguity2basebit(ambig1);
 *	printf("/ * %c * / { ", ambig1);
 *	for (j = 0; j < 14; j++) {
 *	    ambig2 = "ACGTRYMWSKDHVBN"[j];
 *	    bit2 = ambiguity2basebit(ambig2);
 *	    printf("'%c'%c",basebit2ambiguity(bit2 | bit1), " ,"[j < 13]);
 *	}
 *	printf("},\n");
 *  }
 *
 */
char ambiguity_code(char base1, char base2) { 
    char *bases = "ACGTRYMWSKDHVB", *cp;
    char table[15][15] = {
	/*         A   C   G   T   R   Y   M   W   S   K   D   H   V   B   N */
	/* A */ { 'A','M','R','W','R','H','M','W','V','D','D','H','V','N','N'},
	/* C */ { 'M','C','S','Y','V','Y','M','H','S','B','N','H','V','B','N'},
	/* G */ { 'R','S','G','K','R','B','V','D','S','K','D','N','V','B','N'},
	/* T */ { 'W','Y','K','T','D','Y','H','W','B','K','D','H','N','B','N'},
	/* R */ { 'R','V','R','D','R','N','V','D','V','D','D','N','V','N','N'},
	/* Y */ { 'H','Y','B','Y','N','Y','H','H','B','B','N','H','N','B','N'},
	/* M */ { 'M','M','V','H','V','H','M','H','V','N','N','H','V','N','N'},
	/* W */ { 'W','H','D','W','D','H','H','W','N','D','D','H','N','N','N'},
	/* S */ { 'V','S','S','B','V','B','V','N','S','B','N','N','V','B','N'},
	/* K */ { 'D','B','K','K','D','B','N','D','B','K','D','N','N','B','N'},
	/* D */ { 'D','N','D','D','D','N','N','D','N','D','D','N','N','N','N'},
	/* H */ { 'H','H','N','H','N','H','H','H','N','N','N','H','N','N','N'},
	/* V */ { 'V','V','V','N','V','N','V','N','V','N','N','N','V','N','N'},
        /* B */ { 'N','B','B','B','N','B','N','N','B','B','N','N','N','B','N'},
	/* N */ { 'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'}
    };
	
    /* Turn base1 and base2 from A,C,G,T,...,N into 0,1,2,3,...,14 */
    cp = strchr(bases, toupper(base1));
    base1 = cp ? cp - bases : 14;
    cp = strchr(bases, toupper(base2));
    base2 = cp ? cp - bases : 14;

    return table[base1][base2];
}
#endif


