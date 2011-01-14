#include <staden_config.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "renz_utils.h"
#include "getfile.h"
#include "misc.h"                                        /* need for strdup */
#include "dna_utils.h"

#define NAMEDEL "/"
#define SEQDEL "/"
#define CUT '\''
#define PAD 'N'
#define MAXLINE 1024
#define MAXRSEQ 10
#define WHITESPACE " \t\n"
#define COMMENT ';'


/* structure defining ambiguity codes */

struct amb {
    int index;
    int range;
};

static struct amb ambiguity[4];

struct iubc_code {
    int ambiguity;
    char code[2];
    char alternatives[4];
};

struct iubc_code iubc_table[17] =
{         { 1, {'a'}, {'a'} },
	  { 1, {'c'}, {'c'} },
	  { 1, {'g'}, {'g'} },
	  { 1, {'t'}, {'t'} },
	  { 2, {'r'}, {'a','g'} },
	  { 2, {'y'}, {'c','t'} },
	  { 2, {'m'}, {'a','c'} },
	  { 2, {'k'}, {'g','t'} },
	  { 2, {'s'}, {'c','g'} },
	  { 2, {'w'}, {'a','t'} },
	  { 3, {'b'}, {'c','g','t'} },
	  { 3, {'d'}, {'a','g','t'} },
	  { 3, {'h'}, {'a','c','t'} },
	  { 3, {'v'}, {'a','c','g'} },
	  { 4, {'n'}, {'a','c','g','t'} },
	  { 4, {'-'}, {'a','c','g','t'} },
	  { 0, {'?'}, {' '} } };

int hash_word4 ( char *word );

/************************************************************/
int neighbors ( char in_string[], char out_string[256][5] ) {

    int i,i0,i1,i2,i3,out_string_counter;
/*
    for (i=0;i<16;i++) {
	printf("code %s ambibuity %d alternatives %s\n",
	       iubc_table[i].code, iubc_table[i].ambiguity, iubc_table[i].alternatives);
    }
*/
    for (i=0;i<4;i++) {
	ambiguity[i].index = iubc_lookup[(int)in_string[i]];
	ambiguity[i].range = iubc_table[iubc_lookup[(int)in_string[i]]].ambiguity;
    }

/*	now generate all the alternative strings for the ambiguous input */
    out_string_counter = 0;
    for (i0=0;i0<ambiguity[0].range;i0++) {
	for (i1=0;i1<ambiguity[1].range;i1++) {
	    for (i2=0;i2<ambiguity[2].range;i2++) {
		for (i3=0;i3<ambiguity[3].range;i3++) {
			out_string[out_string_counter][0] =
			    iubc_table[ambiguity[0].index].alternatives[i0];	
			out_string[out_string_counter][1] =
			    iubc_table[ambiguity[1].index].alternatives[i1];	
			out_string[out_string_counter][2] =
			    iubc_table[ambiguity[2].index].alternatives[i2];	
			out_string[out_string_counter][3] =
			    iubc_table[ambiguity[3].index].alternatives[i3];	
			out_string_counter++;
		    }
	    }
	}
    }
/*
    for (i=0;i<out_string_counter;i++) {
	printf(" %d %s %d\n",i,out_string[i],hash_word(out_string[i]));
    }
*/
    return out_string_counter;
}

/************************************************************/
int hashed_neighbors ( char in_string[], int word_len, int hash_values[] ) {

    int i,i0,i1,i2,i3,hash_counter;
    char word[4],t_string[4];

/*	Get the ambiguity index and range for each symbol in in_string */

    for (i=0;i<4;i++) t_string[i] = 'n';
    if (word_len > 4) word_len = 4;
    for (i=0;i<word_len;i++) t_string[i] = in_string[i];

    for (i=0;i<4;i++) {
	ambiguity[i].index = iubc_lookup[(int)t_string[i]];
	ambiguity[i].range = iubc_table[iubc_lookup[(int)t_string[i]]].ambiguity;
/*	printf(" %d %d %d\n",i,ambiguity[i].index,ambiguity[i].range);*/
    }

/*	now generate all the alternative strings for the ambiguous input */

    hash_counter = 0;
    for (i0=0;i0<ambiguity[0].range;i0++) {
	word[0] = iubc_table[ambiguity[0].index].alternatives[i0];
	for (i1=0;i1<ambiguity[1].range;i1++) {
	    word[1] = iubc_table[ambiguity[1].index].alternatives[i1];
	    for (i2=0;i2<ambiguity[2].range;i2++) {
		word[2] = iubc_table[ambiguity[2].index].alternatives[i2];
		for (i3=0;i3<ambiguity[3].range;i3++) {
		    word[3] = iubc_table[ambiguity[3].index].alternatives[i3];
		    hash_values[hash_counter] = hash_word4(word);
		    hash_counter++;
		}
	    }
	}
    }
/*    for (i=0;i<hash_counter;i++) {
	printf(" string hash values %d %d\n",i,hash_values[i]);
    }
*/
    return hash_counter;
}

/************************************************************/
int hash_word4 ( char *word ) {

/* given a string word, return its hash value. Assume it is at least 4 chars long */

#define hash_len 4

    int i;
    unsigned char uword;
    for (i=0,uword=0;i<hash_len;i++) {
	uword = ( uword <<2 ) | hash4_lookup[(unsigned int)word[i]];
    }
    return uword;
}

/************************************************************/

int hash_seq4 ( char *seq, int *hash_values, int seq_len) {

/* given a sequence seq, return an array of hash values */
/* note we assume hash_len = 4 */
/* We skip '*' (padding characters) when computing the hash such that the
 * Index into hash_values is correct, but the length of the match in seq may
 * not always be 4 characters.
 */

#define hash_len 4

    int i,j,k;
    unsigned char uword;

    if ( seq_len < hash_len ) return -1;

    for (i=0,uword=0;i<hash_len;i++) {
	uword = ( uword <<2 ) | hash4_lookup[(unsigned int)seq[i]];
    }
    i = 0;
    hash_values[0] = uword;
    k = seq_len - hash_len + 1;

    for (i=1,j=hash_len; i<k; i++,j++) {
	uword = ( uword <<2 ) | hash4_lookup[(unsigned int)seq[j]];
	hash_values[i] = uword;
    }
    return 0;
}

/*
 * A pad-aware version of hash_seq4. It works mainly, but has problems when
 * the sequence starts with pads. Do not use for now...
 */
int hash_seq4_padded ( char *seq, int *hash_values, int seq_len) {
    int i,j;
    unsigned char uword;

    /* Initial word */
    for (i=0,j=0,uword=0;i<hash_len && j < seq_len;j++) {
	if (seq[j] == '*')
	    continue;
	uword = ( uword <<2 ) | hash4_lookup[(unsigned int)seq[j]];
	i++;
    }
    if (j >= seq_len)
	return -1;

    /* Subsequent words */
    i = 0;
    hash_values[0] = uword;
    printf("hash_values[%d] = %x\n", i, uword);

    for (i = 1; i < seq_len && seq[i] == '*'; i++)
	;
    
    for (; j < seq_len; i++, j++) {
	while (seq[j] == '*' && j < seq_len) j++;
	while (seq[i] == '*') hash_values[i++] = 0;
	uword = ( uword <<2 ) | hash4_lookup[(unsigned int)seq[j]];
	hash_values[i] = uword;
	printf("hash_values[%d] = %x\n", i, uword);
    }

    return 0;
}

/************************************************************/

void 
store_hash4 (int *hash_values, /* the hash values for each position in a seq */
	     int seq_len,      /* size of the seq and hash array */
	     int *last_word,   /* last occurrence of this hash value */
	     int *word_count,  /* frequency of each hash value or word */
	     int start_pos,    /* first element of hash_values to process */
	     int size_hash )   /* number of elements in word_count and first_word */
{    
/* 	store the hash values in hash_values: put number of occurrences of
	each hash value in word_count; put the array position of the last 
	occurrence of each hash value in last_word, and previous
	occurrences in hash_values[last_word]. Process array elements
	start_pos onwards (assume start_pos is FORTRAN array element ie 1 not 0.
*/

#define hash_len 4

    int nw;
    register int i,j,n;

/* 	if start_pos is 1 we must zero the word counts */

    if ( 1 == start_pos ) {

	for ( i=0;i<size_hash;i++ ) word_count[i] = 0;
    }

/*	loop for all entries in hash_values	*/

    j = seq_len - hash_len + 1;

    for ( i=start_pos-1;i<j;i++ ) {

	n = hash_values[i];
	nw = word_count[n];

/* 	already an entry for this word ? */

	if ( 0 == nw ) {

/*	no, so put in last_word */

	    last_word[n] = i;
	    word_count[n] += 1;
	}

/*	yes, so put previous last occurrence in hash_values*/

	else {

	    word_count[n] += 1;
	    hash_values[i] = last_word[n];
	    last_word[n] = i;
	}
    }
}

/************************************************************/


int dna_string_search (
		   int *hash_values1, 	/* the hash values for each position in seq1 */
		   int *last_word, 	/* last occurrence of this hash value in seq1 */
		   int *word_count, 	/* frequency of each hash value or word */
		   int *hash_values2,	/* the hash values for seq2 */
		   int num_words,	/* the number of elements in hash_values2 */
		   int *seq1_match,	/* positions of matches in seq */
		   int max_matches,	/* maximum number of matches */
		   char *seq,
		   char *word,		
		   int seq_len, 	/* size of seq */
		   int word_len,	/* the length of the string */
		   char *end_seq,       /* the sequence from around the end of seq */
		   int end_offset,	/* the position of element seq_len in end_seq */
		   int circle	        /* if true we are dealing with a circle */
		   ) {

/*	we have a sequence seq encoded by hash tables:
	word_count[i] contains the number of occurrences of hash value i
	last_word[i] contains the seq position of the last occurrence of hash value i
	hash_values1[i] contains the hash values of previous occurrence of hash value i
	
	We are comparing word, of length word_len, which is composed of iubc ambiguity
	symbols, against seq. Comparison is done by generating all the possible strings
	of a,c,g,t that are represented by the first hash_len characters of word, and
	then generating the hash_values that these correspond to, and sending them here.
	In this routine we use the hashed representation of seq and the list of hash values
	for word, to find possible matches to word. For each possible hit we compare
	word with seq using a routine that understands the iubc codes. At this stage we
	use the correct length of word - until now we have used a fixed length hash_len
	which could be smaller or larger than word.

	algorithm: for each hash value in hash_values2[i] test if it exists in word_count
	if so process all occurrences using last-word and hash_values1.
	we follow matches found this way from start to end to check they really match
	( the hashing routine does not distinguish unknown characters in the sequence, 
	so we must do it here ).


	We must also sort out matches at the 3' end. We have only hashed up to and
	including the last hash_len length words in the sequence. Our word may be shorter
	than hash_len, and/or the sequence may be circular. Deal with this by copying
	the 3' end, and for circles some of the 5' end into a temporary array, and
	search it in a conventional manner. We have to assume that the longest word
	is max_word symbols in length.
	
*/


    int ncw, hash_value, pw;
    register int i,j,n_matches;
    int start_search, stop_search, end_seq_len, end_seq_pos;

    n_matches = 0;

/* 	loop for all num_words words in hash_values2 */

    for (i=0;i<num_words;i++) {

	hash_value = hash_values2[i];

/*	in seq ?	*/

	if ( 0 != (ncw = word_count[hash_value]) ) {

/*	yes, so process all (ncw) of them	*/

/*	first occurrence to process is in last_word 	*/

	    pw = last_word[hash_value];

	    for (j=0;j<ncw;j++) {

		if ( iubc_word_match_padded ( seq,pw,seq_len,
					      word,word_len ) ) {
		    
		    if ( n_matches < max_matches ) {

			seq1_match[n_matches] = pw + 1;
			n_matches += 1;

		    }
		    else {

/*	Error: too many matches */

/*			verror(ERR_WARN, "compare_seq", "too many matches");*/

			return -1;
		    }
		}
/*	get next occurrence of this word */

		pw = hash_values1[pw];

	    }

	}
    }

    end_seq_len = (end_offset+1)*2;

/* 	now deal with searches over and near the end of seq by scanning
	through end_seq in a simple way */

    start_search = ( word_len < hash_len ) ? seq_len - hash_len + 2 : seq_len - word_len + 2;

    stop_search = ( circle ) ? seq_len  : seq_len - word_len + 1;

/* 	for all positions (relative to original seq) from start_search to stop_search
	compare the word with the seq. How does this relate to end_seq?. Well end_offset
	is the element in seq that corresponds to seq_len in seq! */

    for (i=start_search,end_seq_pos=start_search-(seq_len-end_offset);
	 i<=stop_search;i++,end_seq_pos++) {

	if ( iubc_word_match ( end_seq,end_seq_pos,end_seq_len,word,word_len ) ) {

	    if ( n_matches < max_matches ) {

		seq1_match[n_matches] = i;
		n_matches += 1;

	    }
	    else {

/*	Error: too many matches */

/*			verror(ERR_WARN, "compare_seq", "too many matches");*/

		return -1;
	    }
	}
    }


    return n_matches;
}

/************************************************************/
void make_seq_end ( char *seq, int seq_len, char *end_seq, 
		   int max_end_seq, int *end_offset) {

/*	routine to make a short sequence to cover the end of a longer one
	for use during string searches. The search algorithm can miss matches
	either because the search string is shorter than the hash_len, or
	because the sequence is a circle. Here we copy a chunk from the 3' end
	and if a circle the start of the 5' end, to make a linear sequence for
	searching. We take max_end_seq/2 from each end. If the sequence is too
	short we return what we can.
	We also return the element number that corresponds to element seq_len */

    if ( max_end_seq/2 > seq_len ) {

	(void) copy_seq(end_seq,seq,seq_len);
	(void) copy_seq(&end_seq[seq_len],seq,seq_len);
	*end_offset = seq_len - 1;
	return;
    }
    (void) copy_seq(end_seq,&seq[seq_len-max_end_seq/2],max_end_seq/2);
    (void) copy_seq(&end_seq[max_end_seq/2],seq,max_end_seq/2);
    *end_offset = max_end_seq/2 - 1;
    return;
}
    

/* 
 * split "search_dna" into 2 functions: 
 * "hash_dna" which does the sequence hashing and only needs to be done once
 * for searches over the same sequence and
 * "dna_search" which does the actual searching 
 */
void hash_dna(char *seq,                                               /* in */
	      int seq_len,                                             /* in */
	      int *seq_hash_values,                                   /* out */
	      int *last_word,   
	      int *word_count)
{
    int start_pos;

    start_pos = 1;
    hash_seq4(seq, seq_hash_values, seq_len);
    
    store_hash4(seq_hash_values, /* the hash values for each position in a seq */
		 seq_len, 	/* size of the seq and hash array */
		 last_word, 	/* last occurrence of this hash value */
		 word_count, 	/* frequency of each hash value or word */
		 start_pos, 	/* first element of hash_values to process */
		 SIZE_HASH );	/* number of elements in word_count and first_word */
}

int dna_search(char *seq,                                             /* in */
	       int seq_len,                                           /* in */
	       char *word,                                            /* in */
	       int word_len,                                          /* in */
	       int circle,                                            /* in */
	       int *seq_hash_values,                                  /* in */
	       int *last_word,                                        /* in */
	       int *word_count,                                       /* in */
	       int *matches,                                         /* out */
	       int max_matches,                                        /*in */
	       int *num_matches)                                     /* out */
	       
{
#define MAX_END_SEQ 100
    char end_seq[MAX_END_SEQ];
    int end_offset;
    int string_hash_values[SIZE_HASH];
    int num_words;

    /* need to deal with 3' end of seq - make a copy of it for searching */
    (void) make_seq_end(seq, seq_len, end_seq, MAX_END_SEQ, &end_offset);

    if (num_words = hashed_neighbors ( word, word_len, string_hash_values )){

	*num_matches = dna_string_search (
	        seq_hash_values, /* the hash values for each position in seq */
		last_word, 	/* last occurrence of this hash value in seq */
		word_count, 	/* frequency of each hash value or word */
		string_hash_values,	/* the hash values for word */
		num_words,    /* the number of elements in string_hash_values */
		matches,	/* positions of matches in seq */
		max_matches,	/* maximum number of matches */
		seq,		/* the sequence being searched */
		word,	        /* the string being searched for */
		seq_len, 	/* size of seq */
		word_len,	/* the length of the search word */
		end_seq,	/* the bit of seq from the end */
		end_offset,	/* the position of element seq_len in end_seq */
		circle	/* if true seq is a circle */
		   );
/*
	for(i=0;i<*num_matches;i++) {
	    printf(" %d %d %.4s\n",i,matches[i],&seq[matches[i]-1]);
	}
*/
   }
/*   printf("number of matches %d\n",i); */
    return 0;
}

/************************************************************/
int search_dna (char *seq,                                             /* in */
		int seq_len,                                           /* in */
		char *word,                                            /* in */
		int word_len,                                          /* in */
		int circle,                                            /* in */
		int matches[],                                        /* out */
		int max_matches,                                      /* out */
		int *num_matches,                                     /* out */
		int *seq_hash_values)                                 /* out */
{

#define SIZE_HASH 256
#define MAX_END_SEQ 100


    int start_pos,num_words;
    int end_offset;
    char end_seq[MAX_END_SEQ];
    int last_word[SIZE_HASH];
    int word_count[SIZE_HASH];
    int string_hash_values[SIZE_HASH];
    start_pos = 1;

    hash_seq4 ( seq, seq_hash_values, seq_len);

    store_hash4 ( 
	         seq_hash_values, 	/* the hash values for each position in a seq */
		 seq_len, 	/* size of the seq and hash array */
		 last_word, 	/* last occurrence of this hash value */
		 word_count, 	/* frequency of each hash value or word */
		 start_pos, 	/* first element of hash_values to process */
		 SIZE_HASH );	/* number of elements in word_count and first_word */

/* need to deal with 3' end of seq - make a copy of it for searching */

    (void) make_seq_end ( seq, seq_len, end_seq, MAX_END_SEQ, &end_offset);
/*
    for(i=0;i<(end_offset+1)*2;i++){
	printf("%c",end_seq[i]);
    }
    printf("\n");
*/
    if (num_words = hashed_neighbors ( word, word_len, string_hash_values )){

	*num_matches = dna_string_search (
		   seq_hash_values, 	/* the hash values for each position in seq */
		   last_word, 	/* last occurrence of this hash value in seq */
		   word_count, 	/* frequency of each hash value or word */
		   string_hash_values,	/* the hash values for word */
		   num_words,	/* the number of elements in string_hash_values */
		   matches,	/* positions of matches in seq */
		   max_matches,	/* maximum number of matches */
		   seq,		/* the sequenc ebeing searched */
		   word,	/* the string being searched for */
		   seq_len, 	/* size of seq */
		   word_len,	/* the length of the search word */
		   end_seq,	/* the bit of seq from the end */
		   end_offset,	/* the position of element seq_len in end_seq */
		   circle	/* if true seq is a circle */
		   );
/*
	for(i=0;i<*num_matches;i++) {
	    printf(" %d %d %.4s\n",i,matches[i],&seq[matches[i]-1]);
	}
*/
   }
/*   printf("number of matches %d\n",i); */

    return 0;
}

/*
 * extract out only the significant characters in a recognition sequence
 * and find the cut position
 */
void
FindSequence(char *word,                                               /* in */
	     char *seq,                                               /* out */
	     int *cut_pos)                                            /* out */
{
    int i, j;
    int seq_cnt = 0;
    int seq_pos = 0;
    int found_cut = 0;
    int end;

    j = -1;
    /* remove leading pads */
    while ((int)word[++j] == PAD) ;

    end = strlen(word);
    for (i = j; i < end; i++){
	
	/* cut character found */
	if ((int)word[i] == CUT) {
	    *cut_pos = seq_pos;
	    found_cut = 1;
	}
	/* pad character found */
	else if ((int)word[i] == PAD) {
	    /* if found a cut which is before the sequence starts, decrement
	     * the cut position 
	     */
	    if ((seq_cnt == 0) && found_cut) {
		(*cut_pos)--;
	    } else {
		seq[seq_cnt++] = word[i];
	    }
	}
	/* normal char */
	else {
	    seq[seq_cnt++] = word[i];
	}
	seq_pos++;
    }
    /* terminate string */
    seq[seq_cnt] = '\0';

    /* strip off trailing pads */
    end = seq_cnt - 1;
    while ((int)seq[end] == PAD) {
	seq[end] = '\0';
	end--;
    }
}

/*
 * extract the enzyme name, cut position and recognition sequences from a 
 * string
 */
int ParseEnzyme(char *strptr,                                          /* in */
		R_Enz *res_enzyme)                                    /* out */
{
    char *tokptr;
    char *whiteptr;
    char *textptr;
    char word[MAXLINE];
    char name[MAXLINE]; 
    char tmp[MAXLINE]; 
    int num_seq;
    char res_seq[MAXRSEQ][MAXLINE];
    int cut_site[MAXRSEQ];
    int res_seq_len = MAXLINE;
    int i, j;
    int seq_num = 0;
    R_Enz r_enzyme;

    *tmp = 0;           /* MUST initialise tmp for the first strcat to work */

    /* remove white spaces */
    while ( (whiteptr = strtok(strptr, WHITESPACE)) != NULL) {
	strcat(tmp, whiteptr);  /* must be strcat to remove all whitespaces */
	strptr = NULL;    /* need to be NULL for subsequent calls to strtok */
    }

    /* extract name ie characters to the name delimiter */
    if ( ( textptr = strpbrk(tmp, NAMEDEL)) != NULL) {
	strncpy(name, tmp, textptr - tmp);
	name[textptr - tmp] = '\0';
    }

    /* extract sequences */
    while ( (tokptr = strtok (textptr, SEQDEL)) != NULL) {
	strcpy(word, tokptr);
	cut_site[seq_num] = 0;        /* if no cut site found, cut_site = 0 */
	FindSequence(word, res_seq[seq_num], &cut_site[seq_num]);
	textptr = NULL;   /* need to be NULL for subsequent calls to strtok */
	if (++seq_num >= MAXRSEQ) {
	    /* too many recognition sequences */
	    verror(ERR_WARN, "parse enzyme recognition sequences", "Too many recognition sequence");
	    break;
	}
    }

    num_seq = seq_num;

    /* malloc arrays of sequences and cut_positions for each enzyme */ 
    if (NULL == (r_enzyme.name = (char *)xmalloc((strlen(name) + 1) * 
						sizeof(char ))))
	return 0;

    if (NULL == (r_enzyme.seq = (char **)xmalloc((num_seq + 1) * 
						sizeof(char *))))
	return 0;

    if (NULL == (r_enzyme.cut_site = (int *)xmalloc((num_seq + 1) * 
						   sizeof(int))))
	return 0;

    /* malloc character length for sequence array */
    for (j = 0; j < num_seq; j++){

	res_seq_len = strlen(res_seq[j]);
	if (NULL == (r_enzyme.seq[j] = (char *)xmalloc((res_seq_len + 1) * 
						      sizeof(char ))))
	    return 0;
    }

    /* copy temporary variables into r_enzyme structure */
    r_enzyme.num_seq = num_seq;
    strcpy(r_enzyme.name, name);

    for (i = 0; i < num_seq; i++) {
	strcpy(r_enzyme.seq[i], res_seq[i]);
	r_enzyme.cut_site[i] = cut_site[i];

/*
	for (j = 0; j < r_enzyme.num_matches[i]; j++) {
	    printf("%s %s %d\n", r_enzyme.name, r_enzyme.seq[i], r_enzyme.cut_site[i]+r_enzyme.cut_pos[i][j]);

	}
*/
    }
    *res_enzyme = r_enzyme;
    return 1;
}

/*
 * extract only the enzyme name from a string 
 */

int 
GetEnzymeName(char *strptr,                                            /* in */
	      char **names)                                           /* out */
{
    char *whiteptr;
    char *textptr;
    char name[MAXLINE];
    char tmp[MAXLINE]; 
    tmp[0] = '\0';

    /* remove white spaces */
    while ( (whiteptr = strtok(strptr, WHITESPACE)) != NULL) {
	strcat(tmp, whiteptr); 
	strptr = NULL;
    }

    /* extract name ie characters to the name delimiter */
    strptr = tmp;
    if ( ( textptr = strpbrk(strptr, NAMEDEL)) != NULL) {
	strncpy(name, strptr, textptr - strptr);
	name[textptr - strptr] = 0;
    } else {
	return 0;
    }
    if (NULL == ((*names) = (char *)xmalloc((strlen(name) + 1) * 
					sizeof(char ))))
	return 0;
    strcpy((*names), name);
    return 1;

}

/*
 * read a file of enzymes and parse those enzymes that have been selected ie
 * are in 'inlist'
 * store only the enzyme names in array enz_names
 */
int r_enz_file_names(char *file_name,
		     char ***enz_names,
		     int *num_enzymes)
{
    FILE *fp;
    char line[MAXLINE];
    int line_num = 0;
    char exp_filename[FILENAME_MAX];
    char **names;

    expandpath(file_name, exp_filename);

   if (NULL == (fp = fopen(exp_filename, "r"))){
	return 0;
    }
    *num_enzymes = 0;

    /* read in each line of the enzyme file */
    while (fgets(line, MAXLINE, fp) != NULL) {
	if (((int)line[0] != COMMENT) && (strcmp(line, "\n") != 0)) {
	    line_num++;
	}
    }
    rewind(fp);

    if (NULL == (names = (char **)xmalloc(line_num * sizeof(char *))))
	return 0;

    /* for each line in the enzyme file */
    while (fgets(line, MAXLINE, fp) != NULL) {
	/* if the line is not a comment or blank */
	if (((int)line[0] != COMMENT) && (strcmp(line, "\n") != 0)) {
	    /* extract the enzyme name from the line */
	    if (0 == GetEnzymeName(line, &names[(*num_enzymes)++]))
		(*num_enzymes)--;
	}
    }
    fclose(fp);
    *enz_names = names;
    return 1;
}

/*
 * read a file of enzymes and parse those enzymes that have been selected ie
 * are in 'inlist'
 * store the information in array r_enzyme
 */
int open_renz_file (char *file_name,                                   /* in */
		    char *inlist,                                      /* in */
		    int num_items,                                     /* in */
		    R_Enz **res_enzyme,                               /* out */
		    int *num_enzymes)                                 /* out */
{
    FILE *fp;
    char line[MAXLINE];
    R_Enz *r_enzyme;
    int line_num = 0;
    char exp_filename[FILENAME_MAX];
    int line_counter = 0;
    int item;
    char *list_ptr;
    char *list;

    /* expand the name of enzyme environment variable */
    expandpath(file_name, exp_filename);

   if (NULL == (fp = fopen(exp_filename, "r"))){
	return 0;
    }

    /* read in each line of the enzyme file */
/*
    while (fgets(line, MAXLINE, fp) != NULL) {
	num_lines++; 
    }
    rewind(fp);
*/
    /* set inlist to be the active list */
    *num_enzymes = num_items;
    if (NULL == (r_enzyme = (R_Enz *)xmalloc(((*num_enzymes) + 1) * 
						sizeof(R_Enz))))
	return 0;

    /* get the first item on the inlist */
    item = strtol(inlist, &list_ptr, 10);

    /* for each line in the enzyme file */
    while (fgets(line, MAXLINE, fp) != NULL) {

	/* if the line is not a comment */
	if ((int)line[0] != COMMENT) {
	    /* if the line corresponds to a selected enzyme */
	    if (item == line_counter) {
		if (line_num == *num_enzymes)
		    break;
		ParseEnzyme(line, &r_enzyme[line_num]);
		/* r_enzyme[line_num].number = item; */
		line_num++;
		list = list_ptr;
		item = strtol(list, &list_ptr, 10);
	    }
	    line_counter++;
	}
    }


    fclose(fp);
    *res_enzyme = r_enzyme;
    return 1;
}


static int compare_rmatch(const void *p1, const void *p2)
{
    R_Match *r1 = (R_Match *)p1, *r2 = (R_Match *)p2;
    if ((*r1).cut_pos < (*r2).cut_pos) 
	return (-1);
    else if ((*r1).cut_pos == (*r2).cut_pos)
	return (0);
    else
	return (1);
}


static int compareint(const void *p1, const void *p2)
{
    int *i1 = (int *)p1, *i2 = (int *)p2;
    if (*i1 < *i2) 
	return (-1);
    else if (*i1 == *i2)
	return (0);
    else
	return (1);
}

/*
 * find the position of the left hand base of where the enzyme target sequence
 * cuts the dna
 */
int
FindMatches(R_Enz *r_enzyme,                                          /* in */
	    int num_enzymes,                                          /* in */
	    char *sequence,                                           /* in */
	    int sequence_len,                                         /* in */
	    int sequence_type,                                        /* in */
	    R_Match **match,                                         /* out */
	    int *total_matches)                                      /* out */
{
#define SIZE_HASH 256

    int *seq_hash_values;
    int *matches, max_matches=MAXMATCHES;
    int num_matches;
    int i, j, k;
    int cnt = 0;
    int array_size = MAXMATCHES;
    int array_inc = MAXMATCHES;
    int last_word[SIZE_HASH];
    int word_count[SIZE_HASH];

    if ( ! (seq_hash_values = (int *) xmalloc ( sizeof(int)*(sequence_len)))) {
	return -2;
    }

    if ( ! (matches = (int *) xmalloc ( sizeof(int)*(max_matches) ))) {
	return -2;
    }

    /* hash sequence */
    hash_dna(sequence, sequence_len, seq_hash_values, last_word, word_count);

    /* for each enzyme */
    for (i = 0; i < num_enzymes; i++) {
	/* for each recognition sequence */
	for (j = 0; j < r_enzyme[i].num_seq; j++) {
	    /* find the matches between rec seq and dna */
	    dna_search(sequence, sequence_len, r_enzyme[i].seq[j],
		       strlen(r_enzyme[i].seq[j]), sequence_type,
		       seq_hash_values, last_word, word_count, matches, 
		       max_matches, &num_matches);

	    /* error = search_dna (sequence, sequence_len, r_enzyme[i].seq[j], 
				strlen(r_enzyme[i].seq[j]), sequence_type, 
				matches, max_matches, 
				&num_matches, seq_hash_values);
				*/
	    /* store the match results in array 'match' */
	    for (k = 0; k < num_matches; k++) {
		(*match)[cnt].enz_name = i;
		(*match)[cnt].enz_seq = j;
		/* store the cut position in matches as required for 
		 * displaying 
		 */
		if (matches[k] + r_enzyme[i].cut_site[j] == sequence_len) {
		    (*match)[cnt].padded_cut_pos =
			(*match)[cnt].cut_pos = sequence_len;
		} else {
		    /*
		    (*match)[cnt].cut_pos = (matches[k] + 
				      r_enzyme[i].cut_site[j])%sequence_len;
				      */
		    (*match)[cnt].padded_cut_pos =
			(*match)[cnt].cut_pos =
			    matches[k] + r_enzyme[i].cut_site[j];
		}
		   
		cnt++;
		/* cnt is less than MAXMATCHES 
	        * if it is larger, then need to allocate more space to array
		*/
		if (cnt >= array_size) {
	            /* set array_size to be count plus arbitary increment */
	            array_size = cnt + array_inc;

	            if (NULL == ( (*match) = (R_Match *)realloc((*match), 
							     array_size * 
							     sizeof(R_Match))))
			return 0;
		    /* clear the new memory */

		    memset(&(*match)[cnt], 0, sizeof(R_Match) * array_inc); 
		}
	    }
	}
    }

    *total_matches = cnt;
    xfree(seq_hash_values);
    xfree(matches);
    return 1;
}

/*
 * determine the colour for each enzyme
 */
char *
SetREnzColour(int num_enzymes,
	      int enz_num)
{
#define NUM_COL 7
    int r = 0;
    int g = 0;
    int b = 0;
    int max_col = 255;
    int type, cnt;
    int inc;
    static char colour[NUM_COL];

    inc = (int) (max_col /  ((num_enzymes / NUM_COL) + 1));
    type = (enz_num % NUM_COL) + 1;
    cnt = (enz_num / NUM_COL) + 1;
    
    /* printf("inc %d type %d cnt %d\n", inc, type, cnt); */

    if ((type == 1) || (type == 4) || (type == 6) || (type == 7))
	r = inc * cnt;
    if ((type == 2) || (type == 4) || (type == 5) || (type == 7))
	g = inc * cnt;
    if ((type == 3) || (type == 5) || (type == 6) || (type == 7))
	b = inc * cnt;

    sprintf(colour, "#%02x%02x%02x", r, g, b); 
    return colour;

}

/*
 * return the string in the target sequence corresponding to the enzyme
 * sequence - fill in padding chars and add the cut site
 */
void
ExpandRSeq(int match,                                                  /* in */
	   int cut_pos,                                                /* in */
	   char *sequence,                                             /* in */
	   int sequence_len,                                           /* in */
	   int sequence_type,                                          /* in */
	   char *in_seq,                                               /* in */
	   char *out_seq)                                             /* out */
{
    int i, j;
    int cnt = 0;
    int start = 0;
    int length;
    int inseq_len = 0;

    /* Match is the position of the cut site. We want to go back to the first
     * base of the recognition sequence.
     */
    match--; /* starts at offset 1, we want 0 */
    if (cut_pos > 0) {
	int count = cut_pos;
	while (count > 0) {
	    while (--match && match >= 0 && sequence[match] == '*')
		;
	    count--;
	}
    } else {
	match -= cut_pos;
    }

    inseq_len = strlen(in_seq);

    /* only happens if at the end of a circular dna */
    if ((match < 0) && (sequence_type == 1)) { 
	match = sequence_len + match;
    }

    /* padding characters before target string */
    if (cut_pos < 0) {
	length = inseq_len;
	start = cut_pos;
    }
    /* padding characters after target string */
    else if (cut_pos >= inseq_len) {
	length = cut_pos + 1;
    } else {
	length = inseq_len;
    }

    for (j = 0, i = start; i < length; i++) {
	if (i == cut_pos) {
	    out_seq[cnt++] = CUT;
	    if (cut_pos >= inseq_len)
		break;
	}

	if (sequence_type == 0) { /* linear dna */
	    /* at end */
	    while (1) {
		if (match+i+j >= sequence_len || match+i < 0) {
		    out_seq[cnt++] = PAD;
		} else {
		    if (sequence[match+i+j] == '*') {
			j++;
			continue;
		    }
		    out_seq[cnt++] = sequence[match+i+j];
		} 
		break;
	    }

	} else { /* circular dna */
	    while (sequence[(match+i+j+sequence_len)%sequence_len] == '*')
		j++;
	    out_seq[cnt++] = sequence[(match+i+j+sequence_len)%sequence_len];
	}
    }

    /* terminate out_seq */
    out_seq[cnt] = '\0';
}

void
FindFragments(int num_matches,                                         /* in */
	      R_Match *tmp_match,                                      /* in */
	      int sequence_len,                                        /* in */
	      int sequence_type,                                       /* in */
	      int *fragment)                                          /* out */
{
    int i;

    /* if sequence is a circle */
    if (sequence_type == 1) {
	fragment[0] = sequence_len - tmp_match[num_matches-1].cut_pos + 
	    tmp_match[0].cut_pos;
	
	for (i = 1; i < num_matches; i++) {
	    
	    fragment[i] = tmp_match[i].cut_pos - tmp_match[i-1].cut_pos;
	}
    } else {
	fragment[0] = tmp_match[0].cut_pos - 1;
	for (i = 1; i < num_matches; i++) {
	    fragment[i] = tmp_match[i].cut_pos - tmp_match[i-1].cut_pos;
	}
	fragment[num_matches] = sequence_len - tmp_match[num_matches-1].cut_pos + 1;
    }
    

}


int
PrintEnzymeByEnzyme(R_Enz *r_enzyme,                                   /* in */
		    R_Match *match,                                    /* in */
		    int total_matches,                                 /* in */
		    int num_enzymes,                                   /* in */
		    char *sequence,                                    /* in */
		    int sequence_len,                                  /* in */
		    int sequence_type,                                 /* in */
		    int lreg,                                          /* in */
		    int do_all)                                        /* in */
{
    int i,k;
    int cnt =0;
    char r_seq[MAXLINE];
    R_Match *tmp_match;
    int start = 0;
    int num_matches;
    int j = 0;
    int *fragment;
    int *lengths;
    char fbuf[1024], lbuf[1024];
    int fragments_printed = 0;

    if (0 == num_enzymes)
      return 1;

    /* printf("total matches %d\n", total_matches); */
    if (total_matches == 0)
	return 0;
    if (NULL == (tmp_match = (R_Match *)xmalloc(total_matches * sizeof(R_Match))))
	return 0;
    
    /* printf("num enz %d\n", num_enzymes); */

    for (i = 0; i < num_enzymes; i++) {
	while ((j < total_matches) && (i == match[j].enz_name)) {
	    memcpy(&tmp_match[cnt], &match[j], sizeof(R_Match));
	    cnt++;
	    j++;
	}
	if (0 == (num_matches = j - start))
	    continue;

	/* 
	 * malloc total_matches +1 because for linear sequences, have an extra
	 * fragment ie if cut n times, get n+1 fragments
	 */
	if (NULL == (fragment = (int *)xmalloc((num_matches +1)* sizeof(int))))
	    return 0;
	if (NULL == (lengths = (int *)xmalloc((num_matches +1) * sizeof(int))))
	    return 0;
	
	qsort((void *) tmp_match, num_matches, sizeof(R_Match), 
	      compare_rmatch);
	
	vmessage("\n  Matches found= %5d \n", num_matches);
	vmessage("%10s%20s%34s%9s%8s\n", "Name", "Sequence", "Position", 
	       "Fragment", "lengths");

	FindFragments(num_matches, tmp_match, sequence_len, sequence_type,
		      fragment);

	/* linear */
	if (sequence_type == 0) { 
	    memcpy(lengths, fragment, (num_matches + 1) *sizeof(int));
	    qsort((void *) lengths, (num_matches +1), sizeof(int), compareint);
	}
	/* circular */
	else {
	    memcpy(lengths, fragment, (num_matches) *sizeof(int));
	    qsort((void *) lengths, (num_matches), sizeof(int), compareint);
	}

	for (k =0; k < num_matches; k++) {
	    ExpandRSeq(tmp_match[k].cut_pos,
		       r_enzyme[tmp_match[k].enz_name].cut_site[tmp_match[k].enz_seq],
		       sequence, 
		       sequence_len,
		       sequence_type,
		       r_enzyme[tmp_match[k].enz_name].seq[tmp_match[k].enz_seq],
		       r_seq);
	    
	    if (fragment[k] > 0 && fragment[k] <= sequence_len) {
		sprintf(fbuf, "%7d", fragment[k]);
		fragments_printed++;
	    } else
		sprintf(fbuf, "%7s", "-");
	    if (lengths[k] > 0)
		sprintf(lbuf, "%7d", lengths[k]);
	    else
		sprintf(lbuf, "%7s", "-");
	    
	    vmessage("%5d %-15s %-32s%10d%s%s \n",
		     k+1,
		     r_enzyme[tmp_match[k].enz_name].name,
		     r_seq,
		     tmp_match[k].cut_pos + lreg - 1,
		     fbuf, lbuf);
	}

	/* print last fragment if linear sequence */
	if (sequence_type == 0) {
	    if (fragment[num_matches] > 0)
		vmessage("%71d%7d \n",
			 fragment[num_matches],
			 lengths[num_matches]);
	    else if (fragments_printed < 2)
		vmessage("%71d%7d \n",
			 lengths[num_matches],
			 lengths[num_matches]);
	    else
		vmessage("%71s%7d \n",
			 "-",
			 lengths[num_matches]);
	}
	cnt = 0;
	num_matches = 0;
	start = j;

	xfree(fragment);
	xfree(lengths);
    }

    /* 
     * only do zero cutters if displaying information on all restriction
     * enzymes.
     */
    if (do_all) {
	vmessage("Zero cutters:\n");
	cnt = j = start = 0;
	for (i = 0; i < num_enzymes; i++) {
	    while ((j < total_matches) && (i == match[j].enz_name)) {
		cnt++;
		j++;
	    }
	    if (0 == (num_matches = j - start)) {
		vmessage("      %s\n", r_enzyme[i].name);
	    }
	    cnt = 0;
	    num_matches = 0;
	    start = j;
	}
    }

    xfree(tmp_match);
    return 1;
}


int
OrderOnPosition (R_Enz *r_enzyme,
		 R_Match *match,
		 int total_matches,
		 char *sequence,
		 int sequence_len,
		 int sequence_type,
		 int lreg)
{
    int k;
    R_Match *tmp_match;
    char r_seq[MAXLINE];
    int *fragment;
    int *lengths;
    char fbuf[1024], lbuf[1024];
    int fragments_printed = 0;

    /* make temporary copy of match data */
    if (total_matches == 0)
	return 0;
    if (NULL == (tmp_match = (R_Match *)xmalloc(total_matches * sizeof(R_Match))))
	return 0;
    memcpy(tmp_match, match, total_matches * sizeof(R_Match));

    /* sort on position over all the enzymes */
    qsort((void *) tmp_match, total_matches, sizeof(R_Match), compare_rmatch);
 
    vmessage("%10s%20s%34s%9s%8s\n", "Name", "Sequence", "Position", 
	   "Fragment", "lengths");

    /* 
     * malloc total_matches +1 because for linear sequences, have an extra
     * fragment ie if cut n times, get n+1 fragments
     */
    if (NULL == (fragment = (int *)xmalloc((total_matches+1) * sizeof(int))))
	return 0;
    if (NULL == (lengths = (int *)xmalloc((total_matches+1) * sizeof(int))))
	return 0;

    FindFragments(total_matches, tmp_match, sequence_len, sequence_type, 
		  fragment);

    /* linear */
    if (sequence_type == 0) { 
	memcpy(lengths, fragment, (total_matches + 1) *sizeof(int));
	qsort((void *) lengths, (total_matches +1), sizeof(int), compareint);
    }
    /* circular */
    else {
	memcpy(lengths, fragment, (total_matches) *sizeof(int));
	qsort((void *) lengths, (total_matches), sizeof(int), compareint);
    }
    
/*
    memcpy(lengths, fragment, total_matches*sizeof(int));
    qsort((void *) lengths, total_matches, sizeof(int), compareint);
*/
    for (k = 0; k < total_matches; k++) {

	ExpandRSeq(tmp_match[k].cut_pos,
		   r_enzyme[tmp_match[k].enz_name].cut_site[tmp_match[k].enz_seq],
		   sequence, 
		   sequence_len,
		   sequence_type,
		   r_enzyme[tmp_match[k].enz_name].seq[tmp_match[k].enz_seq],
		   r_seq);
	    
	if (fragment[k] >= 0 && fragment[k] <= sequence_len) {
	    sprintf(fbuf, "%7d", fragment[k]);
	    fragments_printed++;
	} else
	    sprintf(fbuf, "%7s", "-");
	if (lengths[k] >= 0)
	    sprintf(lbuf, "%7d", lengths[k]);
	else
	    sprintf(lbuf, "%7s", "-");
	    
	vmessage("%5d %-15s %-32s%10d%s%s \n",
		 k+1,
		 r_enzyme[tmp_match[k].enz_name].name,
		 r_seq,
		 tmp_match[k].cut_pos + lreg - 1,
		 fbuf, lbuf);
    }
    /* print last fragment if linear sequence */
    if (sequence_type == 0) {
	if (fragment[total_matches] >= 0)
	    vmessage("%71d%7d \n",
		     fragment[total_matches],
		     lengths[total_matches]);
	else if (fragments_printed < 2)
	    vmessage("%71d%7d \n",
		     lengths[total_matches],
		     lengths[total_matches]);
	else
	    vmessage("%71s%7d \n",
		     "-",
		     lengths[total_matches]);
    }
    xfree(tmp_match);
    xfree(fragment);
    xfree(lengths);

    return 1;
}

char *
AddCutSites(char *seq,
	    int cut_site)
{
    static char newseq[MAXLINE];
    int i;
    int cnt = 0;

    if (cut_site < 0) {
	newseq[0] = CUT;
	for (i = cut_site; i < 0; i++)
	    newseq[++cnt] = PAD;
	newseq[++cnt] = '\0';
	strcat(newseq, seq);
    } else if (cut_site > (int)strlen(seq)) {
	strcpy(newseq,seq);
	for (i = strlen(seq); i < cut_site; i++)
	     newseq[i] = PAD;
	newseq[cut_site] = CUT;
	newseq[cut_site+1] = '\0';
    } else {
	strncpy(newseq, seq, cut_site);
	newseq[cut_site] = CUT;
	newseq[cut_site+1] = '\0';
	strncat(newseq, &seq[cut_site], strlen(seq) - cut_site);
    }
    return newseq;

}

int find_max_cut_dist(R_Enz * r_enzyme,
		      int num_enzymes)
{
    int i, j;
    int str_len, cut;
    int total_len;
    int max_overlap = 0;

    for (i = 0; i < num_enzymes; i++) {
	for (j = 0; j < r_enzyme[i].num_seq; j++) {
	    str_len = strlen(r_enzyme[i].seq[j]);
	    cut = r_enzyme[i].cut_site[j];
	    if (cut < 0) {
		total_len = abs(cut) + str_len;
	    } else {
		total_len = MAX(cut, str_len);
	    } 
	    if (max_overlap < total_len)
		max_overlap = total_len;
	    
	}
    }

    return max_overlap;
}

