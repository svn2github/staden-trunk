#ifndef _DNA_UTILS_H
#define _DNA_UTILS_H


#ifdef _MSC_VER
#  ifdef BUILDING_SEQ_UTILS_DLL
#    define SEQ_UTILS_EXPORT  extern __declspec(dllexport)
#  else
#    define SEQ_UTILS_EXPORT  extern __declspec(dllimport)
#  endif
#else
#  define SEQ_UTILS_EXPORT extern
#endif


int same_char(char c1, char c2);

int match_len ( char *seq1, int seq1_start, int seq1_len, char *seq2,
	              int seq2_start, int seq2_len );

 /*
 * Find all positions in seq at which string has >= min_match matching
 * characters.
 */
int inexact_match (char *seq, int seq_len, char *string, int string_len,
		   int min_match, int *match, int *score, int max_matches);

int iubc_inexact_match (char *seq, int seq_len, char *string, int string_len,
			int min_match, int use_iub_code, int *match, 
			int *score, int max_matches);

/*
 * Find the best inexact match of 'string' in seq. See inexact_match for
 * more details.
 *
 * Returns the best match value and stores its position in 'match'
 * (if match != NULL).
 */
int best_inexact_match(char *seq, int seq_len, char *string, int string_len,
		       int *match);

/*
 * copy a sequence
 */
void copy_seq ( char *copy, char *original, int seq_len );


/*
 * routine to receive an alignment. It consists of two null terminated
 * strings containing the aligned sequences, their names and the positions
 * of their left ends in the original sequences.
*/
int list_alignment ( char *seq1, char *seq2, char *name1, char *name2,
		     int pos1, int pos2, char *title );

int iubc_list_alignment(char *seq1, char *seq2, char *name1, char *name2,
			int pos1, int pos2, char *title );
/*
 * Complement and reverse a dna sequence
 */
void complement_seq(char *seq, int seq_len);
char complement_base(char base);

/*
 * Complement but not reverse a dna sequence
 */
void complement_dna(char *seq, int seq_len);

/*
 * Reverse a dna sequence
 */
void reverse_dna( char *seq, int seq_len );

/*
 * Does word match seq at seq_pos
 */
int word_match( char *seq, int seq_pos, int seq_len, char *word, int word_len);



SEQ_UTILS_EXPORT int char_set_size, *char_lookup, *char_match;
SEQ_UTILS_EXPORT int dna_lookup[256];
SEQ_UTILS_EXPORT int dna_match[256];
SEQ_UTILS_EXPORT int unknown_char;
SEQ_UTILS_EXPORT int hash4_lookup[256];


/* Useful macros used by alignment routines, they require char_match to be set! */
#define SEQ_MATCH(a,b)    ( ( (char_match[(int)a] < unknown_char) && (char_match[(int)a]) == (char_match[(int)b])) ? 1 : 0 )
#define SEQ_MISMATCH(a,b) ( ( (char_match[(int)a] < unknown_char) && (char_match[(int)a]) == (char_match[(int)b])) ? 0 : 1 )


/*
 * set up table of values for permitted dna characters
 */
void set_dna_lookup(void);

/*
 * set up table of index values for iubc characters
 */
void set_iubc_lookup(void);


/*
 * Does word in iubc codes match seq at seq_pos
 */
int iubc_word_match( char *seq, int seq_pos, int seq_len, char *word,
		    int word_len);

/*
 * Does word in iubc codes match seq at seq_pos, allowing for seq to be
 * padded (but not word).
 */
int iubc_word_match_padded( char *seq, int seq_pos, int seq_len, char *word,
			    int word_len);

SEQ_UTILS_EXPORT int protein_lookup[256];

/*
 * Initialise protein_lookup and protein_match arrays
 */
void set_protein_lookup(void);

/*
 *
 * seq_type == 1 for DNA
 * seq_type != 1 for protein
 *
 */
void set_char_set(int seq_type);

int rotate_seq ( char *seq, int seq_len, int origin );

/*
 * Depad the sequence (length *len) in string str.  Array depad_to_pad is
 * filled with the padded location of each depadded base. The edits to str
 * are made in-situ.
 *
 * Returns: Modified len and str. Fills out depad_to_pad array.
 */
void depad_seq(char *str, int *len, int *depad_to_pad);

/*
 * Given a combination of A, C, G or T, all of which are 0 for not present
 * and 1 for present, this returns an ambiguity code.
 */
char bases2ambiguity(int A, int C, int G, int T);

/*
 * As base2ambiguity, but this time 'bits' encodes A, C, G or T (as bit 3, 2,
 * 1 and 0 respectively)
 */
char basebit2ambiguity(int bits);

/*
 * Given an ambiguity code, this stores in the A, C, G and T pointers either
 * 0 or 1 indicating if this code contains that element. Unknown codes
 * are treated as N.
 */
void ambiguity2bases(char ambig, int *A, int *C, int *G, int *T);

/*
 * As ambiguity2bases, except we return a bit-pattern instead of 4 values.
 */
int ambiguity2basebit(char ambig);

/*
 * Given nucleotides (possibly ambiguity codes themselves) we return
 * the IUB ambiguity codes.
 *
 * Logically speaking, this is equivalent to
 *    return basebit2ambiguity(ambiguity2basebit(b1) | ambiguity2basebit(b2));
 */
char ambiguity_code(char b1, char b2);

#endif
