#ifndef _CODON_CONTENT_H
#define _CODON_CONTENT_H

/* gene search by content */

typedef struct _CodRes {
    int user_start;		/* sequence start position in user units */
    int user_end;		/* sequence end position in user units */
    int num_results;		/* number of results per frame */
    int window_length;		/* length of the window used */
    double error;               /* percentage error */
    double max;	   		/* max observed value */
    double min;			/* min observed value */
    double *frame1;		/* the results for frame 0 */
    double *frame2;		/* the results for frame 1 */
    double *frame3;		/* the results for frame 2 */
    char   *top;
} CodRes;

CodRes *init_CodRes ( int num_results );

void free_CodRes ( CodRes *r );

typedef struct _CodRes1 {
    int user_start;		/* sequence start position in user units */
    int user_end;		/* sequence end position in user units */
    int num_results;		/* number of results per frame */
    int window_length;		/* length of the window used */
    double max;	   		/* max observed value */
    double min;			/* min observed value */
    double *frame1;		/* the results for frame 0 */
} CodRes1;

CodRes1 *init_CodRes1 ( int num_results );

void free_CodRes1 ( CodRes1 *r );

int init_codon_pref (char *file_name, double codon_usage_table[4][4][4],
		     int option);

int init_author_test (char *file_name, char *seq, int seq_length, 
		      double codon_usage_wm[4][4][4], double percent_error,
		      int *window_length );

/* no longer used */
#if 0
int do_pos_base_pref ( char seq[], int seq_length,
		       double codon_usage_table_in[4][4][4], int use_table,
		       CodRes *results );
#endif

int do_codon_pref ( char seq[], int seq_length,
		      double codon_usage_table[4][4][4], CodRes *results );

int do_author_test ( char seq[], int seq_length,
		     double codon_usage_table[4][4][4], CodRes *results );

int do_pos_base_bias ( char seq[], int seq_length, CodRes1 *results );
void init_codon_table(double codon_table[4][4][4]);
void calc_codon_usage(char *seq, int seq_length, double codon_table[4][4][4]);


#endif
