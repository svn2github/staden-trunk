#ifndef _R_ENZ_UTILS_H_
#define _R_ENZ_UTILS_H_

#define SCRIPT_LEN 1024
#define MAXMATCHES 10000
#define SIZE_HASH 256

extern int iubc_lookup[];

/* for each enzyme */
typedef struct r_enz {

    char *name;                                        /* name of the enzyme */
    int num_seq;                          /* number of recognition sequences */
    char **seq;                            /* array of recognition sequences */
    int *cut_site;                                   /* position of cut site */

} R_Enz;

/* for each cut */
typedef struct r_match {

    unsigned short enz_name;         /* position of enzyme in R_Enz array */
    unsigned char enz_seq;           /* refers to particular recognition seq */
    int cut_pos;                     /* cut position in bases */
    int padded_cut_pos;		     /* cut position in padded seq */

} R_Match;

int 
read_r_enz_file(char *sequence,                                        /* in */
		int sequence_len,                                      /* in */
		R_Enz **res_enzyme,                                   /* out */
		int *num_enzymes);                                    /* out */

int
open_renz_file (char *file_name,
		char *inlist,
		int num_items,
		R_Enz **res_enzyme,
		int *num_enzymes);

char *
SetREnzColour(int num_enzymes,
	      int enz_num);

char *
AddCutSites(char *seq,
	    int cut_site);

/*
 * read a file of enzymes and parse those enzymes that have been selected ie
 * are in 'inlist'
 * store only the enzyme names in array enz_names
 */
int
r_enz_file_names(char *file_name,
		 char ***enz_names,
		 int *num_enzymes);

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
	    int *total_matches);                                     /* out */


int search_dna (char *seq,                                             /* in */
		int seq_len,                                           /* in */
		char *word,                                            /* in */
		int word_len,                                          /* in */
		int circle,                                            /* in */
		int matches[],                                        /* out */
		int max_matches,                                      /* out */
		int *num_matches,                                     /* out */
		int *seq_hash_values);                                /* out */

void hash_dna(char *seq,                                               /* in */
	      int seq_len,                                             /* in */
	      int *seq_hash_values,                                   /* out */
	      int *last_word,   
	      int *word_count);

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
	       int *num_matches);                                    /* out */
	       
int
PrintEnzymeByEnzyme(R_Enz *r_enzyme,                                   /* in */
		    R_Match *match,                                    /* in */
		    int total_matches,                                 /* in */
		    int num_enzymes,                                   /* in */
		    char *sequence,                                    /* in */
		    int sequence_len,                                  /* in */
		    int sequence_type,                                 /* in */
		    int lreg,                                          /* in */
		    int do_all);                                       /* in */

int OrderOnPosition (R_Enz *r_enzyme, R_Match *match, int total_matches,
		     char *sequence, int sequence_len, int sequence_type,
		     int lreg);

int find_max_cut_dist(R_Enz * r_enzyme, int num_enzymes);

#endif


