/* screen_seq: a program for finding readings that appear to be known contaminating
   sequence such as E. coli. Readings found to contain a long segment matching a known
   possible contaminant are given a comment line to indicate this and their names
   are not written to the pass file, but go to the fail file.
   The files to screen against are obtained from a file of file names or can be
   a single file.
   It is a very quick screen: the best match is found (say >25) and then all other 8 base
   matches within window/2 either side are added up and the percentage overlap found is
   used as a test.
*/

#include <staden_config.h>

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include "misc.h" 
#include "dna_utils.h"
#include <io_lib/scf.h>
#include <io_lib/expFileIO.h>
#include "getfile.h"

#define NORMAL_MODE 0
#define TEST_ONLY 1

#define MAX_READ 4096 /* the reading */
#define MAX_VECTOR_D 100000 /* the vector sequence */
#define MIN_VECTOR 4096 /* the vector sequence */
#define MAX_READS 10000
#define MAX_VECTORS 5000

char *file_names[MAX_READS];
char *vfile_names[MAX_VECTORS];


int get_text_seq ( char *seq, int max_len, int *seq_len, FILE *fp );


void write_dot(void) {
    fprintf ( stdout, "." );
    (void) fflush ( stdout );
}

int vep_error ( FILE *fp, char *file_name, int error_no ) {
    char *err_mess[] = {
/* 1 */	"Error: could not open experiment file",
/* 2 */	"Error: no sequence in experiment file",
/* 3 */	"Error: sequence too short",
/* 4 */	"Error: could not write to experiment file",
/* 5 */ "Error: hashing problem",
/* 6 */ "Error: invalid sequence for demonstration mode"};

    fprintf ( stderr, "%s\n",err_mess[error_no-1] );
    if ( fp ) fprintf ( fp, "%s %s\n",file_name,err_mess[error_no-1]);
    fprintf ( stdout, "!" );
    (void) fflush ( stdout );
    return 0;
}


int get_filenames ( FILE *fp ) {

  int num_read = -1;
  char file_name[FILENAME_MAX+1];

  while ( fgets ( file_name, FILENAME_MAX, fp )) {

    num_read++;

    if ( num_read == MAX_READS-1) {

      fprintf(stderr, "too many readings\n");
      return -1;
    }

    file_name[strlen(file_name)-1] = '\0';

    if ( NULL == (file_names[num_read] = (char *) malloc ( sizeof(char *)*(strlen(file_name)+1) ))) {
      return -2;
    }

    strcpy ( file_names[num_read], file_name );
  }
  num_read++;

  return num_read;
}

int get_vfilenames ( FILE *fp, char *fofn_name ) {

  int num_read = -1;
  char file_name[FILENAME_MAX+1], tmp[FILENAME_MAX+1];
  char expanded_fn[FILENAME_MAX+1], base_name[FILENAME_MAX+1];
  char *p;

  if (1 != expandpath(fofn_name, expanded_fn)) {
    fprintf(stderr, "Failed to expand input file of vector file names\n");
    return -1;
  }

  if (p = strrchr(expanded_fn, '/')) {
      strncpy(base_name, expanded_fn, p-expanded_fn+1);
      base_name[p-expanded_fn+1] = 0;
  } else {
      base_name[0] = 0;
  }

  while ( fgets ( tmp, FILENAME_MAX, fp )) {
    
    num_read++;

    if ( num_read == MAX_VECTORS-1) {

      fprintf(stderr, "too many sequences\n");
      return -1;
    }

    if (1 != expandpath(tmp, file_name)) {
      fprintf(stderr, "Failed to expand vector file name\n");
      return -1;
    }

    file_name[strlen(file_name)-1] = '\0';
#ifdef _WIN32
    if (file_name[0] != '/' && file_name[1] != ':') {
#else
    if (file_name[0] != '/') {
#endif
	sprintf(tmp, "%s%s", base_name, file_name);
	strcpy(file_name, tmp);
    }

    if ( NULL == (vfile_names[num_read] = (char *) malloc ( sizeof(char *)*(strlen(file_name)+1) ))) {
      return -2;
    }

    strcpy ( vfile_names[num_read], file_name );

  }
  num_read++;


  return num_read;
}

int get_vfilenames_old ( FILE *fp, char *fofn_name ) {

  int num_read = -1;
  char file_name[FILENAME_MAX+1], tmp[FILENAME_MAX+1];
  char base_name[FILENAME_MAX+1];
  char *p;

  if (p = strrchr(fofn_name, '/')) {
      strncpy(base_name, fofn_name, p-fofn_name+1);
      base_name[p-fofn_name+1] = 0;
  } else {
      base_name[0] = 0;
  }

  while ( fgets ( file_name, FILENAME_MAX, fp )) {
    
    num_read++;

    if ( num_read == MAX_VECTORS-1) {

      fprintf(stderr, "too many sequences\n");
      return -1;
    }

    file_name[strlen(file_name)-1] = '\0';
#ifdef _WIN32
    if (file_name[0] != '/' && file_name[1] != ':') {
#else
    if (file_name[0] != '/') {
	strcpy(tmp, file_name);
	sprintf(file_name, "%s%s", base_name, file_name);
#endif
    }

    if ( NULL == (vfile_names[num_read] = (char *) malloc ( sizeof(char *)*(strlen(file_name)+1) ))) {
      return -2;
    }

    strcpy ( vfile_names[num_read], file_name );
  }
  num_read++;


  return num_read;
}

/************************************************************/

void store_hash ( 
	         int *hash_values, 	/* the hash values for each position in a seq */
		 int seq_len, 		/* size of the seq and hash array */
		 int *last_word, 	/* last occurrence of this hash value */
		 int *word_count, 	/* frequency of each hash value or word */
		 int word_length,       /* word length */
		 int size_wc ) {	/* number of elements in word_count and first_word */


/* 	store the hash values in hash_values: put number of occurrences of
	each hash value in word_count; put the array position of the last 
	occurrence of each hash value in last_word, and previous
	occurrences in hash_values[last_word]. 
	Note that words containing unknown characters (like '-') are given
	hash value -1. So we skip them here, and they are ignored.
*/


    int nw;
    register int i,j,n;

/* 	zero the word counts */

    for ( i=0;i<size_wc;i++ ) {
	word_count[i] = 0;
	last_word[i] = 0;
    }

/*	loop for all entries in hash_values	*/

    j = seq_len - word_length + 1;

    for ( i = 0; i < j; i++ ) {

	n = hash_values[i];

/*	is it a good value ? */

	if ( -1 != n ) {

	    nw = word_count[n];

/* 		already an entry for this word ? */

	    if ( 0 == nw ) {

/*		no, so put in last_word */

		last_word[n] = i;
		word_count[n] += 1;
	    }

/*		yes, so put previous last occurrence in hash_values*/

	    else {

		word_count[n] += 1;
		hash_values[i] = last_word[n];
		last_word[n] = i;
	    }
	}
    }
}


int dna_hash8_lookup[256];

void set_hash8_lookup(void) {

/*	hashing values */
/* 	set up table of values for permitted base characters */

    int i;

    for (i=0;i<256;i++) dna_hash8_lookup[i] = 4;

    dna_hash8_lookup['a'] = 0;
    dna_hash8_lookup['c'] = 1;
    dna_hash8_lookup['g'] = 2;
    dna_hash8_lookup['t'] = 3;
    dna_hash8_lookup['A'] = 0;
    dna_hash8_lookup['C'] = 1;
    dna_hash8_lookup['G'] = 2;
    dna_hash8_lookup['T'] = 3;
/*    dna_hash8_lookup['*'] = 0; */

}

int init_hash8 ( int seq1_len, int seq2_len,
	        int **hash_values1, int **last_word,
	        int **word_count, int **hash_values2,
	        int **diag, int **line ) {
    int size_hash, word_length;

    word_length = 8;

    set_hash8_lookup ();

    size_hash = 65536;


    if ( NULL == (*hash_values1 = (int *) xmalloc ( sizeof(int)*(seq1_len) ))) {
	return -2;
    }

    if ( ! (*last_word = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	return -2;
    }

    if ( ! (*word_count = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	return -2;
    }


    if ( ! (*hash_values2 = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	return -2;
    }

    if ( ! (*line = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	return -2;
    }

    if ( ! (*diag = (int *) xmalloc ( sizeof(int)*(seq2_len + seq1_len) ))) {
	return -2;
    }

    return 0;
}

void free_hash8 ( int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2,
	        int *diag) {

    if ( hash_values1 ) xfree ( hash_values1 );
    if ( hash_values2 ) xfree ( hash_values2 );
    if ( word_count )   xfree ( word_count );
    if ( last_word )    xfree ( last_word );
    if ( diag )         xfree ( diag );

}

int hash_word8 ( char *seq, int *start_base, int seq_len,
	      unsigned short *uword) {

    /* 	given a sequence seq, return the hash value for the first word 
     *  after start_base that does not contain an unknown char. Tell 
     *  the caller where this is. If we reach the end of the seq set
     *  start_base and return -1.
     */


    register int i, word_len=8;
    register int end_base,base_index,lstart_base;
    register int unsigned short luword;

    lstart_base = *start_base;
    end_base = lstart_base + word_len;
    if ( seq_len < end_base ) return -1;

    for (i=lstart_base,luword=0,end_base=lstart_base+word_len;i<end_base;i++) {

	base_index = dna_hash8_lookup[(unsigned)seq[i]];
	if ( 4 == base_index ) {

	    /*	weve hit an unknown char, so lets start again */

	    lstart_base = i + 1;
	    end_base = lstart_base + word_len;
	    if ( seq_len < end_base ) {
		*start_base = lstart_base;
		return -1;
	    }
	    luword = 0;
	    i = lstart_base - 1;
	}
	else {
	    luword = ( luword <<2 ) | base_index;
	}
    }
    *start_base = lstart_base;
    *uword = luword;
    return 0;
}


int hash_seq8 ( char *seq, int *hash_values, int seq_len) {

/* given a sequence seq, return an array of hash values
   If we cannot find at least one word to hash on we return -1
   otherwise we return 0.
*/

    register int i,j,k,word_len=8;
    int start_base,prev_start_base,end_base,base_index;
    unsigned short uword;

    if ( seq_len < word_len ) return -1;

    /*	Get the hash value for the first word that contains no unknowns */	
    start_base = 0;
    if (hash_word8 ( seq, &start_base, seq_len, &uword)) return -1;

    for (i=0;i<start_base;i++) hash_values[i] = -1;

    /*	Now do the rest of the sequence */

    hash_values[start_base] = uword;
    k = seq_len - word_len + 1;

    for (i=start_base+1,j=start_base+word_len; i<k; i++,j++) {

	base_index = dna_hash8_lookup[(unsigned)seq[j]];
	if ( 4 == base_index ) {

	    /*	weve hit an unknown char, so lets start again */

	    prev_start_base = i;
	    start_base = j + 1;
	    if (hash_word8 ( seq, &start_base, seq_len, &uword)) {
		for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
		return 0;
	    }

	    for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
	    hash_values[start_base] = uword;
	    end_base = start_base + word_len;
	    i = start_base;
	    j = i + word_len - 1;

	}
	else {
	    uword = ( uword <<2 ) | base_index;
	    hash_values[i] = uword;
	}
    }
    return 0;
}

int do_hash_con (   int seq1_len, int seq2_len,
	        int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2,
		int *diag, int *line,
	        char *seq1, char *seq2, int new_seq1, 
	        int min_match, int half_window, int *x, int *y, int *score ) {

    int nrw, word, pw1, pw2, i, ncw, j, match_length, word_length = 8;
    int size_hash = 65536;
    int diag_pos;
    int top_score, min_line, max_line, line_pos, i1, i2, j1, j2, k;
    int top_score_pos, new_score, xx, yy;
    double n_score;

    /* seq2 is the reading which changes each entry, but seq1 is
       the vector which may be the same for consecutive
       entries, so we use new_seq1 to tell if we need to hash seq1
    */

	
    *score = 0;
    if ( seq1_len < min_match ) return -4; 
    if ( seq2_len < min_match ) return -4; 


    if ( new_seq1 ) {

	if ( hash_seq8 ( seq1, hash_values1, seq1_len )  != 0 ) {
	    return -1;
	}

	(void) store_hash ( hash_values1, seq1_len, last_word, word_count,
			    word_length, size_hash);
    }

    if ( hash_seq8 ( seq2, hash_values2, seq2_len )  != 0 ) {
	return -1;
    }

    j = seq1_len + seq2_len;
    for (i=0;i<j;i++) diag[i] = -word_length;

    nrw = seq2_len - word_length + 1;

/* 	loop for all (nrw) complete words in hash_values2 */


    for (pw2=0;pw2<nrw;pw2++) {

 	word = hash_values2[pw2];

	if ( -1 != word ) {

	    if ( 0 != (ncw = word_count[word]) ) {

		for (j=0,pw1=last_word[word];j<ncw;j++) {

		    diag_pos = seq1_len - pw1 + pw2 - 1;

		    if ( diag[diag_pos] < pw2 ) {

			if ((match_length = match_len ( 
						  seq1, pw1, seq1_len,
						  seq2, pw2, seq2_len))
			    >= min_match ) {

			  *score = match_length;
			  *x = pw1+1;
			  *y = pw2+1;
			  return 1;
			}
			diag[diag_pos] = pw2 + match_length;
		    }
		    pw1 = hash_values1[pw1];
		}
	    }
	}
    }

    if ( *score < min_match ) return 0;

    /* NONE OF THIS USED */


    /* check out this region by adding up the matches
       in all the loca hits to see how much of the
       reading they cover */

    top_score_pos = 0; /* to keep marks compiler happy */
    if ( top_score_pos < seq1_len ) {
	yy = seq1_len - top_score_pos - 1;
	xx = 0;
    }
    else {
	xx = top_score_pos - seq1_len + 1;
	yy = 0;
    }

    top_score = *score;

    min_line = yy;
    max_line = MIN(seq1_len-yy,seq2_len-xx) - (word_length-1);

    for ( i=0;i<max_line;i++) line[i] = 0;
    i1 = MAX(0,xx-half_window);
    i2 = MIN(i1+max_line+half_window,seq2_len-word_length);

    j1 = top_score_pos-half_window;
    j2 = top_score_pos+half_window;

    for (i=i1;i<i2;i++) {
	word = hash_values2[i];
	if ( -1 != word ) {
	    if ( 0 != (ncw = word_count[word]) ) {
		pw1 = last_word[word];
		for (j=0;j<ncw;j++) {
		    k = i - pw1 + seq1_len - 1;
		    if((k>j1)&&(k<j2)) {
			line_pos = pw1 - min_line;

			if ( (line_pos > 0) && (line_pos < max_line) ) line[line_pos]=1;
		    }
		    pw1 = hash_values1[pw1];
		}
	    }
	}
    }
    for (i=0,new_score=0;i<max_line;i++) new_score += line[i];
    n_score = 100.0 * (double)new_score/(double)(max_line-1);
    *score = (int) (n_score + 0.5);
    return 1;
}

int do_it_con ( char *vector_seq, int max_vector,
	       FILE *fp_s, FILE *fp_i, FILE *fp_p, FILE *fp_f,
	       int word_length, int min_match, int percent_cut,
	       int window_size, int tmode, int mode_v, char *fofn_s,
	       int mode_i) {

    char *seq, *expt_file_name;

    Exp_info *e;
    int ql,qr,seq_length,i;

/* for this algorithm */

    int *hash_values1, *hash_values2, *last_word, *word_count;
    int *diag, *line;
    int vector_length, x, y, ret, eret;
    char *vfile_name;
    int num_files, file_num, num_vfiles, vfile_num;
    FILE *vf;
    int sl, sr; /* sequencing vector left and right */
    int new_vector = 0;
    int score, score_f, score_r;
    int lg, rg, xf, xr, yf, yr;
    int match_found, half_window;
    int size_hash = 65536;

    half_window = MAX(1,window_size/2);
    if ( min_match < 8 ) min_match = 8;

    if ( init_hash8 ( max_vector, MAX_READ,
		     &hash_values1, &last_word, &word_count,
		     &hash_values2, &diag, &line ))
	return -1;

    if ( mode_i ) {
	if (( num_files = get_filenames ( fp_i )) < 1 )
	    return -1;
    }
    else {
	num_files = 1;
    }

    if ( mode_v ) {
	if (( num_vfiles = get_vfilenames ( fp_s, fofn_s )) < 1 )
	    return -1;
    }
    else {
	num_vfiles = 1;
    }

    /* for each vector sequence */

    for ( vfile_num = 0; (vfile_num < num_vfiles) && 
	 (vfile_name=vfile_names[vfile_num]); vfile_num++ ) {

      if ( tmode ) {
	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>> %s\n", vfile_names[vfile_num]);
      }

      if ( !(vf = fopen(vfile_names[vfile_num], "r"))) {
	fprintf(stderr, "Error: could not open sequence file %s\n", vfile_name);
	continue;
      }
      ret = get_text_seq ( vector_seq, max_vector, &vector_length, vf);
      fclose(vf);
      if ( ret ) {
	fprintf(stderr, "Error: could not read vector file %s\n", vfile_name);
	continue;
      }

      if ( hash_seq8 ( vector_seq, hash_values1, vector_length ) != 0 ) {
	fprintf(stderr, "Error: could not hash sequence file %s\n", vfile_name);
	continue;
      }
      (void) store_hash ( hash_values1, vector_length, last_word, word_count,
			   word_length, size_hash );

      /* for all the readings that are left */

      for ( file_num = 0; file_num < num_files; file_num++ ) {

	if ( expt_file_name=file_names[file_num]) {


	  if ( tmode ) {
	    printf(">>>>>>>>>>>>>>>>>>>>>> %s\n", expt_file_name );
	  }

	  e = exp_read_info ( expt_file_name );
	  if ( e == NULL ) {
	    eret =  vep_error ( fp_f, expt_file_name, 1 );
	    exp_destroy_info ( e );
	    file_names[file_num] = NULL;
	    continue;
	  }
	  else {
	    if ( exp_Nentries ( e, EFLT_SQ ) < 1 ) {
	      eret =  vep_error ( fp_f, expt_file_name, 2 );
	      exp_destroy_info ( e );
	      file_names[file_num] = NULL;
	      continue;
	    }
	    else {

	      char *expline;

	      seq = exp_get_entry ( e, EFLT_SQ );
	      seq_length = strlen ( seq );
	      ql = 0;
	      qr = seq_length + 1;

	      if ( exp_Nentries ( e, EFLT_QL )) {
		expline = exp_get_entry ( e, EFLT_QL );
		ql = atoi ( expline );
	      }
	      if ( exp_Nentries ( e, EFLT_QR )) {
		expline = exp_get_entry ( e, EFLT_QR );
		qr = atoi ( expline );
	      }

	      sl = 0;
	      sr = seq_length + 1;
	      if ( exp_Nentries ( e, EFLT_SL )) {
		expline = exp_get_entry ( e, EFLT_SL );
		sl = atoi ( expline );
	      }
	      if ( exp_Nentries ( e, EFLT_SR )) {
		expline = exp_get_entry ( e, EFLT_SR );
		sr = atoi ( expline );
	      }

	      /* we have to search both strands so we call do_hash
		 with the read in its original sense, then its
		 complement.
	      */

	      match_found = 0;
	      score_f = score_r = 0;
	      lg = MAX ( ql, sl ) - 1;
	      lg = MAX ( lg, 0 );
	      rg = MIN ( qr, sr ) - 1;
	      rg = MIN ( rg, seq_length-1 );
	      /* 		printf("lg %d rg %d\n",lg,rg); */

	      if ( rg - lg + 1 < min_match ) {
		eret =  vep_error ( fp_f, expt_file_name, 3 );
		exp_destroy_info ( e );
		file_names[file_num] = NULL;
		continue;
	      }
	      ret = do_hash_con ( vector_length, rg-lg+1,
				 hash_values1, last_word,
				 word_count, hash_values2, diag, line,
				 vector_seq, &seq[lg], new_vector,
				 min_match, half_window,
				 &x, &y, &score);
	      if ( ret < 0 ) {
		eret =  vep_error ( fp_f, expt_file_name, 5 );
		exp_destroy_info ( e );
		file_names[file_num] = NULL;
		continue;
	      }

	      if ( ret ) {
		if ( score >= percent_cut ) {
		  xf = x;
		  yf = y + lg;
		  score_f = score;
		  match_found = 1;
		  /* printf("forward score %d at %d %d\n",score_f,xf,yf); */
		  if ( !tmode ) {
		    char mess[2048]; /* twice vfile_name ! */

		    if (exp_put_str(e, EFLT_PS, "contaminated", 
					  strlen("contaminated"))) {
		      eret =  vep_error ( fp_f, expt_file_name, 4 );
		      exp_destroy_info ( e );
		      file_names[file_num] = NULL;
		      continue;
		    }
		    sprintf(mess, "CONT = %d..%d\n%d %d %s",
			    yf,yf+score-1,xf,score,vfile_name);
		    exp_put_str(e, EFLT_TG, mess, strlen(mess));

		    if ( fp_f ) fprintf ( fp_f, "%s\n",expt_file_name);
		    file_names[file_num] = NULL;
		  }
		  else {
		    printf("match %d at %d  %d\n",score,xf,yf);
		  }
		}
	      }

	      if ( !match_found ) {

		complement_seq ( &seq[lg], rg-lg+1);

		ret = do_hash_con ( vector_length, rg-lg+1,
				hash_values1, last_word,
			        word_count, hash_values2, diag, line,
			        vector_seq, &seq[lg], new_vector,
			        min_match, half_window,
			        &x, &y, &score);
		if ( ret < 0 ) {
		  eret =  vep_error ( fp_f, expt_file_name, 5 );
		  exp_destroy_info ( e );
		  file_names[file_num] = NULL;
		  continue;
		}
		if ( ret ) {
		  if ( score >= percent_cut ) {
		    xr = x;
		    yr = rg - lg - y + lg - score + 3;
		    score_r = score;
		    match_found = 1;
		    /*  printf("reverse score %d at %d %d\n",score,xr,yr); */

		    if ( !tmode ) {
		    char mess[2048]; /* twice vfile_name ! */

		      if (exp_put_str(e, EFLT_PS, "contaminated", 
				      strlen("contaminated"))) {
			eret =  vep_error ( fp_f, expt_file_name, 4 );
			exp_destroy_info ( e );
			file_names[file_num] = NULL;
			continue;
		      }
		      sprintf(mess, "CONT = %d..%d\n%d %d %s",
			      yr,yr+score-1,xr,score,vfile_name);
		      exp_put_str(e, EFLT_TG, mess, strlen(mess));

		      if ( fp_f ) fprintf ( fp_f, "%s\n",expt_file_name);
		      file_names[file_num] = NULL;
		    }
		    else {
		      printf("---match %d at %d  %d\n",score,xr,yr);
		    }
		  }
		  match_found = 1;
		}
	      }
	      if ( tmode && !match_found ) printf("no match\n");
	    }
	  }
	  exp_destroy_info ( e );
	  if (!(tmode)) (void) write_dot();
	}
      }
    }
      /*fprintf ( fp_p, "%s\n",expt_file_name); */

      /* screening finished so write out those that have not failed */

    if ( fp_p ) {
      for(i=0;i<num_files;i++) {
	if ( file_names[i] ) fprintf( fp_p, "%s\n", file_names[i]);
      }
    }

    free_hash8 ( hash_values1, last_word,
	         word_count, hash_values2,
	         diag );
    return 0;
}

int get_text_seq ( char *seq, int max_len, int *seq_len, FILE *fp )

/* read in a staden format (yuk) sequence file */

/* Deal with 2 special line types: comments that have ";" in column 0
   and contig consensus sequence headers that have "<----abc.00001---->"
   embedded in them */
{

#define MAX_SEQ_LINE 101

    char line[MAX_SEQ_LINE];
    int j;

    *seq_len = 0;
    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	/* Check for special lines of type ";"*/

	if ( ';' != line[0] ) {

	    for (j = 0;j < MAX_SEQ_LINE && line[j]; j++) {

		if ( '<' == line[j] ) j += 20;
		if (isalpha ( (int) line[j]) || (int) line[j] == '-') {
		    if ( *seq_len >= max_len) return -1;
		     seq[*seq_len] = line[j];
		    *seq_len += 1;
		}
	    }
	}
    }
    return 0;
}


/* johnt 30/6/99 must explicitly import globals from DLLs with Visual C++*/
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif
 
extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;


void usage(int min_match ) {
    fprintf(stderr, 
	    "Usage: screen_seq [options and parameters]\n"
	    "Where options are:\n"
	    "    [-l minimum match (%d)]           [-m Max sequence length (%d)]\n"
	    "    [-i readings to screen fofn]      [-I reading to screen]\n"
	    "    [-s seqs to screen against fofn]  [-S seq to screen against]\n"
	    "    [-t test only]\n"
	    "    [-p passed fofn]                  [-f failed fofn]\n",
	     min_match, MAX_VECTOR_D);
    exit(1);
}

/* program to find readings that are contaminant DNA */

int main(int argc, char **argv) {
    int c;
    int min_match, max_vector, min_match_d, percent_cut;
    int mode_v, mode_i,i,tmode,mr_s,mr_fofn, mv_s,mv_fofn;
    int window_size;
    char *fofn_p, *fofn_f, *fofn_i, *fofn_s = NULL, *vector_seq;
    FILE *fp_p, *fp_f, *fp_i, *fp_s = NULL;

    fofn_p = fofn_f = fofn_i = NULL;
    fp_p = fp_f = fp_i = NULL;

    max_vector = MAX_VECTOR_D;
    min_match_d = 25;
    min_match = -1;
    percent_cut = -1;
    window_size = -1;
    mode_v = mode_i = -1;
    tmode = 0;
    mr_s = mr_fofn = mv_s = mv_fofn = 0;

    while ((c = getopt(argc, argv, "l:m:i:I:p:f:s:S:t")) != -1) {
	switch (c) {
	case 'l':
	    min_match = atoi(optarg);
	    break;
	case 'm':
	    max_vector = atoi(optarg);
	    break;
	case 't':
	    tmode = TEST_ONLY;
	    break;
	case 'i':	/* sequences to screen fofn */
	    mr_fofn = 1;
	    fofn_i = optarg;
	    break;
	case 'I':	/* sequence to screen */
	    mr_s = 1;
	    fofn_i = optarg;
	    break;
	case 's':	/* sequences to screen against fofn */
	    fofn_s = optarg;
	    mv_fofn = 1;
	    break;
	case 'S':	/* sequence to screen against */
	    fofn_s = optarg;
	    mv_s = 1;
	    break;
	case 'p':	/* passed file fofn */
	    fofn_p = optarg;
	    break;
	case 'f':	/* fails file fofn */
	    fofn_f = optarg;
	    break;
	default:
	    usage( min_match_d);
	}
    }

    if ( optind < 2 ) usage( min_match_d );
    if ( mr_s && mr_fofn ) usage( min_match_d );
    if ( mv_s && mv_fofn ) usage( min_match_d);
    if ( min_match < 0 ) min_match = min_match_d;
    if ( max_vector < MIN_VECTOR ) max_vector = MIN_VECTOR;

    if ( mr_fofn ) {
      fp_i = fopen(fofn_i, "r");
      if (fp_i == NULL ) {
	fprintf(stderr, "Failed to open file of file names to screen\n");
	return -1;
      }
      mode_i = 1;
    }

    if ( mr_s ) {
	if ( NULL == (file_names[0] = 
		      (char *) malloc ( sizeof(char *)*(strlen(fofn_i)+1) ))) {
	    fprintf(stderr, "Failed to open single file to screen\n");
	    return -1;
	}
	file_names[0] = fofn_i;
	mode_i = 0;
    }

    if ( mv_fofn ) {
	fp_s = fopen(fofn_s, "r");
	if (fp_s == NULL ) {
	    fprintf(stderr, "Failed to open file of file names to screen against\n");
	    return -1;
	}
	mode_v = 1;
    }

    if ( mv_s ) {
	if ( NULL == (vfile_names[0] = 
		      (char *) malloc ( sizeof(char *)*(strlen(fofn_s)+1) ))) {
	    fprintf(stderr, "Failed to open single file to screen against\n");
	    return -1;
	}
	vfile_names[0] = fofn_s;
	mode_v = 0;
    }

    if ( fofn_p ) {
	fp_p = fopen(fofn_p, "w");
	if (fp_p == NULL ) {
	    fprintf(stderr, "Failed to open file of passed file names\n");
	    return -1;
	}
    }
    if ( fofn_f ) {
	fp_f = fopen(fofn_f, "w");
	if (fp_f == NULL ) {
	    fprintf(stderr, "Failed to open file of failed file names\n");
	    return -1;
	}
    }
    set_dna_lookup();
    set_char_set(1); /* FIXME DNA*/

    if ( ! (vector_seq = (char *) xmalloc ( sizeof(char)*max_vector ))) return -1;

    if ( (mode_v == -1) || (mode_i == -1) ) usage( min_match_d );

    percent_cut = min_match;
    window_size = 7; 

    i = do_it_con ( vector_seq, max_vector, fp_s, fp_i, fp_p, fp_f, 8, min_match,
		   percent_cut, window_size, tmode, mode_v, fofn_s, mode_i );

    fprintf(stdout,"\n");
    xfree ( vector_seq );
    exit(0);
}
