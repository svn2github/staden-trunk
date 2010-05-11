#include <staden_config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <os.h>
#include "seqInfo.h"
#include "dna_utils.h"
#include "xalloc.h"
#include <io_lib/misc.h>

/* johnt 1/6/99 must explicitly import globals from DLLs with Visual C++*/
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif
 
extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;

typedef struct {
  int score; 
  int window_len;
  int verbose;
  int test_mode;
  int min_len;
} params;


/*
 * Scans leftwards in a sequence the percentage A and T drops below
 * a specific threshold
 */
int scan_left(params p, char *seq, int start_pos, int len) {
  int i, left_edge, right_edge;
    int counts[5] = {0,0,0,0,0};
    int win_len = p.window_len;
    int score = p.score;
    int poly_char;

    /*printf("start_pos %d %d %d\n",start_pos,len,win_len);*/

    /* paranoia */

    if ( len - start_pos < win_len ) return len;

    /* initialise the counters */

    for (i=len-win_len;i<len;i++) {
      counts[char_lookup[(unsigned int)seq[i]]]+=1;
    }

    /* if the first window is not polyA or T we are done */

    if ((counts[0]<score) && (counts[3]<score)) return -1;

    if ( counts[0] > counts[3] ) {
      poly_char = 0;
    }
    else {
      poly_char = 3;
    }
    /* keep looking left until not polyA or T */

    for (left_edge=len-win_len-1,right_edge=len-1;
	 left_edge>start_pos;left_edge--,right_edge--) {

      /* last window not polyA or T: clip at rightmost C or G */

      if (counts[poly_char]<score) {
	for (;right_edge>left_edge;right_edge--) {
	  if (char_lookup[(unsigned int)seq[right_edge]] != poly_char) 
	    return right_edge + 2;
	}
	return right_edge;
      }
      counts[char_lookup[(unsigned int)seq[left_edge]]] += 1;
      counts[char_lookup[(unsigned int)seq[right_edge]]] -=1;
    }

    /* all polyA or T! */

    return start_pos + 1;
}


/*
 * Scans rightwards in a sequence the percentage A and T drops below
 * a specific threshold
 */
int scan_right(params p, char *seq, int start_pos, int len) {
  int i, left_edge, right_edge;
    int counts[5] = {0,0,0,0,0};
    int win_len = p.window_len;
    int score = p.score;
    int poly_char;

    /* paranoia */
    if ( len - start_pos < win_len ) return len;

    /* initialise the counters */

    for (i=start_pos;i<start_pos+win_len;i++) {
      counts[char_lookup[(unsigned int)seq[i]]]+=1;
    }

    /* if the first window is not polyA or T we are done */

    if ((counts[0]<score) && (counts[3]<score)) return -1;

    if ( counts[0] > counts[3] ) {
      poly_char = 0;
    }
    else {
      poly_char = 3;
    }
    
    /* keep looking left until not polyA or T */

    for (left_edge=start_pos,right_edge=start_pos+win_len;
	 right_edge<len;left_edge++,right_edge++) {

      /* last window not polyA or T: clip at rightmost C or G */

      if (counts[poly_char]<score) {
	for (;left_edge<right_edge;left_edge++) {

	  if (char_lookup[(unsigned int)seq[left_edge]] != poly_char) 
	    return left_edge;
	}
	return left_edge;
      }
      counts[char_lookup[(unsigned int)seq[left_edge]]] -= 1;
      counts[char_lookup[(unsigned int)seq[right_edge]]] +=1;
    }

    /* all polyA or T! */

    return len + 1;
}

/*
 * Quality clips file 'file', updating the QL and QR records in the process.
 *
 * Returns 0 for success.
 *        -1 for failure.
 */
static int polyA_clip(char *file, params p) {
    SeqInfo *si = NULL;
    int i1, i2, seq_length, right_pos, left_pos;
    char *seq;
    char *expline;
    FILE *fp;

    if (p.verbose)
	printf("Clipping file %s\n", file);

    /* Read the sequence and confidence */
    if (NULL == (si = read_sequence_details(file, 0))) {
	fprintf(stderr, "Failed to read file '%s'\n", file);
	return -1;
    }

    seq = exp_get_entry(si->e, EFLT_SQ);
    seq_length = strlen ( seq );
		
    i1 = i2 = 0;
    if ( exp_Nentries ( si->e, EFLT_QL )) {
      expline = exp_get_entry ( si->e, EFLT_QL );
      i1 = atoi ( expline );
    }
    if ( exp_Nentries ( si->e, EFLT_SL )) {
      expline = exp_get_entry ( si->e, EFLT_SL );
      i2 = atoi ( expline );
    }

    left_pos = MAX(i1,i2);

    i1 = i2 = seq_length - 1;
    if ( exp_Nentries ( si->e, EFLT_QR )) {
      expline = exp_get_entry ( si->e, EFLT_QR );
      i1 = atoi ( expline );
    }
    if ( exp_Nentries ( si->e, EFLT_SR )) {
      expline = exp_get_entry ( si->e, EFLT_SR );
      i2 = atoi ( expline );
    }
    right_pos = MIN(i1,i2);

    i1 = right_pos;
    right_pos = scan_left(p, seq, left_pos, right_pos);
    
    if (right_pos > 0 ) {
      if (right_pos - left_pos < p.min_len) {
	fprintf(stderr, "Sequence too short (length=%d)\n",
		right_pos - left_pos);
	freeSeqInfo(si);
	return -1;
      }

      /* Append details onto the end of the Exp File */
      if (!p.test_mode) {
	if (NULL == (fp = fopen(file, "a"))) {
	  fprintf(stderr, "Failed to write file '%s'\n", file);
	  freeSeqInfo(si);
	  return -1;
	}
	fprintf(fp, "SR   %d\n", right_pos);
	fclose(fp);
      } else {
	printf("%-30s SR %4d\n", file, right_pos);
      }
    }
    else {
      /* first window not polyA or T so do nothing */
      right_pos = i1;
      /*printf("%-30s no PolyA tail\n", file);*/
    }

    left_pos = scan_right(p, seq, left_pos, right_pos);
    if (left_pos > 0 ) {
      if (right_pos - left_pos < p.min_len) {
	fprintf(stderr, "Sequence too short (length=%d)\n",
		right_pos - left_pos);
	freeSeqInfo(si);
	return -1;
      }

      /* Append details onto the end of the Exp File */
      if (!p.test_mode) {
	if (NULL == (fp = fopen(file, "a"))) {
	  fprintf(stderr, "Failed to write file '%s'\n", file);
	  freeSeqInfo(si);
	  return -1;
	}
	fprintf(fp, "SL   %d\n", left_pos);
	fclose(fp);
      } else {
	printf("%-30s SL %4d\n", file, left_pos);
      }
    }
    else {
      /* first window not polyA or T so do nothing */
      /*printf("%-30s no PolyT head\n", file);*/
    }
    freeSeqInfo(si);
    return 0;
}

static void usage(void) {
fprintf(stderr,
	"Usage:\n"
	"polyA_clip [-vt] [-p percent cutoff(95)] [-x min_length(0)]\n"
	"                  [-w window length(50)] file...\n");
    exit(1);
}

int main(int argc, char **argv) {
    int c, i, ret = 0;
double perc;
    params p;

    /* Defaults */
    p.min_len = 0;
    p.verbose = 0;
    p.window_len = 50;
    p.test_mode = 0;
    perc = 95.0;
    set_dna_lookup();
    set_char_set(1);
    while ((c = getopt(argc, argv, "w:x:p:vtx")) != -1) {
	switch (c) {

	case 'v':
	    p.verbose = 1;
	    break;

	case 't':
	    p.test_mode = 1;
	    break;

	case 'x':
	    p.min_len = atoi(optarg);
	    break;

	case 'w':
	    p.window_len = atoi(optarg);
	    break;

	case 'p':
	    perc = atof(optarg);
	    break;

	default:
	    usage();
	}
    }

    if (optind == argc)
	usage();

    p.score = p.window_len * perc/100.0;

    for (i = optind; i < argc; i++) {
	int ret_val;

	ret_val = polyA_clip(argv[i], p);
	if (p.verbose)
	    printf("    polyA_clip() returned %d\n", ret_val);

	ret |= ret_val;
    }

    return ret ? 1 : 0;
}
