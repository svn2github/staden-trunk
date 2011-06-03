#include <string.h>
#include <stdio.h>
#include <math.h>
#include "tkSeqed.h"
#include "xalloc.h"
#include "tkSeqedUtils.h"
#include "text_output.h"
#include "dna_utils.h"
#include "tcl_utils.h"

static int *pos = NULL;
static int *score = NULL;
static int current_match;
static int n_matches;
static int prev_cursor_pos;

void seqed_string_search_free(void)
{
    xfree(pos);
    xfree(score);
    pos = NULL;
    score = NULL;
}

int seqed_next_string(tkSeqed *se,
		      int direction)
{
    int i;

    /* find first match after current cursor position */
    if (direction == 0) {
	/* moved cursor so need to find nearest match */
	if (prev_cursor_pos != se->cursorPos) {
	    current_match = n_matches;
	    for (i = 0; i < n_matches; i++) {
		if (pos[i] >= se->cursorPos) {
		    current_match = i;
		    break;
		}
	    }
	} else {
	    current_match++;
	}
    } else {
	/* reverse search */
	/* moved cursor so need to find nearest match */
	if (prev_cursor_pos != se->cursorPos) {
	    current_match = -1;
	    for (i = n_matches-1; i > 0; i--) {
		if (pos[i] <= se->cursorPos) {
		    current_match = i;
		    break;
		}
	    }
	} else {
	    current_match--;
	}
    }

    if (current_match < 0) {
	bell();
	current_match = 0;
	return -1;
    }

    if (current_match >= n_matches) {
	bell();
	current_match = n_matches - 1;
	return -1;
    }

    seqed_setCursorPos(se, pos[current_match]);
    seqed_showCursor(se, se->cursorSeq, pos[current_match]); 
    prev_cursor_pos = se->cursorPos;
    return 0;
}

int seqed_string_search(char *sequence,
			int seq_len,
			char *seq_name,
			char *string,
			int direction,
			int strand,
			double per_match,
			int cursorPos,
			int use_iub_code)
{
    int max_matches;
    int min_match;
    int i;
    int string_length;
    char *seq_match;
    Tcl_DString input_params;

    vfuncheader("Search string");

    max_matches = seq_len;
    string_length = strlen(string);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    {
      char tmp[10];
      char sstrand[10];
      char sdirection[10];
      if (direction)
	  strcpy(sdirection, "backward");
      else
	  strcpy(sdirection, "forward");

      if (strand)
	  strcpy(sstrand, "reverse");
      else 
	  strcpy(sstrand, "forward");

      if (use_iub_code)
	strcpy(tmp, "iub");
      else
	strcpy(tmp, "literal");
      vTcl_DStringAppend(&input_params, "sequence %s\n"
			 "direction %s\nstrand %s\nuse %s code\nminimum percent match %f\nstring %s\n",
			 seq_name, sdirection, sstrand, tmp, per_match, string);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
    
    /* if previously allocated pos and score, then need to free them here */
    if (pos != NULL) {
	seqed_string_search_free();
    }

    if (NULL == (pos = (int *)xmalloc((max_matches + 1) * sizeof(int))))
	return -1;

    if (NULL == (score = (int *)xmalloc((max_matches + 1) * sizeof(int))))
	return -1;

    if (NULL == (seq_match = (char *)xmalloc((string_length + 1) * sizeof(char))))
	return -1;

    /* convert percentage mis-matches into min matches */
    min_match = ceil(string_length * per_match / 100);

    /* reverse & complement to search from 5' to 3' on complementary strand */
    if (strand == 1) {
	complement_seq(string, string_length);
    }
    
    n_matches = iubc_inexact_match(sequence, seq_len, string, strlen(string),
				   min_match, use_iub_code, pos, score, 
				   max_matches);

    if (n_matches < 0) {
	vmessage("String search: no matches found\n");
	return -1;
    }

    for (i = 0; i < n_matches; i++) {
	vmessage("Position %d score %d", pos[i], score[i]);

	strncpy(seq_match, &sequence[pos[i]-1], string_length);
	seq_match[string_length] = '\0';
	iubc_list_alignment(string, seq_match, "string", seq_name, 1, 
			    pos[i], "");
    }
    current_match = -1;
    prev_cursor_pos = -1;
    xfree(seq_match);
    return 0;
}
