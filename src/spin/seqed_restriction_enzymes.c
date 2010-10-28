#include <string.h>
#include "dna_utils.h"
#include "tkSeqed.h"
#include "renz_utils.h"
#include "seqed_restriction_enzymes.h"
#include "misc.h"
#include "tkSeqedUtils.h"


static R_Enz *global_r_enzyme;
/*
static int num_enzymes;
*/
static int max_lines = 0;
static char **sline = NULL;
static int max_overlap = -1; /* max overlap needed at ends of sequence */
static int max_name_length = -1; /* max renzyme name length */

void free_lines(void)
{
    int i;
    
    for (i = 0; i < max_lines; i++) {
	xfree(sline[i]);
    }
    xfree(sline);
    sline = NULL;
    max_lines = 0;
}

void free_r_enzyme(R_Enz *r_enzyme,
		   int num_enzymes)
{
    int i, j;

    if (r_enzyme) {
	for (i = 0; i < num_enzymes; i++) {
	    xfree(r_enzyme[i].name);
	    for (j = 0; j < r_enzyme[i].num_seq; j++) {
		xfree(r_enzyme[i].seq[j]);
	    }
	    xfree(r_enzyme[i].seq);
	    xfree(r_enzyme[i].cut_site);
	}
	xfree(r_enzyme);
	r_enzyme = NULL;
    }
}

int seqed_get_max_name_length()
{
    return max_name_length;
}

int seqed_add_more_lines(int line_length,
			 char ***lines,
			 int *max_renz_lines)
{
    int incr = 10;
    int prev_lines = *max_renz_lines;
    int i;
    char **alines = *lines;

    *max_renz_lines += incr;

    if (NULL == (alines = (char **)xrealloc(alines, *max_renz_lines * 
					   sizeof(char *))))
	return -1;

    for (i = prev_lines; i < *max_renz_lines; i++) {
	if (NULL == (alines[i] = (char *)xmalloc((line_length+1) 
						* sizeof(char))))
	    return -1;
	memset(alines[i], ' ', line_length);
	alines[i][line_length] = '\0';
    }
    *lines = alines;
#ifdef DEBUG
    printf("seqed_add_more_lines %d line_length %d\n", *max_renz_lines, line_length);
#endif
    return 0;
}

int seqed_redisplay_renzyme(tkSeqed *se, 
			    int pos,
			    int keep_x_pos)
{
    int num_lines;
    int i, j;
    int display_pos;
    int display_width;

    display_width = MIN(se->displayWidth, se->seq_len);
    
    if (-1 == (seqed_write_renzyme(se->sequence, se->sequence_type, 
				   se->r_enzyme, se->num_enzymes, pos, 
				   display_width, 0,
				   &sline, &max_lines, &num_lines)))
	return -1;

    se->renz_lines = num_lines;

    /* need to call set_lines here to give the correct se->lines[RENZ] value */
    set_lines(se, 0, keep_x_pos);
    display_pos = se->lines[SEQ] - se->anchor_pos;
    set_lines(se, display_pos, keep_x_pos);

    j = se->lines[RENZ];
    for (i = num_lines-1; i >= 0; i--) { 
	/* printf("i %d j %d line %s\n", i, j, sline[i]); */
	XawSheetPutText(&se->sw, 0, j++, se->displayWidth, sline[i]); 
    }
    return 0;
}


static int compare_rmatch_name(const void *p1, const void *p2)
{
    R_Match *r1 = (R_Match *) p1, *r2 = (R_Match *) p2;
    return(strcmp(global_r_enzyme[(*r2).enz_name].name, 
		  global_r_enzyme[(*r1).enz_name].name));
}

#ifdef DEBUG
int 
compare_rmatch_name( R_Match *r1, R_Match *r2)
{
    printf("compare \n");
    if (strlen(r_enzyme[(*r1).enz_name].name) < 
	(strlen(r_enzyme[(*r2).enz_name].name)))
	return (-1);
    else if (strlen(r_enzyme[(*r1).enz_name].name) == 
	     (strlen(r_enzyme[(*r2).enz_name].name)))
	return (0);
    else
	return (1);
}
#endif

int compare_rmatch_rev(const void *p1, const void *p2)
{
    R_Match *r1 = (R_Match *) p1, *r2 = (R_Match *) p2;
    if ((*r1).cut_pos > (*r2).cut_pos) {
	return (-1);
    }
    else if ((*r1).cut_pos == (*r2).cut_pos) {
	return (0);
    }
    else {
	return (1);
    }
}

int seqed_write_renzyme(char *seq,
			int sequence_type,
			R_Enz *r_enzyme,
			int num_enzymes,
			int pos, 
			int line_length,
			int name_overlap,
			char ***alines,
			int *max_renz_lines,
			int *num_lines)
{
    int i, j, k;
    int overlap;
    R_Match *match;
    int total_matches;
    int c_pos;
    int name_len;
    int line_num = 0;
    char **lines = *alines;
    int prev_pos, cnt;
    int seq_len = strlen(seq) - 4;
    char sequence[2*MAX_DISPLAY_WIDTH]; /* make large to accomodate overlaps */
    int start;

    *num_lines = 0;

    /* initialise all lines to blank */
    for (i = 0; i < *max_renz_lines; i++) {
	memset(lines[i], ' ', line_length + name_overlap);
	lines[i][line_length+name_overlap] = '\0';
    }

#ifdef DEBUG
    printf("overlap %d pos %d seqlen %d line len %d\n", overlap, pos, seq_len, line_length);
#endif


    /*
     * Initialise sequence by copying the relevant segments of seq[], being
     * careful to handle circular sequences.
     */
    overlap = max_overlap + max_name_length;
    {
	int start_from = pos - overlap;
	int start_to = 0;
	int end_from = pos + line_length + overlap;

	memset(sequence, 'N', line_length + 2*overlap);
	if (start_from < 1) {
	    start_to = overlap - (pos - 1);
	    start_from = 1;
	}
	if (end_from > seq_len) {
	    end_from = seq_len + 1;
	}

	memmove(&sequence[start_to], &seq[start_from-1+2],
		end_from - start_from);

	if (sequence_type == 1) {
	    start_from = pos - overlap;
	    end_from = pos + line_length + overlap;

	    if (start_from < 1) {
		memmove(sequence, &seq[seq_len - (1-start_from) +2],
			1 - start_from);
	    }
	    if (end_from > seq_len) {
		int end_to = overlap + line_length +
		    seq_len - (pos + line_length - 1);
		memmove(&sequence[end_to], &seq[2], end_from - seq_len - 1);
	    }
	}
	sequence[line_length + 2*overlap] = 0;
    }


    seq_len = strlen(sequence);
#ifdef DEBUG
    printf("seq_len %d max_overlap %d \n", seq_len, max_overlap); 
    printf("s1 %s\n", sequence);
#endif

    /* 
     * find the cut sites for each enzymes on the consensus and store all the
     * results in array match of size total_matches
     */
    if (NULL == (match = (R_Match*)xcalloc(MAXMATCHES, sizeof(R_Match))))
	return -1;

    /* 
     * I have already dealt with circular sequences in the preparation of
     * "sequence" which is only line_length + 2*overlap long
     */
    sequence_type = 0;
    if (-2 == FindMatches(r_enzyme, num_enzymes, sequence, 
			  seq_len, sequence_type, &match, &total_matches)) {
	verror(ERR_WARN, "seqed_write_renzyme", "failed FindMatches\n");
	return -1;
    }

    if (total_matches == 0) {
	xfree(match);
	return 0;
    }

    /* remove overlap offset from cut position */
    for (i = 0; i < total_matches; i++) {
#ifdef DEBUG
	printf("before sort %d %s pos %d\n", 
	       i, r_enzyme[match[i].enz_name].name, match[i].cut_pos);
#endif
	/* match[i].cut_pos = match[i].cut_pos - overlap - 2; */
	match[i].cut_pos = match[i].cut_pos - overlap;
    }

    /* order matches on position, highest first ie right to left */
    qsort((void *) match, total_matches, sizeof(R_Match), compare_rmatch_rev);

#ifdef DEBUG
    for (i = 0; i < total_matches; i++) {
	printf("i %d name %s seq %d pos %d\n", i, r_enzyme[match[i].enz_name].name, match[i].enz_seq, match[i].cut_pos);
    }
#endif

    /* sort matches at the same position into alphabetical order */
    prev_pos = match[0].cut_pos;
    cnt = 0;
    for (i = 0; i < total_matches; i++) {

	if (prev_pos == match[i].cut_pos) {
	    cnt++;
	} else {
	    global_r_enzyme = r_enzyme;
	    qsort((void *) &match[i-cnt], cnt, sizeof(R_Match), 
		  compare_rmatch_name);
	    global_r_enzyme = NULL;
	    cnt = 1;
	}
	prev_pos = match[i].cut_pos;
    }
    if (cnt > 1) {
	global_r_enzyme = r_enzyme;
	qsort((void *) &match[i-cnt], cnt, sizeof(R_Match), 
	      compare_rmatch_name);
	global_r_enzyme = NULL;
    }

#ifdef DEBUG
    for (i = 0; i < total_matches; i++) {
	printf("i %d name %s\n", i, r_enzyme[match[i].enz_name].name);
    }
#endif

    for (i = 0; i < total_matches; i++) {
	c_pos = match[i].cut_pos;
	name_len = strlen(r_enzyme[match[i].enz_name].name);
	line_num = 0;

	/*
	 * check that there is enough room to write out the name both one
	 * space before the cut position and and one space after the name
	 * if not, increment line_num and check again
	 */
	j = -1;
	while (j <= name_len) {

	    /* check that c_pos is within the limits of the line length */
	    if (c_pos+j >= line_length) {
		/*
		if (j > -1 && c_pos+j < line_length+name_overlap) {
		} else {
		    break;
		}
		*/
		if (j == -1 || c_pos+j >= line_length+name_overlap) {
		    break;
		}
	    }
#ifdef DEBUG
	    printf("c_pos %d j %d length %d name_len %d overlap %d ", 
		   c_pos, j, line_length, name_len, overlap);

	    printf("lines[%d][%d] lines $%c$\n", line_num, c_pos-1+j, 
		   lines[line_num][c_pos-1+j]);
#endif
      
	    /* if not enough room on current line, move up to next line */
	    if (c_pos + j <= 0) {
		j++;
	    } else if (c_pos-1+j == -1) {
		j++;
	    } else if ((lines[line_num][c_pos-1+j] != ' ')) {

		line_num++;
		j = -1;

		/* need to malloc more space */
		if (line_num >= *max_renz_lines) {
#ifdef DEBUG
		    printf("ALLOC MORE SPACE line_num %d max %d\n", line_num,
			   max_lines);
#endif
		    if (-1 == seqed_add_more_lines(MAX_DISPLAY_WIDTH,
						   &lines,
						   max_renz_lines)) {
			verror(ERR_WARN, "seqed_write_renzyme", 
			       "unable to allocate more space \n");
			return -1;
		    }
		}
	    } else {
		j++;
	    }
	}

	if (line_num > *num_lines)
	    *num_lines = line_num;

#ifdef DEBUG
	printf("name %s line_num %d\n", r_enzyme[match[i].enz_name].name, 
	       line_num);
#endif
	start = 0;
	/* write in name on line line_num */
	for (j = 0; j < name_len; j++) {
	    /*
	    if (c_pos+j > 0 && c_pos-1+j < line_length) {
		lines[line_num][c_pos-1+j] = 
		    r_enzyme[match[i].enz_name].name[j];
	    }
	    */
#ifdef DEBUG
	    printf("pos %d j %d name %c\n", c_pos, j, 
		   r_enzyme[match[i].enz_name].name[j]);
#endif
	    if (c_pos+j > 0 && c_pos-1+j < line_length) {
		/* print overhanging names when save to file */
		start = 1;
		lines[line_num][c_pos-1+j] = 
		    r_enzyme[match[i].enz_name].name[j];
		/* printf("here1 %d %c \n", line_num, lines[line_num][c_pos-1+j]); */
	    }
	    /* write overhanging names on right hand side */
	    if (c_pos-1+j < line_length+name_overlap && j > 0 && start &&
		c_pos-1+j >= 0) {
		lines[line_num][c_pos-1+j] = 
		    r_enzyme[match[i].enz_name].name[j];
	    }
	    /* 
	     * don't want to write letters from previous line on the beginning
	     * of the next line if name_overlap is specified
	     */
	    if (name_overlap && c_pos-1+j < 0) {
		break;
	    }
	}

	/* write in dots down to base line */
	for (k = line_num - 1; k >= 0; k--) {
	    if (c_pos-1 >= 0) {
		if (lines[k][c_pos-1] == ' ') {
		    lines[k][c_pos-1] = '.';
		}
	    }
	}
    }

    (*num_lines)++;
    *alines = lines;

    xfree(match);
    return 0;
}

int seqedREnzyme(tkSeqed *se,
		 char *filename,
		 char *list,
		 int num_items,
		 int width)
{
    int i, j;
    int str_len, cut;
    int total_len;
    int name_len;

    /*
     * open the file of enzyme names and parse the selected enzymes into the
     * r_enzyme array of size num_enzymes
     */
    open_renz_file(filename, list, num_items, &se->r_enzyme, &se->num_enzymes);

    for (i = 0; i < se->num_enzymes; i++) {
#ifdef DEBUG
	printf("name %s seq %d \n", r_enzyme[i].name, r_enzyme[i].num_seq);
#endif
	for (j = 0; j < se->r_enzyme[i].num_seq; j++) {
#ifdef DEBUG
	    printf("  seq %s site %d\n", r_enzyme[i].seq[j], r_enzyme[i].cut_site[j]);
#endif
	    name_len = strlen(se->r_enzyme[i].name);
	    str_len = strlen(se->r_enzyme[i].seq[j]);
	    cut = se->r_enzyme[i].cut_site[j];
	    if (cut < 0) {
		total_len = abs(cut) + str_len;
	    } else {
		total_len = MAX(cut, str_len);
	    } 

	    total_len = MAX(total_len, name_len);
	    if (max_name_length < name_len)
		max_name_length = name_len;

	    if (max_overlap < total_len)
		max_overlap = total_len;
#ifdef DEBUG
	    printf("  strlen %d total_len %d overlap %d\n", str_len, total_len, max_overlap);
#endif
	}
    }
    if (-1 == (seqed_add_more_lines(MAX_DISPLAY_WIDTH, &sline, &max_lines))) {
	verror(ERR_WARN, "seqedREnzyme", "unable to allocate memory\n");
	return -1;
    }
    return 1;
}

void
seqed_delete_renzyme(tkSeqed *se)
{
    free_lines();
    se->renzDisplayed = 0;
    reset_anchor(se);
    seqed_redisplay_seq(se, se->displayPos, 1);
}
