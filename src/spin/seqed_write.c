#include <string.h>
#include "tkSeqed.h"
#include "genetic_code.h"
#include "dna_utils.h"
#include "seqed_restriction_enzymes.h"
#include "seqed_translate.h"
#include "xalloc.h"
#include "seqed_write.h"
/*
 * all the write functions eg write the function to disk
 */

void seqed_write_sequence(char *sequence,
			  int pos, 
			  int line_length,
			  char *line)
{
    if (line_length <= 0) 
        return;

    /* strncpy(line, &sequence[pos], line_length); */
    strncpy(line, sequence, line_length);
    line[line_length] = '\0';
}


void seqed_write_ruler(int pos,
		       int line_length,
		       char *line)
{
    char *k;
    char *tmp;
    int j, lower, times;

    if (line_length <= 0) 
        return;

    lower = (pos - pos%10);
    times = line_length/10 + 3;
    if (NULL == (tmp = (char *)xmalloc((line_length + 31) * sizeof(char))))
	return;
    
    for (j = 0, k = tmp; j < times; j++, k+=10, lower+=10)
	sprintf(k,"%10d",lower);

    strncpy(line, &tmp[9+pos%10], line_length);
    line[line_length] = '\0';
    xfree(tmp);
}

void seqed_write_complement(char *sequence,
			    int pos, 
			    int line_length,
			    char *line)
{
    if (line_length <= 0) 
        return;

    strncpy(line, sequence, line_length);
    line[line_length] = '\0';
    complement_dna(line, line_length);
}

int seqed_init_write_renzyme(int line_length,
			     char ***lines,
			     int max_renz_lines)
{
    char **alines;
    int i;

    if (NULL == (alines = (char **)xmalloc(max_renz_lines * sizeof(char *))))
	return -1;
    for (i = 0; i < max_renz_lines; i++) {
	if (NULL == (alines[i] = (char *)xmalloc((line_length+1) * 
						sizeof(char))))
	    return -1;
	memset(alines[i], ' ', line_length);
	alines[i][line_length] = '\0';
    }
    *lines = alines;
    return 0;
}


int seqed_write(tkSeqed *se,
		FILE *fp, 
		int from,
		int to,
		int line_length) 
{
    int i, j, k;
    char *line;
    char **lines = NULL;
    int num_renz_lines;
    int width;
    int max_renz_lines = 10;
    int overlap = 0;

    if (line_length <= 0) 
      return -1;

    if (NULL == (line = (char *)xmalloc((line_length+4) * sizeof(char))))
	return -1;
    
    if (se->renzDisplayed) {
	overlap = seqed_get_max_name_length()+1;
	seqed_init_write_renzyme(line_length+overlap, &lines, max_renz_lines);
    }
    for (j = from; j < to; j += line_length) {
	
	if (to - j + 1 < line_length) {
	    width = to - j + 1;
	} else {
	    width = line_length;
	}

	/* restriction enzymes */
	if (se->renzDisplayed) {
	    seqed_write_renzyme(se->sequence, se->sequence_type, se->r_enzyme,
				se->num_enzymes, j, width, 
				overlap, &lines, &max_renz_lines, 
				&num_renz_lines);
	    for (k = num_renz_lines-1; k >= 0; k--) {
#ifdef DEBUG
		printf(" %s\n",lines[k]);
#endif
		fprintf(fp, " %s\n",lines[k]);
	    }
	}

	/* forward translations */
	if (se->translationDisplayed) {
	    for (k = 0; k < se->trans_lines; k++) {
		line[0] = ' ';
		if (se->trans[k] < 4) {
		    seqed_write_translation(&se->sequence[j-1], se->trans[k], 
					    se->trans_mode, j, 
					    width, 1, &line[1]);
#ifdef DEBUG
		    printf("%s\n", line);
#endif
		    fprintf(fp, "%s\n",line);
		}
	    }
	}

	/* HACK leave for now */
	if (se->autoDisplayed) {
	}

	/* forward sequence */
	line[0] = ' ';
	seqed_write_sequence(&se->sequence[j+1], j+1, width, &line[1]);
#ifdef DEBUG
	printf("%s\n", line);
#endif
	fprintf(fp, "%s\n",line);

	/* position numbers */
	if (se->rulerDisplayed) {
	    seqed_write_ruler(j, width, &line[1]);
#ifdef DEBUG
	    printf("%s\n", line);
#endif
	    fprintf(fp, "%s\n",line);
	}

	/* complemented sequence */
	if (se->complementDisplayed) {
	    seqed_write_complement(&se->sequence[j+1], j+1, width, &line[1]);
#ifdef DEBUG
	    printf("%s\n", line);
#endif
	    fprintf(fp, "%s\n",line);
	}

	/* complemented translations */
	if (se->translationDisplayed) {
	    for (k = 0; k < se->trans_lines; k++) {
		line[0] = ' ';
		if (se->trans[k] > 3) {
		    seqed_write_translation(&se->sequence[j-1], se->trans[k], 
					    se->trans_mode, j, 
					    width, 1, &line[1]);
#ifdef DEBUG
		    printf("%s\n", line);
#endif
		    fprintf(fp, "%s\n",line);
		}
	    }
	}
#ifdef DEBUG
	printf("\n");
#endif
	fprintf(fp, "\n");
    }
    if (se->renzDisplayed) {
	for (i = 0; i < max_renz_lines; i++) {
	    xfree(lines[i]);
	}
	xfree(lines);
    }
    xfree(line);
    return 0;
}
