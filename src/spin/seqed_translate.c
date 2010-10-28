#include <string.h>
#include "tkSeqed.h"
#include "tkSeqedUtils.h"
#include "genetic_code.h"
#include "seqed_translate.h"

void seqed_write_translation(char *sequence,
			     int frame,
			     int size,
			     int pos,
			     int line_length,
			     int overlap,
			     char *line)
{
    int j, complement, offset;

    complement = (frame > 3);
    offset = (3+frame-pos%3)%3;

    if (size == 3) {
	/* 3 letter code */

	/*char *(*fptr)() = complement ? codon_to_cacid3 : codon_to_acid3;*/
	char *(*fptr)(char *sequence) = complement ? codon_to_cacid3 : codon_to_acid3;

	switch(offset) {
	    char *t;

	case 1:
	    t = fptr(sequence);
	    if (overlap) {
		line[0] = ' ';
	    } else {
		line[0] = t[2];
	    }
	    break;
	    
	case 2:
	    t = fptr(sequence+1);
	    if (overlap) {
		line[-1] = t[0];
	    }
	    line[0] = t[1];
	    line[1] = t[2];
	    break;
	}

	for (j = offset; j < line_length; j+=3) {
	    char *t;
	    t = fptr(&sequence[j+2]);
	    if (j < line_length - overlap) {
		line[j+0] = t[0];
		line[j+1] = t[1];
		line[j+2] = t[2];
	    } else {
		line[j] = '\0';
	    }
	}
	line[j] = '\0';
    } else {
	/* 1 letter code */
	/*char (*fptr)() = complement ? codon_to_cacid1 : codon_to_acid1;*/
	char (*fptr)(char *sequence) = complement ? codon_to_cacid1 : codon_to_acid1;
	memset(line, ' ', line_length);

	if (offset == 2)
	    if ('*' == (line[0] = fptr(&sequence[1])))
		;
	for (j = offset; j < line_length-1; j+=3) {
	    if ('*' == (line[j+1] = fptr(&sequence[j+2]))) {
		;
	    }
	}
	line[line_length] = '\0';
    }
}

void
seqed_translate_frame(tkSeqed *se,
		      char *sequence,
		      int pos,
		      int width,
		      int frame,
		      char *sline,
		      char *name,
		      int size,
		      XawSheetInk *splodge)
{
    int j, complement, offset;

    complement = (frame > 3);

    offset = (3+frame-pos%3)%3;
    
    for (j = 0; j < width; j++)
	splodge[j].sh = sh_default;

    for (j = 0; j < width; j++)
	sline[j] = '.';
	
    if (size == 3) {
	/* 3 letter code */

      /*char *(*fptr)() = complement ? codon_to_cacid3 : codon_to_acid3;*/
	char *(*fptr)(char *sequence) = complement ? codon_to_cacid3 : codon_to_acid3;
	switch(offset) {
	    char *t;

	case 1:
	    t = fptr(sequence);
	    if ('*' == (sline[0] = t[2])) {
#ifdef COLOUR
		splodge[0].fg = se->stop_colour;
		splodge[0].sh |= sh_fg;
#endif
	    }
	    break;
	    
	case 2:
	    t = fptr(sequence+1);
#ifdef COLOUR
	    if ('*' == t[1]) {
		splodge[0].fg = se->stop_colour;
		splodge[1].fg = se->stop_colour;
		splodge[0].sh |= sh_fg;
		splodge[1].sh |= sh_fg;
	    }
#endif

#ifdef DEBUG
	    printf("seq %s \n", sequence+1);
#endif
	    sline[0] = t[1];
	    sline[1] = t[2];
	    break;
	}

	for (j = offset; j < width; j+=3) {
	    char *t;
	    t = fptr(&sequence[j+2]);

#ifdef COLOUR
	    if ('*' == t[0]) {
		splodge[j+0].fg = xx->stop_colour;
		splodge[j+1].fg = xx->stop_colour;
		splodge[j+2].fg = xx->stop_colour;
		splodge[j+0].sh |= sh_fg;
		splodge[j+1].sh |= sh_fg;
		splodge[j+2].sh |= sh_fg;
	    }
#endif
            sline[j+0] = t[0];
            sline[j+1] = t[1];
            sline[j+2] = t[2];

	}
    } else {
	/* 1 letter code */
      /*char (*fptr)() = complement ? codon_to_cacid1 : codon_to_acid1;*/
	char (*fptr)(char *sequence) = complement ? codon_to_cacid1 : codon_to_acid1;

	memset(sline, ' ', width);

	if (offset == 2)
	    if ('*' == (sline[0] = fptr(&sequence[1])))
		;
#ifdef COLOUR
		splodge[0].fg = xx->stop_colour;
#endif
	for (j = offset; j < width-1; j+=3) {
	    if ('*' == (sline[j+1] = fptr(&sequence[j+2]))) {
		;
#ifdef COLOUR
		splodge[j+1].fg = xx->stop_colour;
		splodge[j+1].sh |= sh_fg;
#endif
	    }
	}
    }
    sprintf(name, "Frame %d%c", (frame-1)%3+1, frame>3 ? '-' : '+');

}


/*
 * determine the first codon of an exon using num_char bases from the
 * previous exon
 */
void
first_codon(tkSeqed *se,
	    char *sequence,
	    int num_char,
	    char *codon,
	    region *exon,
	    int exon_num,
	    XawSheetInk *splodge,
	    int seq_pos)
{
    int k, l;
    region prev_exon;

    prev_exon =  exon[exon[exon_num].join];

    /*
     * find num_char bases from the previous exon and put in codon[]
     * NB sequence starts at 0 but is offset by 2 ie the first base is at 
     * position sequence[2] 
     */
    for (k = 0; k < num_char; k++) {
	/* codon[k] = se->sequence[exon[exon_num-1].end + k - num_char + 2]; */
	codon[k] = se->sequence[prev_exon.end + k - num_char + 2];
#ifdef DEBUG
	printf("SPLODGE pos %d k %d num_char %d\n", seq_pos, k, num_char);
#endif
	if (seq_pos+k-num_char >= 0) 
	    splodge[seq_pos+k-num_char].sh |= sh_light;
#ifdef DEBUG
	printf("CODON %d %d %c \n", k, prev_exon.end + 
				k - num_char + 2, codon[k]);
#endif
    }
	
    /* add on the appropriate number of bases from start of current exon */
    for (k = num_char, l = 0; k < 3; k++) {
#ifdef DEBUG
	printf("seq %s\n", sequence); 
	printf("pos %d k %d codon %c\n", seq_pos + l +2, k, sequence[seq_pos + 2 + l]);
#endif
	codon[k] = sequence[seq_pos + 2 + l++];
    }
#ifdef DEBUG
    printf("exon_num %d end %d num_char %d codon %s\n",
	   exon_num, prev_exon.end, num_char, codon);		
#endif
}
	

/*
 * determine the 1st (and 2nd) characters of the automatic translation line
 * depending on offset
 * the four possible conditions:
 * 1)    g | C A  offset = 2 seq_pos = 0  num_char = 1
 * 2)    g C | A  offset = 1 seq_pos = -1 num_char = 1
 * 3)    g | c A  offset = 2 seq_pos = 0  num_char = 2
 * 4)    g c | A  offset = 1 seq_pos = 0  num_char = 2
 */
void
find_line_start3(tkSeqed *se,
		 char *sequence,
		 int pos,
		 int offset,
		 int start,
		 int end,
		 int num_char,
		 int line_num,
		 region *exons,
		 int exon_num,
		 XawSheetInk *splodge,
		 char *(*fptr)(char *sequence),
		 char *sline) 
{
    char *t;
    char codon[3];
    char st_str[3];
    int seq_start = 0;

    st_str[0] = ' ';
    st_str[1] = ' ';

#ifdef DEBUG
    printf("offset %d pos %d start %d\n", offset, pos, start);
#endif
    /* 
     * decide if in translation region (which may be offset by num_char's 
     * carried over from the previous exon)
     */
    if ((pos >= start - num_char) && (pos <= end)) {
	
	/* 
	 * decide if at the start of an exon (not the first), bearing in
	 * mind that you could have up to 2 bases carried over from the
	 * previous exon 
	 */
	if ((exon_num > 0) && (pos - start <= 2 - num_char) && 
	    (exons[exon_num].join > -1)) {
	    
	    /*
	     * HACK to deal with the case when the start of the sequence of
	     * the current exon is one to the left of the display edge
	     * ie g C | A ----
	     * where g is from the previous exon and C is the start of the
	     * current exon and | is the edge of the display
	     */  
	    if (offset == 1 && num_char == 1)
		seq_start = -1;
	    else
		seq_start = 0;

	    first_codon(se, sequence, num_char, codon, exons, 
			exon_num, splodge, seq_start);
	    t = fptr(codon);
#ifdef DEBUG
	    printf("t1 %s codon %s\n", t, codon);
#endif
	    /*
	     * another HACK: to highlight the second base when it is 
	     * displayed without the first
	     */
	    if (offset == 2 && num_char == 2) {
		splodge[0].sh |= sh_light;
	    }
	    strcpy(st_str, t);
	} else {
	    if (offset == 1) 
		t = fptr(sequence);
	    else
		t = fptr(&sequence[1]);
	    strcpy(st_str, t);
#ifdef DEBUG
	    printf("t2 %s st_str %s \n", t, st_str);
#endif
	}
    }

    /* create beginning of display line */
    if (offset == 1) {
	sline[0] = st_str[2];
	splodge[0].fg = exons[exon_num].colour;
	splodge[0].sh |= sh_fg;
    } else if (offset == 2) {
	sline[0] = st_str[1];
	splodge[0].fg = exons[exon_num].colour;
	splodge[0].sh |= sh_fg;
	if (pos + 1 <= end) {
	    sline[1] = st_str[2];
	    splodge[1].fg = exons[exon_num].colour;
	    splodge[1].sh |= sh_fg;
	}
    }
}


void
find_line_start1(tkSeqed *se,
		 char *sequence,
		 int pos,
		 int offset,
		 int start,
		 int end,
		 int num_char,
		 region *exon,
		 int exon_num,
		 XawSheetInk *splodge,
		 char (*fptr)(char *sequence),
		 char *sline) 
{
    char t;
    char codon[3];
    int seq_start = 0;

    t = ' ';
#ifdef DEBUG
    printf("offset %d pos %d start %d end %d\n", offset, pos, start, end);
#endif
    /* 
     * decide if in translation region (which may be offset by num_char's 
     * carried over from the previous exon)
     */
    if ((pos >= start - num_char) && (pos <= end)) {
	
	/* 
	 * decide if at the start of an exon (not the first), bearing in
	 * mind that you could have up to 2 bases carried over from the
	 * previous exon 
	 */
	if ((exon_num > 0) && (pos - start <= 2 - num_char)) {
	    
	    /*
	     * HACK to deal with the case when the start of the sequence of
	     * the current exon is one to the left of the display edge
	     * ie g C | A ----
	     * where g is from the previous exon and C is the start of the
	     * current exon and | is the edge of the display
	     */  
	    if (offset == 1 && num_char == 1)
		seq_start = -1;
	    else
		seq_start = 0;

	    first_codon(se, sequence, num_char, codon, exon, 
			exon_num, splodge, seq_start);
	    t = fptr(codon);
#ifdef DEBUG
	    printf("t1 %c codon %s\n", t, codon);
#endif
	} else {
	    t = fptr(&sequence[1]);
#ifdef DEBUG
	    printf("t2 %c st_str %s \n", t, st_str);
#endif
	}
    }

    /* create beginning of display line */
    sline[0] = t;
}

/*
 * not very efficient because most of the time the exon hasn't
 * changed
 */
int
check_in_exon(region *exon,
	      int num_exon,
	      int pos,
	      int *exon_num,
	      int *start,
	      int *end)
{
    int i;
    
    for (i = 0; i < num_exon; i++) {
	if ((pos >= exon[i].start) && (pos <= exon[i].end)) {
	   
	    /*
	       printf("pos %d start %d end %d i %d\n", pos,
	       exon[i].start, exon[i].end, i);
	       */
	    *exon_num = i;
	    *start = exon[i].start;
	    *end = exon[i].end;
	    return 1;
	}
   }
    return 0;
}

int
check_overlap(region exon1,
	      region exon2)
{
    if (exon1.end >= exon2.start) {
#ifdef DEBUG
	printf("FOUND overlap %d %d\n", exon1.end, exon2.start);
#endif
	return 1;
    } else {
	return 0;
    }
}

int
find_auto_lines(region **exons,
		int num_exons,
		int width,
		int pos,
		int complement)
{	 
    int i, k;
    int num_lines = 0;
    int max_lines = 0;

#ifdef DEBUG
    printf("           find_auto_lines\n");
#endif
    for (i = 0; i < num_exons; i++) {
	for (k = 0; k < width; k++) {

	    /* check within display region */
	    if ((pos+k >= (*exons)[i].start) && (pos+k <= (*exons)[i].end) &&
		((*exons)[i].complement == complement)) {
		
		/* printf("FOUND exon %d\n", i); */
		/* if overlapping with previous exon */
		if (i > 0 &&  (*exons)[i-1].end >= (*exons)[i].start) {

		    /* check within display region */
		    if ((pos+k >= (*exons)[i-1].start) && 
			(pos+k <= (*exons)[i-1].end) &&
			((*exons)[i-1].complement == complement)) {

			/* printf("OVERLAPPING exon %d %d\n", i, i-1); */

			(*exons)[i].line_num = num_lines;
			num_lines++;
		    } else {
			(*exons)[i].line_num = 0;
			num_lines = 1;
		    }
		} else {
		    (*exons)[i].line_num = 0;
		    num_lines = 1;
		}
		if (max_lines < num_lines)
		    max_lines = num_lines;
		break;
	    }
	}
    }

    /* printf("MAXLINES %d\n", max_lines); */
    return max_lines;
}

void
seqed_auto_translate(tkSeqed *se,
		     char *sequence,
		     int pos,
		     int width,
		     char *sline,
		     char *name,
		     XawSheetInk *splodge,
		     int size,
		     region *exons,
		     int exon_num,
		     int start,
		     int end,
		     int line_num,
		     int complement)
{
    int i, j, k, offset;
    char codon[3];
    static int frame = 1;
    int num_char;
    region prev_exon;
    
#ifdef DEBUG
    printf("                 AUTO TRANSLATE pos %d exon %d\n", pos, exon_num);
#endif
    for (j = 0; j < width; j++) 
	splodge[j].sh = sh_default;

    for (j = 0; j < width; j++)
	sline[j] = ' ';

    if (complement != exons[exon_num].complement)
	return;

    if (exon_num > 0 && exons[exon_num].join > -1) {
	prev_exon = exons[exons[exon_num].join];
	frame = (exons[exon_num].start - prev_exon.num_char) %3;
	num_char = prev_exon.num_char;

    } else {
	frame = exons[exon_num].start % 3;
	num_char = 0;
    }
    offset = (3+frame-pos%3)%3;

    if (size == 3) {
	char *(*fptr)(char *sequence) = complement ? codon_to_cacid3 : codon_to_acid3;

	if (offset < 3) 
	    find_line_start3(se, sequence, pos, offset, start, end, num_char,
			     line_num, exons, exon_num, splodge, fptr, 
			     sline);

	/*
	 * loop through the visible sequence in codons
	 */
	for (j = offset; j < width; j+=3) {
	    char *t;
	    t = fptr(&sequence[j+2]);

	    for (i = 0; i < 3; i++) {

		/* check if within translation region */
		if ((j+pos+i >= start) && (j+pos+i <= end)) {

		    splodge[j+i].fg = exons[exon_num].colour;
		    splodge[j+i].sh |= sh_fg;

		    if ((exon_num > 0) && (j+pos+i == start)) {

			first_codon(se, sequence, num_char, codon, exons, 
				    exon_num, splodge, i+j);
			t = fptr(codon);
			for (k = 0; k < 3; k++) {
			    sline[j+k+i-num_char] = t[k];
			}
		    } else {
			sline[j+i] = t[i];
			/* printf("sline2 %d %s\n", line_num, &(sline[line_num][3])); */
		    }
		} else {
		    sline[j+i] = ' ';
		}
	    }
	}
			    
    } else {
	/* 1 letter code */
	char (*fptr)(char *sequence) = complement ? codon_to_cacid1 : codon_to_acid1;

	if (offset == 2) 
	    find_line_start1(se, sequence, pos, offset, start, end, num_char,
			     exons, exon_num, splodge, fptr, sline);

	for (j = offset; j < width-1; j+=3) {
	    char t;
	    t = fptr(&sequence[j+2]);

	    for (i = 0; i < 3; i++) {

		/* check if within translation region */
		if ((j+pos+i >= start) && (j+pos+i <= end)) {
		    
		    if ((exon_num > 0) && (j+pos+i == start)) {

			first_codon(se, sequence, num_char, codon, exons, 
				    exon_num, splodge, i+j);
			t = fptr(codon);
#ifdef DEBUG
			printf("t %c codon %s pos %d\n", t, codon, j-num_char+1);
#endif
			sline[j+1+i-num_char] = t;
			break;
		    } else {
			sline[j+1] = t;
			/* printf("sline2 %s\n", &(sline[3])); */
		    }
		} else {
		    sline[j+1] = ' ';
		}
	    }
	}
	
    }
    /* sprintf(name, "Frame %d%c", (frame-1)%3+1, frame>3 ? '-' : '+'); */
}
