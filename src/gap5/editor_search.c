/*
 * editor_search.c:
 * Contains functions for the Contig Editor "search" button.
 */

#include <tg_gio.h>

#include "dna_utils.h"
#include "search_utils.h"

#include "editor_view.h"
#include "consensus.h"

#define WIN_WIDTH 65536

int edview_search_position(edview *xx, int dir, int strand, char *value) { 
    puts("position search");
    return 0;
}

int edview_search_uposition(edview *xx, int dir, int strand, char *value) {
    puts("uposition search");
    return 0;
}

int edview_search_sequence(edview *xx, int dir, int strand, char *value) {
    int mismatches = 0; /* exact match */
    int where = 2;      /* consensus */
    char *p;
    int start, end;
    char cons[WIN_WIDTH+1];
    int patlen;
    char *uppert, *upperb;
    int found = 0, at_end = 0;
    int fseq, fpos, i, j;

    /*
     * Parse value search string. It optionally includes two extra params
     * separated by #. Ie:
     *     <string>#<N.mismatches>#<where>.
     * <where> is 1 for readings, 2 for consensus, 3 for both.
     */
    if (p = strchr(value, '#')) {
	mismatches = atoi(p+1);
	*p = 0;
	if (p = strchr(p+1, '#'))
	    where = atoi(p+1);
    }


    /* uppercase search string, remove pads, and store fwd/rev copies */
    patlen = strlen(value);
    depad_seq(value, &patlen, NULL);
    if (NULL == (uppert = (char *)xmalloc(patlen + 1)))
	return 0;
    if (NULL == (upperb = (char *)xmalloc(patlen + 1)))
	return 0;

    uppert[patlen] = upperb[patlen] = 0;
    for (i = patlen-1; i >= 0; i--) {
	upperb[i] = uppert[i] = toupper(value[i]);
    }
    complement_seq(upperb, patlen);


    /* Loop */
    if (dir) {
	start = xx->cursor_apos + (dir ? 1 : -1);
	end   = start + (WIN_WIDTH-1);
    } else {
	end   = xx->cursor_apos + (dir ? 1 : -1);
	start = end - (WIN_WIDTH-1);
    }
    fpos = xx->cursor_apos;

    do {
	char *ind, *indt = NULL, *indb = NULL;

	calculate_consensus_simple(xx->io, xx->cnum, start, end, cons, NULL);
	cons[WIN_WIDTH] = 0;

	if (dir) {
	    if (strand == '+' || strand == '=')
		indt = pstrstr_inexact(cons, uppert, mismatches, NULL);
	    if (strand == '-' || strand == '=')
		indb = pstrstr_inexact(cons, upperb, mismatches, NULL);
	} else {
	    if (strand == '+' || strand == '=')
		indt = prstrstr_inexact(cons, uppert, mismatches, NULL);
	    if (strand == '-' || strand == '=')
		indb = prstrstr_inexact(cons, upperb, mismatches, NULL);
	}

	if (indt && indb)
	    ind = MIN(indt, indb);
	else if (indt)
	    ind = indt;
	else if (indb)
	    ind = indb;
	else
	    ind = NULL;

	if (ind != NULL) {
	    found = 1;
	    if (dir) {
		if (fpos <= start + ind-cons) {
		    fpos = start + ind-cons;
		    fseq = xx->contig->rec;
		}
	    } else {
		if (fpos >= start + ind-cons) {
		    fpos = start + ind-cons;
		    fseq = xx->contig->rec;
		}
	    }
	    break;
	}

	/* Next search region - overlapping by patlen+pads */
	if (dir) {
	    for (i = WIN_WIDTH-1, j = patlen; j && i; i--) {
		if (cons[i] != '*')
		    j--;
	    }
	    if (i == 0)
		break;
	    start += i;
	    end   += i;

	    if (start > xx->contig->end)
		at_end = 1;
	} else {
	    for (i = 0, j = patlen && i < WIN_WIDTH; j; i++) {
		if (cons[i] != '*')
		    j--;
	    }
	    if (i == WIN_WIDTH)
		break;

	    start -= WIN_WIDTH-i;
	    end   -= WIN_WIDTH-i;

	    if (end < xx->contig->start)
		at_end = 1;
	}
    } while (!at_end);

    if (found) {
	edSetCursorPos(xx, fseq == xx->contig->rec ? GT_Contig : GT_Seq,
		       fseq, fpos, 1);
    }

    free(uppert);
    free(upperb);

    return found ? 0 : -1;
}

int edview_search_consquality(edview *xx, int dir, int strand, char *value) {
    int start, end;
    float qual[WIN_WIDTH+1];
    int found = 0, at_end = 0;
    int fpos, i, qval = atoi(value);

    /* Set initial start positions */
    if (dir) {
	start = xx->cursor_apos + (dir ? 1 : -1);
	end   = start + (WIN_WIDTH-1);
    } else {
	end   = xx->cursor_apos + (dir ? 1 : -1);
	start = end - (WIN_WIDTH-1);
    }
    fpos = xx->cursor_apos;

    /* Loop WIN_WIDTH block at a time */
    do {
	calculate_consensus_simple(xx->io, xx->cnum, start, end, NULL, qual);

	if (dir) {
	    for (i = 0; i < WIN_WIDTH; i++) {
		if (qual[i] < qval) {
		    found = 1;
		    break;
		}
	    }
	} else {
	    for (i = WIN_WIDTH-1; i; i--) {
		if (qual[i] < qval) {
		    found = 1;
		    break;
		}
	    }
	}

	if (found) {
	    fpos = start + i;
	    break;
	}

	/* Next search region - overlapping by patlen+pads */
	if (dir) {
	    start += WIN_WIDTH;
	    end   += WIN_WIDTH;

	    if (start > xx->contig->end)
		at_end = 1;
	} else {
	    start -= WIN_WIDTH;
	    end   -= WIN_WIDTH;

	    if (end < xx->contig->start)
		at_end = 1;
	}
    } while (!at_end);

    if (found) {
	edSetCursorPos(xx, GT_Contig, xx->contig->rec, fpos, 1);
	return 0;
    }

    return -1;
}

int edview_search_name(edview *xx, int dir, int strand, char *value)
{
    int rec;
    
    rec = sequence_index_query(xx->io, value);
    if (rec <= 0)
	return -1;

    edSetCursorPos(xx, GT_Seq, rec, 0, 1);
    return 0;
}

int edview_search(edview *xx, int dir, int strand,
		  char *type, char *value) {
    struct {
	char *name;
	int (*func)(edview *xx, int dir, int strand, char *value);
    } types[] = {
	{"position",     edview_search_position},
	{"uposition",    edview_search_uposition},
	{"sequence",     edview_search_sequence},
	{"consquality",  edview_search_consquality},
	{"name",         edview_search_name},
    };
    int i;

    for (i = 0; i < sizeof(types)/sizeof(*types); i++) {
	if (0 == strcmp(types[i].name, type))
	    return types[i].func(xx, dir, strand, value);
    }

    return -1;
}
