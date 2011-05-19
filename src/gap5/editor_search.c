/*
 * editor_search.c:
 * Contains functions for the Contig Editor "search" button.
 */

#include <tg_gio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "dna_utils.h"
#include "search_utils.h"

#include "editor_view.h"
#include "consensus.h"
#include "reg_exp.h"

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
    tg_rec fseq;
    int fpos, i, j;

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
	    if (dir) {
		if (fpos <= start + ind-cons) {
		    found = 1;
		    fpos = start + ind-cons;
		    fseq = xx->contig->rec;
		}
	    } else {
		if (fpos >= start + ind-cons) {
		    found = 1;
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

int edview_search_cons_het(edview *xx, int dir, int strand, char *value) {
    int start, end;
    int found = 0, at_end = 0;
    int fpos, i, qval = atoi(value);
    consensus_t cons[WIN_WIDTH+1];

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
	calculate_consensus(xx->io, xx->cnum, start, end, cons);

	if (dir) {
	    for (i = 0; i < WIN_WIDTH; i++) {
		if (cons[i].scores[6] >= qval) {
		    found = 1;
		    break;
		}
	    }
	} else {
	    for (i = WIN_WIDTH-1; i; i--) {
		if (cons[i].scores[6] >= qval) {
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

int edview_search_cons_discrep(edview *xx, int dir, int strand, char *value) {
    int start, end;
    int found = 0, at_end = 0;
    int fpos, i;
    double qval = atof(value);
    consensus_t cons[WIN_WIDTH+1];

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
	calculate_consensus(xx->io, xx->cnum, start, end, cons);

	if (dir) {
	    for (i = 0; i < WIN_WIDTH; i++) {
		if (cons[i].discrep >= qval) {
		    found = 1;
		    break;
		}
	    }
	} else {
	    for (i = WIN_WIDTH-1; i; i--) {
		if (cons[i].discrep >= qval) {
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
    tg_rec rec, *rp, cnum = -1, best_rec;
    int best_pos;
    int nr, i;
    rangec_t *(*ifunc)(GapIO *io, contig_iterator *ci);
    int start, end;
    contig_iterator *iter;

    /* Check for #num where num is a sequence record in this contig */
    if (*value == '#') {
	char *endp;
	int64_t v = strtol64(value+1, &endp, 10);
	rec = v;

	if (*endp == '\0' && cache_exists(xx->io, GT_Seq, rec)) {
	    sequence_get_position(xx->io, rec, &cnum, NULL, NULL, NULL);
	    if (cnum == xx->cnum) {
		edSetCursorPos(xx, GT_Seq, rec, 0, 1);
		return 0;
	    }
	}
    }

    /* Find all hits matching this name */
    rp = sequence_index_query_all(xx->io, value, 1, &nr);

    /* Also get an position-based iterator */
    if (dir) {
	start    = xx->cursor_apos + 1;
	end      = xx->contig->end;
	ifunc    = contig_iter_next;
	best_pos = end + 1;
    } else {
	start    = xx->contig->start;
	end      = xx->cursor_apos - 1;
	ifunc    = contig_iter_prev;
	best_pos = start - 1;
    }

    iter = contig_iter_new_by_type(xx->io, xx->cnum, 1,
				   dir == 1 ? CITER_FIRST : CITER_LAST,
				   start-1, end+1, GRANGE_FLAG_ISSEQ);
    if (!iter)
	return -1;

    /*
     * The iterator also finds overlapping objects, not just ones beyond this
     * point. That's fine if we're on the consensus as we probably want to
     * jump to the first seq-name overlapping this point.
     *
     * However if we're on a sequence already, we want the first one
     * after or before that sequence. So we skip along iterator until we're
     * at the current record.
     */
    if (xx->cursor_type == GT_Seq) {
	rangec_t *r;
	while ((r = ifunc(xx->io, iter))) {
	    if (r->rec == xx->cursor_rec)
		break;
	}
    }
				   

    /* Alternate between the by-name and by-position scan */
    best_rec = -1;
    for (i = 0; i < nr; i++) {
	int start, end;
	rangec_t *r;

	/* From name index */
	rec = rp[i++];
	sequence_get_position(xx->io, rec, &cnum, &start, &end, NULL);
	if (cnum == xx->cnum) {
	    if ((dir  && best_pos > start && start > xx->cursor_apos) ||
		(!dir && best_pos < start && start < xx->cursor_apos)) {
		best_pos = start;
		best_rec = rec;
	    }
	}

	/* From iterator */
	if ((r = ifunc(xx->io, iter))) {
	    seq_t *s;
	    if (NULL == (s = cache_search(xx->io, GT_Seq, r->rec))) {
		/* No match */
		best_rec = -1;
		break;
	    }
	    
	    if (strncmp(s->name, value, strlen(value)) == 0) {
		/* prefix match */
		puts("Found by pos iterator");
		best_rec = r->rec;
		break;
	    }
	} else {
	    /* End of contig - bail out early */
	    best_rec = -1;
	    break;
	}
    }

    contig_iter_del(iter);
    if (rp)
	free(rp);
    
    if (best_rec != -1) {
	edSetCursorPos(xx, GT_Seq, best_rec, 0, 1);
	return 0;
    }

    return -1;
}

int edview_search_tag_type(edview *xx, int dir, int strand, char *value) {
    contig_iterator *iter;
    int start, end;
    rangec_t *r;
    rangec_t *(*ifunc)(GapIO *io, contig_iterator *ci);
    int type = str2type(value);

    if (dir) {
	start = xx->cursor_apos + (dir ? 1 : -1);
	end   = xx->contig->end;
	ifunc = contig_iter_next;
    } else {
	start = xx->contig->start;
	end   = xx->cursor_apos + (dir ? 1 : -1);
	ifunc = contig_iter_prev;
    }

    iter = contig_iter_new_by_type(xx->io, xx->cnum, 1,
				   dir == 1 ? CITER_FIRST : CITER_LAST,
				   start, end, GRANGE_FLAG_ISANNO);
    if (!iter)
	/* Can happen legitimately when we're already at the end of contig */
	return -1;

    while (r = ifunc(xx->io, iter)) {
	if ((dir  && r->start < start) ||
	    (!dir && r->start > end))
	    continue;

	if (r->mqual == type) 
	    break;
    }

    if (r) {
	if (r->flags & GRANGE_FLAG_TAG_SEQ) {
	    int pos;
	    sequence_get_position(xx->io, r->pair_rec, NULL, &pos, NULL, NULL);
	    pos = r->start - pos;
	    edSetCursorPos(xx, GT_Seq, r->pair_rec, pos, 1);
	} else {
	    edSetCursorPos(xx, GT_Contig, xx->cnum, r->start, 1);
	}
	contig_iter_del(iter);
	return 0;
    }

    contig_iter_del(iter);
    return -1;
}

int edview_search_tag_anno(edview *xx, int dir, int strand, char *value) {
    contig_iterator *iter;
    int start, end;
    rangec_t *r;
    rangec_t *(*ifunc)(GapIO *io, contig_iterator *ci);
    char *r_exp = NULL;

    if (value) {
	if (NULL == (r_exp = REGCMP(xx->interp, value))) {
	    verror(ERR_WARN, "Search by anno", "invalid regular expression");
	    return -1;
	}
    }

    if (dir) {
	start = xx->cursor_apos + (dir ? 1 : -1);
	end   = xx->contig->end;
	ifunc = contig_iter_next;
    } else {
	start = xx->contig->start;
	end   = xx->cursor_apos + (dir ? 1 : -1);
	ifunc = contig_iter_prev;
    }

    iter = contig_iter_new_by_type(xx->io, xx->cnum, 1,
				   dir == 1 ? CITER_FIRST : CITER_LAST,
				   start, end, GRANGE_FLAG_ISANNO);
    if (!iter)
	return -1;

    while (r = ifunc(xx->io, iter)) {
	anno_ele_t *ae;

	if ((dir  && r->start < start) ||
	    (!dir && r->start > end))
	    continue;

	if (!r_exp)
	    break; /* blank expr => match all */

	ae = cache_search(xx->io, GT_AnnoEle, r->rec);
	if (!ae->comment)
	    continue;

	if (REGEX(xx->interp, ae->comment, r_exp))
	    break;
    }

    REGFREE(xx->interp, r_exp);

    if (r) {
	if (r->flags & GRANGE_FLAG_TAG_SEQ) {
	    int pos;
	    sequence_get_position(xx->io, r->pair_rec, NULL, &pos, NULL, NULL);
	    pos = r->start - pos;
	    edSetCursorPos(xx, GT_Seq, r->pair_rec, pos, 1);
	} else {
	    edSetCursorPos(xx, GT_Contig, xx->cnum, r->start, 1);
	}
	contig_iter_del(iter);
	return 0;
    }

    contig_iter_del(iter);
    return -1;
}

int edview_search_tag_indel(edview *xx, int dir, int strand, char *value) {
    contig_iterator *iter;
    int start, end;
    rangec_t *r;
    rangec_t *(*ifunc)(GapIO *io, contig_iterator *ci);

    if (dir) {
	start = xx->cursor_apos + (dir ? 1 : -1);
	end   = xx->contig->end;
	ifunc = contig_iter_next;
    } else {
	start = xx->contig->start;
	end   = xx->cursor_apos + (dir ? 1 : -1);
	ifunc = contig_iter_prev;
    }

    iter = contig_iter_new_by_type(xx->io, xx->cnum, 1,
				   dir == 1 ? CITER_FIRST : CITER_LAST,
				   start, end, GRANGE_FLAG_ISREFPOS);
    if (!iter)
	return -1;

    while (r = ifunc(xx->io, iter)) {
	anno_ele_t *ae;

	if ((dir  && r->start < start) ||
	    (!dir && r->start > end))
	    continue;

	break;
    }

    if (r) {
	edSetCursorPos(xx, GT_Contig, xx->cnum, r->start, 1);
	contig_iter_del(iter);
	return 0;
    }

    contig_iter_del(iter);
    return -1;
}

/*
 * Performs a search within the editor.
 * type   is a string indicating the search type - name, tag, sequence, ...
 * dir    is 1 for forward, 0 for reverse.
 * strand is '+', '-' or '=' but only applicable for some search types.
 *
 * Returns 0 on success (and moves the editor cursor)
 *        -1 on failure / not found.
 */
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
	{"tag",          edview_search_tag_type},
	{"annotation",   edview_search_tag_anno},
	{"indel",        edview_search_tag_indel},
	{"conshet",      edview_search_cons_het},
	{"consdiscrep",  edview_search_cons_discrep},
    };
    int i;

    for (i = 0; i < sizeof(types)/sizeof(*types); i++) {
	if (0 == strcmp(types[i].name, type))
	    return types[i].func(xx, dir, strand, value);
    }

    fprintf(stderr, "Unrecognised search type '%s'\n", type);

    return -1;
}
