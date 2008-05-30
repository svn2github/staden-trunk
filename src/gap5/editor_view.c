#include <xalloc.h>
#include <math.h>
#include <ctype.h>

#include "editor_view.h"
#include "tkSheet.h"
#include "tman_interface.h"
#include "gap_globals.h"
#include "qualIO.h"
#include "qual.h"
#include "tg_gio.h"
#include "misc.h"
#include "consensus.h"

/*
 * Allocates and initialises a new edview
 */
edview *edview_new(GapIO *io, int contig,
		   Editor *ed, edNames *names,
		   void (*dispFunc)(void *, int, int, int, void *))
{
    edview *xx;
    static int editor_id = 1;

    xx = (edview *)xcalloc(1, sizeof(*xx));
    if (!xx)
	return NULL;
    
    xx->editor_id = editor_id++;

    xx->interp = NULL; /* filled out elsewhere */
    xx->io = io; /* model */
    xx->cnum = contig;
    xx->contig = (contig_t *)cache_search(io, GT_Contig, xx->cnum);
    cache_incr(xx->io, xx->contig);

    xx->ed = ed;
    xx->displayPos = 1;
    xx->displayYPos = 0;
    xx->displayWidth = xx->ed->sw.columns;
    xx->displayHeight = xx->ed->sw.rows;
    xx->dispFunc = dispFunc;
    xx->editorState = StateUp;

    xx->y_cons = 0;
    xx->y_numbers = 1;
    xx->y_seq_start = 2;
    xx->y_seq_end = 0;

    xx->names = names;
    xx->names_xPos = 0;

    xx->cursor_type = GT_Contig;
    xx->cursor_rec  = contig;
    xx->cursor_pos  = xx->displayPos;
    xx->cursor_apos = xx->displayPos;

    xx->trace_lock = 1;

    return xx;
}

/*
 * Deallocates an edview
 */
static void edview_destroy(edview *xx) {
    xfree(xx);
}

/* ACGT to 0123 conversion */
static unsigned char lookup[256], lookup_done = 0;

static seq_t *get_seq(GapIO *io, int rec) {
    return (seq_t *)cache_search(io, GT_Seq, rec);
}


/* ----- 'brief' line manipulation ----- */

static void add_number(char *buf, int *j, int l1, int l2, int val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*d", l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*d", l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*d", l2, val);
	else
	    *j += sprintf(buf + *j, "%d", val);
}

static void add_double(char *buf, int *j, int l1, int l2, double val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*f", l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*f", l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*f", l2, val);
	else
	    *j += sprintf(buf + *j, "%f", val);
}

static void add_string(char *buf, int *j, int l1, int l2, char *str) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*s", l1, l2, str);
	else
	    *j += sprintf(buf + *j, "%*s", l1, str);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*s", l2, str);
	else
	    *j += sprintf(buf + *j, "%s", str);
}


/*
 * Formats reading information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %n	Reading name (Raw: number)
 * %#	Reading number
 * %t	Trace name
 * %p	Position
 * %l	Clipped length
 * %L	Total length
 * %s	Start of clip
 * %e	End of clip
 * %m   Mapping quality
 * %S   Sense (whether complemented, +/-, Raw: 0/1)
 * %a	Chemistry (primer/terminator, Raw: integer)
 * %d	Strand (+/-, Raw 0/1)
 * %i	Reading freetext 'info' comment
 * %P	Primer (unknown/forward universal/reverse universal/forward custom/
 *              reverse custom,  Raw: 0/1/2/3/4)
 * %t   Trace name
 * %Tn	Template name (Raw: template number)
 * %T#	Template number
 * %Tv	Template vector (Raw: template vector number)
 * %Tc	Template consistency (Raw: as a number)
 * %Ti	Template insert size
 * %Cn	Clone name (Raw: clone number)
 * %C#	Clone number
 * %Cv	Clone vector (Raw: clone vector number)
 * %b   Base call
 * %c   Base confidence
 *
 * Additionally specifying %<number><format> forces AT MOST that many
 * characters to be displayed.
 * Specifying %R<format> (or %<number>R<format>) requests the raw data to
 * be displayed. This only works for some fields. Eg %Rp displays 0 to 4, but
 * %p displays, for instance, "forward universal"
 */
char *edGetBriefSeq(edview *xx, int seq, int pos, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    GapIO *io = xx->io;
    seq_t *s = get_seq(io, seq);
    
    for (i = j = 0; format[i]; i++) {
	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case '#':
	    add_number(status_buf, &j, l1, l2, seq);
	    break;

	case 'n':
	    if (raw)
		add_number(status_buf, &j, l1, l2, seq);
	    else
		add_string(status_buf, &j, l1, l2, s->name);
	    break;

	case 'p': {
	    int cnum, cpos;
	    sequence_get_position(xx->io, seq, &cnum, &cpos);
	    add_number(status_buf, &j, l1, l2, cpos);
	    break;
	}

	case 'l':
	    add_number(status_buf, &j, l1, l2, ABS(s->len));
	    break;

	case 'L':
	    add_number(status_buf, &j, l1, l2, s->right - s->left + 1);
	    break;

	case 's':
	    add_number(status_buf, &j, l1, l2, s->left);
	    break;

	case 'e':
	    add_number(status_buf, &j, l1, l2, s->right);
	    break;

	case 'S':
	    if (raw)
		add_number(status_buf, &j, l1, l2, s->len < 0);
	    else
		add_string(status_buf, &j, l1, l2, s->len < 0 ? "-" : "+");
	    break;

	case 'd':
	    {
		int strand = sequence_get_len(&s) < 0 ? 1 : 0;

		if (raw)
		    add_number(status_buf, &j, l1, l2, strand);
		else {
		    char *str;
		    if      (strand == 0) str = "+";
		    else if (strand == 1) str = "-";
		    else                  str = "?";
		    add_string(status_buf, &j, l1, l2, str);
		}
	    }
	    break;

	case 'b':
	    if (pos >= 0 && pos < ABS(s->len)) {
		char base[2];
		base[0] = s->seq[pos];
		base[1] = 0;
		add_string(status_buf, &j, l1, l2, base);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'c':
	    if (pos >= 0 && pos < ABS(s->len)) {
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       1 - pow(10, s->conf[pos]/-10.0));
		} else {
		    add_number(status_buf, &j, l1, l2, s->conf[pos]);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'm':
	    if (raw) {
		add_double(status_buf, &j, l1, l2,
			   1 - pow(10, s->mapping_qual/-10.0));
	    } else {
		add_number(status_buf, &j, l1, l2, s->mapping_qual);
	    }
	    break;
	    

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return status_buf;
}

/*
 * Formats consensus information for the status line.
 * This is done using a format string where certain % rules are replaced by
 * appropriate components.
 *
 * %%	Single % sign
 * %n	Contig name
 * %#	Contig number
 * %p	Position
 * %l	Length
 * %s	Start of clip
 * %e	End of clip
 * %b   Base call
 * %c   Base confidence log-odds (raw for probability value)
 * %A   A confidence log-odds (raw for probability value)
 * %C   C confidence log-odds (raw for probability value)
 * %G   G confidence log-odds (raw for probability value)
 * %T   T confidence log-odds (raw for probability value)
 * %*   * (gap) confidence
 *
 * Additionally specifying %<number><format> forces AT MOST that many
 * characters to be displayed.
 * Specifying %R<format> (or %<number>R<format>) requests the raw data to
 * be displayed. This only works for some fields. Eg %Rp displays 0 to 4, but
 * %p displays, for instance, "forward universal"
 */
char *edGetBriefCon(edview *xx, int crec, int pos, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    
    for (i = j = 0; format[i]; i++) {
	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case '#':
	    add_number(status_buf, &j, l1, l2, crec);
	    break;

	case 'n':
	    if (raw)
		add_number(status_buf, &j, l1, l2, crec);
	    else
		add_string(status_buf, &j, l1, l2,
			   contig_get_name(&xx->contig));
	    break;

	case 'p': {
	    add_number(status_buf, &j, l1, l2, pos);
	    break;
	}

	case 'l':
	    add_number(status_buf, &j, l1, l2, contig_get_length(&xx->contig));
	    break;

	case 's':
	    add_number(status_buf, &j, l1, l2, contig_get_start(&xx->contig));
	    break;

	case 'e':
	    add_number(status_buf, &j, l1, l2, contig_get_end(&xx->contig));
	    break;

	case 'b':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		char base[2];
		base[0] =  xx->displayedConsensus[pos - xx->displayPos];
		base[1] = 0;
		add_string(status_buf, &j, l1, l2, base);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'c':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q =
		    xx->cachedConsensus[p].scores[xx->cachedConsensus[p].call];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'A':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[0];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'C':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[1];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'G':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[2];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'T':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[3];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case '*':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[4];
		if (raw) {
		    add_double(status_buf, &j, l1, l2,
			       pow(10, q/10.0) / (1 + pow(10, q/10.0)));
		} else {
		    add_double(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return status_buf;
}

/*
 * Compute a basic non-weighted consensus. We simply pick the basecall
 * most frequently used.
 *
 * FIXME: use a weighted sum based on confidence values instead?
 */
int calc_cons(GapIO *io, rangec_t *r, int nr, int xpos, int wid,
	      char *cons) {
    int i, j;
    int (*cvec)[6] = (int (*)[6])calloc(wid, 6 * sizeof(int));

    if (!lookup_done) {
	memset(lookup, 5, 256);
	lookup_done = 1;
	lookup['A'] = lookup['a'] = 0;
	lookup['C'] = lookup['c'] = 1;
	lookup['G'] = lookup['g'] = 2;
	lookup['T'] = lookup['t'] = 3;
	lookup['*'] = lookup[','] = 4;
    }

    /* Accumulate */
    for (i = 0; i < nr; i++) {
	int sp = r[i].start;
	seq_t *s = get_seq(io, r[i].rec);
	seq_t *sorig = s;
	int l = s->len > 0 ? s->len : -s->len;
	unsigned char *seq;
	int left, right;

	/* Complement data on-the-fly */
	if ((s->len < 0) ^ r[i].comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	if (s->len < 0) {
	    sp += s->len+1;
	}

	seq = (unsigned char *)s->seq;
	left = s->left;
	right = s->right;

	if (sp < xpos) {
	    seq   += xpos - sp;
	    l     -= xpos - sp;
	    left  -= xpos - sp;
	    right -= xpos - sp;
	    sp = xpos;
	}
	if (l > wid - (sp-xpos))
	    l = wid - (sp-xpos);
	if (left < 1)
	    left = 1;

	for (j = left-1; j < right; j++) {
	    if (sp-xpos+j < wid)
		cvec[sp-xpos+j][lookup[seq[j]]]++;
	}
	cache_decr(io, sorig);

	if (s != sorig)
	    free(s);
    }

    memset(cons, ' ', wid);

    /* and speculate :-) */
    for (i = 0; i < wid; i++) {
	int max, max_base = 5;
	for (max = j = 0; j < 6; j++) {
	    if (max < cvec[i][j]) {
		max = cvec[i][j];
		max_base = j;
	    }
	}
	cons[i] = "ACGT*N"[max_base];
    }

    free(cvec);

    return 0;

}

/*
 * Given an X,Y coordinate return the reading id under this position.
 *
 * Returns record number on success
 *         -1 on failure.
 */
int edGetGelNumber(edview *xx, int x, int y) {
    if (xx->editorState == StateDown)
	return -1;

    if (y < 0 || y >= xx->displayHeight ||
	x < 0 || x >= xx->displayWidth)
	return -1;

    puts("edGetGelNumber unimplemented");
    return 0;
}

/* Sort comparison function for range_t; sort by ascending position */
static int sort_range(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    return r1->start - r2->start;
}

static int ed_set_xslider_pos(edview *xx, int offset) {
    char buf[100];
    double len = contig_get_length(&xx->contig);

    sprintf(buf, " %.20f %.20f",
	    offset / len,
	    (offset + xx->displayWidth) / len);

    if (Tcl_VarEval(xx->interp, xx->ed->xScrollCmd, buf, NULL)
	!= TCL_OK) {
        Tcl_AddErrorInfo(xx->interp, "\n(xscrollcommand executed by Editor)");
        Tcl_BackgroundError(xx->interp);
	return -1;
    }

    return 0;
}

static int ed_set_yslider_pos(edview *xx, int offset, int size, int total) {
    char buf[100];
    double len = total;

    sprintf(buf, " %.20f %.20f",
	    offset / len,
	    (offset + size) / len);

    if (Tcl_VarEval(xx->interp, xx->ed->yScrollCmd, buf, NULL)
	!= TCL_OK) {
        Tcl_AddErrorInfo(xx->interp, "\n(yscrollcommand executed by Editor)");
        Tcl_BackgroundError(xx->interp);
	return -1;
    }

    return 0;
}

/* Update X scrollbar of names display */
void ed_set_nslider_pos(edview *xx, int pos) {
    edNames *en = xx->names;
    char buf[1024];

    if (!en || xx->editorState == StateDown)
	return;

    if (en->xScrollCmd) {
	double fract1, fract2;
	fract1 = pos / (double)MAX_NAME_LEN;
	fract2 = (pos + en->sw.columns) / (double)MAX_NAME_LEN;
	sprintf(buf, " %.20f %.20f", fract1, fract2);
	if (Tcl_VarEval(EDINTERP(en), en->xScrollCmd, buf, NULL) != TCL_OK) {
	    printf("Error in editor names scroll: %s\n", EDINTERP(en)->result);
	}
    }
}

static void tk_redisplaySeqSequences(edview *xx, rangec_t *r, int nr) {
    int i, j;

    /*
    sheet_clear(&xx->ed->sw);
    sheet_clear(&xx->names->sw);
    */

    for (j = xx->y_seq_start, i = xx->displayYPos;
	 j < xx->displayHeight - xx->y_seq_end && i < nr;
	 i++, j++) {
	seq_t *s = get_seq(xx->io, r[i].rec);
	seq_t *sorig = s;
	int sp = r[i].start;
	int l = s->len > 0 ? s->len : -s->len;
	unsigned char seq_a[MAX_SEQ_LEN+1], *seq = seq_a;
	XawSheetInk ink[MAX_DISPLAY_WIDTH];
	int dir = '+';
	int left, right;
	char *conf;
	int seq_p = 0;

	/* Optimisation for single sequence only */
	if (xx->refresh_flags & ED_DISP_SEQ &&
	    !(xx->refresh_flags & ED_DISP_SEQS)) {
	    if (xx->refresh_seq != r[i].rec) {
		continue;
	    }
	}
	
	/* Complement data on-the-fly */
	if ((s->len < 0) ^ r[i].comp) {
	    dir = '-';
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	if (s->len < 0) {
	    sp += s->len+1;
	}

	left = s->left;
	right = s->right;

	memcpy(seq, s->seq, l);
	conf = s->conf;

	if (sp < xx->displayPos) {
	    seq_p += xx->displayPos - sp;
	    //seq   += xx->displayPos - sp;
	    //conf  += xx->displayPos - sp;
	    l     -= xx->displayPos - sp;
	    left  -= xx->displayPos - sp;
	    right -= xx->displayPos - sp;
	    sp = xx->displayPos;
	}
	if (l > xx->displayWidth - (sp-xx->displayPos))
	    l = xx->displayWidth - (sp-xx->displayPos);

	/* Sequence */
	if (xx->refresh_flags & (ED_DISP_READS | ED_DISP_SEQ)) {
	    char line[MAX_DISPLAY_WIDTH+1];
	    int p, p2;

	    memset(ink, 0, MAX_DISPLAY_WIDTH * sizeof(*ink));
	    memset(line, ' ', MAX_DISPLAY_WIDTH);

	    for (p2 = sp - xx->displayPos, p = 0; p < l; p++, p2++) {
		char base = seq[p+seq_p];
		int qual = conf[p+seq_p];
		line[p2] = base;

		/* Cutoffs */
		if (p < left-1 || p > right-1) {
		    if (xx->ed->display_cutoffs) {
			ink[p2].sh |= sh_light;
			//seq[p+seq_p] = tolower(seq[p+seq_p]);
		    } else {
			line[p2] = ' ';
		    }
		}

		/* Quality values */
		if (xx->ed->display_cutoffs || (p >= left-1 && p <= right-1)) {
		    if (xx->ed->display_quality) {
			int qbin = qual / 10;
			if (qbin < 0) qbin = 0;
			if (qbin > 9) qbin = 9;
			ink[p2].sh |= sh_bg;
			ink[p2].bg = xx->ed->qual_bg[qbin]->pixel;
		    }
		}

		/* Highlight disagreements */
		if (xx->ed->display_differences && line[p2] != ' ') {
		    char ubase = xx->ed->display_differences_case
			? base
			: toupper(base);
		    switch (xx->ed->display_differences) {
		    case 1:
			if (ubase == xx->displayedConsensus[p2])
			    line[p2] = '.';
			else if (qual < xx->ed->display_differences_qual)
			    line[p2] = ':';
			break;

		    case 2:
			if (ubase != xx->displayedConsensus[p2]) {
			    ink[p2].sh |= sh_fg;
			    if (qual >= xx->ed->display_differences_qual)
				ink[p2].fg = xx->ed->diff2_fg->pixel;
			    else
				ink[p2].fg = xx->ed->diff1_fg->pixel;
			}
			break;

		    case 3:
			if (ubase != xx->displayedConsensus[p2]) {
			    ink[p2].sh |= sh_bg;
			    if (qual >= xx->ed->display_differences_qual)
				ink[p2].bg = xx->ed->diff2_bg->pixel;
			    else
				ink[p2].bg = xx->ed->diff1_bg->pixel;
			}
			break;
		    }
		}
	    }

	    XawSheetPutJazzyText(&xx->ed->sw, 0, j, xx->displayWidth,
				 line, ink);
	}

	/* Name */
	if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME)) {
	    char name[1024];
	    XawSheetInk ink[1024];
	    int nl = s->name_len - xx->names_xPos;
	    int ncol = xx->names->sw.columns;
	    
	    memset(name, ' ', ncol);
	    name[0] = dir;
	    if (nl > 0)
		memcpy(&name[1], s->name + xx->names_xPos, nl);
	    if (xx->ed->display_mapping_quality) {
		int i, qbin;
		qbin = s->mapping_qual / 10;
		if (qbin < 0) qbin = 0;
		if (qbin > 9) qbin = 9;
		for (i = 0; i < ncol && i < 1024; i++) {
		    ink[i].sh = sh_bg;
		    ink[i].bg = xx->ed->qual_bg[qbin]->pixel;
		}
		XawSheetPutJazzyText(&xx->names->sw, 0, j, ncol, name, ink);
	    } else {
		XawSheetPutText(&xx->names->sw, 0, j, ncol, name);
	    }
	}

	cache_decr(xx->io, sorig);

	if (s != sorig)
	    free(s);
    }

    /*
     * Clear any blank region too.
     */
    if (xx->refresh_flags & ED_DISP_SEQS) {
	char line[MAX_DISPLAY_WIDTH];
	memset(line, ' ', MAX_DISPLAY_WIDTH);
	for (; j < xx->displayHeight - xx->y_seq_end; j++) {
	    XawSheetPutText(&xx->ed->sw, 0, j, xx->displayWidth, line);
	    XawSheetPutText(&xx->names->sw, 0, j, xx->names->sw.columns, line);
	}
    }
 
    //ed_set_xslider_pos(xx, xx->displayPos);
    //ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight, nr);

    return;
}


static void tk_redisplaySeqConsensus(edview *xx, rangec_t *r, int nr) {
    int pos = xx->displayPos;
    int wid =  xx->displayWidth;
    char name[] = " Consensus";
    int i;

    /* Names panel */
    XawSheetPutText(&xx->names->sw, 0, xx->y_cons, strlen(name), name);

    /* Editor panel */
    //calc_cons(xx->io, r, nr, pos, wid, xx->displayedConsensus);
    /*
    calculate_consensus_simple(xx->io, xx->cnum, pos, pos+wid,
			       xx->displayedConsensus, xx->displayedQual);
    */

    calculate_consensus(xx->io, xx->cnum, pos, pos+wid, xx->cachedConsensus);
    for (i = 0; i < wid; i++) {
	xx->displayedConsensus[i] = "ACGT*N"[xx->cachedConsensus[i].call];
    }

    if (xx->ed->display_quality) {
	int i, qbin;
	XawSheetInk ink[MAX_DISPLAY_WIDTH];
	memset(ink, 0, MAX_DISPLAY_WIDTH * sizeof(*ink));

	for (i = 0; i < wid; i++) {
	    qbin = xx->cachedConsensus[i].phred/10;
	    if (qbin < 0) qbin = 0;
	    if (qbin > 9) qbin = 9;
	    ink[i].sh |= sh_bg;
	    ink[i].bg = xx->ed->qual_bg[qbin]->pixel;
	}
	XawSheetPutJazzyText(&xx->ed->sw, 0, xx->y_cons, wid,
			     xx->displayedConsensus, ink);
    } else {
	XawSheetPutText(&xx->ed->sw, 0, xx->y_cons, wid,
			xx->displayedConsensus);
    }
}

/*
 * Calculates the numbers for the contig editor ruler line.
 * The return value is the index into this ruler buffer to plot.
 */
static int generate_ruler(edview *xx, char *ruler, int pos, int width) {
    char *k = ruler;
    int j;

    //int padded_pos[MAX_DISPLAY_WIDTH+21];

    memset(ruler, ' ', MAX_DISPLAY_WIDTH+21);
#if 0
    if (DBI(xx)->reference_seq) {
	/* Number relative to a specific sequence number */
	char *bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
	int reflen = DB_Length(xx, DBI(xx)->reference_seq);
	int rp = DB_RelPos(xx, DBI(xx)->reference_seq);

	/* 
	 * Compute unpadded base positions for this window.
	 * WARNING: This code is full of warts. I carefully worked it out, but
	 * it was not quite perfect. Working out again gave a slightly
	 * different, but equally imperfect answer. So I applied the next
	 * logical step - I hope you like fudge!
	 *
	 * If the change the code - check for these cases and combinations:
	 *   Position from left end when uncomplemented
	 *   Position from right end when complemented
	 *   Negative offsets
	 *   When relpos of ref is 1 and when it is not, both when refseq
	 *     is complemented and when it is not.
	 *   Effect of pads in ref seq
	 *   That base numbers are never positioned above the pads
	 *   Circular sequences (both orientations)
	 *   Offset base numbers
	 */
	if ((DB_Comp(xx, DBI(xx)->reference_seq) == UNCOMPLEMENTED)) {
	    int unpadded = MIN(0, pos);

	    for (j = unpadded; j < pos + width + 9; j++) {
		if (j >= pos)
		    padded_pos[j-pos] = unpadded-rp;
		if (j-rp < 0 || j-rp >= reflen || bases[j-rp] != '*')
		    unpadded++;
	    }
	} else {
	    int unpadded = reflen - MAX(pos + width + 8, reflen);

	    for (j = MAX(pos + width + 8, reflen); j >= pos; j--) {
		if (j <= pos + width + 8)
		    padded_pos[j-pos] = unpadded+(rp-1);
		if (j-rp-1 < 0 || j-rp-1 >= reflen || bases[j-rp-1] != '*')
		    unpadded++;
	    }
	}

	for (j = pos; j < pos + width + 9; j++, k++) {
	    int unpadded2 = (padded_pos[j-pos] + DBI(xx)->reference_offset);

	    if (DBI(xx)->reference_len) {
		unpadded2 %= DBI(xx) -> reference_len;
		while (unpadded2 < 0) {
		    unpadded2 += DBI(xx)->reference_len;
		}
		if (unpadded2 == 0)
		    continue; /* Don't display 0 for circular seqs */
	    }
	    if (!(unpadded2 % 10)) {
		sprintf(k, "%10d", unpadded2);
	    }
	}

	/*
	 * Pads in the consensus can sometimes leave nulls in the buffer,
	 * which turn into square boxes on some fonts. Replace these with
	 * spaces.
	 */
	for (j = 0; j < MAX_DISPLAY_WIDTH+21; j++)
	    if (ruler[j] == '\0')
		ruler[j] = ' ';

	return 9;

    } else if (xx->unpadded_ruler) {
	/* Number by unpadded base calls, using the consensus */
	int unpadded;
	
	edUnpaddedBaseNumber(xx, pos, width + 9);
	for (j = pos; j < pos + width + 9; j++, k++) {
	    unpadded = edUnpaddedBaseNumber(xx, j, 0);
	    if (unpadded%10)
		continue;
	    sprintf(k, "%10d", unpadded);
	    if (k+10-ruler < MAX_DISPLAY_WIDTH+21)
		k[10]=' ';
	}
	edUnpaddedBaseNumber(xx, pos, -1);

	return 9;
    } /* else */
#endif
    {
	/* Basic numbering */
	int lower,times;
	lower = (pos - pos%10);
	times = width/10 + 3;
	for (j=0;j<times;j++,k+=10,lower+=10)
	    sprintf(k,"%10d",lower);
	return 9+pos%10;
	
    }
}

static void tk_redisplaySeqNumbers(edview *xx) {
    char ruler[MAX_DISPLAY_WIDTH+21];
    int off;

    off = generate_ruler(xx, ruler, xx->displayPos, xx->displayWidth);
    XawSheetPutText(&xx->ed->sw, 0, xx->y_numbers, xx->displayWidth,
		    &ruler[off]);
}


/* Handle scrolling changes */
static void tk_redisplaySeqScroll(edview *xx, rangec_t *r, int nr) {
    if (xx->refresh_flags & ED_DISP_YSCROLL) {
	ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight, nr);

	/* No need to redraw the consenus/numbers */
	xx->refresh_flags |= ED_DISP_NAMES | ED_DISP_SEQS | ED_DISP_CURSOR;
    }

    if (xx->refresh_flags & ED_DISP_XSCROLL) {
	ed_set_xslider_pos(xx, xx->displayPos);

	/* Changing X may also change height */
	if (!(xx->refresh_flags & ED_DISP_YSCROLL))
	    ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight, nr);

	xx->refresh_flags |= ED_DISP_ALL;
    }
}

/*
 * Forces the sheet widget to redisplay the cursor.
 */
static void tk_redisplayCursor(edview *xx, rangec_t *r, int nr) {
    int x, y;

    /* If visible, find the screen x/y coord */
    if (xx->cursor_rec == xx->cnum) {
	y = xx->y_cons;
    } else {
	int i, j;
	y = -1;
	for (j = xx->y_seq_start, i = xx->displayYPos;
	     j < xx->displayHeight - xx->y_seq_end && i < nr;
	     i++, j++) {
	    if (r[i].rec == xx->cursor_rec) {
		y = j;
		break;
	    }
	}
	
	if (-1 == y) {
	    XawSheetDisplayCursor(&xx->ed->sw, False);
	    return; /* not visible */
	}
    }

    x = xx->cursor_apos - xx->displayPos;
    XawSheetDisplayCursor(&xx->ed->sw, True);
    XawSheetPositionCursor(&xx->ed->sw, x, y);
}


/*
 * Force the cursor to be visible. If x_safe or y_safe are true then
 * we omit some of the searching and assume there is no reason to check
 * that x or y is still visible.
 *
 * Returns 1 if redraw has taken place
 *         0 if not
 */
static int showCursor(edview *xx, int x_safe, int y_safe) {
    int y_pos = 0;
    int do_x = 0;
    int do_y = 0;
    int nr;

    /* X position */
    if (!x_safe) {
	int w = xx->displayWidth > 10 ? 10 : xx->displayWidth;
	if (xx->cursor_apos < xx->displayPos) {
	    xx->displayPos = xx->cursor_apos + 1 - w;
	    do_x = 1;
	}
	if (xx->cursor_apos >= xx->displayPos + xx->displayWidth) {
	    xx->displayPos = xx->cursor_apos - xx->displayWidth + w;
	    do_x = 1;
	}
    }

    if (do_x)
	y_safe = 0;

    /* Y position */
    if (!y_safe && xx->cursor_type != GT_Contig) {
	rangec_t *r;
	int i;
	int sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;
	
	/* Find out what's visible - cache this */
	r = contig_seqs_in_range(xx->io, &xx->contig, xx->displayPos,
				 xx->displayPos + xx->displayWidth, &nr);
	qsort(r, nr, sizeof(*r), sort_range);
	
	for (i = 0; i < nr; i++) {
	    if (r[i].rec == xx->cursor_rec) {
		y_pos = i;
		break;
	    }
	}
	free(r);

	/* If above, scroll so this is the first row */
	if (y_pos < xx->displayYPos) {
	    xx->displayYPos = y_pos;
	    do_y = 1;
	}

	/* If below, scroll so this is the last row */
	if (y_pos >= xx->displayYPos + sheight) {
	    xx->displayYPos = y_pos - sheight + 1;
	    do_y = 1;
	}
    }

    tman_reposition_traces(xx, xx->cursor_apos, 0);

    if (do_x)
	ed_set_xslider_pos(xx, xx->displayPos);

    if (do_y)
	ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight, nr);

    if (do_x || do_y) {
	xx->refresh_flags = ED_DISP_ALL;
	edview_redraw(xx);
	return 1;
    }

    return 0;
}

int set_displayPos(edview *xx, int pos) {
    char buf[100];
    xx->displayPos = pos;

    sprintf(buf, "%d", pos);
    Tcl_SetVar2(xx->interp, xx->edname, "displayPos", buf, TCL_GLOBAL_ONLY);

    xx->refresh_flags = ED_DISP_XSCROLL;
    return edview_redraw(xx);
}

/*
 * Resets the absolute position in the contig based on the contents of the
 * cursor_pos and cursor_rec fields.
 */
void edSetApos(edview *xx) {
    if (xx->cursor_type == GT_Contig) {
	xx->cursor_apos = xx->cursor_pos;
    } else {
	int cnum, cpos;
	sequence_get_position(xx->io, xx->cursor_rec, &cnum, &cpos);
	xx->cursor_apos = cpos + xx->cursor_pos;
    }
}

int edSetCursorPos(edview *xx, int type, int rec, int pos) {
    if (type == GT_Seq) {
	seq_t *s = get_seq(xx->io, rec);

	if (pos < 0)
	    pos = 0;
	if (pos > ABS(s->len))
	    pos = ABS(s->len);
    }

    xx->cursor_type = type;
    xx->cursor_rec  = rec;
    xx->cursor_pos  = pos;

    edSetApos(xx);

    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}


int edCursorUp(edview *xx) {
    int j, nr;
    int cpos = xx->cursor_apos;
    rangec_t *r;

    /* Find out what's visible - cache this*/
    r = contig_seqs_in_range(xx->io, &xx->contig, xx->displayPos,
			     xx->displayPos + xx->displayWidth, &nr);
    qsort(r, nr, sizeof(*r), sort_range);

    if (nr == 0)
	return 0;

    /* Find the current sequence number */
    if (xx->cursor_type == GT_Contig) {
	j = nr;
    } else {
	for (j = 0; j < nr; j++) {
	    if (r[j].rec == xx->cursor_rec) {
		break;
	    }
	}
    }

    /* Step up until we find something overlapping */
    for (j--; j >= 0; j--) {
	if (r[j].start <= cpos && r[j].end >= cpos) {
	    xx->cursor_type = GT_Seq;
	    xx->cursor_pos = cpos - r[j].start;
	    xx->cursor_rec = r[j].rec;
	    break;
	}
    }

    /* Otherwise we've hit the consensus */
    if (j < 0) {
	xx->cursor_type = GT_Contig;
	xx->cursor_rec = xx->cnum;
	xx->cursor_pos = cpos;
    }

    free(r);

    if (!showCursor(xx, 1, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorDown(edview *xx) {
    int j, nr;
    int cpos = xx->cursor_apos;
    rangec_t *r;

    /* Find out what's visible - cache this*/
    r = contig_seqs_in_range(xx->io, &xx->contig, xx->displayPos,
			     xx->displayPos + xx->displayWidth, &nr);
    qsort(r, nr, sizeof(*r), sort_range);

    if (nr == 0)
	return 0;

    /* Find the current sequence number */
    if (xx->cursor_type == GT_Contig) {
	cpos = xx->cursor_pos;
	j = -1;
    } else {
	for (j = 0; j < nr; j++) {
	    if (r[j].rec == xx->cursor_rec) {
		cpos = r[j].start + xx->cursor_pos;
		break;
	    }
	}
    }

    /* Step up until we find something overlapping */
    for (j++; j < nr; j++) {
	if (r[j].start <= cpos && r[j].end >= cpos) {
	    xx->cursor_type = GT_Seq;
	    xx->cursor_pos = cpos - r[j].start;
	    xx->cursor_rec = r[j].rec;
	    break;
	}
    }

    /* Otherwise we've hit the consensus */
    if (j >= nr) {
	xx->cursor_type = GT_Contig;
	xx->cursor_rec = xx->cnum;
	xx->cursor_pos = cpos;
    }

    free(r);

    if (!showCursor(xx, 1, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorLeft(edview *xx) {
    if (xx->cursor_type == GT_Seq) {
	seq_t *s = get_seq(xx->io, xx->cursor_rec);
	seq_t *sorig = s;

	if (s->len < 0) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	if (xx->ed->display_cutoffs) {
	    if (xx->cursor_pos > 0) {
		xx->cursor_pos--;
		xx->cursor_apos--;
	    }
	} else {
	    if (xx->cursor_pos >= s->left) {
		xx->cursor_pos--;
		xx->cursor_apos--;
	    }

	}

	if (s != sorig)
	    free(s);
    } else {
	xx->cursor_pos--;
	xx->cursor_apos--;
    }

    if (!showCursor(xx, 0, 1)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorRight(edview *xx) {
    if (xx->cursor_type == GT_Seq) {
	seq_t *s = get_seq(xx->io, xx->cursor_rec);
	seq_t *sorig = s;

	if (s->len < 0) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	if (xx->ed->display_cutoffs) {
	    if (xx->cursor_pos < ABS(s->len)) {
		xx->cursor_pos++;
		xx->cursor_apos++;
	    }
	} else {
	    if (xx->cursor_pos < s->right) {
		xx->cursor_pos++;
		xx->cursor_apos++;
	    }
	}

	if (s != sorig)
	    free(s);
    } else {
	xx->cursor_pos++;
	xx->cursor_apos++;
    }

    if (!showCursor(xx, 0, 1)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadStart(edview *xx) {
    if (xx->cursor_type == GT_Seq) {
	xx->cursor_pos = 0;
    } else {
	xx->cursor_pos = xx->contig->start;
    }

    edSetApos(xx);

    if (!showCursor(xx, 0, 1)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadStart2(edview *xx) {
    return edReadStart(xx);
}

int edReadEnd(edview *xx) {
    if (xx->cursor_type == GT_Seq) {
	seq_t *s = get_seq(xx->io, xx->cursor_rec);
	xx->cursor_pos = ABS(s->len);
    } else {
	xx->cursor_pos = xx->contig->end;
    }

    edSetApos(xx);

    if (!showCursor(xx, 0, 1)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadEnd2(edview *xx) {
    return edReadEnd(xx);
}

int edContigStart(edview *xx) {
    xx->cursor_pos = xx->contig->start;
    xx->cursor_type = GT_Contig;
    xx->cursor_rec = xx->cnum;
    xx->cursor_apos = xx->cursor_pos;

    if (!showCursor(xx, 0, 1)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edContigEnd(edview *xx) {
    xx->cursor_pos = xx->contig->end;
    xx->cursor_type = GT_Contig;
    xx->cursor_rec = xx->cnum;
    xx->cursor_apos = xx->cursor_pos;

    if (!showCursor(xx, 0, 1)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

/*
 * The main editor redraw function
 */
int edview_redraw(edview *xx) {
    rangec_t *r;
    int nr;

    if (!xx->ed || xx->editorState == StateDown)
	return -1;

    if (xx->refresh_flags & ED_DISP_SEQ && xx->refresh_seq == xx->cnum)
	xx->refresh_flags |= ED_DISP_CONS;

    if (xx->displayWidth > MAX_DISPLAY_WIDTH)
	xx->displayWidth = MAX_DISPLAY_WIDTH;

#if 0
    /* Work out the status line; this may control window height */
    cur_depth = xx->status_depth;
    if (xx->refresh_flags & (ED_DISP_STATUS | ED_DISP_SCROLL)) {
	tk_redisplaySeqStatusCompute(xx, xx->displayPos, xx->displayWidth);
    }

#endif

    /* Find out what's visible */
    /* FIXME: Cache this */
    r = contig_seqs_in_range(xx->io, &xx->contig, xx->displayPos,
			     xx->displayPos + xx->displayWidth, &nr);
    qsort(r, nr, sizeof(*r), sort_range);

    /* Deal with ED_DISP_XSCROLL and ED_DISP_YSCROLL events */
    if (xx->refresh_flags & (ED_DISP_XSCROLL | ED_DISP_YSCROLL))
	tk_redisplaySeqScroll(xx, r, nr);

    /* Consensus emacs-style edit-status ----, -%%-, -**- */
    //tk_redisplaySeqEditStatus(xx);

    /* Redraw the consensus and/or numbers */
    if (xx->refresh_flags & ED_DISP_CONS) {
	tk_redisplaySeqConsensus(xx, r, nr);
    }
    if (xx->refresh_flags & ED_DISP_RULER) {
	tk_redisplaySeqNumbers(xx);
    }

    /* Redraw the main sequences or names section */
    if (xx->refresh_flags & (ED_DISP_SEQS  | ED_DISP_SEQ |
			     ED_DISP_NAMES | ED_DISP_NAME)) {
	tk_redisplaySeqSequences(xx, r, nr);
    }

    /* Editor cursor position */
    if (xx->refresh_flags & ED_DISP_CURSOR) {
	tk_redisplayCursor(xx, r, nr);
    }

#if 0
    /* We've already computed them, but now we actually display them */
    if (xx->refresh_flags & ED_DISP_STATUS) {
	tk_redisplaySeqStatusDisplay(xx);
    }

    /* Underlining for current selection */
    if (xx->refresh_flags & ED_DISP_SELECTION) {
	redisplaySelection(xx);
    }
#endif

    free(r);

#if 0
    /* FIXME: only need to redraw here if major change => scrolling etc */
    if (xx->refresh_flags & (ED_DISP_SEQS))
	sheet_display(&xx->ed->sw);
    if (xx->refresh_flags & (ED_DISP_NAMES))
	sheet_display(&xx->names->sw);
#endif

    xx->refresh_flags = 0;
    xx->refresh_seq = 0;

    return 0;
}

/*
 * Identifies the type of object underneath a specific row and column.
 *
 * Returns the item type GT_* on success and the record/pos in *rec, *pos
 *         -1 on failure (eg numbers, off screen, etc)
 */
int edview_item_at_pos(edview *xx, int row, int col, int *rec, int *pos) {
    rangec_t *r;
    int nr, i, j;
    int type = -1;

    /* Special case - the reserve row numbers */
    if (row == xx->y_cons) {
	*rec = xx->cnum;
	*pos = col + xx->displayPos;
	return GT_Contig;
    }

    if (! (row >= xx->y_seq_start && row < xx->displayHeight - xx->y_seq_end))
	return -1;

    /* A sequence, so find out what's visible */
    r = contig_seqs_in_range(xx->io, &xx->contig, xx->displayPos,
			     xx->displayPos + xx->displayWidth, &nr);
    qsort(r, nr, sizeof(*r), sort_range);

    /* Inefficient, but just a copy from tk_redisplaySeqSequences() */
    for (j = xx->y_seq_start, i = xx->displayYPos;
	 j < xx->displayHeight - xx->y_seq_end && i < nr;
	 i++, j++) {
	if (j == row) {
	    *rec = r[i].rec;
	    *pos = col + xx->displayPos - r[i].start;
	    type = GT_Seq;
	    break;
	}
    }

    free(r);
    return type;
}


int inJoinMode(edview *xx) {return 0;}

void edDisplayTrace(edview *xx) {
    seq_t *s;

    if (xx->cursor_type == GT_Seq) {
	/* Single sequence */
	s = get_seq(xx->io, xx->cursor_rec);
	tman_manage_trace("ANY", sequence_get_name(&s), xx->cursor_pos,
			  0, 0, /* left/right clips */
			  s->len < 0, /* complemented */
			  1, /* base spacing */
			  sequence_get_name(&s),
			  xx, xx->cursor_rec, 0, 0);
    } else if (xx->cursor_type == GT_Contig) {
	/* Consensus click */
	rangec_t *r;
	int nr, i;

	/* Shut down existing traces */
	tman_shutdown_traces(xx, 2);

	/* And display the new ones */
	r = contig_seqs_in_range(xx->io, &xx->contig,
				 xx->cursor_apos, xx->cursor_apos, &nr);
	qsort(r, nr, sizeof(*r), sort_range);

	for (i = 0; i < nr; i++) {
	    s = get_seq(xx->io, r[i].rec);
	    tman_manage_trace("ANY", sequence_get_name(&s), xx->cursor_pos,
			      0, 0, /* left/right clips */
			      s->len < 0, /* complemented */
			      1, /* base spacing */
			      sequence_get_name(&s),
			      xx, r[i].rec, 0, 0);
	}
    }

    tman_reposition_traces(xx, xx->cursor_apos, 0);
}

/*
 * Given a sequence record number this identifies all other sequence
 * records from the same template. The returned array is malloced and should
 * be freed by the caller once finished with.
 *
 * Returns pointer to array of records of size *nrec on success
 *         NULL on failure (or zero found)
 */
int *edGetTemplateReads(edview *xx, int seqrec, int *nrec) {
    seq_t *s = get_seq(xx->io, seqrec);
    int *r = NULL;

    if (!s)
	return NULL;

    /* FIXME: support s->parent_rec and s->parent_type */

    /* Solexa data is easy: we have just one other end */
    if (s->other_end) {
	*nrec = 1;
	r = (int *)malloc(sizeof(*r));
	*r = s->other_end;
    } else {
	*nrec = 0;
    }

    return r;
}
