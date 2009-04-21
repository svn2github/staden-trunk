#include <xalloc.h>
#include <math.h>
#include <ctype.h>
#include <X11/Xatom.h> /* XA_PRIMARY - included in Tk distribrution */

#include "editor_view.h"
#include "tkSheet.h"
#include "tman_interface.h"
#include "gap_globals.h"
#include "qualIO.h"
#include "qual.h"
#include "tg_gio.h"
#include "misc.h"
#include "consensus.h"

static void redisplaySelection(edview *xx);

/*
 * A C interface to the edit_contig and join_contig Tcl functions.
 */
int edit_contig(GapIO *io, int cnum, int rnum, int pos) {
    char cmd[1024];

    sprintf(cmd, "edit_contig -io %s -contig %d -reading %d -pos %d\n",
	    io_obj_as_string(io), cnum, rnum, pos);
    return Tcl_Eval(GetInterp(), cmd);
}

int join_contig(GapIO *io, int cnum[2], int rnum[2], int pos[2]) {
    char cmd[1024];
    int ret;

    sprintf(cmd, "join_contig -io %s -contig %d -reading %d -pos %d "
	    "-contig2 %d -reading2 %d -pos2 %d",
	    io_obj_as_string(io),
	    cnum[0], rnum[0], pos[0],
	    cnum[1], rnum[1], pos[1]);
    ret = Tcl_Eval(GetInterp(), cmd);
    if (ret != TCL_OK) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(GetInterp()));
    }
    return ret;
}

/*
 * Allocates and initialises a new edview
 */
edview *edview_new(GapIO *io, int contig, int crec, int cpos,
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

    xx->cursor_pos  = cpos;
    xx->cursor_type = (xx->cursor_rec == 0 || xx->cursor_rec == contig)
	? GT_Contig : GT_Seq;
    xx->cursor_rec  = crec ? crec : contig;
    edSetApos(xx);
    xx->displayPos = xx->cursor_apos;

    xx->trace_lock = 1;
    
    if (!xx->ed->consensus_at_top) {
	xx->ed->sw.yflip = 1;
	xx->names->sw.yflip = 1;
    }

    xx->r = NULL;
    xx->anno_hash = NULL;
    
    return xx;
}

/*
 * Deallocates an edview
 */
static void edview_destroy(edview *xx) {
    if (xx->r)
	free(xx->r);
    xfree(xx);
}

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
 * %A   A confidence log-odds (raw for probability value)
 * %C   C confidence log-odds (raw for probability value)
 * %G   G confidence log-odds (raw for probability value)
 * %T   T confidence log-odds (raw for probability value)
 *
 * Additionally specifying %<number><format> forces AT MOST that many
 * characters to be displayed.
 * Specifying %R<format> (or %<number>R<format>) requests the raw data to
 * be displayed. This only works for some fields. Eg %Rp displays 0 to 4, but
 * %p displays, for instance, "forward universal"
 * Specying %*<format> indicates that the above formats should be applied to
 * the other end of a read-pair instead of this read.
 * The special format %** is used to terminate decoding of the format if
 * the sequence is single-ended.
 */
char *edGetBriefSeq(edview *xx, int seq, int pos, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    GapIO *io = xx->io;
    seq_t *s1 = get_seq(io, seq), *s2 = NULL, *s;
    int pair = 0;
    
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

	if (format[i] == '*') {
	    if (pair == 0)
		pair = sequence_get_pair(io, s1);
	    if (pair > 0 && !s2)
		s2 = get_seq(io, pair);
	    s = s2 ? s2 : s1;
	    i++;
	} else {
	    s = s1;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case '*':
	    if (!s2)
		goto bail_out;
	    break;

	case '#':
	    add_number(status_buf, &j, l1, l2, s->rec);
	    break;

	case 'n':
	    if (raw)
		add_number(status_buf, &j, l1, l2, s->rec);
	    else
		add_string(status_buf, &j, l1, l2, s->name);
	    break;

	case 'p': {
	    int cnum, cpos;
	    if (0 == sequence_get_position(xx->io, s->rec, &cnum, &cpos, NULL,
					   NULL))
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
		sequence_get_base(xx->io, &s, pos, &base[0], NULL, 1);
		base[1] = 0;
		add_string(status_buf, &j, l1, l2, base);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'c':
	    if (pos >= 0 && pos < ABS(s->len)) {
		int q;
		sequence_get_base(xx->io, &s, pos, NULL, &q, 1);
		if (raw) {
		    add_double(status_buf, &j, l1, l2, 1 - pow(10, q/-10.0));
		} else {
		    add_number(status_buf, &j, l1, l2, q);
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'A':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[0]));
		else {
		    double p = exp(q[0]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'C':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[1]));
		else {
		    double p = exp(q[1]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'G':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[2]));
		else {
		    double p = exp(q[2]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
		}
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'T':
	    if (pos >= 0 && pos < ABS(s->len)) {
		double q[4];
		sequence_get_base4(xx->io, &s, pos, NULL, q, 1);
		if (raw)
		    add_double(status_buf, &j, l1, l2, exp(q[3]));
		else {
		    double p = exp(q[3]);
		    add_double(status_buf, &j, l1, l2, 4.342945*log(p/(1-p)));
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
 bail_out:
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
 * Populates the cache of visible items in xx->r and xx->nr.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int edview_visible_items(edview *xx, int start, int end) {
    int i;
    int mode = xx->ed->stack_mode
	? CSIR_ALLOCATE_Y_MULTIPLE
	: CSIR_ALLOCATE_Y_SINGLE;

    /* Always reload for now as we can't spot edits yet */
    if (0 && xx->r && xx->r_start == start && xx->r_end == end)
	return 0;

    /* Query sequences */
    if (xx->r)
	free(xx->r);
    xx->r_start = start;
    xx->r_end = end;
    xx->r = contig_items_in_range(xx->io, &xx->contig, start, end,
				  CSIR_SORT_BY_Y | mode, &xx->nr);
    if (!xx->r)
	return -1;
    
    /* Work out Y dimension */
    xx->max_height = 0;
    for (i = 0; i < xx->nr; i++) {
	if (xx->max_height < xx->r[i].y)
	    xx->max_height = xx->r[i].y;
    }
    xx->max_height += 3; /* +1 for from 0, +2 for consensus+ruler */

    /* Fast map of annotations to sequences */
    if (xx->anno_hash)
	HacheTableDestroy(xx->anno_hash, 0);

    xx->anno_hash = HacheTableCreate(8192, HASH_DYNAMIC_SIZE |
				     HASH_ALLOW_DUP_KEYS);
    for (i = 0; i < xx->nr; i++) {
	int key = xx->r[i].pair_rec; /* aka obj_rec */
	HacheData hd;
	anno_ele_t *a;

	if (!(xx->r[i].flags & GRANGE_FLAG_ISANNO))
	    continue;

	hd.i = i;
	HacheTableAdd(xx->anno_hash, &key, sizeof(key), hd, NULL);

	/*
	a = (anno_ele_t *)cache_search(xx->io, GT_AnnoEle, xx->r[i].rec);
	printf("Add anno %d/%d/%s to seq %d,%d,%d\n",
	       i, a->rec, a->comment,
	       xx->r[i].pair_rec,
	       xx->r[i].mqual,
	       xx->r[i].flags);
	*/
    }

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

static int ed_set_xslider_pos(edview *xx, int offset) {
    char buf[100];
    double len = contig_get_length(&xx->contig);

    offset -= contig_get_start(&xx->contig);

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

    if (xx->ed->consensus_at_top) {
	sprintf(buf, " %.20f %.20f",
		offset / len,
		(offset + size) / len);
    } else {
	sprintf(buf, " %.20f %.20f",
		(len - offset - size) / len,
		(len - offset) / len);
    }

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

	if (xx->ed->stack_mode) {
	    fract1 = 0;
	    fract2 = 1;
	} else {
	    fract1 = pos / (double)MAX_NAME_LEN;
	    fract2 = (pos + en->sw.columns) / (double)MAX_NAME_LEN;
	}
	sprintf(buf, " %.20f %.20f", fract1, fract2);
	if (Tcl_VarEval(EDINTERP(en), en->xScrollCmd, buf, NULL) != TCL_OK) {
	    printf("Error in editor names scroll: %s\n", EDINTERP(en)->result);
	}
    }
}

static void tk_redisplaySeqSequences(edview *xx, rangec_t *r, int nr) {
    int i, j, k;

    /*
    sheet_clear(&xx->ed->sw);
    sheet_clear(&xx->names->sw);
    */

    /* Work down the screen line by line */
    for (j = xx->y_seq_start, i = 0;
	 j < xx->displayHeight - xx->y_seq_end && i < nr;
	 j++) {
	int sp, l;
	unsigned char seq_a[MAX_SEQ_LEN+1], *seq = seq_a;
	XawSheetInk ink[MAX_DISPLAY_WIDTH], nink[MAX_NAME_WIDTH];
	char line[MAX_DISPLAY_WIDTH+1], nline[MAX_NAME_WIDTH];
	int dir;
	int left, right;
	int seq_p;

	if (xx->refresh_flags & (ED_DISP_READS | ED_DISP_SEQ)) {
	    memset(ink, 0, MAX_DISPLAY_WIDTH * sizeof(*ink));
	    memset(line, ' ', MAX_DISPLAY_WIDTH);
	}
	if (xx->ed->stripe_mode) {
	    int n = xx->ed->stripe_mode;
	    for (k = n-xx->displayPos%n; k < xx->displayWidth; k+=n) {
		ink[k].sh |= sh_bg;
		ink[k].bg = xx->ed->stripe_bg->pixel;
	    }
	}

	if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME)) {
	    memset(nink, 0, MAX_NAME_WIDTH * sizeof(*nink));
	    memset(nline, ' ', MAX_NAME_WIDTH);
	}

	/* Iterate through all sequences on this line */
	if (i >= nr || xx->r[i].y - xx->displayYPos > j - xx->y_seq_start)
	    continue;

	while (i < nr && xx->r[i].y - xx->displayYPos <= j - xx->y_seq_start) {
	    seq_t *s, *sorig;

	    if (r[i].flags & GRANGE_FLAG_ISANNO) {
		i++;
		continue;
	    }

	    if (xx->r[i].y - xx->displayYPos < j - xx->y_seq_start) {
		i++;
		continue;
	    }

	    s = sorig = get_seq(xx->io, r[i].rec);
	    sp = r[i].start;
	    l = s->len > 0 ? s->len : -s->len;
	    seq_p = 0;
	    dir = '>';

	    /* Optimisation for single sequence only */
	    if (xx->refresh_flags & ED_DISP_SEQ &&
		!(xx->refresh_flags & ED_DISP_SEQS)) {
		if (xx->refresh_seq != r[i].rec) {
		    continue;
		}
	    }
	
	    /* Complement data on-the-fly */
	    if ((s->len < 0) ^ r[i].comp) {
		dir = '<';
		s = dup_seq(s);
		complement_seq_t(s);
	    }

	    left = s->left;
	    right = s->right;

	    memcpy(seq, s->seq, l);

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
		int p, p2;

		for (p2 = sp - xx->displayPos, p = 0; p < l; p++, p2++) {
		    char base;
		    int qual;

		    if (p+seq_p >= 0 && p+seq_p < ABS(s->len)) {
			sequence_get_base(xx->io, &s, p+seq_p, &base, &qual,0);
			//ink[p2].sh &= ~sh_bg;
		    } else {
			base = ' ';
			qual = 0;
		    }
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
		    if (xx->ed->display_cutoffs ||
			(p >= left-1 && p <= right-1)) {
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

		/* Annotations */
		if (xx->anno_hash) {
		    HacheItem *hi;
		    int key = s->rec;
		    for (hi = HacheTableSearch(xx->anno_hash, &key,
					       sizeof(key));
			 hi; hi = HacheTableNext(hi, &key, sizeof(key))) {
			int ai = hi->data.i;
			printf("Seq %d anno %d\n", key, ai);

			for (p2 = sp - xx->displayPos, p = 0; p<l; p++,p2++) {
			    if (p2 + xx->displayPos >= xx->r[ai].start &&
				p2 + xx->displayPos <= xx->r[ai].end)
				ink[p2].sh |= sh_inverse;
			}
		    }
		}
	    }

	    /* Name */
	    if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME)) {
		int nl = s->name_len - xx->names_xPos;
		int ncol = xx->names->sw.columns;
	    
		if (xx->ed->stack_mode) {
		    int p  = r[i].start - xx->displayPos;
		    int p2 = r[i].end   - xx->displayPos;
		    int bg = -1;
		    double nc = xx->names->sw.columns;
		    if (p < 0) p = 0;
		    p = p * (nc / xx->displayWidth);
		    if (p2 < 0) p2 = 0;
		    p2 = p2 * (nc / xx->displayWidth);
		    if (p2 > xx->names->sw.columns)
			p2 = xx->names->sw.columns;
		    while (nline[p] != ' ')
			p++;

		    if (xx->ed->display_mapping_quality) {
			int qbin;
			qbin = s->mapping_qual / 10;
			if (qbin < 0) qbin = 0;
			if (qbin > 9) qbin = 9;
			bg = xx->ed->qual_bg[qbin]->pixel;
		    }

		    nline[p] = dir;
		    if (bg != -1) {
			nink[p].sh |= sh_bg;
			nink[p].bg = bg;
		    }
		    
		    for (++p; p < p2; p++) {
			nline[p] = '.';
			if (bg != -1) {
			    nink[p].sh |= sh_bg;
			    nink[p].bg = bg;
			}
		    }

		} else {
		    nline[0] = dir;
		    if (nl > 0)
			memcpy(&nline[1], s->name + xx->names_xPos, nl);
		    nink[0].sh = sh_bg;
		    if (r[i].pair_rec) {
			nink[0].bg = xx->ed->qual_bg[9]->pixel;
		    } else {
			nink[0].bg = xx->ed->qual_bg[0]->pixel;
		    }
		    if (xx->ed->display_mapping_quality) {
			int i, qbin;
			qbin = s->mapping_qual / 10;
			if (qbin < 0) qbin = 0;
			if (qbin > 9) qbin = 9;
			for (i = 1; i < ncol && i < MAX_NAME_WIDTH; i++) {
			    nink[i].sh = sh_bg;
			    nink[i].bg = xx->ed->qual_bg[qbin]->pixel;
			}
		    } else {
			int i;
			for (i = 1; i < ncol && i < 1024; i++) {
			    nink[i].sh = sh_default;
			}
		    }
		}
	    }

	    cache_decr(xx->io, sorig);

	    if (s != sorig)
		free(s);

	    i++;
	}

	if (xx->refresh_flags & (ED_DISP_READS | ED_DISP_SEQ))
	    XawSheetPutJazzyText(&xx->ed->sw, 0, j, xx->displayWidth,
				 line, ink);

	if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME))
	    XawSheetPutJazzyText(&xx->names->sw, 0, j, xx->names->sw.columns,
				 nline, nink);
    }

    /*
     * Clear any blank lines too.
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
    XawSheetInk ink[MAX_DISPLAY_WIDTH];

    /* Names panel */
    XawSheetPutText(&xx->names->sw, 0, xx->y_cons, strlen(name), name);

    /* Editor panel */

    calculate_consensus(xx->io, xx->cnum, pos, pos+wid, xx->cachedConsensus);
    for (i = 0; i < wid; i++) {
	xx->displayedConsensus[i] = "ACGT*N"[xx->cachedConsensus[i].call];
    }

    memset(ink, 0, MAX_DISPLAY_WIDTH * sizeof(*ink));
    if (xx->ed->display_quality) {
	int i, qbin;

	for (i = 0; i < wid; i++) {
	    qbin = xx->cachedConsensus[i].phred/10;
	    if (qbin < 0) qbin = 0;
	    if (qbin > 9) qbin = 9;
	    ink[i].sh |= sh_bg;
	    ink[i].bg = xx->ed->qual_bg[qbin]->pixel;
	}
    }

    /* Consensus annotations */
    for (i = 0; i < xx->nr; i++) {
	int sp  = xx->r[i].start;
	int j;

	if (!(xx->r[i].flags & GRANGE_FLAG_ISANNO))
	    continue;

	if (xx->r[i].mqual /* obj_type */ != GT_Contig)
	    continue;

	if (sp < xx->displayPos)
	    sp = xx->displayPos;

	for (j = sp; j <= xx->r[i].end; j++) {
	    if (j >= xx->displayPos + xx->displayWidth)
		break;
	    ink[j-xx->displayPos].sh |= sh_inverse;
	}
    }

    XawSheetPutJazzyText(&xx->ed->sw, 0, xx->y_cons, wid,
			 xx->displayedConsensus, ink);
}

/*
 * Returns the other editor in a join-editor pair.
 *         NULL if not joined.
 */
edview *linked_editor(edview *xx) {
    if (!xx->link)
	return NULL;
    return xx->link->xx[0] == xx ? xx->link->xx[1] : xx->link->xx[0];
}

static int tk_redisplaySeqDiffs(edview *xx) {
    char diff[MAX_DISPLAY_WIDTH+1];
    edview *xx0, *xx1;
    int i;

    xx0 = xx->link->xx[0];
    xx1 = xx->link->xx[1];

    for (i = 0; i < xx->displayWidth; i++) {
	diff[i] = xx0->displayedConsensus[i] == xx1->displayedConsensus[i]
	    ? ' ' : '!';
    }
    XawSheetPutText(&xx->link->diffs->sw, 0, 0,
		    xx->displayWidth, diff);
    return 0;
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
	ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight,
			   xx->max_height);

	/* No need to redraw the consenus/numbers */
	xx->refresh_flags |= ED_DISP_NAMES | ED_DISP_SEQS | ED_DISP_CURSOR |
	    ED_DISP_SELECTION;
    }

    if (xx->refresh_flags & ED_DISP_XSCROLL) {
	ed_set_xslider_pos(xx, xx->displayPos);

	/* Changing X may also change height */
	if (!(xx->refresh_flags & ED_DISP_YSCROLL))
	    ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight,
			       xx->max_height);

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
	int i;
	y = -1;
	for (i = 0; i < nr; i++) {
	    if (r[i].rec == xx->cursor_rec) {
		y = r[i].y + xx->y_seq_start - xx->displayYPos;
		break;
	    }
	}
	
	if (y < xx->y_seq_start || y >= xx->displayHeight) {
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

    /* X position */
    if (!x_safe) {
	int w = xx->displayWidth > 10 ? 10 : xx->displayWidth;
	if (xx->cursor_apos < xx->displayPos) {
	    set_displayPos(xx, xx->cursor_apos + 1 - w);
	    do_x = 1;
	}
	if (xx->cursor_apos >= xx->displayPos + xx->displayWidth) {
	    set_displayPos(xx, xx->cursor_apos - xx->displayWidth + w);
	    do_x = 1;
	}
    }

    if (do_x)
	y_safe = 0;

    /* Y position */
    if (!y_safe && xx->cursor_type != GT_Contig) {
	int i;
	int sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;
	
	edview_visible_items(xx, xx->displayPos,
			     xx->displayPos + xx->displayWidth);

	for (i = 0; i < xx->nr; i++) {
	    if (xx->r[i].rec == xx->cursor_rec) {
		y_pos = xx->r[i].y;
		break;
	    }
	}

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
	ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight,
			   xx->max_height);

    if (do_x || do_y) {
	xx->refresh_flags = ED_DISP_ALL;
	edview_redraw(xx);
	return 1;
    }

    return 0;
}

int set_displayPos(edview *xx, int pos) {
    char buf[100];
    int i, ret = 0;
    int delta = pos - xx->displayPos;
    edview *xx2[2];

    if (xx->link && xx->link->locked)
	xx = xx->link->xx[0];

    for (i = 0; i < 2; i++) {
	xx2[i] = xx;

	if (!xx)
	    break;

	xx->displayPos += delta;

	sprintf(buf, "%d", pos);
	Tcl_SetVar2(xx->interp, xx->edname, "displayPos", buf,
		    TCL_GLOBAL_ONLY);

	xx->refresh_flags = ED_DISP_XSCROLL;
	if (i == 1)
	    xx->refresh_flags |= ED_DISP_NO_DIFFS;

	xx = (xx->link && xx->link->locked)
	    ? xx->link->xx[1] : NULL;
    }

    if (xx2[0]->link)
	xx2[0]->link->lockOffset =
	    xx2[0]->link->xx[1]->displayPos - xx2[0]->link->xx[0]->displayPos;


    if (xx2[1]) ret |= edview_redraw(xx2[1]);
    if (xx2[0]) ret |= edview_redraw(xx2[0]);

    return ret;
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
	sequence_get_position(xx->io, xx->cursor_rec, &cnum, &cpos, NULL,
			      NULL);
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
    int j;
    int cpos = xx->cursor_apos;

    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (xx->nr == 0)
	return 0;

    /* Find the current sequence number */
    if (xx->cursor_type == GT_Contig) {
	j = xx->nr;
    } else {
	for (j = 0; j < xx->nr; j++) {
	    if (xx->r[j].rec == xx->cursor_rec) {
		break;
	    }
	}
    }

    /* Step up until we find something overlapping */
    for (j--; j >= 0; j--) {
	if (xx->r[j].start <= cpos && xx->r[j].end >= cpos) {
	    xx->cursor_type = GT_Seq;
	    xx->cursor_pos = cpos - xx->r[j].start;
	    xx->cursor_rec = xx->r[j].rec;
	    break;
	}
    }

    /* Otherwise we've hit the consensus */
    if (j < 0) {
	xx->cursor_type = GT_Contig;
	xx->cursor_rec = xx->cnum;
	xx->cursor_pos = cpos;
    }

    if (!showCursor(xx, 1, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorDown(edview *xx) {
    int j;
    int cpos = xx->cursor_apos;

    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (xx->nr == 0)
	return 0;

    /* Find the current sequence number */
    if (xx->cursor_type == GT_Contig) {
	cpos = xx->cursor_pos;
	j = -1;
    } else {
	for (j = 0; j < xx->nr; j++) {
	    if (xx->r[j].rec == xx->cursor_rec) {
		cpos = xx->r[j].start + xx->cursor_pos;
		break;
	    }
	}
    }

    /* Step up until we find something overlapping */
    for (j++; j < xx->nr; j++) {
	if (xx->r[j].start <= cpos && xx->r[j].end >= cpos) {
	    xx->cursor_type = GT_Seq;
	    xx->cursor_pos = cpos - xx->r[j].start;
	    xx->cursor_rec = xx->r[j].rec;
	    break;
	}
    }

    /* Otherwise we've hit the consensus */
    if (j >= xx->nr) {
	xx->cursor_type = GT_Contig;
	xx->cursor_rec = xx->cnum;
	xx->cursor_pos = cpos;
    }

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
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* Deal with ED_DISP_XSCROLL and ED_DISP_YSCROLL events */
    if (xx->refresh_flags & (ED_DISP_XSCROLL | ED_DISP_YSCROLL))
	tk_redisplaySeqScroll(xx, xx->r, xx->nr);

    /* Consensus emacs-style edit-status ----, -%%-, -**- */
    //tk_redisplaySeqEditStatus(xx);

    /* Redraw the consensus and/or numbers */
    if (xx->refresh_flags & ED_DISP_CONS) {
	tk_redisplaySeqConsensus(xx, xx->r, xx->nr);
    }
    if (xx->refresh_flags & ED_DISP_RULER) {
	tk_redisplaySeqNumbers(xx);
    }

    /* Redraw the main sequences or names section */
    if (xx->refresh_flags & (ED_DISP_SEQS  | ED_DISP_SEQ |
			     ED_DISP_NAMES | ED_DISP_NAME)) {
	tk_redisplaySeqSequences(xx, xx->r, xx->nr);
    }

    /* Editor cursor position */
    if (xx->refresh_flags & ED_DISP_CURSOR) {
	tk_redisplayCursor(xx, xx->r, xx->nr);
    }

#if 0
    /* We've already computed them, but now we actually display them */
    if (xx->refresh_flags & ED_DISP_STATUS) {
	tk_redisplaySeqStatusDisplay(xx);
    }
#endif

    /* Underlining for current selection */
    if (xx->refresh_flags & ED_DISP_SELECTION) {
	redisplaySelection(xx);
    }

    if (inJoinMode(xx) && !(xx->refresh_flags & ED_DISP_NO_DIFFS))
	tk_redisplaySeqDiffs(xx);

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
 * 'name' is a boolean which when true indicates the row,col are in the
 * names panel instead of the sequence panel.
 *
 * Returns the item type GT_* on success and the record/pos in *rec, *pos
 *         -1 on failure (eg numbers, off screen, etc)
 */
int edview_item_at_pos(edview *xx, int row, int col, int name, int exact,
		       int *rec, int *pos) {
    int i;
    int type = -1;
    int best_delta = INT_MAX;
    //    int exact = (name && xx->ed->stack_mode) || !name;

    /* Special case - the reserve row numbers */
    if (row == xx->y_cons) {
	*rec = xx->cnum;
	*pos = col + xx->displayPos;
	return GT_Contig;
    }

    if (row < xx->y_seq_start)
	return -1;
    
    /* A sequence, so find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* Inefficient, but just a copy from tk_redisplaySeqSequences() */
    for (i = 0; i < xx->nr; i++) {
	if (xx->r[i].y + xx->y_seq_start - xx->displayYPos == row) {
	    int delta;

	    /* Find distance from object, in X */
	    if (xx->ed->stack_mode && name) {
		/* In names display during stacking mode */
		int p1 = xx->r[i].start - xx->displayPos;
		int p2 = xx->r[i].end   - xx->displayPos;
		double nc = xx->names->sw.columns;

		if (p1 < 0) p1 = 0;
		p1 = p1 * (nc / xx->displayWidth);
		if (p2 < 0) p2 = 0;
		p2 = p2 * (nc / xx->displayWidth);
		if (col >= p1 && col < p2)
		    delta = 0;
		else
		    delta = INT_MAX;
	    } else {
		/* In sequence display, or only 1 seq per line */
		if (col + xx->displayPos < xx->r[i].start)
		    delta = xx->r[i].start - (col + xx->displayPos);
		else if (col + xx->displayPos > xx->r[i].end)
		    delta = col + xx->displayPos - xx->r[i].end;
		else {
		    delta = 0;
		}
	    }

	    /* And if this is closest match, use it */
	    if (best_delta > delta) {
		best_delta = delta;
		*rec = xx->r[i].rec;
		*pos = col + xx->displayPos - xx->r[i].start;
		type = GT_Seq;
	    }

	    if (!delta)
		/* Can't get better than perfect hit */
		break;
	}
    }

    return !exact || best_delta == 0 ? type : -1;
}

/*
 * More or less the opposite of the above function. Given a record and type
 * this returns the row number on screen. If non-NULL xmin and xmax are the
 * X coordinates of the extends of this record. These may be beyond the bounds
 * of the window.
 *
 * Returns Y coordinate (and optionally min/max X coordinates) if found
 *        -1 if not (with xmin/xmax unset).
 */
int edview_row_for_item(edview *xx, int rec, int *xmin, int *xmax) {
    int i, r = -1;

    if (rec == xx->cnum) {
	if (xmin) *xmin = -xx->displayPos;
	if (xmax) *xmax = -xx->displayPos;
	return 0;
    }

    /* A sequence, so find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* And search for rec in this list */
    for (i = 0; i < xx->nr; i++) {
	if (xx->r[i].rec == rec) {
	    if (xmin) *xmin = xx->r[i].start - xx->displayPos;
	    if (xmax) *xmax = xx->r[i].end   - xx->displayPos;
	    r = xx->r[i].y + xx->y_seq_start - xx->displayYPos;
	    break;
	    
	}
    }
    
    return r >= xx->y_seq_start ? r : -1;
}


int inJoinMode(edview *xx) {
    return xx->link ? 1 : 0;
}

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
	puts("FIXME: reuse existing cache of items");
	r = contig_seqs_in_range(xx->io, &xx->contig,
				 xx->cursor_apos, xx->cursor_apos,
				 CSIR_SORT_BY_X, &nr);

	for (i = 0; i < nr; i++) {
	    s = get_seq(xx->io, r[i].rec);
	    tman_manage_trace("ANY", sequence_get_name(&s), xx->cursor_pos,
			      0, 0, /* left/right clips */
			      s->len < 0, /* complemented */
			      1, /* base spacing */
			      sequence_get_name(&s),
			      xx, r[i].rec, 0, 0);
	}
	free(r);
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
    int *r = NULL, p;

    if (!s)
	return NULL;

    /* FIXME: support s->template_rec */

    /* Solexa data is easy: we have just one other end */
    p = sequence_get_pair(xx->io, s);
    if (p > 0) {
	*nrec = 1;
	r = (int *)malloc(sizeof(*r));
	*r = p;
    } else {
	*nrec = 0;
    }

    cache_decr(xx->io, s);

    return r;
}


/* ---------------------------------------------------------------------- */
/* Selection aka cut'n'paste handling code */

/*
 * (Un)draws the selection - toggles it on or off
 */
static void toggle_select(edview *xx, int seq, int from_pos, int to_pos) {
    int row, xmin, xmax;

    if (from_pos > to_pos) {
	int temp = from_pos;
	from_pos = to_pos;
	to_pos = temp;
    }

    /* Find out the X and Y coord, and exit now if it's not visible */
    if (-1 == (row = edview_row_for_item(xx, seq, &xmin, NULL)))
	return;
    
    /* Convert xmin/xmax to the region we wish to view */
    xmin += from_pos;
    xmax = xmin + to_pos - from_pos;

    /* clip to screen */
    if (xmin < 0) xmin = 0;
    if (xmax >= xx->displayWidth) xmax = xx->displayWidth-1;

    if (xmin > xmax)
	return;

    /* Toggle the line */
    XawSheetOpHilightText(&xx->ed->sw, xmin, row, xmax-xmin+1,
			  sh_select, HOP_TOG);
}

static void redisplaySelection(edview *xx) {
    toggle_select(xx, xx->select_seq, xx->select_start, xx->select_end);
}

/*
 * Callback from Tk_OwnSelection().
 */
static void edSelectionLost(ClientData cd) {
    edview *xx = (edview *)cd;

    /* Undisplay the selection */
    redisplaySelection(xx);

    xx->select_made = 0;
    xx->select_seq = 0;
    xx->select_start = 0;
    xx->select_end = 0;
}

int edSelectClear(edview *xx) {
    if (xx->select_made && EDTKWIN(xx->ed))
	Tk_ClearSelection(EDTKWIN(xx->ed), XA_PRIMARY);
    edSelectionLost((ClientData)xx);

    return 0;
}

void edSelectFrom(edview *xx, int pos) {
    /* Undisplay an old selection */
    if (xx->select_made)
	redisplaySelection(xx);
    else
	xx->select_made = 1;

    /* Set start/end */
    xx->select_seq = xx->cursor_rec;
    pos += xx->displayPos;
    if (xx->select_seq != xx->cnum) {
	int cnum, cpos;
	seq_t *s = get_seq(xx->io, xx->select_seq);

	sequence_get_position(xx->io, xx->select_seq,
			      &cnum, &cpos, NULL, NULL);
	pos -= cpos;
	if (pos < 0)
	    pos = 0;
	if (pos >= ABS(s->len))
	    pos = ABS(s->len)-1;
    }
    xx->select_start = xx->select_end = pos;

    Tk_OwnSelection(EDTKWIN(xx->ed), XA_PRIMARY, edSelectionLost,
		    (ClientData)xx);

    /* Display new selection */
    redisplaySelection(xx);
}

void edSelectTo(edview *xx, int pos) {
    if (!xx->select_made) {
	edSelectFrom(xx, pos);
    }

    /* Undisplay old selection */
    redisplaySelection(xx);

    /* Set start/end */
    pos += xx->displayPos;
    if (xx->select_seq != xx->cnum) {
	int cnum, cpos;
	seq_t *s = get_seq(xx->io, xx->select_seq);

	sequence_get_position(xx->io, xx->select_seq,
			      &cnum, &cpos, NULL, NULL);
	pos -= cpos;
	if (pos < 0)
	    pos = 0;
	if (pos >= ABS(s->len))
	    pos = ABS(s->len)-1;
    }
    xx->select_end = pos;

    /* Display new selection */
    redisplaySelection(xx);
}

/*
 * Automatically called when X wishes to obtain a selection. We register
 * this procedure in the initialise code of tkEditor.c.
 *
 * Return codes expected:
 *    -1  Failure
 *    >0  Number of bytes
 */
int edGetSelection(ClientData clientData, int offset, char *buffer,
		   int bufsize) {
    Editor *ed = (Editor *)clientData;
    int start, end, len;
    edview *xx = ed->xx;

    /* Do we have a selection? */
    if (!xx->select_made)
	return -1;

    start = xx->select_start + offset;
    end   = xx->select_end;

    if (start > end) {
	len = start;
	start = end;
	end = len;
    }

    len = end - start+1 > bufsize ? bufsize : end - start + 1;

    if (len && xx->select_seq) {
	if (xx->select_seq != xx->cnum) {
	    seq_t *s, *sorig;
	    sorig = s = get_seq(xx->io, xx->select_seq);
	    if (s->len < 0) {
		s = dup_seq(s);
		complement_seq_t(s);
	    }
	    memcpy(buffer, s->seq+start, len);
	    if (s != sorig)
		free(s);
	} else {
	    calculate_consensus_simple(xx->io, xx->cnum, start, start+len-1,
				       buffer, NULL);
	}
    }

    return len;
}
