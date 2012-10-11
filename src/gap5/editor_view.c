#include <xalloc.h>
#include <math.h>
#include <ctype.h>
#include <X11/Xatom.h> /* XA_PRIMARY - included in Tk distribrution */
#include <assert.h>

#include "editor_view.h"
#include "tkSheet.h"
#include "tman_interface.h"
#include "gap_globals.h"
#include "qualIO.h"
#include "qual.h"
#include "tg_gio.h"
#include "misc.h"
#include "consensus.h"
#include "tagdb.h"
#include "active_tags.h"
#include "io_utils.h"
#include "depad_seq_tree.h"
#include "align.h"
#include "align_lib.h"

/*
 * Use this if you wish to make the cached consensus visible as its own
 * sequence, for debugging only.
 *
 * #define CACHED_CONS_VISIBLE
 */

/*
 * Similarly for debugging, define this if you want the REFPOS markers to be
 * visible. They appear on line 1 as square boxes. Not ideal and it doesn't
 * reveal the type (I/D) or size of deletion, but it's still useful debugging
 * data.
 *
 * #define REFPOS_VISIBLE
 */
#define REFPOS_VISIBLE

static void redisplaySelection(edview *xx);

/*
 * Global hash of active edview widgets, used for finding existing editors
 * for remote control.
 */
static HacheTable *edview_hash = NULL;

/*
 * A C interface to the edit_contig and join_contig Tcl functions.
 */
int edit_contig(GapIO *io, tg_rec cnum, tg_rec rnum, int pos) {
    char cmd[1024];

    sprintf(cmd, "edit_contig -io %s -contig %"PRIrec
	    " -reading #%"PRIrec" -pos %d\n",
	    io_obj_as_string(io), cnum, rnum, pos);
    return Tcl_Eval(GetInterp(), cmd);
}

int join_contig(GapIO *io, tg_rec cnum[2], tg_rec rnum[2], int pos[2]) {
    char cmd[1024];
    int ret;

    sprintf(cmd, "join_contig -io %s -contig %"PRIrec" -reading #%"PRIrec
	    " -pos %d -contig2 %"PRIrec" -reading2 #%"PRIrec" -pos2 %d",
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
edview *edview_new(GapIO *io, tg_rec contig, tg_rec crec, int cpos,
		   Editor *ed, edNames *names,
		   void (*dispFunc)(void *, int, int, int, void *),
		   Tcl_Interp *interp)
{
    edview *xx;
    static int editor_id = 1;
    char *cp;

    xx = (edview *)xcalloc(1, sizeof(*xx));
    if (!xx)
	return NULL;
    
    xx->editor_id = editor_id++;

    xx->interp = interp;
    xx->io = io; /* model */
    xx->cnum = contig;

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
    xx->cursor_rec  = crec ? crec : contig;
    xx->cursor_type = (xx->cursor_rec == 0 || xx->cursor_rec == contig)
	? GT_Contig : GT_Seq;

    xx->trace_lock = 1;
    
    if (!xx->ed->consensus_at_top) {
	xx->ed->sw.yflip = 1;
	xx->names->sw.yflip = 1;
    }

    xx->r = NULL;
    xx->anno_hash = NULL;
    xx->rec_hash = NULL;

    /* Private cursor */
    cp = Tcl_GetVar2(xx->interp, Tk_PathName(xx->ed->sw.tkwin), "reg",
		     TCL_GLOBAL_ONLY);
    xx->reg_id = cp ? atoi(cp) : 0;
    if (io->base)
	xx->cursor = create_contig_cursor(gio_base(io), contig, 1, xx->reg_id);
    edSetApos(xx);
    xx->displayPos = xx->cursor_apos;
    
    edview_set_sort_order(xx);

    /* Add to our global hash table */
    {
	HacheData hd;

	if (!edview_hash)
	    edview_hash = HacheTableCreate(16, HASH_DYNAMIC_SIZE);
	
	hd.p = xx;
	HacheTableAdd(edview_hash, (char *)&contig, sizeof(tg_rec), hd, NULL);
    }

    xx->trace_hash = HacheTableCreate(256, HASH_DYNAMIC_SIZE);
    
    return xx;
}

/*
 * Deallocates an edview
 */
void edview_destroy(edview *xx) {
    xx->editorState = StateDown;
    if (xx->link) {
	xx->link->xx[0]->editorState = StateDown;
	xx->link->xx[1]->editorState = StateDown;

	/* Set other_xx->link = NULL */
	xx->link->xx[xx == xx->link->xx[0]]->link = NULL;
	free(xx->link);
	xx->link = NULL;
    }

    if (xx->cursor)
	delete_contig_cursor(gio_base(xx->io), xx->cnum, xx->cursor->id, 1);

    if (xx->r)
	free(xx->r);

    if (xx->anno_hash)
	HacheTableDestroy(xx->anno_hash, 0);
    
    if (xx->rec_hash)
	HacheTableDestroy(xx->rec_hash, 0);

    if (xx->trace_hash) {
	HacheIter *iter = HacheTableIterCreate();
	HacheItem *hi;

	while (hi = HacheTableIterNext(xx->trace_hash, iter)) {
	    if (hi->data.p)
		read_deallocate(hi->data.p);
	}
	HacheTableDestroy(xx->trace_hash, 0);
	HacheTableIterDestroy(iter);
    }

    HacheTableRemove(edview_hash, (char *)&xx->cnum, sizeof(tg_rec), 0);

    xfree(xx);
}

/*
 * Finds an existing editor widget for a specific contig.
 * Returns edview on success
 *         NULL on failure
 */
edview *edview_find(GapIO *io, tg_rec contig) {
    HacheIter *iter;
    HacheItem *hi;

    if (!edview_hash)
	return NULL;

    iter = HacheTableIterCreate();
    while (hi = HacheTableIterNext(edview_hash, iter)) {
	edview *xx = (edview *)hi->data.p;
	if (!xx->link && xx->cnum == contig)
	    return xx;
    }
    HacheTableIterDestroy(iter);

    return NULL;
}

static seq_t *get_seq(GapIO *io, tg_rec rec) {
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

static void add_number64(char *buf, int *j, int l1, int l2, int64_t val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*"PRId64, l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*"PRId64, l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*"PRId64, l2, val);
	else
	    *j += sprintf(buf + *j, "%"PRId64, val);
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

static void add_char(char *buf, int *j, int l1, int l2, char chr) {
    if (l1)
	*j += sprintf(buf + *j, "%*c", l1, chr);
    else
	*j += sprintf(buf + *j, "%c", chr);
}

/*
 * Formats tag information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %p	Tag position
 * %t	Tag type (always 4 characters)
 * %l	Tag length
 * %#	Tag number (0 if unknown)
 * %c	Tag comment
 * %d   Tag direction
 *
 * Additionally, some formats (p, l, n and c) can be specified as
 * %<number><format> (eg %80c) to allow AT MOST that many characters.
 */
char *edGetBriefTag(edview *xx, tg_rec anno_ele, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    GapIO *io = xx->io;
    anno_ele_t *e;

    if (!anno_ele)
	return "";

    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, anno_ele);
    
    for (i = j = 0; format[i]; i++) {
	char type[5];

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

	case 't': /* Type */
	    (void)type2str(e->tag_type, type);
	    status_buf[j++] = type[0];
	    status_buf[j++] = type[1];
	    status_buf[j++] = type[2];
	    status_buf[j++] = type[3];
	    break;

	case 'p': { /* Position */
	    range_t *r = anno_get_range(io, anno_ele, NULL, 0);
	    add_number(status_buf, &j, l1, l2, r->start);
	    break;
	}

	case 'l': { /* Length */
	    range_t *r = anno_get_range(io, anno_ele, NULL, 0);
	    add_number(status_buf, &j, l1, l2, r->end - r->start + 1);
	    break;
	}

	case '#': /* Number */
	    add_number64(status_buf, &j, l1, l2, e->rec);
	    break;

	case 'd':
	    add_char(status_buf, &j, l1, l2, e->direction);
	    break;

	case 'c': /* Comment */
	    add_string(status_buf, &j, l1, l2,
		       e->comment ? e->comment : "(no comment)");
	    break;

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return status_buf;
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
 * %V   Vendor/platform
 * %B   Bin record number
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
char *edGetBriefSeq(edview *xx, tg_rec seq, int pos, char *format) {
    static char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    GapIO *io = xx->io;
    seq_t *s1 = get_seq(io, seq), *s2 = NULL, *s = s1;
    tg_rec pair = 0;

    cache_incr(io, s1);
    
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
	    if (pair > 0 && !s2) {
		if ((s2 = get_seq(io, pair)) != NULL) {
		    cache_incr(io, s2);
		    cache_decr(io, s1);
		}
	    }
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
	    add_number64(status_buf, &j, l1, l2, s->rec);
	    break;

	case 'n':
	    if (raw)
		add_number64(status_buf, &j, l1, l2, s->rec);
	    else
		add_string(status_buf, &j, l1, l2, s->name);
	    break;

	case 'B':
	    add_number64(status_buf, &j, l1, l2, s->bin);
	    break;

	case 'p': {
	    tg_rec cnum;
	    int cpos;
	    if (0 == sequence_get_position(xx->io, s->rec, &cnum, &cpos, NULL,
					   NULL)) {
		
		if (raw || cnum == xx->cnum) {
		    add_number(status_buf, &j, l1, l2, cpos);
		} else {
		    char buf[1024];
		    sprintf(buf, "%d@%s", cpos, get_contig_name(io, cnum));
		    add_string(status_buf, &j, l1, l2, buf);
		}
	    }
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

	case 'S': {
	    int orient = sequence_get_orient(xx->io, s->rec);
	    if (raw)
		add_number(status_buf, &j, l1, l2, orient);
	    else
		add_string(status_buf, &j, l1, l2, orient ? "<<" : ">>");
	    break;
	}

	case 'd': {
	    range_t *r = sequence_get_range(xx->io, s);
	    int end = r
		? ((r->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD
		   ? 0 : 1)
		: 0;

	    if (raw)
		add_number(status_buf, &j, l1, l2, end);
	    else
		add_string(status_buf, &j, l1, l2, end ? "-" : "+");
	    break;
	}

	case 'b':
	    if (pos >= 0 && pos < ABS(s->len)) {
		char base[2];
		int cut;
		sequence_get_base(xx->io, &s, pos, &base[0], NULL, &cut, 1);
		base[1] = 0;
		if (cut)
		    base[0] = tolower(base[0]);
		else
		    base[0] = toupper(base[0]);
		add_string(status_buf, &j, l1, l2, base);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'c':
	    if (pos >= 0 && pos < ABS(s->len)) {
		int q;
		sequence_get_base(xx->io, &s, pos, NULL, &q, NULL, 1);
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
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
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
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
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
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
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
		sequence_get_base4(xx->io, &s, pos, NULL, q, NULL, 1);
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

	case 'V':
	    if (raw) {
		add_number(status_buf, &j, l1, l2, s->seq_tech);
	    } else {
		switch(s->seq_tech) {
		case STECH_SANGER:
		    add_string(status_buf, &j, l1, l2, "Sanger");
		    break;
		case STECH_SOLEXA:
		    add_string(status_buf, &j, l1, l2, "Illumina");
		    break;
		case STECH_SOLID:
		    add_string(status_buf, &j, l1, l2, "SOLiD");
		    break;
		case STECH_454:
		    add_string(status_buf, &j, l1, l2, "454");
		    break;
		default:
		    add_string(status_buf, &j, l1, l2, "unknown");
		    break;
		}
	    }
	    break;
	    
	default:
	    status_buf[j++] = format[i];
	}
    }
 bail_out:
    status_buf[j] = 0;

    cache_decr(io, s);

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
 * %p	padded position
 * %r	reference position
 * %l	Length
 * %s	Start of clip
 * %e	End of clip
 * %b   Base call
 * %c   Base confidence log-odds (raw for probability value)
 * %A   A confidence log-odds (raw for probability value)
 * %C   C confidence log-odds (raw for probability value)
 * %G   G confidence log-odds (raw for probability value)
 * %T   T confidence log-odds (raw for probability value)
 * %H   Het confidence log-odds
 * %h	Het call
 * %d   Depth
 * %*   * (gap) confidence
 *
 * Additionally specifying %<number><format> forces AT MOST that many
 * characters to be displayed.
 * Specifying %R<format> (or %<number>R<format>) requests the raw data to
 * be displayed. This only works for some fields. Eg %Rp displays 0 to 4, but
 * %p displays, for instance, "forward universal"
 */
char *edGetBriefCon(edview *xx, tg_rec crec, int pos, char *format) {
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
	    add_number64(status_buf, &j, l1, l2, crec);
	    break;

	case 'n':
	    if (raw) {
		add_number64(status_buf, &j, l1, l2, crec);
	    } else {
		contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
		add_string(status_buf, &j, l1, l2,
			   contig_get_name(&c));
	    }
	    break;

	case 'p': {
	    add_number(status_buf, &j, l1, l2, pos);
	    break;
	}

	case 'r': {
	    add_number(status_buf, &j, l1, l2,
		       padded_to_reference_pos(xx->io, xx->cnum, pos,
					       NULL, NULL));
	    break;
	}

	case 'l': {
	    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
	    add_number(status_buf, &j, l1, l2, contig_get_length(&c));
	    break;
	}

	case 's': {
	    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
	    add_number(status_buf, &j, l1, l2, contig_get_start(&c));
	    break;
	}

	case 'e': {
	    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
	    add_number(status_buf, &j, l1, l2, contig_get_end(&c));
	    break;
	}

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

	case 'H':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		double q = xx->cachedConsensus[p].scores[6];
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

	case 'h':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		int h = xx->cachedConsensus[p].het_call;
		char str[3];
		str[0] = "ACGT*"[h/5];
		str[1] = "ACGT*"[h%5];
		str[2] = 0;
		add_string(status_buf, &j, l1, l2, str);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'd':
	    if (pos >= xx->displayPos &&
		pos < xx->displayPos + xx->displayWidth) {
		int p = pos - xx->displayPos;
		int d = xx->cachedConsensus[p].depth;
		//d = get_uniqueness_pos(xx->displayedConsensus,
		//		       xx->displayWidth,
		//		       p);
		add_number(status_buf, &j, l1, l2, d);
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
    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
    int i;
    int mode = xx->ed->stack_mode
	? CSIR_ALLOCATE_Y_MULTIPLE
	: CSIR_ALLOCATE_Y_SINGLE;

    if (NULL == c) return -1;

    /* sort... */
    mode |= CSIR_DEFAULT;
    
    /* Always reload for now as we can't spot edits yet */
    if (xx->r && xx->r_start == start && xx->r_end == end)
	return 0;
	
    /* Query sequences */
    if (xx->r)
	free(xx->r);
	
    xx->r_start = start;
    xx->r_end = end;
    xx->r = contig_items_in_range(xx->io, &c, start, end,
				  CSIR_SORT_BY_Y | mode, CSIR_DEFAULT,
				  &xx->nr);
    if (!xx->r) {
	xx->nr = 0;
	return -1;
    }

    if (xx->rec_hash) {
	HacheTableDestroy(xx->rec_hash, 0);
    }

    xx->rec_hash = HacheTableCreate(8192, HASH_DYNAMIC_SIZE);
    if (NULL == xx->rec_hash) return -1;
    xx->rec_hash->name = "rec_hash";

    /* Work out Y dimension */
    xx->max_height = 0;
    for (i = 0; i < xx->nr; i++) {
	HacheData hd;
	tg_rec key = xx->r[i].rec;

	if (xx->max_height < xx->r[i].y)
	    xx->max_height = xx->r[i].y;

	hd.i = i;
	if (!HacheTableAdd(xx->rec_hash, (char *)&key, sizeof(key), hd, NULL))
	    return -1;
    }
    xx->max_height += 3; /* +1 for from 0, +2 for consensus+ruler */

    /* Fast map of annotations to sequences */
    if (xx->anno_hash)
	HacheTableDestroy(xx->anno_hash, 0);

    xx->anno_hash = HacheTableCreate(8192, HASH_DYNAMIC_SIZE |
				     HASH_ALLOW_DUP_KEYS);
    if (!xx->anno_hash) return -1;
    xx->anno_hash->name = "anno_hash";
    for (i = 0; i < xx->nr; i++) {
	tg_rec key = xx->r[i].pair_rec; /* aka obj_rec */
	HacheData hd;

	if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)
	    continue;

	/*
	 * Work around a bug (or design flaw?) in break_contig.
	 * When breaking a contig we need to reparent our tags to a new
	 * contig record. This isn't feasible as tags shouldn't "know"
	 * which contigs they're in. Instead we use the flag.
	 */
	if (!(xx->r[i].flags & GRANGE_FLAG_TAG_SEQ))
	    key = xx->cnum;

	hd.i = i;
	if (!HacheTableAdd(xx->anno_hash, (char *)&key, sizeof(key), hd, NULL))
	    return -1;
    }

    /* Reverse the order of annotations in the anno_hash. */
    HacheTableReverse(xx->anno_hash);

#if 0
    puts("");
    for (i = 0; i < xx->nr; i++) {
	tg_rec rec;

	if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    anno_ele_t *a = (anno_ele_t *)cache_search(xx->io, GT_AnnoEle,
						       xx->r[i].rec);
	    rec = a->rec;
	} else {
	    seq_t *s = (seq_t *)cache_search(xx->io, GT_Seq,
					     xx->r[i].rec);
	    rec = s->rec;
	}
	printf("%d\t%d\t%s%d/%d\t%d\t%s\t%d..%d\n",
	       i, xx->r[i].y,
	       xx->r[i].rec == rec ? "" : "*", xx->r[i].rec, rec,
	       ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
	           ? xx->r[i].pair_rec : 0,
	       ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
	           ? "tag" : "seq",
	       xx->r[i].start, xx->r[i].end);
    }
#endif


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

/*
 * Find the largest i such that r[i].y <= y.
 */
static int edview_binary_search_y(rangec_t *r, int nr, int y) {
    int i_start, i_end, i_mid;

    if (nr <= 0)
	return 0;

    i_start = 0;
    i_end = nr;

    while (i_start < i_end) {
	i_mid = (i_end - i_start) / 2 + i_start;

	if (r[i_mid].y < y) {
	    i_start = i_mid+1;
	} else {
	    i_end = i_mid;
	}
    }

    return i_mid;
}

static int ed_set_xslider_pos(edview *xx, int offset) {
    char buf[100];
    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
    double len = contig_get_length(&c);

    offset -= contig_get_start(&c);

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
	    printf("Error in editor names scroll: %s\n", Tcl_GetStringResult(EDINTERP(en)));
	}
    }
}

static void tk_redisplaySeqTags(edview *xx, XawSheetInk *ink, seq_t *s,
				int sp, int left, int right) {
    char type[5];
    HacheItem *hi;

    if (xx->ed->hide_annos)
	return;

    if (xx->nr == 0)
	return;

    /* Identify overlapping tags and mark using boxes */
#define sh_tmp (1L<<30)

    if (s) {
	/* A sequence */
	tg_rec key = s->rec;
	int p, p2, l = s->len >= 0 ? s->len : -s->len;

	if (l > xx->displayWidth - (sp-xx->displayPos))
	    l = xx->displayWidth - (sp-xx->displayPos);

	for (hi = HacheTableSearch(xx->anno_hash,
			       (char *)&key, sizeof(key));
	     hi; hi = HacheTableNext(hi, (char *)&key, sizeof(key))) {
	    int ai = hi->data.i;
	    int db = idToIndex(type2str(xx->r[ai].mqual, type));

	    /* FIXME: Inefficient! */
	    for (p2 = sp - xx->displayPos, p = 0; p<l; p++,p2++) {
		if (p2 + xx->displayPos >= xx->r[ai].start &&
		    p2 + xx->displayPos <= xx->r[ai].end) {
		    if (xx->ed->display_cutoffs ||
			(p >= left-1 && p <= right-1)) {
			if (ink[p2].sh & sh_tmp)
			    ink[p2].sh |= sh_box;
			ink[p2].sh |= sh_tmp;
			if (tag_db[db].fg_colour!=NULL) {
			    ink[p2].sh|=sh_fg;
			    ink[p2].fg=tag_db[db].fg_pixel;
			}
			if (tag_db[db].bg_colour!=NULL) {
			    ink[p2].sh|=sh_bg;
			    ink[p2].bg=tag_db[db].bg_pixel;
			}
		    }
		}
	    }
	}
    } else {
	/* Consensus */
	int j, p2;
	tg_rec key = xx->cnum;

	for (hi = HacheTableSearch(xx->anno_hash,
			       (char *)&key, sizeof(key));
	     hi; hi = HacheTableNext(hi, (char *)&key, sizeof(key))) {
	    int ai = hi->data.i;
	    int db = idToIndex(type2str(xx->r[ai].mqual, type));

	    sp = xx->r[ai].start;
	    if (sp < xx->displayPos)
		sp = xx->displayPos;

	    for (j = sp; j <= xx->r[ai].end; j++) {
		if (j >= xx->displayPos + xx->displayWidth)
		    break;

		p2 = j-xx->displayPos;

		if (ink[p2].sh & sh_tmp)
		    ink[p2].sh |= sh_box;
		ink[p2].sh |= sh_tmp;

		if (tag_db[db].fg_colour!=NULL) {
		    ink[p2].sh|=sh_fg;
		    ink[p2].fg=tag_db[db].fg_pixel;
		}
		if (tag_db[db].bg_colour!=NULL) {
		    ink[p2].sh|=sh_bg;
		    ink[p2].bg=tag_db[db].bg_pixel;
		}
	    }
	}
    }
}

/* Returns 1 if rec is in the global "readings" list, 0 if not */
static int seq_in_readings_list(edview *xx, tg_rec rec) {
    char srec[20], list[1024];

    sprintf(list, "NGList_read_hash_%s", xx->ed->output_list);
    sprintf(srec, "#%"PRIrec, rec);
    return Tcl_GetVar2(xx->interp, list, srec, TCL_GLOBAL_ONLY) ? 1 : 0;
}

static void tk_redisplaySeqSequences(edview *xx, rangec_t *r, int nr) {
    int i, j, k, box_alt;

    /*
    sheet_clear(&xx->ed->sw);
    sheet_clear(&xx->names->sw);
    */

    i = edview_binary_search_y(r, nr, xx->displayYPos);

    /* Work down the screen line by line */
    for (j = xx->y_seq_start;
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
	box_alt = 0; /* alternating 0/1 */
	while (i < nr && xx->r[i].y - xx->displayYPos <= j - xx->y_seq_start) {
	    seq_t *s, *sorig;

	    if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO
#ifndef CACHED_CONS_VISIBLE
		|| (r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISCONS
#endif
		|| (r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS
		) {
		i++;
		continue;
	    }

	    if (xx->r[i].y - xx->displayYPos < j - xx->y_seq_start) {
		i++;
		continue;
	    }

	    if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
		int p2 = r[i].start - xx->displayPos;
		ink[p2].sh |= sh_box;
		i++;
		continue;
	    }

	    s = sorig = get_seq(xx->io, r[i].rec);
	    cache_incr(xx->io, sorig);
	    sp = r[i].start;
	    l = s->len > 0 ? s->len : -s->len;
	    seq_p = 0;
	    dir = '>';

	    /* Optimisation for single sequence only */
	    if (xx->refresh_flags & ED_DISP_SEQ &&
		!(xx->refresh_flags & ED_DISP_SEQS)) {
		if (xx->refresh_seq != r[i].rec) {
		    cache_decr(xx->io, sorig);
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

	    if (l > MAX_SEQ_LEN)
		seq = malloc(l);
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
			sequence_get_base(xx->io, &s, p+seq_p, &base, &qual,
					  NULL, 0);
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
			    int qbin = qual / 7;
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
		    tk_redisplaySeqTags(xx, ink, s, sp, left, right);
		}
	    }

	    /* Name */
	    if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME)) {
		int nl = s->name_len - xx->names_xPos;
		int ncol = xx->names->sw.columns;
		XColor **qual_bg = xx->ed->qual_bg;

		box_alt ^= 1;

		if (xx->ed->stack_mode) {
		    int p  = r[i].start - xx->displayPos;
		    int p2 = r[i].end   - xx->displayPos;
		    int bg = -1, t;
		    double nc = xx->names->sw.columns;
		    if (p < 0) p = 0;
		    p = p * (nc / xx->displayWidth);
		    if (p2 < 0) p2 = 0;
		    p2 = p2 * (nc / xx->displayWidth);
		    if (p2 > xx->names->sw.columns)
			p2 = xx->names->sw.columns;
		    while (nline[p] != ' ')
			p++;

		    if (seq_in_readings_list(xx, s->rec)) {
			int ptmp = p;
			do {
			    nink[ptmp].sh |= sh_bold;
			    nink[ptmp++].sh |= box_alt ? sh_box : sh_box_alt;
			} while (ptmp < p2);
			qual_bg = xx->ed->qual_bg2;
			bg = xx->ed->qual_bg2[9]->pixel;
		    }

		    if (xx->ed->display_mapping_quality) {
			int qbin = s->mapping_qual / 10;
			if (qbin < 0) qbin = 0;
			if (qbin > 9) qbin = 9;
			bg = qual_bg[qbin]->pixel;
		    } else {
			t = sequence_get_template_info(xx->io, sorig,
						       NULL, NULL);
			bg = xx->ed->tmpl_bg[t+1]->pixel;
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
		    XColor **qual_bg = xx->ed->qual_bg;
		    int t;

		    nline[0] = dir;
		    if (nl > 0)
			memcpy(&nline[1], s->name + xx->names_xPos, nl);

		    t = sequence_get_template_info(xx->io, sorig, NULL, NULL);
		    nink[0].sh = sh_bg;
		    nink[0].bg = xx->ed->tmpl_bg[t+1]->pixel;

		    if (seq_in_readings_list(xx, s->rec)) {
			qual_bg = xx->ed->qual_bg2;
			for (k = 1; k < ncol && k < MAX_NAME_WIDTH; k++) {
			    nink[k].sh |= sh_bg |
				(box_alt ? sh_box : sh_box_alt);
			    nink[k].bg = qual_bg[9]->pixel;
			}
		    }

		    /*if (xx->ed->display_mapping_quality)*/ {
			int qbin = s->mapping_qual / 10;
			if (qbin < 0) qbin = 0;
			if (qbin > 9) qbin = 9;

			for (k = 1; k < ncol && k < MAX_NAME_WIDTH; k++) {
			    nink[k].sh |= sh_bg;
			    nink[k].bg = qual_bg[qbin]->pixel;
			}
		    }
		}
	    }

	    cache_decr(xx->io, sorig);

	    if (s != sorig)
		free(s);

	    i++;

	    if (seq != seq_a) {
		free(seq);
		seq = seq_a;
	    }
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

    calculate_consensus(xx->io, xx->cnum, pos, pos+wid-1, xx->cachedConsensus);
    for (i = 0; i < wid; i++) {
	xx->displayedConsensus[i] = "ACGT*N "[xx->cachedConsensus[i].call];
    }

    memset(ink, 0, MAX_DISPLAY_WIDTH * sizeof(*ink));
    if (xx->ed->display_quality) {
	int qbin;

	for (i = 0; i < wid; i++) {
	    qbin = xx->cachedConsensus[i].phred/10;
//	    qbin = get_uniqueness_pos(xx->displayedConsensus,
//				      xx->displayWidth,
//				      i);
//	    qbin = 9-2*log(qbin>1?qbin:1);
	    if (qbin < 0) qbin = 0;
	    if (qbin > 9) qbin = 9;
	    ink[i].sh |= sh_bg;
	    ink[i].bg = xx->ed->qual_bg[qbin]->pixel;
	}
    }

    /* Consensus annotations */
    if (xx->anno_hash)
	tk_redisplaySeqTags(xx, ink, NULL, 0, 0, 0);

#ifdef REFPOS_VISIBLE
    for (i = 0; i < nr; i++) {
	int p2;

	if ((r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS)
	    continue;

	p2 = r[i].start - xx->displayPos;
	if (p2 < 0 || p2 >= wid)
	    continue;

	if ((r[i].flags & GRANGE_FLAG_REFPOS_INDEL) == GRANGE_FLAG_REFPOS_INS)
	    ink[p2].sh |= sh_underline | sh_indel;
	else {
	    ink[p2].sh |= sh_caret_l | sh_indel;
	    if (p2 > 0)
		ink[p2-1].sh |= sh_caret_r | sh_indel;
	}
    }
#endif

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
static int generate_ruler(edview *xx, char *ruler, XawSheetInk *ink,
			  int pos, int width) {
    char *k = ruler;
    XawSheetInk *K = ink;
    int j;

    //int padded_pos[MAX_DISPLAY_WIDTH+21];

    memset(ruler, ' ', MAX_DISPLAY_WIDTH+21);
    memset(ink, 0, (MAX_DISPLAY_WIDTH+21) * sizeof(*ink));
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
    if (xx->ed->pos_type == 'P') {
	/* Basic padded coordinate numbering */
	int lower,times;
	lower = (pos - pos%10);
	times = width/10 + 3;
	for (j=0;j<times;j++,k+=10,K+=10,lower+=10) {
	    sprintf(k,"%10d",lower);
	    K[9].sh |= sh_underline;
	}
	return 9+pos%10;
    } else {
	int last_x = -100;
	/* Reference based coordinates, where known */
	int rpos[MAX_DISPLAY_WIDTH+11], rid[MAX_DISPLAY_WIDTH+11], i;

	padded_to_reference_array(xx->io, xx->cnum, xx->displayPos,
				  xx->displayPos+xx->displayWidth+10-1,
				  rpos, rid);

	k += 10;
	K += 10;
	for (i = 0; i < xx->displayWidth+10; i++) {
	    int len = log(ABS(rpos[i] ? rpos[i] : 1)) * 0.4342945;
	    len++;

	    if (rpos[i] % 10 == 0 && rid[i] != -1) {
		if (i - last_x > len) {
		    sprintf(&k[i-(len-1)], "%.*d", len, rpos[i]);
		    k[i+1+(rpos[i]<0)] = ' ';
		    K[i].sh |= sh_underline;
		} else {
		    k[i] = '|';
		    K[i].sh |= sh_underline;
		}
		last_x = i;
	    }
	}
	return 10;
    }
}

static void tk_redisplaySeqNumbers(edview *xx) {
    char ruler[MAX_DISPLAY_WIDTH+21];
    XawSheetInk ink[MAX_DISPLAY_WIDTH+21];
    int off;

    off = generate_ruler(xx, ruler, ink, xx->displayPos, xx->displayWidth);
    XawSheetPutJazzyText(&xx->ed->sw, 0, xx->y_numbers, xx->displayWidth,
    			 &ruler[off], &ink[off]);
}


/* Handle scrolling changes */
static void tk_redisplaySeqScroll(edview *xx, rangec_t *r, int nr) {
    if (xx->refresh_flags & ED_DISP_YSCROLL) {
	ed_set_yslider_pos(xx, xx->displayYPos, xx->displayHeight,
			   xx->max_height);

	/* No need to redraw the consenus/numbers */
	xx->refresh_flags |= ED_DISP_NAMES | ED_DISP_READS | ED_DISP_CURSOR |
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
	tg_rec key;
	HacheItem *hi;
	
	key = xx->cursor_rec;
	if (!xx->rec_hash || !r)
	    return;

	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	y = hi && hi->data.i < nr
	    ? r[hi->data.i].y + xx->y_seq_start - xx->displayYPos
	    : -1;

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
int showCursor(edview *xx, int x_safe, int y_safe) {
    int y_pos = 0;
    int do_x = 0;
    int do_y = 0;

    /* X position */
    if (!x_safe) {
	//int w = xx->displayWidth > 50 ? 50 : xx->displayWidth;
	int w = xx->displayWidth >= 20
	    ? xx->displayWidth-10
	    : xx->displayWidth/2;
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
	tg_rec key;
	int sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;
	HacheItem *hi;
	
	edview_visible_items(xx, xx->displayPos,
			     xx->displayPos + xx->displayWidth);

	key = xx->cursor_rec;
	if (!xx->rec_hash)
	    return 0;
	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	if (!hi)
	    return 0;
	i = hi->data.i;

	if (!xx->r)
	    return 0;

	y_pos = xx->r[i].y;
	if (y_pos == -1) {
	    y_pos = 0; /* tag on consensus */
	    xx->cursor_rec = xx->cnum;
	    xx->cursor_type = GT_Contig;
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
	if (!do_x)
	    xx->refresh_flags &= ~(ED_DISP_CONS | ED_DISP_XSCROLL);
	edview_redraw(xx);
	return 1;
    }

    return 0;
}

/*
 * Sends out a notification of our cursor movement
 */
static void cursor_notify(edview *xx) {
    reg_cursor_notify cn;

    if (!xx->cursor)
	return;

    xx->cursor->seq = xx->cursor_rec;
    xx->cursor->pos = xx->cursor_pos;
    xx->cursor->abspos = xx->cursor_apos;
    xx->cursor->job = CURSOR_MOVE;
    xx->cursor->sent_by = xx->reg_id;
    cn.job = REG_CURSOR_NOTIFY;
    cn.cursor = xx->cursor;
    contig_notify(gio_base(xx->io), xx->cnum, (reg_data *)&cn);
}

/*
 * Returns  1 if seq is visible in the current editor position.
 *          0 if not.
 *
 * It may also fill out the 'new_y' (if non-null) value to indicate
 * the preferred new displayYPos parameter. This is set to -1 when
 * there is no point in changing it as 'seq' will never be visible
 * (it's not overlapping this X coord).
 */
int edview_seq_visible(edview *xx, tg_rec seq, int *new_y) {
    int i, y_pos, vis = 0;
    int sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;
    HacheItem *hi;
	
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (new_y)
	*new_y = xx->displayYPos;

    if (!xx->rec_hash)
	return 0;
    hi = HacheTableSearch(xx->rec_hash, (char *)&seq, sizeof(seq));
    if (!hi || !xx->r)
	return 0;
    i = hi->data.i;

    y_pos = xx->r[i].y;
    vis = 1;
    if (y_pos == -1) {
	/* tag? */
	return 1;
    }

    if (!vis) {
	if (new_y)
	    *new_y = -1;
	return 0;
    }

    /* seq is above, scroll so this is the first row */
    if (y_pos < xx->displayYPos) {
	if (new_y)
	    *new_y = y_pos;
	return 0;
    }

    /* seq is below, scroll so this is the last row */
    if (y_pos >= xx->displayYPos + sheight) {
	if (new_y)
	    *new_y = y_pos - sheight + 1;;
	return 0;
    }

    /* Otherwise it's already on-screen */
    if (new_y)
	*new_y = y_pos;
    return 1;
}

/*
 * This is called by an editor xview method - ie the X scrollbar. It also
 * gets called programmatically on other conditions (such as making the
 * editing cursor visibile).
 *
 * We attempt to make sure that sequences that were previously visible on
 * screen are now also visible on screen (where possible). To achieve
 * this we may have to adjust the Y scrollbar position too.
 */
int set_displayPos(edview *xx, int pos) {
    char buf[100];
    int i, ret = 0;
    int delta = pos - xx->displayPos;
    edview *xx2[2];

    if (xx->link && xx->link->locked)
	xx = xx->link->xx[0];

    for (i = 0; i < 2; i++) {
	int new_y = -1, vis_pos;
	tg_rec vis_rec1, vis_rec2;
	int sheight;
	tg_rec vis_cur;

	xx2[i] = xx;

	if (!xx)
	    break;

	sheight = xx->displayHeight - xx->y_seq_end - xx->y_seq_start;

	edview_visible_items(xx, xx->displayPos,
			     xx->displayPos + xx->displayWidth);
	/*
	 * Pick an appropriate sequence to try and keep track of.
	 * This is one we'll attempt to keep visible before and after
	 * scrolling.
	 */
	vis_cur = edview_seq_visible(xx, xx->cursor_rec, NULL);

	edview_item_at_pos(xx, xx->y_seq_start,
			   0, 0, 0, 1, &vis_rec1, &vis_pos);
	edview_item_at_pos(xx, xx->displayHeight - xx->y_seq_end - 1,
			   0, 0, 0, 1, &vis_rec2, &vis_pos);
	    
	xx->displayPos += delta;

	sprintf(buf, "%d", pos);
	Tcl_SetVar2(xx->interp, xx->seq_win, "displayPos", buf,
		    TCL_GLOBAL_ONLY);

	xx->refresh_flags = ED_DISP_XSCROLL;
	if (i == 1)
	    xx->refresh_flags |= ED_DISP_NO_DIFFS;

	/* Try and ensure 'vis_rec' is still visible */
	if (vis_rec1 == -1 || !edview_seq_visible(xx, vis_rec1, &new_y)) {
	    if (new_y == -1 && vis_rec2 != -1) {
		if (edview_seq_visible(xx, vis_rec2, &new_y)) {
		    /* Already visible, but new_y is bottom loc */
		    new_y -= sheight-1;
		}
	    }

	    if (new_y != -1) {
		xx->displayYPos = new_y;
		xx->refresh_flags |= ED_DISP_YSCROLL;
	    }
	} else {
	    /* Still visible, but potentially changed Y */
	    if (new_y != -1 && new_y != xx->displayYPos) {
		xx->displayYPos = new_y;
		xx->refresh_flags |= ED_DISP_YSCROLL;
	    }
	}

	/* If editing cursor was visible, ensure it still is too */
	if (vis_cur) {
	    if (!edview_seq_visible(xx, xx->cursor_rec, &new_y)) {
		xx->displayYPos = new_y;
		xx->refresh_flags |= ED_DISP_YSCROLL;
	    }
	}

	if (xx->displayYPos + sheight > xx->nr) {
	    xx->displayYPos = xx->nr - sheight;
	    xx->refresh_flags |= ED_DISP_YSCROLL;
	}

	if (xx->displayYPos < 0) {
	    xx->displayYPos = 0;
	    xx->refresh_flags |= ED_DISP_YSCROLL;
	}

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
    switch (xx->cursor_type) {
    case GT_Contig:
	xx->cursor_apos = xx->cursor_pos;
	break;

    case GT_Seq: {
	tg_rec cnum;
	int cpos;
	sequence_get_position(xx->io, xx->cursor_rec, &cnum, &cpos, NULL,
			      NULL);
	xx->cursor_apos = cpos + xx->cursor_pos;
	break;
    }

    case GT_AnnoEle: {
	tg_rec cnum;
	range_t *r = anno_get_range(xx->io, xx->cursor_rec, &cnum, 0);
	xx->cursor_apos = r->start + xx->cursor_pos;
	break;
    }

    default:
	fprintf(stderr, "Unknown item type in edSetApos(): %d\n",
		xx->cursor_type);
    }

    /* Send a notification of cursor movement */
    cursor_notify(xx);
}

int edSetCursorPos(edview *xx, int type, tg_rec rec, int pos, int visible) {
    if (!xx)
	return 0;

    if (type == GT_Seq) {
	seq_t *s = get_seq(xx->io, rec);
	int left = s->left-1;
	int right = s->right;

	if (xx->ed->display_cutoffs) {
	    left = 0;
	    right = ABS(s->len);
	} else {
	    if (sequence_get_orient(xx->io, rec)) {
		s = get_seq(xx->io, rec);
		left  = ABS(s->len) - (s->right-1) -1;
		right = ABS(s->len) - (s->left-1);
	    }
	}

	/* If out of bounds, punt it to the consensus */
	if (pos < left || pos > right) {
	    if (visible) {
		/*
		int cpos;
		sequence_get_position(xx->io, rec, NULL, &cpos, NULL, NULL);
		type = GT_Contig;
		pos += cpos;
		rec = xx->cnum;
		*/
		if (pos < 0 || pos > ABS(s->len))
		    return 0;
		xx->ed->display_cutoffs = 1;
		Tcl_SetVar2(xx->interp, xx->edname, "Cutoffs", "1",
			    TCL_GLOBAL_ONLY);
	    } else {
		return 0;
	    }
	}
    }

    if (type != GT_Seq) {
	int ustart, uend;

	if (xx->ed->display_cutoffs) {
	    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
	    ustart = c->start;
	    uend   = c->end;
	} else {
	    char con;
	    calculate_consensus_simple(xx->io, xx->cnum, pos, pos, &con, NULL);
	    if (con != 'N') {
		/* Must be valid, so fake boundaries */
		ustart = pos;
		uend = pos;
	    } else {
		consensus_valid_range(xx->io, xx->cnum, &ustart, &uend);
	    }
	}

	uend++;
	if (pos < ustart)
	    pos = ustart;
	if (pos > uend)
	    pos = uend;
    }

    xx->cursor_type = type;
    xx->cursor_rec  = rec;
    xx->cursor_pos  = pos;

    edSetApos(xx);

    if (visible && !showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    } else {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

/*
 * Convert from a record and position to a window X,Y coordinate in sheet
 * units.
 *
 * Returns 0 on success and stores via x and y pointers.
 *        -1 on failure (rec/pos not visible).
 */
int edGetXY(edview *xx, int rec_type, tg_rec rec, int pos, int *x, int *y) {
    int i;

    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    if (xx->nr == 0)
	return -1;

    if (rec == xx->cnum) {
	int col = pos - xx->displayPos;

	if (col < 0 || col > xx->displayWidth)
	    return -1;
	
	*x = col;
	*y = 0;
	return 0;
    }

    for (i = 0; i < xx->nr; i++) {
	if (xx->r[i].rec == rec) {
	    int row, col;

	    row = xx->r[i].y + xx->y_seq_start - xx->displayYPos;
	    col = xx->r[i].start - xx->displayPos + pos;

	    if (col < 0 || col >= xx->displayWidth)
		return -1;

	    if (row < xx->y_seq_start ||
		row >= xx->displayHeight - xx->y_seq_end)
		return -1;

	    *x = col;
	    *y = row;
	    return 0;
	}
    }

    return -1;
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
	HacheItem *hi;
	tg_rec key = xx->cursor_rec;

	if (!xx->rec_hash)
	    return 0;
	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	if (!hi)
	    return 0;
	j = hi->data.i;
    }

    /* Step up until we find something overlapping */
    for (j--; j >= 0; j--) {
	if (xx->r[j].start <= cpos && xx->r[j].end+1 >= cpos
#ifndef CACHED_CONS_VISIBLE
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISCONS)
#endif
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS)
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)) {
	    if (!xx->ed->display_cutoffs) {
		seq_t *s = get_seq(xx->io, xx->r[j].rec);
		int left = s->left;
		int right = s->right;
		if (sequence_get_orient(xx->io, xx->r[j].rec)) {
		    s = get_seq(xx->io, xx->r[j].rec);
		    left  = ABS(s->len) - (s->right-1);
		    right = ABS(s->len) - (s->left-1);
		}
		if (cpos - xx->r[j].start < left-1 ||
		    cpos - xx->r[j].start > right)
		    continue; /* Sequence present, but hidden */
	    }
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

    cursor_notify(xx);
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
	HacheItem *hi;
	tg_rec key = xx->cursor_rec;

	if (!xx->rec_hash)
	    return 0;
	hi = HacheTableSearch(xx->rec_hash, (char *)&key, sizeof(key));
	if (!hi)
	    return 0;
	j = hi->data.i;
	cpos = xx->r[j].start + xx->cursor_pos;
    }

    /* Step up until we find something overlapping */
    for (j++; j < xx->nr; j++) {
	if (xx->r[j].start <= cpos && xx->r[j].end+1 >= cpos
#ifndef CACHED_CONS_VISIBLE
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISCONS)
#endif
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS)
	    && ((xx->r[j].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)) {
	    if (!xx->ed->display_cutoffs) {
		seq_t *s = get_seq(xx->io, xx->r[j].rec);
		int left = s->left;
		int right = s->right;
		if (sequence_get_orient(xx->io, xx->r[j].rec)) {
		    s = get_seq(xx->io, xx->r[j].rec);
		    left  = ABS(s->len) - (s->right-1);
		    right = ABS(s->len) - (s->left-1);
		}
		if (cpos - xx->r[j].start < left-1 ||
		    cpos - xx->r[j].start > right)
		    continue; /* Sequence present, but hidden */
	    }
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

    cursor_notify(xx);
    if (!showCursor(xx, 1, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorLeft(edview *xx) {
    if (xx->cursor_type == GT_Seq) {
	if (xx->ed->display_cutoffs) {
	    if (xx->cursor_pos > 0) {
		xx->cursor_pos--;
		xx->cursor_apos--;
	    }
	} else {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);
	    int left = s->left;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		left = ABS(s->len) - (s->right-1);
	    }

	    if (xx->cursor_pos >= left) {
		xx->cursor_pos--;
		xx->cursor_apos--;
	    }

	}
    } else {
	xx->cursor_pos--;
	xx->cursor_apos--;
    }

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edCursorRight(edview *xx) {
    if (xx->cursor_type == GT_Seq) {
	seq_t *s = get_seq(xx->io, xx->cursor_rec);

	if (xx->ed->display_cutoffs) {
	    if (xx->cursor_pos < ABS(s->len)) {
		xx->cursor_pos++;
		xx->cursor_apos++;
	    }
	} else {
	    int right = s->right;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		right = ABS(s->len) - (s->left-1);
	    }

	    if (xx->cursor_pos < right) {
		xx->cursor_pos++;
		xx->cursor_apos++;
	    }
	}
    } else {
	xx->cursor_pos++;
	xx->cursor_apos++;
    }

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadStart(edview *xx) {
    if (xx->ed->display_cutoffs) {
	if (xx->cursor_type == GT_Seq) {
	    xx->cursor_pos = 0;
	} else {
	    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
	    xx->cursor_pos = c->start;
	}
    } else {
	if (xx->cursor_type == GT_Seq) {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);

	    xx->cursor_pos = s->left-1;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		xx->cursor_pos = ABS(s->len) - (s->right-1) - 1;
	    }
	} else {
	    int ustart, uend;
	    consensus_valid_range(xx->io, xx->cursor_rec, &ustart, &uend);

	    xx->cursor_pos = ustart;
	}
    }

    edSetApos(xx);

    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadStart2(edview *xx) {
    return edReadStart(xx);
}

int edReadEnd(edview *xx) {
    if (xx->ed->display_cutoffs) {
	if (xx->cursor_type == GT_Seq) {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);
	    xx->cursor_pos = ABS(s->len);
	} else {
	    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);
	    xx->cursor_pos = c->end+1;
	}
    } else {
	if (xx->cursor_type == GT_Seq) {
	    seq_t *s = get_seq(xx->io, xx->cursor_rec);

	    xx->cursor_pos = s->right;

	    if (sequence_get_orient(xx->io, xx->cursor_rec)) {
		s = get_seq(xx->io, xx->cursor_rec);
		xx->cursor_pos = ABS(s->len) - (s->left-1);
	    }
	} else {
	    int ustart, uend;
	    consensus_valid_range(xx->io, xx->cursor_rec, &ustart, &uend);

	    xx->cursor_pos = uend+1;
	}
    } 

    edSetApos(xx);

    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edReadEnd2(edview *xx) {
    return edReadEnd(xx);
}

int edContigStart(edview *xx) {
    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

    xx->cursor_pos = c->start;
    xx->cursor_type = GT_Contig;
    xx->cursor_rec = xx->cnum;
    xx->cursor_apos = xx->cursor_pos;

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
	xx->refresh_flags = ED_DISP_CURSOR;
	edview_redraw(xx);
    }

    return 0;
}

int edContigEnd(edview *xx) {
    contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

    xx->cursor_pos = c->end;
    xx->cursor_type = GT_Contig;
    xx->cursor_rec = xx->cnum;
    xx->cursor_apos = xx->cursor_pos;

    cursor_notify(xx);
    if (!showCursor(xx, 0, 0)) {
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
 * 'seq_only' forces the item to be a sequence or consensus, and not
 * an object on them (eg annotation).
 *
 * Returns the item type GT_* on success and the record/pos in *rec, *pos
 *         -1 on failure (eg numbers, off screen, etc)
 */
int edview_item_at_pos(edview *xx, int row, int col, int name, int exact,
		       int seq_only, tg_rec *rec, int *pos) {
    int i;
    int type = -1;
    int best_delta = INT_MAX;
    char nline[MAX_NAME_WIDTH];
    //    int exact = (name && xx->ed->stack_mode) || !name;

    assert(rec);
    assert(pos);

    *rec = -1;
    *pos =  0;

    if (!xx->r)
	return -1;

    /* Special case - the reserve row numbers */
    if (row == xx->y_cons) {
	*rec = xx->cnum;
	*pos = col + xx->displayPos;
	type = GT_Contig;

	if (xx->ed->hide_annos || seq_only)
	    return type;

	/* Look for consensus tags */
	for (i = 0; i < xx->nr; i++) {
	    if (xx->r[i].y != -1)
		break;

	    if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)
		continue;

	    if (col + xx->displayPos >= xx->r[i].start &&
		col + xx->displayPos <= xx->r[i].end) {
		*rec = xx->r[i].rec;
		*pos = col + xx->displayPos - xx->r[i].start;
		type = GT_AnnoEle;
	    }
	}

	return type;
    }

    if (row < xx->y_seq_start)
	return -1;
    
    /* A sequence, so find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* Inefficient, but just a copy from tk_redisplaySeqSequences() */
    i = edview_binary_search_y(xx->r, xx->nr, xx->displayYPos);
    memset(nline, ' ', MAX_NAME_WIDTH);
    for (; i < xx->nr; i++) {
	if ((xx->ed->hide_annos || seq_only || name) &&
	    ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO))
	    continue;

#ifndef CACHED_CONS_VISIBLE
	if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISCONS)
	    continue;
#endif

	if ((xx->r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
	    continue;

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

		while (p1 < nc && nline[p1] != ' ')
		    p1++;

		if (col >= p1 && (col < p2 || col == p1))
		    delta = 0;
		else
		    delta = INT_MAX;

		if (p2 > nc)
		    p2 = nc;

		do {
		    nline[p1++] = '.';
		} while (p1 < p2);

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
	    if (best_delta >= delta) {
		best_delta =  delta;
		*rec = xx->r[i].rec;
		*pos = col + xx->displayPos - xx->r[i].start;
                type = (xx->r[i].flags & GRANGE_FLAG_ISMASK)
                       == GRANGE_FLAG_ISANNO
		    ? GT_AnnoEle
		    : GT_Seq;
	    }
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
int edview_row_for_item(edview *xx, tg_rec rec, int *xmin, int *xmax) {
    int i, r = -1;
    HacheItem *hi;

    if (rec == 0)
	return -1;

    if (rec == xx->cnum) {
	if (xmin) *xmin = -xx->displayPos;
	if (xmax) *xmax = -xx->displayPos;
	return 0;
    }

    if (xx->nr <= 0 || !xx->r)
	return -1;

    /* A sequence, so find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* And search for rec in this list */
    if (!xx->rec_hash)
	return -1;
    hi = HacheTableSearch(xx->rec_hash, (char *)&rec, sizeof(rec));
    if (!hi)
	return -1;

    i = hi->data.i;
    if (xmin) *xmin = xx->r[i].start - xx->displayPos;
    if (xmax) *xmax = xx->r[i].end   - xx->displayPos;
    r = xx->r[i].y + xx->y_seq_start - xx->displayYPos;
    
    return r >= xx->y_seq_start ? r : -1;
}


/*
 * As above, but an absolute Y location ignoring the scroll bar. This means
 * we can query the location of items off screen, as needed by the
 * edview_items_between function.
 *
 * Returns Y coordinate (and optionally min/max X coordinates) if found
 *        -1 if not (with xmin/xmax unset).
 */
int edview_abs_row_for_item(edview *xx, tg_rec rec, int *xmin, int *xmax) {
    int i, r = -1;
    HacheItem *hi;

    if (rec == 0)
	return -1;

    if (rec == xx->cnum) {
	if (xmin) *xmin = -xx->displayPos;
	if (xmax) *xmax = -xx->displayPos;
	return 0;
    }

    if (xx->nr <= 0 || !xx->r)
	return -1;

    /* A sequence, so find out what's visible */
    edview_visible_items(xx, xx->displayPos,
			 xx->displayPos + xx->displayWidth);

    /* And search for rec in this list */
    if (!xx->rec_hash)
	return -1;
    hi = HacheTableSearch(xx->rec_hash, (char *)&rec, sizeof(rec));
    if (!hi)
	return -1;

    i = hi->data.i;
    if (xmin) *xmin = xx->r[i].start - xx->displayPos;
    if (xmax) *xmax = xx->r[i].end   - xx->displayPos;
    return xx->r[i].y;
}


/*
 * Returns an array of record numbers after from_rec and before to_rec.
 * The "between" function uses the Y coordinates. It is designed for use
 * with multi-selection methods; eg click followed by shift+click.
 *
 * When in stack_mode where multiple sequences reside on the same line, the
 * between function will take heed of the X location of sequences in the
 * names panel too.
 *
 * Returns an array on success (caller to free)
 *         NULL on failure
 */
Array edview_items_between(edview *xx, tg_rec from_rec, tg_rec to_rec) {
    int from_x, from_y;
    int to_x, to_y;
    int i;
    Array a = ArrayCreate(sizeof(tg_rec), 0);
    tg_rec *r;

    if (xx->nr <= 0 || !xx->r)
	return NULL;

    from_y = edview_abs_row_for_item(xx, from_rec, &from_x, NULL);
    to_y   = edview_abs_row_for_item(xx, to_rec,   &to_x, NULL);

    if (xx->ed->stack_mode) {
	if (from_x < 0) from_x = 0;
	from_x = from_x * ((double)xx->names->sw.columns / xx->displayWidth);
	if (from_x > xx->names->sw.columns)
	    from_x = xx->names->sw.columns;
	if (to_x < 0) to_x = 0;
	to_x = to_x * ((double)xx->names->sw.columns / xx->displayWidth);
	if (to_x > xx->names->sw.columns)
	    to_x = xx->names->sw.columns;
    }

    if (from_y > to_y) {
	int tmp;
	tmp = from_y; from_y = to_y; to_y = tmp;
	tmp = from_x; from_x = to_x; to_x = tmp;
    }


    for (i = 0; i < xx->nr; i++) {
	int y = xx->r[i].y;
	
	if (y < from_y || y > to_y)
	    continue;

	if ((y == from_y || y == to_y) && xx->ed->stack_mode) {
	    int p1 = xx->r[i].start - xx->displayPos;
	    if (p1 < 0) p1 = 0;
	    p1 = p1 * ((double)xx->names->sw.columns / xx->displayWidth);

	    if (y == from_y && p1 <= from_x)
		continue;
	    if (y == to_y   && p1 >= to_x)
		continue;
	}

	r = (tg_rec *)ArrayRef(a, ArrayMax(a));
	if (!r)
	    return NULL;
	*r = xx->r[i].rec;
    }
    
    return a;
}


int inJoinMode(edview *xx) {
    return xx->link ? 1 : 0;
}

void edDisplayTrace(edview *xx) {
    seq_t *s;

    if (xx->cursor_type == GT_Seq) {
	/* Single sequence */
	s = get_seq(xx->io, xx->cursor_rec);
	if (NULL == s) return;
	tman_manage_trace("ANY", sequence_get_name(&s), xx->cursor_pos,
			  0, 0, /* left/right clips */
			  sequence_get_orient(xx->io, xx->cursor_rec),
			  1, /* base spacing */
			  sequence_get_name(&s),
			  xx, xx->cursor_rec, 0, 0);
    } else if (xx->cursor_type == GT_Contig) {
	/* Consensus click */
	rangec_t *r;
	int nr, i;
	contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

	if (NULL == c) return;

	/* Shut down existing traces */
	tman_shutdown_traces(xx, 2);

	/* And display the new ones */
	puts("FIXME: reuse existing cache of items");
	r = contig_seqs_in_range(xx->io, &c,
				 xx->cursor_apos, xx->cursor_apos,
				 CSIR_SORT_BY_X, &nr);
	if (NULL == r) return;

	for (i = 0; i < nr; i++) {
	    s = get_seq(xx->io, r[i].rec);
	    /* For now don't try to bring up mass-sequencing data from cons */
	    if (NULL == s ||
		s->seq_tech == STECH_SOLEXA ||
		s->seq_tech == STECH_SOLID)
		continue;

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
tg_rec *edGetTemplateReads(edview *xx, tg_rec seqrec, int *nrec) {
    seq_t *s = get_seq(xx->io, seqrec);
    tg_rec *r = NULL, p;

    if (!s)
	return NULL;

    /* FIXME: support s->template_rec */

    /* Solexa data is easy: we have just one other end */
    p = sequence_get_pair(xx->io, s);
    if (p > 0) {
	*nrec = 1;
	r = malloc(sizeof(*r));
	*r = p;
    } else {
	*nrec = 0;
    }

    return r;
}


/* ---------------------------------------------------------------------- */
/* Selection aka cut'n'paste handling code */

/*
 * (Un)draws the selection - toggles it on or off
 */
static void toggle_select(edview *xx, tg_rec seq, int from_pos, int to_pos) {
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

void redisplaySelection(edview *xx) {
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
	tg_rec cnum;
	int cpos, left, right, orient;
	seq_t *s = get_seq(xx->io, xx->select_seq);

	cache_incr(xx->io, s);
	sequence_get_position(xx->io, xx->select_seq,
			      &cnum, &cpos, NULL, &orient);
	pos -= cpos;

	if (xx->ed->display_cutoffs) {
	    left  = 0;
	    right = ABS(s->len);
	} else {
	    if ((s->len < 0) ^ orient) {
		left  = ABS(s->len) - (s->right-1) - 1;
		right = ABS(s->len) - (s->left-1);
	    } else {
		left  = s->left - 1;
		right = s->right;
	    }
	}

	if (pos < left)
	    pos = left;
	if (pos > right-1)
	    pos = right-1;

	cache_decr(xx->io, s);
    } else {
	contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

	/* Clip to contig extents */
	if (pos < c->start)
	    pos = c->start;
	if (pos > c->end)
	    pos = c->end;
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
 	tg_rec cnum;
	int cpos, left, right, orient;
	seq_t *s = get_seq(xx->io, xx->select_seq);

	cache_incr(xx->io, s);
	sequence_get_position(xx->io, xx->select_seq,
			      &cnum, &cpos, NULL, &orient);
	pos -= cpos;

	if (xx->ed->display_cutoffs) {
	    left  = 0;
	    right = ABS(s->len);
	} else {
	    if ((s->len < 0) ^ orient) {
		left  = ABS(s->len) - (s->right-1) - 1;
		right = ABS(s->len) - (s->left-1);
	    } else {
		left  = s->left - 1;
		right = s->right;
	    }
	}

	if (pos < left)
	    pos = left;
	if (pos > right-1)
	    pos = right-1;

	cache_decr(xx->io, s);
    } else {
	contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

	/* Clip to contig extents */
	if (pos < c->start)
	    pos = c->start;
	if (pos > c->end)
	    pos = c->end;
    }
    xx->select_end = pos;

    /* Display new selection */
    redisplaySelection(xx);
}

void edSelectSet(edview *xx, tg_rec rec, int start, int end) {
    int do_x = 0;

    /* Undisplay an old selection */
    if (xx->select_made)
	redisplaySelection(xx);

    xx->select_made = 0;


    xx->select_seq   = rec;
    xx->select_start = start;
    xx->select_end   = end;
    xx->select_made  = 1;

    /* Scroll and redraw if appropriate */
    if (xx->select_end+2 >= xx->displayPos + xx->displayWidth) {
	set_displayPos(xx, xx->select_end+2 - xx->displayWidth);
	do_x = 1;
    }
    if (xx->select_start-1 <= xx->displayPos) {
	set_displayPos(xx, xx->select_start-1);
	do_x = 1;
    }

    if (do_x) {
	xx->refresh_flags = ED_DISP_ALL;
	ed_set_xslider_pos(xx, xx->displayPos);
    }

    xx->refresh_flags |= ED_DISP_SELECTION;
    edview_redraw(xx);
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

    start = xx->select_start;
    end   = xx->select_end;

    if (start > end) {
	len = start;
	start = end;
	end = len;
    }

    start += offset;

    len = end - start+1 > bufsize ? bufsize : end - start + 1;

    if (len && xx->select_seq) {
	if (xx->select_seq != xx->cnum) {
	    seq_t *s, *sorig;

	    sorig = s = get_seq(xx->io, xx->select_seq);
	    if (sequence_get_orient(xx->io, xx->select_seq)) {
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

/*
 * Finds the next/previous difference between a pair of join editors.
 *
 * Returns 0 on sucess
 *        -1 on failure
 */
int edNextDifference(edview *xx) {
    int pos0, pos1;
    contig_t *c0, *c1;

    if (!xx->link)
	return -1;

    c0 = cache_search(xx->link->xx[0]->io, GT_Contig, xx->link->xx[0]->cnum);
    cache_incr(xx->link->xx[0]->io, c0);

    c1 = cache_search(xx->link->xx[1]->io, GT_Contig, xx->link->xx[1]->cnum);
    cache_incr(xx->link->xx[1]->io, c1);

    /* Find the positions, anchored from top contig incase they differ */
    pos1 = xx->link->xx[1]->cursor_apos + 1;
    pos0 = pos1 - xx->link->lockOffset;

    while (pos0 <= c0->end && 
	   pos1 <= c1->end) {
	int len = 1023, i;
	char cons0[1024], cons1[1024];

	if (pos0 + len > c0->end)
	    len = c0->end - pos0 + 1;
	if (pos1 + len > c1->end)
	    len = c1->end - pos1 + 1;

	/* Compute consensus fragments and compare */
	calculate_consensus_simple(xx->link->xx[0]->io, c0->rec,
				   pos0, pos0+len-1, cons0, NULL);
	calculate_consensus_simple(xx->link->xx[1]->io, c1->rec,
				   pos1, pos1+len-1, cons1, NULL);

	for (i = 0; i < len; i++) {
	    if (cons0[i] != cons1[i])
		break;
	}

	pos0 += i;
	pos1 += i;

	if (i != len)
	    break;
    }

    /* Found a difference, or at the end of the contig */
    edSetCursorPos(xx->link->xx[0], GT_Contig, c0->rec, pos0, 1);
    edSetCursorPos(xx->link->xx[1], GT_Contig, c1->rec, pos1, 1);

    cache_decr(xx->link->xx[0]->io, c0);
    cache_decr(xx->link->xx[1]->io, c1);

    return 0;
}

int edPrevDifference(edview *xx) {
    int pos0, pos1;
    contig_t *c0, *c1;

    if (!xx->link)
	return -1;

    /* Find the positions, anchored from top contig incase they differ */
    pos1 = xx->link->xx[1]->cursor_apos - 1;
    pos0 = pos1 - xx->link->lockOffset;

    c0 = cache_search(xx->link->xx[0]->io, GT_Contig, xx->link->xx[0]->cnum);
    cache_incr(xx->link->xx[0]->io, c0);

    c1 = cache_search(xx->link->xx[1]->io, GT_Contig, xx->link->xx[1]->cnum);
    cache_incr(xx->link->xx[1]->io, c1);

    while (pos0 >= c0->start && 
	   pos1 >= c1->start) {
	int len = 1023, i;
	char cons0[1024], cons1[1024];

	if (pos0 - len < c0->start)
	    len = pos0 - c0->start + 1;
	if (pos1 - len < c1->start)
	    len = pos1 - c1->start + 1;

	/* Compute consensus fragments and compare */
	calculate_consensus_simple(xx->link->xx[0]->io, c0->rec,
				   pos0-(len-1), pos0, cons0, NULL);
	calculate_consensus_simple(xx->link->xx[1]->io, c1->rec,
				   pos1-(len-1), pos1, cons1, NULL);

	for (i = len-1; i >= 0; i--) {
	    if (cons0[i] != cons1[i])
		break;
	}

	pos0 -= (len-i-1);
	pos1 -= (len-i-1);

	if (i != -1)
	    break;
    }

    /* Found a difference, or at the end of the contig */
    edSetCursorPos(xx->link->xx[0], GT_Contig, c0->rec, pos0, 1);
    edSetCursorPos(xx->link->xx[1], GT_Contig, c1->rec, pos1, 1);

    cache_decr(xx->link->xx[0]->io, c0);
    cache_decr(xx->link->xx[1]->io, c1);

    return 0;
}

void edview_set_sort_order(edview *xx) {
    contig_set_default_sort(xx->ed->group_primary, xx->ed->group_secondary);
    
    if (xx->r) xx->r_start = xx->r_end; // force re-calc in edview_visible_items
}

int depad_and_opos(char *str, int len, char *depad, int *opos) {
    int i, j;
    for (i = j = 0; i < len; i++) {
	opos[i] = j;
	if (str[i] != '*') {
	    depad[j++] = str[i];
	}
    }
    
    return j;
}

/* Compute original positions array via alignments */
int origpos(edview *xx, tg_rec srec, int pos) {
    seq_t *s = cache_search(xx->io, GT_Seq, srec);
    int tpos;
    OVERLAP *overlap = NULL;
    ALIGN_PARAMS *params = NULL;
    int *opos = NULL, op;
    Read *r;
    char *fn;
    align_int *S;
    char *depadded;
    int ulen, dlen, i, j, k, mis = 0;
    char *seq_out = NULL, *trace_out = NULL;
    int seq_out_len, trace_out_len, new_i = 1;
    HacheData hd;
    HacheItem *hi;
    int seq_comp;

    if (!s)
	return pos+1;

    seq_comp = sequence_get_orient(xx->io, srec);

    /* Fetch trace */
    hi = HacheTableSearch(xx->trace_hash, (char *)&srec, sizeof(srec));
    if (hi) {
	if (hi->data.i == -1) {
	    /* Previously tried and failed */
	    return pos+1;
	} else {
	    r = hi->data.p;
	}
    } else {
	fn = s->trace_name && *s->trace_name
	    ? s->trace_name
	    : s->name;
	
	if (!(r = read_reading(fn, TT_ANYTR))) {
	    hd.i = -1;
	    HacheTableAdd(xx->trace_hash, (char *)&srec, sizeof(srec), hd,
			  NULL);

	    return pos+1;
	}

	hd.p = r;
	HacheTableAdd(xx->trace_hash, (char *)&srec, sizeof(srec), hd, NULL);
    }

    /* Depad seq and do simple match to trace */
    ulen = ABS(s->len);
    depadded = malloc(ulen);
    opos = malloc(ulen*sizeof(*opos));

    dlen = depad_and_opos(s->seq, ulen, depadded, opos);

    if (dlen == r->NBases) {
	mis = 0;
	for (i = 0; i < dlen; i++) {
	    if (r->base[i] == '-')
		r->base[i] = 'N';
	    if (r->base[i] != depadded[i])
		mis = 1;
	}
    } else {
	mis = 1;
	for (i = 0; i < r->NBases; i++) {
	    if (r->base[i] == '-')
		r->base[i] = 'N';
	}
    }

    if (!mis) {
	op = seq_comp
	    ? opos[ABS(s->len) - pos]+1
	    : opos[pos]+1;
	free(depadded);
	free(opos);
	return op;
    }

    /* Not identical, so align it fully */
    //printf("Mismatch %d so align\n", mis);
    //printf("%.*s\n", dlen, depadded);
    //printf("%.*s\n", r->NBases, r->base);

    /* Align it within a tight band */
    overlap = create_overlap();
    init_overlap(overlap, depadded, r->base, dlen, r->NBases);

    params = create_align_params();
    set_align_params(params,
		     10,               // band
		     0,                // gap_open
		     0,                // gap_extend
		     EDGE_GAPS_COUNT,  // edge_mode
		     0,                // job
		     0,                // seq1_start
		     0,                // seq2_start
		     '.',              // new_pad_sym
		     '*',              // old_pad_sym
		     0);               // set_job

    affine_align(overlap, params);
    destroy_alignment_params (params);

    //print_overlap(overlap,stdout);
	
    seq_out_len = overlap->seq_out_len;
    seq_out = overlap->seq1_out;
    trace_out = overlap->seq2_out;


    /*
     * Rewrite the opos[] array to map to trace coords now.
     *
     * i = offset in orig padded seq (unaligned, but with pads in editor)
     * j = offset in aligned seq (orig vs trace, aligned to match)
     * k = corresponding trace base (unaligned).
     */
    i = j = k = 0;
    while (i < ulen && j < seq_out_len && k < r->NBases) {
	opos[i] = k;

	//printf("%d:%c %d:%c %d:%c %d:%c\n",
	//       i, s->seq[i],
	//       j, seq_out[j],
	//       j, trace_out[j],
	//       k, r->base[k]);

	if (s->seq[i] == '*') {
	    if (seq_out[j] == '.') {
		if (trace_out[j] != '.')
		    k++;
		i++; j++;
	    } else {
		i++;
	    }
	} else {
	    if (seq_out[j] == '.') {
		if (trace_out[j] != '.')
		    k++;
		j++;
	    } else {
		if (trace_out[j] != '.')
		    k++;
		i++; j++;
	    }
	}
    }
    while (i < ulen)
	opos[i++] = k;

    if (depadded)
	free(depadded);

    if (overlap)
	destroy_overlap(overlap);

    //printf("2 Origpos %d=%d\n", pos, s->opos[pos]+1);
    op = seq_comp
	? opos[ABS(s->len) - pos]+1
	: opos[pos]+1;
    free(opos);

    return op;
}

void ed_set_base_sort_point(edview *xx) {
    contig_set_base_sort_point(xx->cursor_apos);
}

void ed_set_sequence_sort(edview *xx) {
    contig_set_sequence_sort(xx->select_seq == xx->cnum ? GT_Contig : GT_Seq,
    	    	    	     xx->select_seq, xx->select_start, xx->select_end);
}
    
