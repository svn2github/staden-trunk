/*
 * Code to suggest various reading reactions to perform based on problems
 * found.
 *
 * NASTY ASSUMPTIONS: That sizeof(int) <= sizeof(void *). See add_item() calls
 * for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "io_utils.h"
#include "gap-dbstruct.h"
#include "IO.h"
#include "fort.h"
#include "qual.h"
#include "qualP.h"
#include "tagUtils.h"
#include "list.h"
#include "edUtils.h"
#include "xalloc.h"
#include "template.h"
#include "misc.h"
#include "gap_globals.h"

#define abs(a) ((a) > 0 ? (a) : -(a))

#define AVG_READ_LEN 300

struct taq_solution {
    int gel;
    char name[DB_NAMELEN+1];
    int position;
    int distance;
};

/*
 * ---------------------------------------------------------------------------
 * LONG GEL CODE
 * ---------------------------------------------------------------------------
 */

/*
 * Outputs a long gel selection result.
 */
void report_long(GapIO *io, int gel, int position, int coverage, int size,
		 int end, template_c **tarr) {
    if (gel) {
	GTemplates t;
	GReadings r;
	char *name;
	template_c *tc;
	
	gel_read(io, gel, r);
	name = io_rname(io, gel);
	
	GT_Read(io, arr(GCardinal, io->templates, r.template-1),
		&t, sizeof(t), GT_Templates);

	tc = tarr[r.template];
	
	vmessage("%c%c  Long %*s %5d. T_pos=%3d, T_size=%d..%d (%d), cov %d%s"
		 ,
		 (tc->consistency & TEMP_CONSIST_UNKNOWN) ? '?' : ' ',
		 (tc->consistency & TEMP_CONSIST_DIST)    ? 'D' : ' ',
		 DB_NAMELEN, name, position,
		 last_template_base(io, tc, gel),
		 t.insert_length_min, t.insert_length_max,
		 TEMP_LENGTH(tc),
		 coverage, (!end && coverage >= size) ? "*\n" : "\n");
    } else {
	vmessage("    No solution.\n");
    }
}

/*
 * Picks readings for rerunning as a long gel.
 * 'read' is a gel reading number guaranteed to be leftwards (or the
 * read itself) of 'pos'.
 *
 * Selection criteria is based on primarily read position in the contig and
 * its template. We output all the possible solutions.
 */
void pick_long(GapIO *io, int read, int pos, int size, int strand,
	       int end, int avg_long_read_len, template_c **tarr) {
    GReadings r;
    int e, fail = 1;

    for(; read; read = r.right) {
	gel_read(io, read, r);

	/* Skip 'majorly' inconsistent templates */
	if (tarr[r.template] &&
	    (tarr[r.template]->consistency & ~(TEMP_CONSIST_DIST |
					       TEMP_CONSIST_UNKNOWN)))
	    continue;

	if (strand) {
	    /*
	     * assume left end is clipped for a reason - eg vector,
	     * was: e = r.position + r.length - r.start;
	     */

	    e = r.position + r.end - r.start - 2;

	    /*
	    if (r.position - r.start + r.length - io->db.max_gel_len > pos)
		break;
	    */
	    
	    if (r.sense == 1 && e > pos && e - avg_long_read_len < pos) {
		int coverage;

		coverage = pos -
		    (e - (avg_long_read_len > TEMP_LENGTH(tarr[r.template])
			  ? TEMP_LENGTH(tarr[r.template])
			  : avg_long_read_len));

		if (coverage > 0) {
		    report_long(io, read, e, coverage, size, end, tarr);
		    fail = 0;
		}
	    }

	} else {
	    /*
	     * assume left end is clipped for a reason - eg vector,
	     * so we ignore the r.start value
	     */
	    e = r.position /* - r.start */;
	    if (e >= pos)
		break;
	    else if (r.sense == 0 && e + avg_long_read_len > pos) {
		int coverage;

		coverage = e +
		    (avg_long_read_len > TEMP_LENGTH(tarr[r.template])
		     ? TEMP_LENGTH(tarr[r.template])
		     : avg_long_read_len) - pos;

		if (coverage > 0) {
		    report_long(io, read, e, coverage, size, end, tarr);
		    fail = 0;
		}
	    }
	}
    }

    if (fail)
	report_long(io, 0, 0, 0, 0, 0, NULL);

    return;
}

/*
 * Finds locations suitable for longs gels.
 *
 * These are:
 *  i) All holes requiring double stranding
 * ii) The ends of contigs.
 *
 * Note that in case ii) we may have already dealt with this problem in the
 * guise of double stranding. Only when both strands end at exactly the same
 * base will this not be the case.
 * We can detect these cases for each end. For the left end, our first problem
 * should be at position 1.
 * For the right end, if pos+len != contig length, we've also missed the case.
 *
 * In either case, the region is defined by lreg and rreg. Setting these to
 * zero will use the entire contig instead of a specified region.
 *
 * FIXME: what about long walks? Is that in this option or the walk option?
 *        Also we always scan from the leftmost read and not the first read
 *        nearer than max_gel_len leftwards.
 */
int find_long_gels_single(GapIO *io, int contig, int lreg, int rreg,
			  int avg_len) {
    int pos, len, last = -1, first = -1, left_read, i;
    char *reason;
    GContigs c;
    template_c **tarr;

    (void)GT_Read(io, arr(GCardinal, io->contigs, contig-1),
		  &c, sizeof(GContigs), GT_Contigs);

    left_read = c.left;

    if (!lreg)
	lreg = 1;
    if (!rreg)
	rreg = c.length;

    /* Initialise templates to the single contig */
    if (NULL == (tarr = init_template_checks(io, 1, &contig, 1)))
	return -1;
    check_all_templates(io, tarr);
    for (i = 0; i <= Ntemplates(io); i++) {
	if (tarr[i]) {
	    if (tarr[i]->flags & TEMP_FLAG_SPANNING) {
		get_template_positions(io, tarr[i], contig);
	    }
	}
    }
    
    /* init hole detection */
    if (-1 == next_hole(contig, lreg, rreg,
			consensus_cutoff, quality_cutoff, NULL, NULL,
			database_info, (void *)io))
	return -1;

    pos = lreg - 1;
    while (pos <= rreg &&
	   (pos = next_hole(0, pos+1, 0, 0.0, quality_cutoff, &reason, &len,
			    NULL, NULL)) > 0) {
	    /*
	     * Our first problem is a +ve strand one, we'll do our extend
	     * left end case here then. Or possibly it is a -ve one, but
	     * not at the start (ie both strands start together on the
	     * contig).
	     */
	if (first && lreg == 1 &&
	    (*reason == R_NONE_GOOD ||
	     *reason == R_NONE_BAD ||
	     pos != -1)) {
	    vmessage("Prob 1..1:\tExtend contig start for joining.\n");
	    pick_long(io, left_read, 1, INT_MAX, 1, 1, avg_len, tarr);
	    vmessage("\n");
	    
	    first = 0;
	}

	vmessage("Prob %d..%d:\t", pos, pos + len - 1);

	if (*reason == R_GOOD_NONE || *reason == R_BAD_NONE) {
	    if (first) {
		vmessage("Extend contig start for joining.\n");
		pick_long(io, left_read, pos + len, len, 1, 1, avg_len, tarr);
		first = 0;
	    } else {
		vmessage("No -ve strand data.\n");
		pick_long(io, left_read, pos + len, len, 1, 0, avg_len, tarr);
	    }

	} else if (*reason == R_NONE_GOOD || *reason == R_NONE_BAD) {
	    if (pos == c.length) {
		vmessage("Extend contig end for joining.\n");
		pick_long(io, left_read, pos, len, 0, 1, avg_len, tarr);
	    } else {
		vmessage("No +ve strand data.\n");
		pick_long(io, left_read, pos, len, 0, 0, avg_len, tarr);
	    }
	    
	    last = pos;
	    
	} else {
	    vmessage("No data available!\n");
	}

	vmessage("\n");

	pos += len-1;
    }
	
    /* special case - no problems? */
    if (first && lreg == 1) {
	vmessage("Prob 1..1:\tExtend contig start for joining.\n");
	pick_long(io, left_read, 1, INT_MAX, 1, 1, avg_len, tarr);
    }

    if (rreg == c.length && last != c.length) {
	vmessage("Prob %d..%d:\tExtend contig end for joining.\n",
	       c.length, c.length);
	pick_long(io, left_read, c.length, INT_MAX, 0, 1, avg_len, tarr);
    }

    uninit_template_checks(io, tarr);
    return 0;
}

int find_long_gels(GapIO *io, int num_contigs, contig_list_t *contigs,
		   int avg_len) {
    int i, e = 0;

    for (i = 0; i < num_contigs; i++) {
	vmessage("\n-- Searching contig %s --\n\n",
		 get_contig_name(io, contigs[i].contig));
	e |= find_long_gels_single(io, contigs[i].contig, contigs[i].start,
				   contigs[i].end, avg_len);
    }

    return e;
}



/*
 * ---------------------------------------------------------------------------
 * TAQ TERMINATOR CODE
 * ---------------------------------------------------------------------------
 */

/*
 * Reports taq solutions. We pass a 'list' of taq_solution structs.
 */
void report_taq(list_t *lp) {
    item_t *i;
    struct taq_solution *ts;

    if (!head(lp)) {
	vmessage("   No solution.\n");
	return;
    }

    for (i = head(lp); i; i = i->next) {
	ts = (struct taq_solution *)i->data;

	vmessage("   Taq for %-20s %6d %3d\n",
	       ts->name, ts->position, ts->distance);
    }
}

/*
 * Picks a template for a taq reaction. We'll only display readings less
 * than the average taq length from the problem and that are on the same
 * strand as the compression/stop.
 *
 * The 'from' parameter is used internally to decide where to start searching
 * from. It should be set to the first gel in the contig initially and then
 * left alone between subsequent calls. It's used for efficiencies sake.
 */
list_t *pick_taq(GapIO *io, int pos, int len, int strand,
	       int avg_len, int *from) {
    int gel, updated = 0;
    GReadings r;
    list_t *lp;
    struct taq_solution *ts;

    /* init */
    lp = new_list();
    gel = *from;

    /* find readings within avg_len of pos */
    while (gel) {
	int pe, ps;

	gel_read(io, gel, r);

	/*
	if (r.position + r.end - r.start - io->db.max_gel_len > pos)
	    break;
	*/

	ps = r.sense       ? r.position + r.end - r.start - 2: r.position;
	pe = ps + (r.sense ? -avg_len                        : avg_len);
	
	if ( ( (r.sense == GAP_SENSE_ORIGINAL &&
 	        pe > pos &&
	        ps < pos)
	      ||
	       (r.sense == GAP_SENSE_REVERSE &&
	        pe < pos + len - 1 &&
	        ps > pos + len - 1))
	    &&
	     r.sense == strand) {

	    ts = (struct taq_solution *)xmalloc(sizeof(struct taq_solution));
	    if (ts == NULL)
		return NULL;

	    ts->gel = gel;
	    ts->position = ps;
	    strcpy(ts->name, io_rname(io, gel));
	    ts->distance = abs(ps - pos);

	    add_item(lp, (void *)ts);

	    if (!updated) {
		updated = 1;
		*from = gel;
	    }
	}

	gel = r.right;
    }

    return lp;
}

/*
 * Creates an array of compression and stop tags within a given region.
 * Internally we use the list mechanism to build up our list of tags, and
 * then sort this data into an array.
 */
int sort_tags(const void *v1, const void *v2) {
    return (*(GAnnotations **)v1)->position - (*(GAnnotations **)v2)->position;
}

GAnnotations **list_comps(GapIO *io, int contig, int lreg, int rreg,
			  int *nitems) {
    GContigs c;
    GReadings r;
    GAnnotations *ap, *ap2, **tags = NULL;
    int gel, j;
    list_t *l;
    item_t *i;
    char *types[] = {"COMP", "STOP"};
    
    /* init */
    GT_Read(io, arr(GCardinal, io->contigs, contig-1),
	    &c, sizeof(c), GT_Contigs);
    gel = c.left;
    l = new_list();
    *nitems = 0;

    /* Loop around finding tags */
    while (gel) {
	gel_read(io, gel, r);

	if (r.position > rreg)
	    break;

	if (r.position >= lreg) {
	    ap = vtagget(io, gel, sizeof(types)/sizeof(*types), types);
	    
	    while (ap && ap != (GAnnotations *)-1) {
		int t_pos = r.sense
		    ? (r.position + (r.length - r.start) -
		       (ap->position + ap->length - 1))
			: (r.position - r.start + ap->position - 1);

		if (t_pos > rreg || t_pos + ap->length < lreg)
		    break;
		
		/* found tag - so add it */
		ap2 = (GAnnotations *)xmalloc(sizeof(GAnnotations));
		memcpy(ap2, ap, sizeof(GAnnotations));
		/* 'fix' a few fields to be more useful */
		ap2->position = t_pos;
		ap2->strand = r.sense;
		if (-1 == add_item(l, ap2)) {
		    verror(ERR_FATAL, "list_comps",
			   "Failed to add item to tag list");
		    return NULL;
		}
		(*nitems)++;
		
		ap = vtagget(io, 0, sizeof(types)/sizeof(*types), types);
	    }
	}

	gel = r.right;
    }

    if (*nitems) {
	/* Sort tags. Need to create an array of item->data pointers first. */
	tags = (GAnnotations **)xmalloc(*nitems * sizeof(GAnnotations *));
	for (i=head(l), j = 0; i; i = i->next, j++)
	    tags[j] = (GAnnotations *)i->data;

	qsort ((void *)tags, *nitems, sizeof(GAnnotations *), sort_tags);
    }

    /* remove list (only interested in array now) */
    free_list(l, NULL);

    return tags;
}

/*
 * Finds suitable locations for running taq terminators.
 *
 * These are:
 * i) All locations where the compression tag is used.
 *ii) All locations where the stop tag is used.
 *
 * To avoid suggestions of running a taq terminator even when a terminator has
 * been ran and the compression has been solved, the COMP tag should be
 * changed to a resolved compression (RCMP) tag. This is currently a task
 * performed by the user.
 */
int find_taq_terms_single(GapIO *io, int contig, int lreg, int rreg,
			  int avg_len) {
    int ntags, i, pos;
    GAnnotations **tags;
    GContigs c;
    char t[5];
    list_t *sols = NULL;
#if 0
    int last_prob = -avg_len;
    char solution[100];
#endif

    /* Initialise tags */
    tags = list_comps(io, contig, lreg, rreg, &ntags);

    GT_Read(io, arr(GCardinal, io->contigs, contig-1),
	    &c, sizeof(c), GT_Contigs);

    /*
     * Loops through our tags. Once we've "solved" a problem, we've also
     * solved all other problems within 'avg_len' bases. So we output the
     * same solution for these.
     */
    for (i=0; i<ntags; i++) {
	pos = tags[i]->position;

	vmessage("\nProb %d..%d: %s tag on strand %d (%s)\n",
	       pos, pos + tags[i]->length - 1,
	       type2str(tags[i]->type, t),
	       tags[i]->strand,
	       tags[i]->strand ? "reverse" : "forward");

#if 0
	if (pos >= last_prob + avg_len) {
	    if (sols) free_list(sols, xfree);
	    sols = pick_taq(io, pos, tags[i]->length, tags[i]->strand,
			    avg_len, &c.left);
	    last_prob = pos;
	}
#else
	if (sols) free_list(sols, free);
	sols = pick_taq(io, pos, tags[i]->length, tags[i]->strand,
			avg_len, &c.left);
#endif

	report_taq(sols);
    }
	
    /* tidyup */
    for (i=0; i<ntags; i++)
	xfree(tags[i]);
    xfree(tags);

    if (sols)
	free_list(sols, free);

    return 0;
}

int find_taq_terms(GapIO *io, int num_contigs, contig_list_t *contigs,
		   int avg_len) {
    int i, e = 0;

    for (i = 0; i < num_contigs; i++) {
	vmessage("\n-- Searching contig %s --\n\n",
		 get_contig_name(io, contigs[i].contig));
	e |= find_taq_terms_single(io, contigs[i].contig, contigs[i].start,
				   contigs[i].end, avg_len);
    }

    return e;
}
