#include "misc.h"
#include "gap-dbstruct.h"
#include "IO.h"
#include "template.h"
#include "tagUtils.h"

/*
 * Lookup tables to calculate the primer and strand values. We have two
 * primer tables, one where we always make a best guess and never return
 * UNKNOWN.
 * In case of conflicting strand and primer information, the primer takes
 * precedence (as the strand info will vanish in time).
 * The exception to this rule is the case where primer is CUSTFOR, but
 * strand is REVERSE, as this corresponds to old databases where we didn't
 * distinguish between forward and reverse custom primers.
 */
int primer_type_arr[5][2] = {
    {GAP_PRIMER_UNKNOWN, GAP_PRIMER_UNKNOWN},
    {GAP_PRIMER_FORWARD, GAP_PRIMER_FORWARD},
    {GAP_PRIMER_REVERSE, GAP_PRIMER_REVERSE},
    {GAP_PRIMER_CUSTFOR, GAP_PRIMER_CUSTREV},
    {GAP_PRIMER_CUSTREV, GAP_PRIMER_CUSTREV}
};

int primer_type_guess_arr[5][2] = {
    {GAP_PRIMER_CUSTFOR, GAP_PRIMER_CUSTREV},
    {GAP_PRIMER_FORWARD, GAP_PRIMER_FORWARD},
    {GAP_PRIMER_REVERSE, GAP_PRIMER_REVERSE},
    {GAP_PRIMER_CUSTFOR, GAP_PRIMER_CUSTREV},
    {GAP_PRIMER_CUSTREV, GAP_PRIMER_CUSTREV}
};

int strand_arr[5][2] = {
    {GAP_STRAND_FORWARD, GAP_STRAND_REVERSE},
    {GAP_STRAND_FORWARD, GAP_STRAND_FORWARD},
    {GAP_STRAND_REVERSE, GAP_STRAND_REVERSE},
    {GAP_STRAND_FORWARD, GAP_STRAND_FORWARD},
    {GAP_STRAND_REVERSE, GAP_STRAND_REVERSE}
};

static int check_universal_primers = 1;
static int standard_primer_tol = 100;

/*
 * Allocate and initialise a new gel_cont_t structure.
 */
gel_cont_t *new_gel_cont() {
    gel_cont_t *gc;

    if (NULL == (gc = (gel_cont_t *)xcalloc(sizeof(gel_cont_t), 1)))
	return NULL;

    return gc;
}

/*
 * Deallocate a template_c structure.
 */
void free_gel_cont(gel_cont_t *gc) {
    xfree(gc);
}

/*
 * Deallocate a template_c structure.
 */
static void free_template_c(template_c *t) {
    if (t->gel_cont)
	free_list(t->gel_cont, (void (*)(void *))free_gel_cont);
    xfree(t);
}

/*
 * Allocate and initialise a new template_c structure.
 */
static template_c *new_template_c(void) {
    template_c *t;

    if (NULL == (t = (template_c *)xcalloc(1, sizeof(template_c))))
	return NULL;

    t->consistency = TEMP_CONSIST_OK;
    t->score = 1; /* perfect score */
    if (NULL == (t->gel_cont = new_list())) {
	free_template_c(t);
	return NULL;
    }

    return t;
}

static int init_template_contig(GapIO *io, int contig,
				template_c **template_check,
				int connected) {
    int i, old_contig;
    GReadings r;
    gel_cont_t *gc;
    item_t *item;

    for (i = io_clnbr(io, contig); i; i = io_rnbr(io, i)) {
	gel_read(io, i, r);
	
	/*
	 * Skip over this if we're only looking to add readings to templates
	 * that we already know about.
	 */
	if (connected && !template_check[r.template])
	    continue;

	if (!template_check[r.template]) {
	    if (!(template_check[r.template] = new_template_c()))
		return -1;
	    template_check[r.template]->num = r.template;
	}
	if (!(gc = new_gel_cont()))
	    return -1;
	gc->contig = contig;
	gc->read = i;
	if ((item = tail(template_check[r.template]->gel_cont))) {
	    old_contig = ((gel_cont_t *)item->data)->contig;
	    if (gc->contig != old_contig)
		template_check[r.template]->flags |=
		    TEMP_FLAG_SPANNING;
	}
	add_item(template_check[r.template]->gel_cont, (void *)gc);
    }

    return 0;
}

/*
 * Initialise memory used for template checks.
 * num_contigs and contig_array specifies which contigs to include templates
 * from. If num_contigs==0 then use all.
 *
 * If 'connected' is set then we then scan all contigs again adding reading
 * information to only those templates that we already have information on.
 *
 * Returns "template_c **" for success, NULL for error.
 */
template_c **init_template_checks(GapIO *io, int num_contigs,
				  int *contig_array, int connected) {
    template_c **template_check = NULL;
    int c;

    /* NB: We'll start the array from 1 for simplicity */
    template_check = (template_c **)xcalloc(Ntemplates(io) + 1,
					    sizeof(template_c *));

    if (NULL == template_check)
	return NULL;

    /*
     * foreach contig
     *     foreach gel
     *         Add contig/gel to template_c structure
     */
    if (num_contigs) {
	for (c = 0; c < num_contigs; c++) {
	    if (-1 == init_template_contig(io, contig_array[c],
					   template_check, 0))
		return NULL;
	}

	/*
	 * Note that by initialising the requested contigs first, and 
	 * any connected contigs afterwards, we make sure that the
	 * first contig in the t->gel_cont structure is the first one in the
	 * specified contig_array. This has (correct) implications for the 
	 * t->direction computed. See the check_template_strand() function.
	 */
	if (connected) {
	    int i, j;

	    for (i = 1; i <= NumContigs(io); i++) {
		/* Ignore those we've already done */
		for (j = 0; j < num_contigs; j++) {
		    if (contig_array[j] == i) {
			j = -1;
			break;
		    }
		}
		if (j == -1)
		    continue;

		if (-1 == init_template_contig(io, i, template_check, 1))
		    return NULL;
	    }
	}
    } else {
	for (c = 1; c <= NumContigs(io); c++) {
	    if (-1 == init_template_contig(io, c, template_check, 0))
		return NULL;
	}
    }

    return template_check;
}

/*
 * Free memory used for template checks.
 */
void uninit_template_checks(GapIO *io, template_c **template_check) {
    int i;

    for (i = 1; i <= Ntemplates(io); i++) {
	if (template_check[i])
	    free_template_c(template_check[i]);
    }

    xfree(template_check);
}


/*
 * Check the primer and sense (whether complemented) information
 * for each reading is consistent with other readings.
 *
 * This spots problems such as:
 *
 *  <----     ---->  (forward & reverse reads on same strand)
 *  ---------------
 *
 *  ---->	     (two forward reads, but on different strands)
 *  ---------------
 *            ---->
 */
static void check_template_strand(GapIO *io, template_c *t) {
    item_t *item, *at;
    int contig;
    gel_cont_t *gc;
    GReadings r;
    int first_contig;
    
    at = head(t->gel_cont);

    /*
     * Loop around readings in list checking template consistency for
     * readings within this template with others within this template that
     * are also on this template.
     *
     * When finished, if we've spotted some on another contig then start
     * again.
     *
     * With readings spanning multiple contigs this will iterate around
     * multiple times. The template direction returned is from the first
     * contig observed. Note that when init_template_checks() is called with
     * a single contig and including connected contigs then the single contig
     * requested is always the first in the list, hence this is the one from
     * which the direction is computed.
     */
    first_contig = 1;
    while (at) {
	int dir_plus = 0;
	int dir_minus = 0;
	item = at;
	at = NULL;
	contig = 0;

	/*
	 * Foreach reading in list
	 *     Compute overall direction of the template from this one read.
	 *     Accumulate direction counts in dir_plus and dir_minus
	 *     Negate the contig number for this read/contig pair to
	 *       signify that it's been checked.
	 * Set direction based on dir_plus/dir_minus
	 * if dir_plus and dir_minus are both non zero => strand inconsistency
	 */
	for (; item; item = item->next) {
	    int sense, strand;
	    gc = (gel_cont_t *)(item->data);
	    
	    if (contig == 0)
		contig = gc->contig;
	    else if (gc->contig != contig && gc->contig > 0) {
		at = item;
		continue;
	    }

	    gel_read(io, gc->read, r);
	    
	    sense = r.sense == GAP_SENSE_ORIGINAL ? 0 : 1;
	    strand = (PRIMER_TYPE_GUESS(r) == GAP_PRIMER_FORWARD ||
		      PRIMER_TYPE_GUESS(r) == GAP_PRIMER_CUSTFOR) ? 0 : 1;
	    
	    if (sense == strand)
		dir_plus++;
	    else
		dir_minus++;

	    gc->contig = -gc->contig;
	}

	if (dir_plus && dir_minus)
	    t->consistency |= TEMP_CONSIST_STRAND;

	if (first_contig) {
	    if (dir_plus > dir_minus)
		t->direction = 0;
	    else if (dir_plus < dir_minus)
		t->direction = 1;
	    else
		t->direction = -1;

	    first_contig = 0;
	}
    }

    /*
     * Now make all contigs in the read/contig pairs positive again.
     */
    for (item = head(t->gel_cont); item; item = item->next) {
	gc = (gel_cont_t *)(item->data);

	gc->contig = -gc->contig;
    }

    /* ... and finally, set template_c flags */
    t->flags |= TEMP_FLAG_CHECKED_STRAND;
}


/*
 * Looks for SVEC tags on a reading in an attempt to better identify the true
 * start offset of this template.
 *
 * Sets 'start' and 'end' to contain computed start and ends of templates,
 * or UNKNOWN_POS if not known.
 */
static void get_template_tag(GapIO *io, int gel, int primer,
			     int *start_p, int *end_p, GReadings *r,
			     int *tag_start, int *tag_end, int use_cvec) {
    GAnnotations *a;
    char *types[] = {"SVEC", "CVEC"};
    int start = UNKNOWN_POS;
    int end = UNKNOWN_POS;
    int ntypes;
    int cosmid_end;

    ntypes = use_cvec ? 2 : 1;

    a = vtagget(io, gel, ntypes, types);

    *tag_start = 0;
    *tag_end = 0;

    while (a && a != (GAnnotations *)-1) {
	cosmid_end = 0;
	if (ntypes == 2 && a->type == str2type("CVEC")) {
	    /* Cosmid tag, treat as bound to nearest end */
	    if (a->position < r->length - (a->position + a->length - 1))
		cosmid_end = 1;
	    else
		cosmid_end = 2;
	}

	if (a->position == 1 || cosmid_end == 1) { /* template start */

	    if (primer == GAP_PRIMER_FORWARD ||
		primer == GAP_PRIMER_CUSTFOR) {
		*tag_start = 1;

		start = r->position - r->start +
		    (r->sense
		     ? r->length - 1 - (a->length + a->position - 1)
		     : a->length + a->position - 1);

	    } else if (primer == GAP_PRIMER_REVERSE ||
		       primer == GAP_PRIMER_CUSTREV) {
		*tag_end = 1;

		end = r->position - r->start +
		    (r->sense
		     ? r->length - 1 - (a->length + a->position - 1)
		     : a->length + a->position - 1);
	    } /* else primer is unknown */
	    
	} else if (a->position + a->length == r->length + 1 ||
		   cosmid_end == 2) {
	    
	    if (primer == GAP_PRIMER_FORWARD ||
		primer == GAP_PRIMER_CUSTFOR) {
		*tag_end = 1;

		end = r->position - r->start +
		    (r->sense
		     ? r->length - 1 - a->position
		     : a->position);

	    } else if (primer == GAP_PRIMER_REVERSE ||
		       primer == GAP_PRIMER_CUSTREV) {
		*tag_start = 1;

		start = r->position - r->start +
		    (r->sense
		     ? r->length - 1 - a->position
		     : a->position);
	    } /* else primer is unknown */
	}

	a = vtagget(io, 0, ntypes, types);
    }

    /*
     * If start and end aren't known, then default to start and end of reading
     */
    if (start != UNKNOWN_POS) {
	*start_p = start;
    } else if (primer != GAP_PRIMER_UNKNOWN) {
	*start_p =
	    ((r->sense == 0 && (primer == GAP_PRIMER_FORWARD ||
				primer == GAP_PRIMER_CUSTFOR)) ||
	     (r->sense == 1 && (primer == GAP_PRIMER_REVERSE ||
				primer == GAP_PRIMER_CUSTREV)))
	    ? r->position - r->start + 1
	    : r->position - r->start + 1 + r->length - 1;
    } else {
	*start_p = UNKNOWN_POS;
    }

    if (end != UNKNOWN_POS) {
	*end_p = end;
    } else if (primer != GAP_PRIMER_UNKNOWN) {
	*end_p =
	    ((r->sense == 0 && (primer == GAP_PRIMER_FORWARD ||
				primer == GAP_PRIMER_CUSTFOR)) ||
	     (r->sense == 1 && (primer == GAP_PRIMER_REVERSE ||
				primer == GAP_PRIMER_CUSTREV)))
	    ? r->position - r->start + 1 + r->length - 1
	    : r->position - r->start + 1;
    } else {
	*end_p = UNKNOWN_POS;
    }
}


/*
 * Sets the template_c start, end, min and max fields for the template in
 * a specified contig. This function also checks the validity of the
 * primer values.
 *
 * Returns 'orient' - a guess of template orientation based upon
 * available, primer information. This is not reliable when the primer
 * information is inconsistent, but then again, what is.
 */
static int get_template_positions2(GapIO *io, template_c *t, int contig,
				   int *found_forp, int *found_revp) {
    item_t *item;
    gel_cont_t *gc;
    GReadings r;
    int st, end, found_for = 0, found_rev = 0;
    int temp_min   = UNKNOWN_POS; /* Min,   from reading positions */
    int temp_max   = UNKNOWN_POS; /* Max,   from reading positions */
    int temp_start = UNKNOWN_POS; /* Start, from PRIMER_FORWARD */
    int temp_end   = UNKNOWN_POS; /* End,   from PRIMER_REVERSE */
    int orient = 0;
    int tag_start, tag_end;

    /*
     * Scan through gel+cont list looking for gels on this contig.
     */
    for (item = head(t->gel_cont); item; item = item->next) {
	gc = (gel_cont_t *)(item->data);

	if (gc->contig != contig)
	    continue;

	gel_read(io, gc->read, r);

	/* Find minimum and maximum extents of readings on the template */
	if (r.position < temp_min || temp_min == UNKNOWN_POS)
	    temp_min = r.position;

	if (r.position + r.sequence_length - 1 > temp_max)
	    temp_max = r.position + r.sequence_length - 1;
	
	/*
	 * Get estimated start and end coordinates for this template.
	 */
	get_template_tag(io, gc->read, PRIMER_TYPE_GUESS(r), &st, &end, &r,
			 &tag_start, &tag_end,
			 t->oflags & TEMP_OFLAG_CVEC);
	
	/* Check primer type, and derive start/end accordingly */
	switch (PRIMER_TYPE_GUESS(r)) {
	case GAP_PRIMER_FORWARD:
	    orient = r.sense == GAP_SENSE_ORIGINAL ? 1 : -1;
	    found_for = 1;
	    if (check_universal_primers) {
		if (temp_start != UNKNOWN_POS && tag_start) {
		    if (st != UNKNOWN_POS) {
			/* Already exists - are they near each other? */
			if (ABS(temp_start - st) > standard_primer_tol)
			    t->consistency |= TEMP_CONSIST_PRIMER;
		    
			temp_start = MIN(temp_start, st);
		    }
		} else {
		    temp_start = st;
		}
	    } else {
		temp_start = (temp_start == UNKNOWN_POS)
		    ? st : MIN(temp_start, st);
	    }

	    if (tag_end) {
		if (temp_end != UNKNOWN_POS) {
		    if (end != UNKNOWN_POS) {
			/* redfining end - does it match? */
			if (ABS(temp_end - end) > standard_primer_tol)
			    t->consistency |= TEMP_CONSIST_PRIMER;
		    }

		} else {
		    /* FIXME: Not strictly true, but good enough? */
		    found_rev = 1;
		    temp_end = end;
		}
	    }
	    break;
	    
	case GAP_PRIMER_CUSTFOR:
	    orient = r.sense == GAP_SENSE_ORIGINAL ? 1 : -1;
	    found_for = 1;
	    if (temp_start != UNKNOWN_POS && st != UNKNOWN_POS) {
		/* Already exists - is custom leftwards of univeral? */
		if (r.sense == GAP_SENSE_ORIGINAL) {
		    if (st + standard_primer_tol < temp_start)
			t->consistency |= TEMP_CONSIST_PRIMER;
		} else {
		    if (st - standard_primer_tol > temp_start)
			t->consistency |= TEMP_CONSIST_PRIMER;
		}
	    }
	    if (temp_end == UNKNOWN_POS && tag_end) {
		found_rev = 1; /* FIXME: Not strictly true, but good enough? */
		temp_end = end;
	    }
	    break;
	    
	case GAP_PRIMER_REVERSE:
	    orient = r.sense == GAP_SENSE_REVERSE ? 1 : -1;
	    found_rev = 1;
	    if (check_universal_primers) {
		if (temp_end != UNKNOWN_POS && tag_end) {
		    if (end != UNKNOWN_POS) {
			/* Already exists - are they near each other? */
			if (ABS(temp_end - end) > standard_primer_tol)
			    t->consistency |= TEMP_CONSIST_PRIMER;
		    
			temp_end = MAX(temp_end, end);
		    }
		} else {
		    temp_end = end;
		}
	    } else {
		temp_end = (temp_end == UNKNOWN_POS)
		    ? end : MAX(temp_end, end);
	    }

	    if (tag_start) {
		if (temp_start != UNKNOWN_POS) {
		    if (st != UNKNOWN_POS) {
			/* redfining start - does it match? */
			if (ABS(temp_start - st) > standard_primer_tol)
			    t->consistency |= TEMP_CONSIST_PRIMER;
		    }
		} else {
		    /* FIXME: Not strictly true, but good enough? */
		    found_for = 1;
		    temp_start = st;
		}
	    }
	    break;
	    
	case GAP_PRIMER_CUSTREV:
	    orient = r.sense == GAP_SENSE_REVERSE ? 1 : -1;
	    found_rev = 1;
	    if (temp_end != UNKNOWN_POS && end != UNKNOWN_POS) {
		/* Already exists - is custom rightwards of universal? */
		if (r.sense == GAP_SENSE_REVERSE) {
		    if (end - standard_primer_tol > temp_end)
			t->consistency |= TEMP_CONSIST_PRIMER;
		} else {
		    if (end + standard_primer_tol < temp_end)
			t->consistency |= TEMP_CONSIST_PRIMER;
		}
	    }
	    if (temp_start == UNKNOWN_POS && tag_start) {
		found_for = 1; /* FIXME: Not strictly true, but good enough? */
		temp_start = st;
	    }
	    break;

	case GAP_PRIMER_UNKNOWN:
	    if (temp_end == UNKNOWN_POS && tag_end) {
		found_rev = 1; /* FIXME: Not strictly true, but good enough? */
		temp_end = end;
	    }
	    if (temp_start == UNKNOWN_POS && tag_start) {
		found_for = 1; /* FIXME: Not strictly true, but good enough? */
		temp_start = st;
	    }
	    t->consistency |= TEMP_CONSIST_UNKNOWN;
	}
    }

    t->start = temp_start;
    t->end   = temp_end;
    t->min   = temp_min;
    t->max   = temp_max;

    *found_forp = found_for;
    *found_revp = found_rev;

    return orient;
}


/*
 * Sets the template_c start, end, min and max fields as performed in
 * get_template_positions2(). The function then also attempts to guess
 * the start and end coordinates when they are unknown.
 */
void get_template_positions(GapIO *io, template_c *t, int contig) {
    int found_for, found_rev, orient;
    GTemplates te;

    t->consistency &= ~(TEMP_CONSIST_DIST | TEMP_CONSIST_PRIMER |
			TEMP_CONSIST_UNKNOWN);

    /*
     * Initialise template_c start, end, min and max fields
     */
    orient = get_template_positions2(io, t, contig, &found_for,
				     &found_rev);
    /*
     * Adjust template_c start and end fields when unknown
     */
    template_read(io, t->num, te);
    if (t->start == UNKNOWN_POS) {
	if (found_for)
	    t->start = t->min;
	t->flags |= TEMP_FLAG_GUESSED_START;
    } else {
	t->flags &= ~TEMP_FLAG_GUESSED_START;
    }

    if (t->end == UNKNOWN_POS) {
	if (found_rev)
	    t->end = t->max;
	t->flags |= TEMP_FLAG_GUESSED_END;
    } else {
	t->flags &= ~TEMP_FLAG_GUESSED_END;
    }

    /* Set min/max sizes to be the same to start with */
    t->start2 = t->start;
    t->end2   = t->end;

    if (!found_for && !found_rev) {
	/*
	 * No forward or reverse primer information at all!?
	 * In this case guestimate by using min and max
	 */
	t->start = t->start2 = t->min;
	t->end   = t->end2   = t->max;

	t->flags |= TEMP_FLAG_GUESSED_START | TEMP_FLAG_GUESSED_END;

    } else if (t->consistency & ~TEMP_CONSIST_UNKNOWN) {
	/*
	 * If it's inconsistent then just pick any old template to
	 * draw, so that the template line drawn (eg. in the
	 * template display) is at least clearly visible (and hence
	 * highlights the problem).
	 */
	if (!found_for) {
	    t->start = t->start2 = t->min;
	    t->flags |= TEMP_FLAG_GUESSED_START;
	} else {
	    t->end = t->end2 = t->max;
	    t->flags |= TEMP_FLAG_GUESSED_END;
	}

    } else {
	/*
	 * Here we need to guess the other end of the template by looking
	 * at the expected template insert size.
	 * If TEMP_OFLAG_MINMAX_SIZE is set then we compute both the min and
	 * max possible insert sizes, otherwise we compute the average
	 * (expected) size.
	 */
	int tsize1, tsize2;

	tsize1 = (t->oflags & TEMP_OFLAG_MINMAX_SIZE)
	    ? te.insert_length_min
	    : (te.insert_length_min + te.insert_length_max) / 2;

	tsize2 = (t->oflags & TEMP_OFLAG_MINMAX_SIZE)
	    ? te.insert_length_max
	    : (te.insert_length_min + te.insert_length_max) / 2;

	if (!found_for) {
	    t->start  = t->end  - (orient ? orient : 1) * tsize1;
	    t->start2 = t->end2 - (orient ? orient : 1) * tsize2;
	    t->flags |= TEMP_FLAG_GUESSED_START | TEMP_FLAG_EXPECTED;
	}
	
	if (!found_rev) {
	    t->end  = t->start  + (orient ? orient : 1) * tsize1;
	    t->end2 = t->start2 + (orient ? orient : 1) * tsize2;
	    t->flags |= TEMP_FLAG_GUESSED_END | TEMP_FLAG_EXPECTED;
	}
    }

    t->flags |= TEMP_FLAG_DONE_POSITIONS;
}


/*
 * Check the relative positions of the start and end primers on a template.
 * This spots problems like:
 *
 *  ----> .. big gap ..         (too large, or small, template length)
 *  --------------------------
 *                       <----
 *
 *  <----			(readings pointing in opposite directions)
 *  --------------------------
 *                       ---->
 */
static void check_template_position(GapIO *io, template_c *t) {
    item_t *item, *new_contig_l;
    gel_cont_t *gc;
    GTemplates te;
    int contig = 0, tmp;

    new_contig_l = head(t->gel_cont);

    /*
     * Loop around the list starting at new_contig_l looking for
     * all readings in this contig (which are implicitly on this
     * template).
     */
    while (new_contig_l) {
	gc = (gel_cont_t *)(new_contig_l->data);

	/*
	 * Call get_template_pos to check this template/contig pairing.
	 */
	tmp = t->consistency;
	get_template_positions(io, t, gc->contig);
	t->consistency |= tmp;

	/*
	 * FIXME: fill out an array indexed by contig containing start, end
	 * direction, etc, so we can work out the NxN analysis of the maximum
	 * expected contig overlap for spanning templates.
	 */

	/*
	 * Check the primer size is within range.
	 */
	template_read(io, t->num, te);
	if ((ABS(t->end - t->start) < te.insert_length_min) ||
	    (ABS(t->end - t->start) > te.insert_length_max))
	    t->consistency |= TEMP_CONSIST_DIST;
	
	if (t->max - t->min > te.insert_length_max)
	    t->consistency |= TEMP_CONSIST_DIST;
	
	/*
	 * Now we step through the list starting from new_contig_l looking
	 * for a new contig, and one that we haven't visited yet. We mark the
	 * ones that we've visited by temporarily negating the contig number.
	 */
	contig = gc->contig;
	item = new_contig_l;
	new_contig_l = NULL;
	for (; item; item = item->next) {
	    gc = (gel_cont_t *)(item->data);
	    if (gc->contig == contig)
		gc->contig = -gc->contig;
	    else if (gc->contig > 0 && !new_contig_l)
		new_contig_l = item;
	}
    }

    /*
     * Now make all contigs in the read/contig pairs positive again.
     */
    for (item = head(t->gel_cont); item; item = item->next) {
	gc = (gel_cont_t *)(item->data);
	if (gc->contig < 0)
	    gc->contig = -gc->contig;
    }

    /* ... and finally, set some more template_c flags */
    t->flags |= TEMP_FLAG_CHECKED_PRIMER;
    t->flags |= TEMP_FLAG_CHECKED_DIST;
}

/*
 * Given a template this checks for a valid length based on the start/end,
 * min/max and for spanning templates distance from contig ends.
 * The overlap parameter is positive when the contig sequences overlap
 * and negative when they do not. If you do not know any overlap information
 * specify this as zero.
 *
 * The start/end values filled up will be wrong for inconsistent templates.
 * Eg        fwd              rev
 *         ------>          ------->
 *         ^                ^      ^
 *       start/            end    max
 *       min
 *
 * The min/max values here make more sense though. So the computed size will
 * be the max of the two.
 *
 * Also, if we have multiple contigs, compute the distance from the ends.
 *
 * eg:                       overlap size (-ve in this eg)
 *                              |<-->|
 *            --fwd-->          |    |           <--rev--
 *     -------------------------      ---------------------------
 *            |                                         |
 *            |<----------- template size ------------->|
 *
 * This allows us to mark a template as inconsistent, due to probable
 * misassemble, based on the assumption that the contigs should be joined if
 * there is a real overlap.
 */
void check_template_length(GapIO *io, template_c *t, int overlap) {
    item_t *item;
    int distf, distr, c1;
    GReadings r;
    GTemplates te;
    
    template_read(io, t->num, te);

    if (t->start < t->end) {
	if (t->start > t->min)
	    t->start = t->min;
	if (t->start2 > t->min)
	    t->start2 = t->min;
	if (t->end < t->max)
	    t->end = t->max;
	if (t->end2 < t->max)
	    t->end2 = t->max;
    } else {
	if (t->start < t->max)
	    t->start = t->max;
	if (t->start2 < t->max)
	    t->start2 = t->max;
	if (t->end > t->min)
	    t->end = t->min;
	if (t->end2 > t->min)
	    t->end2 = t->min;
    }
    t->computed_length = ABS(MAX(t->end, t->end2) - MIN(t->start, t->start2));

    if (t->computed_length > te.insert_length_max)
	t->consistency |= TEMP_CONSIST_DIST;

    if (!(t->flags & TEMP_FLAG_SPANNING))
	return;

    c1 = 0;
    distf = distr = 0;

    /*
     * Find first contig.
     * The compare this against each successive contig.
     */
    for (item = head(t->gel_cont); item; item = item->next) {
	gel_cont_t *gc = (gel_cont_t *)item->data;
	int tmp_distf, tmp_distr;

	if (!c1)
	    c1 = gc->contig;
	else if (gc->contig == c1)
	    continue;

	gel_read(io, gc->read, r);

	switch(PRIMER_TYPE(r)) {
	case GAP_PRIMER_UNKNOWN:
	    continue;

	case GAP_PRIMER_FORWARD:
	case GAP_PRIMER_CUSTFOR:
	    /* Distance forward to 'end' of contig */
	    tmp_distf = r.sense
		? r.position + r.sequence_length - 1
		: io_clength(io, gc->contig) - r.position + 1;
	    distf = MAX(distf, tmp_distf);
	    break;

	case GAP_PRIMER_REVERSE:
	case GAP_PRIMER_CUSTREV:
	    /* Distance backwards to 'start' of contig */
	    tmp_distr = r.sense
		? r.position + r.sequence_length - 1
		: io_clength(io, gc->contig) - r.position + 1;
	    distr = MAX(distr, tmp_distr);
	}
    }

    if (distf && distr) {
	int length;
	length = distf + distr - overlap;

	t->computed_length = length;

	if (length > te.insert_length_max) {
	    t->consistency |= TEMP_CONSIST_DIST;
	}
	if (overlap > 0 && length < te.insert_length_max) {
	    t->consistency |= TEMP_CONSIST_DIST;
	}
    }
}

/*
 * As check_template_length, but taking into account a specific ordering
 * of two contigs (c1 left, c2 right) and an expected overlap size.
 *
 * Computed length is set to the template length, but if the template does
 * not span c1 and c2 then it is set to zero.
 */
void check_template_length_overlap(GapIO *io, template_c *t,
				   int c1, int c2, int offset) {
    int c1_start, c1_end;
    int c2_start, c2_end;
    int f1, f2;
    int both_guessed = TEMP_FLAG_GUESSED_START | TEMP_FLAG_GUESSED_END;
    
    if (!(t->flags & TEMP_FLAG_SPANNING)) {
	t->computed_length = 0;
	return;
    }


    /* Compute positions in each contig, adjusting for offset in c2 */
    get_template_positions(io, t, c1);
    c1_start = t->start;
    c1_end = t->end;
    f1 = t->flags;
    int length;
	
    get_template_positions(io, t, c2);
    c2_start = t->start + offset;
    c2_end = t->end + offset;
    f2 = t->flags;

    t->consistency = 0; /* Mark as valid, until we determine otherwise */

    /* Check if it's spanning at all... */
    if ((f1 & both_guessed) == both_guessed ||
	(f2 & both_guessed) == both_guessed) {
	t->computed_length = 0;
	return;
    }
	
    /* Check length */
    length = -1;
    if (!(f1 & TEMP_FLAG_GUESSED_START) &&
	!(f2 & TEMP_FLAG_GUESSED_END)) {
	length = ABS(c2_end - c1_start);
    } else if (!(f1 & TEMP_FLAG_GUESSED_END) &&
	       !(f2 & TEMP_FLAG_GUESSED_START)) {
	length = ABS(c1_end - c2_start);
    }

    /* Check validity of length */
    if (length != -1) {
	GTemplates te;

	template_read(io, t->num, te);
	t->computed_length = length;
	if (length > te.insert_length_max ||
	    length < te.insert_length_min) {
	    t->consistency |= TEMP_CONSIST_DIST;
	}
    }

    /*
     * If we have fwd reads in both contigs, or rev reads in both contigs,
     * then the 'primer positions' should align.
     */
    if (!(f1 & TEMP_FLAG_GUESSED_START) &&
	!(f2 & TEMP_FLAG_GUESSED_START)) {
	if (ABS(c1_start - c2_start) > standard_primer_tol)
	    t->consistency |= TEMP_CONSIST_PRIMER;
    }
    if (!(f1 & TEMP_FLAG_GUESSED_END) &&
	!(f2 & TEMP_FLAG_GUESSED_END)) {
	if (ABS(c1_end - c2_end) > standard_primer_tol)
	    t->consistency |= TEMP_CONSIST_PRIMER;
    }
}

/*
 * Computes a template score based on consistency, flags and size.
 */
void score_template(GapIO *io, template_c *t) {
    if (t->consistency & TEMP_CONSIST_STRAND) {
	/* Most probable cause is a naming error */
	t->score *= 0.5;
    }

    if (t->consistency & TEMP_CONSIST_PRIMER) {
	t->score *= 0.7;
    }

    if (t->consistency & TEMP_CONSIST_DIST) {
	GTemplates te;
	int size = ABS(t->end - t->start);
	double scale;

	template_read(io, t->num, te);
	if (size < te.insert_length_min)
	    scale = (double)size/te.insert_length_min;
	else if (size > te.insert_length_max)
	    scale = (double)te.insert_length_max/size;
	else
	    scale = 1;
	if (t->computed_length > te.insert_length_max)
	    scale = MIN(scale,
			(double)te.insert_length_max/t->computed_length);

	/* _too_ far apart, probably naming error instead */
	if (scale < 0.5)
	    scale = 0.5;

	t->score *= scale;
    }

    if (t->flags & TEMP_FLAG_GUESSED_START) {
	t->score *= 0.9;
    }

    if (t->flags & TEMP_FLAG_GUESSED_END) {
	t->score *= 0.9;
    }
}

/*
 * Checks a single template_c structure for consistency.
 *
 * Returns the new consistency value (-1 for error).
 */
int check_template_c(GapIO *io, template_c *t) {
    t->consistency = 0;
    check_template_strand(io, t);
    check_template_position(io, t);
    check_template_length(io, t, 0);
    score_template(io, t);

    return t->consistency;
}


void check_all_templates(GapIO *io, template_c **template_check) {
    int i;

    for (i = 1; i <= Ntemplates(io); i++) {
	if (template_check[i])
	    (void)check_template_c(io, template_check[i]);
    }
}

/*
 * Removes any templates from the template_check[] that don't have more
 * than one reading on them.
 */
void remove_single_templates(GapIO *io, template_c **template_check) {
    int i;
    item_t *item;

    for (i = 1; i <= Ntemplates(io); i++) {
	if (template_check[i]) {
	    item = head(template_check[i]->gel_cont);
	    if (!item->next) {
		free_template_c(template_check[i]);
		template_check[i] = NULL;
	    }
	}
    }
}


/*
 * Remove any templates from template_check[] that don't span multiple
 * contigs.
 */
void contig_spanning_templates(GapIO *io, template_c **template_check) {
    int i;

    for (i = 1; i <= Ntemplates(io); i++) {
	if (template_check[i] &&
	    !(template_check[i]->flags & TEMP_FLAG_SPANNING)) {
	    free_template_c(template_check[i]);
	    template_check[i] = NULL;
	}
    }
}


/*
 * Sorts templates based on consistency. Best templates are listed first.
 * The template_check[] array is left intact, but the template_order[]
 * array is initialised and sortted. It contains indexes into the
 * template_check[] array.
 *
 * Returns order for success, NULL for error.
 *
 * The order of precedence is:
 *   No problems
 *   Missing information
 *   Distance problem
 *   Primer problem
 *   Strand problem
 */
static template_c **tarr_tmp;

int *sort_templates(GapIO *io, template_c **template_check) {
    int *template_order;
    int i, j;
    int sort_template_func(const void *, const void *);

    if (NULL == (template_order = (int *)xmalloc((Ntemplates(io)+1) *
						 sizeof(int))))
	return NULL;

    /* Initialise order to currently available templates */
    for (i=1, j=0; i <= Ntemplates(io); i++)
	if (template_check[i])
	    template_order[j++] = i;
    template_order[j] = 0;

    /* Sort */
    tarr_tmp = template_check;
    qsort(template_order, j, sizeof(int), sort_template_func);

    return template_order;
}

int sort_template_func(const void *vi, const void *vj) {
    int c1, c2;
    double d;
    
    c1 = tarr_tmp[*((int *)vi)]->consistency;
    c2 = tarr_tmp[*((int *)vj)]->consistency;

    if (c1 == c2) {
	d = tarr_tmp[*((int *)vj)]->score - tarr_tmp[*((int *)vi)]->score;
	return d > 0 ? 1 : (d == 0 ? 0 : -1);
    }

    if ((c1 & TEMP_CONSIST_STRAND) && !(c2 & TEMP_CONSIST_STRAND))
	return 1;
    else if ((c2 & TEMP_CONSIST_STRAND) && !(c1 & TEMP_CONSIST_STRAND))
	return -1;

    if ((c1 & TEMP_CONSIST_PRIMER) && !(c2 & TEMP_CONSIST_PRIMER))
	return 1;
    else if ((c2 & TEMP_CONSIST_PRIMER) && !(c1 & TEMP_CONSIST_PRIMER))
	return -1;

    if ((c1 & TEMP_CONSIST_DIST) && !(c2 & TEMP_CONSIST_DIST))
	return 1;
    else if ((c2 & TEMP_CONSIST_DIST) && !(c1 & TEMP_CONSIST_DIST))
	return -1;

    if ((c1 & TEMP_CONSIST_UNKNOWN) && !(c2 & TEMP_CONSIST_UNKNOWN))
	return 1;
    else if ((c2 & TEMP_CONSIST_UNKNOWN) && !(c1 & TEMP_CONSIST_UNKNOWN))
	return -1;

    return 0;
}


void dump_template(template_c *t) {
    item_t *it;

    printf("%3d: %04d-%04d, %04d-%04d, 0x%02x, 0x%x:", 
	   t->num, t->start, t->end, t->min, t->max,
	   t->consistency, t->flags);
    
    for (it = head(t->gel_cont); it; it=it->next) {
	gel_cont_t *gc = (gel_cont_t *)(it->data);
	
	printf(" %02d.%03d", gc->contig, gc->read);
    }
    
    putchar('\n');
}

void dump_templates(GapIO *io, template_c **template_check, int *order) {
    int i;

    if (order) {
	for (i = 0; order[i]; i++)
	    dump_template(template_check[order[i]]);
    } else {
	for (i = 1; i <= Ntemplates(io); i++) {
	    if (template_check[i])
		dump_template(template_check[i]);
	}
    }
}


/*
 * Find the last used base at the 3' end of the specified reading for this
 * template.
 * When used in conjunction with the template length this gives an estimate
 * for how much undetermined sequence we have left on this template.
 *
 * Returns the position (relative to the start of the template, whichever
 * orientation it may be).
 */
int last_template_base(GapIO *io, template_c *t, int g) {
    if (t->start <= t->end) {
	if (io_length(io, g) > 0)
	    return io_length(io, g);
	else
	    return t->end - io_relpos(io, g);
    } else {
	if (io_length(io, g) > 0)
	    return t->start - (io_relpos(io, g) + io_length(io, g));
	else
	    return -io_length(io, g);
    }
}

/*
 * Checks whether any of the sequences on template 't' for contig 'c' are 
 * in the range of range_start to range_end inclusive.
 *
 * Returns number of bases covered.
 *	   ie. "0" when not covered and
 *	   "range_end-range_start+1" when totally covered.
 */
int template_covered(GapIO *io, template_c *t, int c,
		     int range_start, int range_end) {
    int covered = 0;
    item_t *item;
    gel_cont_t *gc;
    GReadings r;

    for (item = head(t->gel_cont); item; item = item->next) {
	gc = (gel_cont_t *)(item->data);

	if (gc->contig != c)
	    continue;

	gel_read(io, gc->read, r);

	if (r.position <= range_start &&
	    r.position + r.sequence_length >= range_end) {
	    covered = range_end - range_start + 1; /* all */
	    break;
	}

	if (r.position <= range_start &&
	    r.position + r.sequence_length >= range_start)
	    covered += r.position + r.sequence_length - range_start + 1;

	if (r.position <= range_end &&
	    r.position + r.sequence_length >= range_end)
	    covered += range_end - r.position + 1;
    }

    return covered;
}
