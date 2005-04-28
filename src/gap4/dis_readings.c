#include "misc.h"
#include "IO.h"
#include "IO2.h"
#include "fort.h"
#include "list_proc.h"
#include "dis_readings.h"
#include "tagUtils.h"
#include "notes.h"
#include "dbcheck.h"

/**
 * Unlinks an individual read by linking its left and right neighbours to
 * one another.
 *
 * Only updates certain in-memory caches...
 */
static void unlink_reading(GapIO *io, int rnum, int cnum, int *rnum_changed) {
    GReadings *r, *rl, *rr;

    r = arrp(GReadings, io->reading, rnum-1);
    if (r->left) {
	rl = arrp(GReadings, io->reading, r->left-1);
	rl->right = r->right;
	io_rnbr(io, r->left) = r->right;
	rnum_changed[r->left] = -1;
    } else {
	io_clnbr(io, cnum) = io_rnbr(io, rnum);
    }

    if (r->right) {
	rr = arrp(GReadings, io->reading, r->right-1);
	rr->left = r->left;
	io_lnbr(io, r->right) = r->left;
	rnum_changed[io_rnbr(io, rnum)] = -1;
    } else {
	io_crnbr(io, cnum) = io_lnbr(io, rnum);
    }

    update_rnumtocnum(io, rnum, 0);
}


/**
 * Removes and deallocates a reading by swapping it with the last reading
 * numbers (relinking as needed) and decrementing the number of readings.
 *
 * Only updates certain in-memory caches...
 */
static void remove_and_swap_reading(GapIO *io, int rnum, int *rnum2cnum,
				   int *rnum_changed) {
    GReadings *rdel, *rmov;
    int cnum;

    if (rnum != NumReadings(io)) {
	/* Swap the last used reading structures with this one. */
	swap_read(io, rnum, NumReadings(io));


	/*
	 * And now relink the neighbours too for the reading we shuffled
	 * down.
	 * rmov = the one just moved. rdel = the one to delete.
	 */
	rmov = arrp(GReadings, io->reading, rnum-1);
	rdel = arrp(GReadings, io->reading, NumReadings(io)-1);
	rnum_changed[rnum] = -1;

	cnum = rnum2cnum[NumReadings(io)];
	rnum2cnum[NumReadings(io)] = rnum2cnum[rnum];
	update_rnumtocnum(io, NumReadings(io), rnum2cnum[rnum]);
	rnum2cnum[rnum] = cnum;
	update_rnumtocnum(io, rnum, cnum);

	if (rmov->left) {
	    GReadings *rl = arrp(GReadings, io->reading, rmov->left-1);
	    rl->right = rnum;
	    rnum_changed[rmov->left] = -1;
	} else {
	    io_clnbr(io, cnum) = rnum;
	}

	if (rmov->right) {
	    GReadings *rr = arrp(GReadings, io->reading, rmov->right-1);
	    rr->left = rnum;
	    rnum_changed[rmov->right] = -1;
	} else {
	    io_crnbr(io, cnum) = rnum;
	}
    } else {
	rdel = arrp(GReadings, io->reading, rnum-1);
    }

    /*
     * Deallocate reading and decrement database count.
     */
    remove_gel_tags(io, NumReadings(io), 0, 0);
    delete_note_list(io, rdel->notes);
    io_deallocate_reading(io, NumReadings(io));
    NumReadings(io)--;
}


/**
 * Removes contigs with left and right neighbours set to zero (ie empty).
 */
static void delete_empty_contigs(GapIO *io) {
    int i;

    for (i = 1; i <= NumContigs(io); i++) {
	if (io_clnbr(io, i) != 0 || io_crnbr(io, i) != 0)
	    continue;

	io_delete_contig(io, i);
	i--;
    }
}

/**
 * Moves an unlinked reading into a specific contig.
 * This inserts the reading at the appropriate position.
 * FIXME: we could optimise this by inserting all at the start, and
 * sorting at the end.
 */
static void move_read_to_contig(GapIO *io, int rnum, int cnum) {
    int rnbr;
    int pos = io_relpos(io, rnum);
    int len = ABS(io_length(io, rnum));
    GReadings *r, *rr;

    /* Find read to right of this one, or zero if none */
    rnbr = io_clnbr(io, cnum);
    while (rnbr && io_relpos(io, rnbr) <= pos)
	rnbr = io_rnbr(io, rnbr);

    if (rnbr) {
	io_rnbr(io, rnum) = rnbr;
	if ((io_lnbr(io, rnum) = io_lnbr(io, rnbr))) {
	    rr = arrp(GReadings, io->reading, io_lnbr(io, rnum)-1);
	    rr->right = rnum;
	    io_rnbr(io, io_lnbr(io, rnum)) = rnum;
	} else {
	    io_clnbr(io, cnum) = rnum;
	}
	rr = arrp(GReadings, io->reading, rnbr-1);
	rr->left = rnum;
	io_lnbr(io, rnbr) = rnum;
    } else {
	if ((io_lnbr(io, rnum) = io_crnbr(io, cnum))) {
	    rr = arrp(GReadings, io->reading, io_lnbr(io, rnum)-1);
	    rr->right = rnum;
	    io_rnbr(io, io_lnbr(io, rnum)) = rnum;
	} else {
	    io_clnbr(io, cnum) = rnum;
	}
	io_rnbr(io, rnum) = 0;
	io_crnbr(io, cnum) = rnum;
    }

    r = arrp(GReadings, io->reading, rnum-1);
    r->right = io_rnbr(io, rnum);
    r->left  = io_lnbr(io, rnum);

    if (io_clength(io, cnum) < pos + len-1)
	io_clength(io, cnum) = pos + len-1;

    update_rnumtocnum(io, rnum, cnum);
}

static int rsort_int(const void *pi1, const void *pi2) {
    int i1 = *(int *)pi1;
    int i2 = *(int *)pi2;

    return i2 - i1;
}

/**
 * Copies all annotations from contig cfrom to contig cto where a reading
 * in cto exists. cto may have holes, so we skip copying annotations in these
 * holes. Annotations overlapping a hole are copied (without clipping).
 * The expectation is that all of these inconsistencies (holes, clipped
 * annotations, etc) will be resolved by a call to remove_contig_holes().
 */
static void copy_consensus_annotations(GapIO *io, int cfrom, int cto) {
    GContigs cf, ct;
    GAnnotations a, last;
    int anno, lastnum;
    int rnum;

    /* Find start of annotation list. Bail out early if none exist */
    contig_read(io, cfrom, cf);
    contig_read(io, cto,   ct);
    if (!(anno = cf.annotations))
	return;
    tag_read(io, anno, a);

    /*
     * Walk along reading list in cto. For each reading, scan from
     * current 'anno' to beyond the end of this reading to see whether or
     * not this annotation overlaps (and hence requires copying).
     */
    lastnum = 0;
    for (rnum = io_clnbr(io, cto); rnum; rnum = io_rnbr(io, rnum)) {
	int start = io_relpos(io, rnum);
	int end = start + ABS(io_length(io, rnum))-1;

	while (anno && a.position <= end) {
	    if (a.position + a.length >= start) {
		int newanno;
		newanno = get_free_tag(io);
		if (lastnum) {
		    last.next = newanno;
		    tag_write(io, lastnum, last);
		} else {
		    ct.annotations = newanno;
		    contig_write(io, cto, ct);
		}
		lastnum = newanno;
		tag_read(io, lastnum, last);
		last.position = a.position;
		last.length = a.length;
		last.type = a.type;
		last.strand = a.strand;
		last.next = 0;
		if (a.annotation) {
		    char *comment = get_comment(io, a.annotation);
		    last.annotation = put_comment(io, comment);
		    xfree(comment);
		}
	    }
	    anno = a.next;
	    if (anno)
		tag_read(io, anno, a);
	}
    }

    if (lastnum)
	tag_write(io, lastnum, last);
}


/**
 * Removes a set of readings from either the contig or the database.
 *
 * When removing from the database we need to delete everything related to
 * that reading (annotations, template if not used elsewhere, etc). When
 * moving to a new contig, we prefer to keep all readings clustered together.
 *
 * Ie if we remove A & B (overlapping) from one contig and C from another
 * then we create two new contigs containing A & B in one and C in the other.
 *
 * When creating new contigs, we have the option of copying over any
 * overlapping consensus tags to the new contigs. This choice only refers to
 * consensus tags; reading tags are always copied.
 *
 * move == 0   => remove
 * move == 1   => split to new single-read contigs
 * move == 2   => move to new still-joined contigs
 */
int disassemble_readings(GapIO *io, int *rnums, int nreads, int move,
			 int duplicate_tags)
{
    int i, j;
    int *rnum2cnum = NULL;
    int *rnum_changed = NULL; /* 0 for unchanged, new cnum for changed */
    int *new_cnum = NULL;
    int cn, rn;
    int last_read;

    /*
     * To implement this we firstly take the readings out of the contigs
     * (moving/removing as appropriate) and then tidy up the holes in the
     * existing contig, splitting the contig if desired.
     *
     * So the plan is:
     *
     * 0. Reverse sort all reading numbers (so as to never renumber a
     *    reading we wish to delete).
     * 1. Produce a table of which contig each reading number is in.
     * 2. Unlink readings from their neighbours
     * 3. If removing, delete the unlinked readings (and shuffle numbers
     *    around)
     * 4. If keeping, form new contigs based on the contigs the readings
     *    were in. No need to worry about positions or overlaps.
     * 5. Remove the holes in the new contigs, splitting and shifting
     *    as necessary.
     */


    /*
     * Reverse sort by number. We remove the highest one first so that
     * we never get in hot water by swapping the reading number of a reading
     * that we still need to delete.
     */
    qsort(rnums, nreads, sizeof(rnums[0]), rsort_int);

    /*
     * Filter out any duplicates.
     */
    last_read = 0;
    for (i = j = 0; i < nreads; i++) {
	if (rnums[i] == last_read) {
	    continue;
	}
	last_read = rnums[j++] = rnums[i];
    }
    nreads = j;

    /*
     * Produce a table of which contig each reading number is in.
     */
    rnum2cnum = (int *)xmalloc((NumReadings(io)+1) * sizeof(int));
    rnum_changed = (int *)xmalloc((NumReadings(io)+1) * sizeof(int));
    for (cn = 1; cn <= NumContigs(io); cn++) {
	for (rn = io_clnbr(io, cn); rn; rn = io_rnbr(io, rn)) {
	    rnum2cnum[rn] = cn;
	    rnum_changed[rn] = 0;
	}
    } 

    
    /*
     * A table of new contig numbers to be used for moving all reads from
     * contig X into the new contig Y, preserving overlaps where possible.
     */
    new_cnum = (int *)xcalloc(NumContigs(io)+1, sizeof(int));


    /*
     * Unlink readings from their neighbours.
     * Need to consider the boundary cases where the "neighbour" is a contig
     * rather than another reading.
     * Do this entirely using in-memory data structures for now...
     */
    for (i = 0; i < nreads; i++) {
	unlink_reading(io, rnums[i], rnum2cnum[rnums[i]], rnum_changed);
    }


    /*
     * If removing, delete the unlinked readings (and shuffle numbers around).
     * We know how many we are removing, so we swap 
     */
    if (move == 0) {
	for (i = 0; i < nreads; i++) {
	    remove_and_swap_reading(io, rnums[i], rnum2cnum, rnum_changed);
	}
	/* Update global GDatabase record to update num_readings */
	GT_Write(io,GR_Database,&io->db,sizeof(io->db),GT_Database);
    } else {
	/* split into single-read contigs or move as a whole to new contig */
	for (i = 0; i < nreads; i++) {
	    int rnum = rnums[i];
	    int cnum = new_cnum[rnum2cnum[rnum]];

	    /* Create a new contig if desired */
	    if (move == 1 || !cnum) {
		io_init_contig(io, NumContigs(io)+1);
		cnum = NumContigs(io);
		io_clnbr(io, cnum) = 0;
		io_crnbr(io, cnum) = 0;
		io_clength(io, cnum) = 0;
		new_cnum[rnum2cnum[rnum]] = cnum;
		vmessage("New contig created containing reading %s\n",
			 io_rname(io, rnum));
	    }

	    /* Move reading to the contig */
	    move_read_to_contig(io, rnum, cnum);

	    rnum_changed[rnum] = cnum;
	}
    }

    /*
     * Write back all contig records.
     * The io_c[lr]nbr records have been updated, but not the GContigs
     * structures.
     */
    for (i = 1; i <= NumContigs(io); i++) {
	GContigs c;
	contig_read(io, i, c);
	c.left  = io_clnbr(io, i);
	c.right = io_crnbr(io, i);
	c.length = io_clength(io, i);
	contig_write(io, i, c);
    }


    /* Copy consensus annotations */
    for (i = 1; i <= NumReadings(io); i++) {
	int j;
	if (rnum_changed[i] <= 0)
	    continue;

	/* New contig detected, so copy over the annotations */
	if (duplicate_tags)
	    copy_consensus_annotations(io, rnum2cnum[i], rnum_changed[i]);
	
	/* Prevent copy of annotations for other reads in this new contig */
	for (j = i+1; j <= NumReadings(io); j++)
	    if (rnum_changed[j] == rnum_changed[i])
		rnum_changed[j] = -1;
    }

    /* Delete contigs that have entirely vanished */
    delete_empty_contigs(io);

    /*
     * Write back all edited reading records.
     * The cached GReadings structures are up to date, but not the
     * io_* arrays. So we recalculate those too.
     */
    for (i = 1; i <= NumReadings(io); i++) {
	GReadings r;

	if (!rnum_changed[i])
	    continue;
	
	gel_read(io, i, r);
	io_lnbr(io, i) = r.left;
	io_rnbr(io, i) = r.right;
	io_length(io, i) = (r.sense ? -1 : 1) * r.sequence_length;
	io_relpos(io, i) = r.position;
	gel_write(io, i, r); /* Force write to sync mem & disk */
    }
    
    /* flush2t(io); */

    /* Remove contig holes. This may break contigs too */
    remove_contig_holes_all(io);

    flush2t(io);

    xfree(rnum2cnum);
    xfree(rnum_changed);
    xfree(new_cnum);

    return 0;
}


/*
 * ---------------------------------------------------------------------------
 * Contig tidying functions.
 *
 * If we break the logical consistency of the databse then these will try to
 * fix it. This is useful to do for simplicities sake. For example to enter
 * directed assembly data we can just slurp up all the data and then scan
 * through afterwards to make sure we have no contigs and that the contig
 * starts at base 1.
 */

/**
 * remove_contig_holes - checks for gaps in a contig at the start, end or
 * internal. Internal holes require splitting the contig in two, while start
 * and end gaps simply require adjustment of the length and shuffling sequence
 * positions.
 *
 * This function obtains all information from the disk (or cached) database 
 * structures rather than the internal arrays incase these are not
 * consistent (yet). It does however assume that the reading positions
 * are sorted left to right.
 *
 * Returns 0 for success
 *        -1 for error
 */
int remove_contig_holes(GapIO *io, int cnum) {
    int rnum, prev_rnum, new_contig;
    int furthest;
    int cstart;
    int shift;
    int first_loop = 1;
    GContigs c;
    GReadings r;
    int rm_contig, rm_left, rm_right;
    
    if (contig_read(io, cnum, c))
	return -1;

    /*
     * Firstly remove any annotations overhanging the left and right ends.
     * This can happen when we remove a reading at the very start or end of
     * contig.
     */
    if (c.annotations)
	c.annotations = rmanno(io, c.annotations, -INT_MAX, 0);
    if (c.annotations)
	c.annotations = rmanno(io, c.annotations, c.length+1, INT_MAX);

    do {
	new_contig = 0;
	
	rnum = c.left;
	prev_rnum = 0;
	cstart = 1;
	furthest = 1;
	shift = 0;

	while (rnum) {
	    if (gel_read(io, rnum, r))
		return -1;

	    /* Ensure the reading to contig number cache is up to date */
	    update_rnumtocnum(io, rnum, cnum);
	    
	    /* First read in contig => clip & shift consensus tags */
	    if (cstart) {
		shift = r.position - 1;
		if (shift && first_loop) {
		    /*
		     * This may be further sped up. Shifting tags may be
		     * slow (it writes back the results to disk), and if
		     * this contig has many holes then some tags will be 
		     * shifted multiple times.
		     */
		    c.annotations = rmanno(io, c.annotations, 1, shift+1);
		    contig_write(io, cnum, c);
		    shift_contig_tags(io, cnum, 1, -shift);
		}
	    }
	    
	    r.position -= shift;
	    io_relpos(io, rnum) -= shift;

	    /* If there's a gap - start a new contig */
	    if (!cstart && r.position > furthest) {
		new_contig = 1;
		break;
	    }

	    /* keep track of rightmost sequence end, to spot gaps */
	    if (furthest < r.position + r.sequence_length - 1)
		furthest = r.position + r.sequence_length - 1;

	    if (shift) {
		gel_write(io, rnum, r);
		io_relpos(io, rnum) = r.position;
	    }

	    prev_rnum = rnum;
	    rnum = r.right;
	    cstart = 0;
	}

	/*
	 * Last read => remove any annotations off contig end.
	 * Cannot do this here as we need to do the "if (new_contig)" and
	 * split_contig_tags() code below first, so we just record the 
	 * fact the tags should be removed.
	 */
	if (furthest < c.length) {
	    rm_left = furthest+1;
	    rm_right = c.length+1;
	    rm_contig = cnum;
	} else {
	    rm_contig = 0;
	    rm_left = rm_right = 0; /* avoids compiler warnings */
	}

	/* Update contig size etc */
	c.length = furthest;
	c.right = prev_rnum;
	contig_write(io, cnum, c);
	io_crnbr(io, cnum) = c.right;
	io_clength(io, cnum) = c.length;

	if (new_contig) {
	    int left_cnum = cnum;

	    vmessage("Breaking contig %s at reading %s\n",
		     io_rname(io, io_clnbr(io, cnum)),
		     io_rname(io, rnum));

	    if (-1 == io_init_contig(io, cnum = NumContigs(io)+1))
		return -1;

	    split_contig_tags(io, left_cnum, cnum, r.position, furthest);

	    contig_read(io, cnum, c);
	    c.left = rnum;
	    io_clnbr(io, cnum) = c.left;

	    /* Terminate read link list for existing contig */
	    gel_read(io, prev_rnum, r);
	    r.right = 0;
	    io_rnbr(io, prev_rnum) = 0;
	    gel_write(io, prev_rnum, r);

	    /* Start read link list for new contig */
	    gel_read(io, rnum, r);
	    r.left = 0;
	    io_lnbr(io, rnum) = 0;
	    gel_write(io, rnum, r);

	    update_rnumtocnum(io, rnum, cnum);
	}

	/* Now perform the tag removal, if detected as necessary earlier */
	if (rm_contig) {
	    GContigs ct;
	    contig_read(io, rm_contig, ct);
	    ct.annotations = rmanno(io, ct.annotations, rm_left, rm_right);
	    contig_write(io, rm_contig, ct);
	}

	first_loop = 0;
    } while (new_contig);

    return 0;
}

/**
 * Calls remove_contig_holes for all contigs in the DB.
 * Returns the same codes.
 *
 * Returns 0 for success
 *        -1 for error
 */
int remove_contig_holes_all(GapIO *io) {
    int i, ret = 0;
    for (i = 1; i <= NumContigs(io); i++) {
	ret |= remove_contig_holes(io, i);
    }

    return ret;
}


/**
 * Deletes an entire contig icluding all the readings, annotations, notes, etc
 * This is a wrapper around disassemble readings.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int delete_contig(GapIO *io, int contig) {
    int ret;
    int *reads, nreads, rnum;

    /*
     * Create an array to hold all reads in the contig.
     * We don't know how large this is up front, but as it's temporary
     * memory we just use the worst case of all readings in the
     * database.
     */
    if (NULL == (reads = (int *)xmalloc(NumReadings(io) * sizeof(int))))
	return -1;

    nreads = 0;
    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum))
	reads[nreads++] = rnum;

    /* Remove them all */
    ret = disassemble_readings(io, reads, nreads, 0, 0);

    xfree(reads);

    return ret;
}
