#include <math.h>
#include <stdio.h>

#include "gap_globals.h"
#include "clip.h"
#include "misc.h"
#include "qual.h"
#include "dna_utils.h"
#include "tagUtils.h"
#include "dis_readings.h"

/*
 * Perform quality clipping on an individual contig.
 * This identifies low quality regions at the ends of contigs and clips them
 *
 * Note that this can produce holes, which need to be fixed by the fix_holes()
 * function, and can also reorder the readings.
 */
static int quality_clip_contig(GapIO *io, int contig, int start, int end,
			       int qual_avg, int *old_start, int *old_end)
{
    int rnum;
    int win_len = 31; /* odd num */
    int i, i_end;
    int left, right;
    GReadings r;
    int1 *conf = NULL;
    int conf_alloc = 10000; /* realloced as needed */
    int total;
    int qual_tot, win_len2;
    
    qual_tot = win_len * qual_avg;
    win_len2 = win_len / 2;

    if (NULL == (conf = (int1 *)xmalloc(conf_alloc)))
	return -1;

    /* Find start reading */
    for (rnum = io_clnbr(io, contig);
	 io_relpos(io, rnum) < start;
	 rnum = io_rnbr(io, rnum))
	;

    /*
     * Loop through readings up to 'end' position.
     * The leftmost and rightmost sequences need to be treated specially
     * as we do not want to clip these. Doing so will change the length of the
     * consensus. Doing that means that we need to deal with consensus tags
     * (shifting, clipping and deleting them).
     */
    for (; rnum && io_relpos(io, rnum) <= end; rnum = io_rnbr(io, rnum)) {
	/* Find left and right clip points */
	gel_read(io, rnum, r);
	if (r.length < win_len)
	    continue;

	/* Realloc conf if necessary */
	if (r.length > conf_alloc) {
	    int1 *conf2;
	    conf_alloc = r.length + 100;
	    conf2 = (int1 *)xrealloc(conf, conf_alloc);
	    if (!conf2) {
		xfree(conf);
		return -1;
	    }
	    conf = conf2;
	}

	if (DataRead(io, r.confidence, conf, r.length * sizeof(*conf),
		     sizeof(*conf)))
	    continue;

	/* Clip left by quality */
	if (rnum != io_clnbr(io, contig)) {
	    total = 0;
	    for (i = 0; i < win_len; i++)
		total += conf[i];

	    if (total > qual_tot) {
		i = r.start;
	    } else {
		i_end = r.length - win_len2 - 1;
		i = win_len2 + 1;
		do {
		    total += conf[i + win_len2] - conf[i - win_len2 - 1];
		    i++;
		} while (total < qual_tot && i < i_end);
		i--;
	    }
	    left = MAX(i, r.start);
	} else {
	    left = r.start;
	}

	/* Clip right by quality */
	if (r.position + r.sequence_length <= io_clength(io, contig)) {
	    total = 0;
	    for (i = 0; i < win_len; i++)
		total += conf[r.length - i - 1];
	    
	    if (total > qual_tot) {
		i = r.end;
	    } else {
		i = r.length - win_len2 - 2;
		i_end = win_len2 + 1;
		do {
		    total += conf[i - win_len2] - conf[i + win_len2 + 1];
		    i--;
		} while (total < qual_tot && i > i_end);
		i++;
	    }
	    right = MIN(i, r.end);
	} else {
	    right = r.end;
	}

	if (left >= r.end - 1)
	    left = r.end - 2;
	if (right <= r.start + 1)
	    right = r.start + 2;

	if (right < left + 2)
	    right = left + 2;

	/* printf("Gel %d: L %d->%d    R %d->%d\n",
	       rnum, r.start, left, r.end, right); */

	r.position += left - r.start;
	old_start[rnum] = r.start;
	old_end[rnum] = r.end;
	r.start = left;
	r.end = right;
	r.sequence_length = r.end - r.start - 1;

	gel_write(io, rnum, r);
	io_relpos(io, rnum) = r.position;
	io_length(io, rnum) = r.sense ? -r.sequence_length : r.sequence_length;
    }

    xfree(conf);
    return 0;
}

/*
 * Perform quality clipping on an individual contig based on runs of Ns.
 * The idea is to clip of sections that are left in by phrap due to flukey
 * matches (an N match scores zero, so will not be rejected unless there are
 * other mismatches).
 *
 * Note that this can produce holes, which need to be fixed by the fix_holes()
 * function, and can also reorder the readings.
 */
static int N_clip_contig(GapIO *io, int contig, int start, int end,
			 int *old_start, int *old_end)
{
    int rnum;
    int i;
    int left, right;
    GReadings r;
    int scoretab[256];
    
    /* Find start reading */
    for (rnum = io_clnbr(io, contig);
	 io_relpos(io, rnum) < start;
	 rnum = io_rnbr(io, rnum))
	;

    for (i = 0; i < 256; i++) {
	scoretab[i] = 0; /* K, R, Y, etc */
    }
    scoretab['N'] = scoretab['n'] = scoretab['-'] = 1;
    scoretab['A'] = scoretab['a'] = -1;
    scoretab['C'] = scoretab['c'] = -1;
    scoretab['G'] = scoretab['g'] = -1;
    scoretab['T'] = scoretab['t'] = -1;

    /*
     * Loop through readings up to 'end' position.
     * The leftmost and rightmost sequences need to be treated specially
     * as we do not want to clip these. Doing so will change the length of the
     * consensus. Doing that means that we need to deal with consensus tags
     * (shifting, clipping and deleting them).
     */
    for (; rnum && io_relpos(io, rnum) <= end; rnum = io_rnbr(io, rnum)) {
	char *seq;

	gel_read(io, rnum, r);
	io_aread_seq(io, rnum, NULL, NULL, NULL,
			 &seq, NULL, NULL, 0);

	/* Clip left */
	if (rnum != io_clnbr(io, contig)) {
	    int score = 0, best_score = 0, best_pos = -1;
	    
	    for (i = r.start; i < r.end-1 && score >= -10; i++) {
		score += scoretab[(unsigned char)seq[i]];
		if (best_score <= score) {
		    best_score = score;
		    best_pos = i;
		}
	    }
	    left = (best_pos != -1) ? best_pos+1 : r.start;
	} else {
	    left = r.start;
	}

	/* Clip right */
	if (rnum != io_crnbr(io, contig)) {
	    int score = 0, best_score = 0, best_pos = -1;
	    
	    for (i = r.end-2; i >= r.start && score >= -10; i--) {
		score += scoretab[(unsigned char)seq[i]];
		if (best_score <= score) {
		    best_score = score;
		    best_pos = i;
		}
	    }
	    right = (best_pos != -1) ? best_pos+1 : r.end;
	} else {
	    right = r.end;
	}

	if (left >= r.end - 1)
	    left = r.end - 2;
	if (right <= r.start + 1)
	    right = r.start + 2;

	if (right < left + 2)
	    right = left + 2;

	if (r.start < left)
	    vmessage("Read #%d: clipping %d base%s from left end\n",
		     rnum, left-r.start, left-r.start ? "s" : "");

	if (r.end > right)
	    vmessage("Read #%d: clipping %d base%s from right end\n",
		     rnum, r.end - right, r.end - right  ? "s" : "");

	/*
	printf("Gel %d: L %d->%d    R %d->%d\n",
	       rnum, r.start, left, r.end, right);
	*/

	r.position += left - r.start;
	old_start[rnum] = r.start;
	old_end[rnum] = r.end;
	r.start = left;
	r.end = right;
	r.sequence_length = r.end - r.start - 1;

	gel_write(io, rnum, r);
	io_relpos(io, rnum) = r.position;
	io_length(io, rnum) = r.sense ? -r.sequence_length : r.sequence_length;

	xfree(seq);
    }

    return 0;
}

/*
 * Perform difference clipping on an individual contig.
 * This identifies reading regions that differs with the consensus and
 * tags and clips them accordingly.
 *
 * Note that this can produce holes, which need to be fixed by the fix_holes()
 * function, and can also reorder the readings.
 *
 * Returns the number of bases clipped
 */
#define SCORE_MATCH 1
#define SCORE_MISMATCH -2
static int difference_clip_contig(GapIO *io, int contig, int start, int end,
				int *old_start, int *old_end, int add_tags)
{
    int rnum;
    int win_len = 31;
    double disfract = 0.25;
    char *seq;
    int win_len2;
    int start2, end2, length2;
    int i, i_end, j;
    char *con;
    int count, worst_count;
    int lowest_dis, lowest_pos;
    int score, highest_score, highest_pos;
    int new_left, new_right;
    GReadings r;
    int nbases = 0;
    
    win_len2 = win_len / 2;

    /* Get consensus sequence */
    if (NULL == (con = (char *)xmalloc(io_clength(io, contig)+1)))
	return -1;
    calc_consensus(contig, 1, io_clength(io, contig), CON_SUM,
		   con, NULL,  NULL, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);
    

    /* Find start reading */
    for (rnum = io_clnbr(io, contig);
	 io_relpos(io, rnum) < start;
	 rnum = io_rnbr(io, rnum))
	;

    /* Loop through readings up to 'end' position */
    for (; rnum && io_relpos(io, rnum) <= end; rnum = io_rnbr(io, rnum)) {
	/*
	 * Find segment with the lowest difference count, averaged
	 * over 'win_len' bases.
	 */
	io_aread_seq(io, rnum, &length2, &start2, &end2, &seq, NULL, NULL, 0);
	if (length2 < win_len) {
	    if (seq) xfree(seq);
	    continue;
	}

	count = 0;
	worst_count = 0;

	if (win_len > end2 - start2) {
	    if (seq) xfree(seq);
	    continue; /* too short to compare over win_len */
	}

	for (j = io_relpos(io, rnum) - 1, i = start2;
	     i < start2 + win_len && i < end2 - 1; i++, j++) {
	    if (!same_char(seq[i], con[j]))
		count++;
	}
	i_end = end2 - 2 - win_len2;
	lowest_dis = count;
	lowest_pos = start2 + win_len2;
	for (i = start2 + win_len2, j = io_relpos(io, rnum) - 1 + win_len2;
	     i < i_end; i++, j++) {
	    if (count < lowest_dis) {
		lowest_dis = count;
		lowest_pos = i;
	    }
	    if (count > worst_count)
		worst_count = count;
	    count -= !same_char(seq[i-win_len2], con[j-win_len2]);
	    count += !same_char(seq[i+win_len2+1], con[j+win_len2+1]);
	}

	if (worst_count < disfract * win_len) {
	    if (seq) xfree(seq);
	    continue;
	}

	gel_read(io, rnum, r);

	/*
	 * Work outwards from lowest_pos in both directions finding the
	 * location where the alignment score is highest. We used a fixed
	 * score/penalty system.
	 */
	score = highest_score = 0;
	highest_pos = lowest_pos;
	for (j = io_relpos(io, rnum) - 1 + lowest_pos - start2, i = lowest_pos;
	     i < r.end - 1; i++, j++) {
	    if (same_char(seq[i], con[j]))
		score += SCORE_MATCH;
	    else
		score += SCORE_MISMATCH;
	    if (score > highest_score) {
		highest_score = score;
		highest_pos = i;
	    }
	}
	/* printf("Read %d: Best match at %d, extend R to %d (score %d)\n",
	       rnum, lowest_pos, highest_pos, highest_score); */
	new_right = highest_pos + 2;

	score = highest_score = 0;
	highest_pos = lowest_pos;
	for (j = io_relpos(io, rnum) - 1 + lowest_pos - start2, i = lowest_pos;
	     i >= r.start; i--, j--) {
	    if (same_char(seq[i], con[j]))
		score += SCORE_MATCH;
	    else
		score += SCORE_MISMATCH;
	    if (score > highest_score) {
		highest_score = score;
		highest_pos = i;
	    }
	}
	/* printf("Read %d: Best match at %d, extend L to %d (score %d)\n",
	       rnum, lowest_pos, highest_pos, highest_score); */
	new_left = highest_pos;

	if (new_right < new_left + 2)
	    new_right = new_left + 2;
	
	/* Tag the change */
	if (add_tags && new_left != r.start) {
	    char buf[100];
	    int tag_st;

	    if (r.sense) {
		tag_st = r.length - new_left + 1;
	    } else {
		tag_st = r.start + 1;
	    }
	    sprintf(buf, "Difference clipped from old start at %d\n", r.start);
	    insert_NEW_tag(io, rnum, tag_st, new_left - r.start, "DIFF",
			   buf, 2);
	}

	if (add_tags && new_right != r.end) {
	    char buf[100];
	    int tag_st;

	    if (r.sense) {
		tag_st = r.length - r.end + 2;
	    } else {
		tag_st = new_right;
	    }
	    sprintf(buf, "Difference clipped from old end at %d\n", r.end);
	    insert_NEW_tag(io, rnum, tag_st, r.end - new_right, "DIFF",
			   buf, 2);
	}

	/* Update reading positions/length */
	gel_read(io, rnum, r);
	r.position += new_left - r.start;
	old_start[rnum] = r.start;
	old_end[rnum] = r.end;
	r.start = new_left;
	r.end = new_right;
	nbases += r.sequence_length - (r.end - r.start - 1);
	r.sequence_length = r.end - r.start - 1;
	gel_write(io, rnum, r);
	io_relpos(io, rnum) = r.position;
	io_length(io, rnum) = r.sense ? -r.sequence_length : r.sequence_length;

	if (seq) xfree(seq);
    }

    xfree(con);
    return nbases;
}


typedef struct s_reads_t {
    int position;
    int num;
} reads_t;

/*
 * Our local reading sort function. Sorts on position.
 */
static int sort_reads(const void *r1, const void *r2)
{
    return ((reads_t *)r1)->position - ((reads_t *)r2)->position;
}


/*
 * Ensures that the order of readings is correct. This also updates the
 * contig length and shifts readings if the leftmost reading no longer
 * starts at base 1.
 */
static void reorder_readings(GapIO *io, int contig)
{
    reads_t *reads;
    int rnum;
    int i;
    int count;
    int left;
    int length;
    int shift;
    GReadings r;
    GContigs c;
    
    if (NULL == (reads = (reads_t *)xmalloc(NumReadings(io)*sizeof(reads_t))))
	return;

    for (i = 0, rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	reads[i].position = io_relpos(io, rnum);
	reads[i].num = rnum;
	i++;
    }
    count = i;

    qsort(reads, count, sizeof(*reads), sort_reads);

    left = length = 0;
    shift = reads[0].position - 1;
    for (i = 0; i < count; i++) {
	int old_left, old_right;

	rnum = reads[i].num;
	gel_read(io, rnum, r);

	old_left  = r.left;  r.left = left;
	old_right = r.right; r.right = i < count-1 ? reads[i+1].num : 0;
	r.position -= shift;

	io_lnbr(io, rnum) = r.left;
	io_rnbr(io, rnum) = r.right;
	io_relpos(io, rnum) = r.position;

	if (shift || old_left != r.left || old_right != r.right) {
	    gel_write(io, rnum, r);
	}

	if (r.position + r.sequence_length - 1 > length)
	    length = r.position + r.sequence_length - 1;

	left = rnum;
    }

    contig_read(io, contig, c);
    c.length = length;
    c.left = reads[0].num;
    c.right = reads[count-1].num;
    io_clength(io, contig) = c.length;
    io_clnbr(io, contig) = c.left;
    io_crnbr(io, contig) = c.right;
    contig_write(io, contig, c);

    xfree(reads);
}

/*
 * Fills 'errs' with the accumulated number of errors since the reading's
 * cutoff end/start. 'errs' is allocated to be size of hole_end-hole_start+1.
 */
static void get_acc_errs_l(GapIO *io, double *errs, int rnum, int old_end,
			   int hole_start, int hole_end)
{
    int1 *conf = NULL;
    int i, i_end, j;
    double accum;
    GReadings r;

    /* get confidence values */
    io_aread_seq(io, rnum, NULL, NULL, NULL,  NULL, &conf, NULL, 0);

    /* Find number of errors before the hole */
    gel_read(io, rnum, r);
    accum = 0.0;
    i_end = hole_start - r.position + r.start;
    for (i = r.end - 1; i < i_end; i++) {
	accum += pow(10.0, conf[i] / -10.0);
    }

    /* Fill errs */
    i_end = hole_end - r.position + r.start + 1;
    for (j = 0; i < i_end && i < old_end - 1; i++, j++) {
	accum += pow(10.0, conf[i] / -10.0);
	errs[j] = accum;
    }

    /* Fill any remaining bits of errs with a very high value */
    for (; j < hole_end - hole_start + 1; j++)
	errs[j] = 99999.0;

    if (conf)
	xfree(conf);
}

static void get_acc_errs_r(GapIO *io, double *errs, int rnum, int old_start,
			   int hole_start, int hole_end)
{
    int1 *conf = NULL;
    int i, i_end, j;
    double accum;
    GReadings r;

    /* get confidence values */
    io_aread_seq(io, rnum, NULL, NULL, NULL,  NULL, &conf, NULL, 0);

    /* Find number of errors after the hole */
    gel_read(io, rnum, r);
    accum = 0.0;
    i_end = hole_end - r.position + r.start;
    for (i = r.start - 1; i > i_end; i--) {
	accum += pow(10.0, conf[i] / -10.0);
    }

    /* Fill errs */
    i_end = hole_start - r.position + r.start - 1;
    for (j = hole_end - hole_start; i > i_end && i > old_start - 1; i--, j--) {
	accum += pow(10.0, conf[i] / -10.0);;
	errs[j] = accum;
    }

    /* Fill any remaining bits of errs with a very high value */
    for (; j >= 0; j--)
	errs[j] = 99999.0;

    if (conf)
	xfree(conf);
}

static void adjust_extensions(GapIO *io, int hole_start, int hole_end,
			      double **left_errs, int *left_rnums,
			      int left_count,
			      double **right_errs, int *right_rnums,
			      int right_count)
{
    int i;
    int h;
    int best_left, best_right;
    int best_hole_left = 0, best_hole_right = 0, best_hole_pos;
    double best_left_val, best_right_val, best_hole_val;
    GReadings r;
    int hole_size = hole_end - hole_start + 1;
    int ext_len;

    /*
     * Brute force search through all possible 'joint' positions in the hole
     * (unless either left_count or right_count is 0)
     */
    if (left_count && right_count) {
	best_hole_val = 1e20;
	best_hole_left = best_hole_right = best_hole_pos = 0;
	for (h = 0; h < hole_size; h++) {
	    best_left_val = best_right_val = 1e10;
	    best_left = best_right = 0;
	    
	    /* Pick left and right extension with fewest errors */
	    for (i = 0; i < left_count; i++) {
		if (left_errs[i][h] < best_left_val) {
		    best_left_val = left_errs[i][h];
		    best_left = i;
		}
	    }
	    for (i = 0; i < right_count; i++) {
		if (right_errs[i][h] < best_right_val) {
		    best_right_val = right_errs[i][h];
		    best_right = i;
		}
	    }
	    
	    if (best_left_val + best_right_val < best_hole_val) {
		best_hole_val = best_left_val + best_right_val;
		best_hole_left = best_left;
		best_hole_right = best_right;
		best_hole_pos = h;
	    }
	}

    } else if (left_count == 0) {
	/* Pick best right only */
	h = 0;
	best_hole_right = 0;
	best_right_val = 1e10;
	
	if (hole_size != 0) {
	    for (i = 0; i < right_count; i++) {
		if (right_errs[i][h] <  best_right_val) {
		    best_right_val = right_errs[i][h];
		    best_right = i;
		}
	    }
	} else {
	    best_right = 0;
	}

	best_hole_val = best_right_val;
	best_hole_pos = -1;

    } else {
	/* Pick best left only */
	h = hole_end - hole_start;
	best_hole_left = 0;
	best_left_val = 1e10;

	if (hole_size != 0) {
	    for (i = 0; i < left_count; i++) {
		if (left_errs[i][h] < best_left_val) {
		    best_left_val = left_errs[i][h];
		    best_left = i;
		}
	    }
	} else {
	    best_left = 0;
	}

	best_hole_val = best_left_val;
	best_hole_pos = hole_end - hole_start + 1;
    }

    if (left_count && right_count)
	vmessage(" extend #%d and #%d with %f expected errors\n",
		 left_rnums[best_hole_left], right_rnums[best_hole_right],
		 best_hole_val);
    else if (left_count)
	vmessage(" extend #%d with %f expected errors\n",
		 left_rnums[best_hole_left], best_hole_val);
    else
	vmessage(" extend #%d with %f expected errors\n",
		 right_rnums[best_hole_right], best_hole_val);

    if (left_count) {
	gel_read(io, left_rnums[best_hole_left], r);
	ext_len = best_hole_pos + hole_start -
	    (r.position + r.sequence_length) + 1;
	r.end += ext_len;
	r.sequence_length += ext_len;
	gel_write(io, left_rnums[best_hole_left], r);
	io_length(io, left_rnums[best_hole_left]) =
	    r.sense ? -r.sequence_length : r.sequence_length;
    }

    if (right_count) {
	gel_read(io, right_rnums[best_hole_right], r);
	ext_len = hole_size - best_hole_pos + r.position - hole_end - 1;
	r.start -= ext_len;
	r.position -= ext_len;
	r.sequence_length += ext_len;
	gel_write(io, right_rnums[best_hole_right], r);
	io_relpos(io, right_rnums[best_hole_right]) = r.position;
	io_length(io, right_rnums[best_hole_right]) =
	    r.sense ? -r.sequence_length : r.sequence_length;
    }
}

/*
 * Finds holes in contigs (gaps where no sequence covers a region) and extends
 * the highest quality segment over that hole. old_start and old_end contain
 * the previous start and end quality clips. We make sure that we never extend
 * beyond these, even if it gives good quality. The reason for this is two
 * fold. Firstly it means we no the data is aligned that far (and hence save
 * ourselves a lot of work), and secondly it avoids worrying about vector tags.
 *
 * Our strategy, once we've found a hole, is to find all the readings that
 * extend over this hole from the left side and the right side. For these
 * readings we produce an array of length `size of hole' of the accumulated
 * expected number of errors to extend to that position in the hole.
 * Next, for each position in the hole, we pick the lowest error total from
 * readings from the left side of the hole and the lowest from the right side
 * and add these together. This gives us an estimated number of errors in the
 * consensus that we would create by extending these two readings to fill the
 * hole. The optimal solution is then the lowest summation; we have one sum
 * per position in the hole - ie `size of hole').
 */
static void fix_holes(GapIO *io, int contig, int *old_start, int *old_end)
{
    int rnum;
    int hole_start = 2, hole_end, hole_size;
    int left_count, right_count;
    int *left_rnums, *right_rnums, left_rnums_size, right_rnums_size;
    double **left_errs, **right_errs, *left_errs2, *right_errs2;
    int i;
    int max_len = find_max_gel_len(io, contig, 0)+10;
    GReadings r;

    left_rnums_size = right_rnums_size = 16;
    if (NULL == (left_rnums = xmalloc(left_rnums_size * sizeof(int))))
	return;
    if (NULL == (right_rnums = xmalloc(right_rnums_size * sizeof(int))))
	return;

    /* Find holes */
    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	if (io_relpos(io, rnum) >= hole_start) {
	    hole_end = io_relpos(io, rnum) - 1;
	    hole_size = hole_end - hole_start + 1;
	    vmessage("Hole from %d to %d:", hole_start, hole_end);

	    /* Find overlapping left reads */
	    left_count = right_count = 0;
	    for (i = io_lnbr(io, rnum);
		 i && io_relpos(io, i) > hole_start - max_len;
		 i = io_lnbr(io, i)) {
		gel_read(io, i, r);
		if (r.position - r.start + old_end[i] - 1 > hole_start) {
		    left_rnums[left_count++] = i;
		    if (left_count >= left_rnums_size) {
			left_rnums_size *= 2;
			left_rnums = xrealloc(left_rnums,
					      left_rnums_size * sizeof(int));
			if (NULL == left_rnums)
			    return;
		    }
		}
	    }

	    /* Find overlapping right reads */
	    for (i = rnum;
		 i && io_relpos(io, i) < hole_start + max_len;
		 i = io_rnbr(io, i)) {
		gel_read(io, i, r);
		if (r.position - r.start + old_start[i] <= hole_end) {
		    right_rnums[right_count++] = i;
		    if (right_count >= right_rnums_size) {
			right_rnums_size *= 2;
			right_rnums = xrealloc(right_rnums,
					       right_rnums_size * sizeof(int));
			if (NULL == right_rnums)
			    return;
		    }
		}
	    }

	    if (left_count == 0 && right_count == 0) {
		verror(ERR_WARN, "quality_clip", "Hole in contig");
		continue;
	    }

	    /* Alloc accumulated error rate arrays - one per overlapping seq */
	    left_errs = (double **)xmalloc(left_count*sizeof(*left_errs)+1);
	    right_errs = (double **)xmalloc(right_count*sizeof(*right_errs)+1);
	    left_errs2 = (double *)xmalloc(left_count * hole_size *
					   sizeof(*left_errs2) + 1);
	    right_errs2 = (double *)xmalloc(right_count * hole_size *
					    sizeof(*right_errs2) + 1);
	    if (NULL == left_errs || NULL == right_errs ||
		NULL == left_errs2 || NULL == right_errs2)
		return;

	    for (i = 0; i < left_count; i++)
		left_errs[i] = left_errs2 + i * hole_size;
	    for (i = 0; i < right_count; i++)
		right_errs[i] = right_errs2 + i * hole_size;

	    /* Fill the accumulated error rate arrays */
	    for (i = 0; i < left_count; i++) {
		get_acc_errs_l(io, left_errs[i], left_rnums[i],
			       old_end[left_rnums[i]],
			       hole_start, hole_end);
	    }
	    for (i = 0; i < right_count; i++) {
		get_acc_errs_r(io, right_errs[i], right_rnums[i],
			       old_start[right_rnums[i]],
			       hole_start, hole_end);
	    }

	    /* Find the best left and right extensions and adjust them */
	    adjust_extensions(io, hole_start, hole_end,
			      left_errs, left_rnums, left_count,
			      right_errs, right_rnums, right_count);

	    xfree(left_errs);
	    xfree(right_errs);
	    xfree(left_errs2);
	    xfree(right_errs2);
	}

	if (io_relpos(io, rnum) + ABS(io_length(io, rnum)) > hole_start)
	    hole_start = io_relpos(io, rnum) + ABS(io_length(io, rnum));
    }

    xfree(left_rnums);
    xfree(right_rnums);
}

void quality_clip(GapIO *io, int num_contigs, contig_list_t *cl, int qual_avg)
{
    int i;
    int *old_start, *old_end;

    if (NULL == (old_start = (int *)xcalloc(NumReadings(io)+1, sizeof(int))))
	return;
    if (NULL == (old_end = (int *)xcalloc(NumReadings(io)+1, sizeof(int))))
	return;

    for (i = 0; i < num_contigs; i++) {
	quality_clip_contig(io, cl[i].contig, cl[i].start, cl[i].end, qual_avg,
			    old_start, old_end);
	reorder_readings(io, cl[i].contig);
	fix_holes(io, cl[i].contig, old_start, old_end);
	reorder_readings(io, cl[i].contig);
	flush2t(io);
    }

    xfree(old_start);
    xfree(old_end);
}


void N_clip(GapIO *io, int num_contigs, contig_list_t *cl)
{
    int i;
    int *old_start, *old_end;

    if (NULL == (old_start = (int *)xcalloc(NumReadings(io)+1, sizeof(int))))
	return;
    if (NULL == (old_end = (int *)xcalloc(NumReadings(io)+1, sizeof(int))))
	return;

    for (i = 0; i < num_contigs; i++) {
	N_clip_contig(io, cl[i].contig, cl[i].start, cl[i].end,
		      old_start, old_end);
	reorder_readings(io, cl[i].contig);
	fix_holes(io, cl[i].contig, old_start, old_end);
	reorder_readings(io, cl[i].contig);
	flush2t(io);
    }

    xfree(old_start);
    xfree(old_end);
}


void difference_clip(GapIO *io, int num_contigs, contig_list_t *cl,
		     int add_tags)
{
    int i, nbases;
    int *old_start, *old_end;
    int all;

    if (NULL == (old_start = (int *)xcalloc(NumReadings(io)+1, sizeof(int))))
	return;
    if (NULL == (old_end = (int *)xcalloc(NumReadings(io)+1, sizeof(int))))
	return;

    for (i = 0; i < num_contigs; i++) {
	vmessage("--Contig %s--\n", io_rname(io, io_clnbr(io, cl[i].contig)));
	all = (cl[i].start == 1 && cl[i].end == io_clength(io, cl[i].contig))
	    ? 1 : 0;
	nbases = difference_clip_contig(io, cl[i].contig, cl[i].start,
					cl[i].end, old_start, old_end,
					add_tags);
	/*
	 * if (all)
	 *   nbases += quality_clip_ends(io, cl[i].contig, add_tags);
	 */
	vmessage("  Clipped %d bases\n", nbases);
	reorder_readings(io, cl[i].contig);
	fix_holes(io, cl[i].contig, old_start, old_end);
	reorder_readings(io, cl[i].contig);
	flush2t(io);
	vmessage("\n");
    }

    xfree(old_start);
    xfree(old_end);
}

/*
 * Starting from the left end, this identifies the point at which the
 * average quality (over win_len bases) is >= avg_qual. This base number is
 * returned.
 * Returns -1 for failure.
 */
static int avg_clip(GapIO *io, int1 *conf, int conf_len,
		    int win_len, int avg_qual) {
    int total, qual_tot;
    int win_len2;
    int i, i_end;

    /* Win_len should be an odd number */
    if ((win_len & 1) == 0)
	win_len |= 1;
    win_len2 = win_len/2;

    if (conf_len <= win_len)
	return -1;

    total = 0;
    for (i = 0; i < win_len; i++)
	total += conf[i];
    qual_tot = avg_qual * win_len;

    if (total >= qual_tot)
	return -1;

    i_end = conf_len - win_len2 - 1;
    i = win_len2 + 1;
    do {
	total += conf[i + win_len2] - conf[i - win_len2 - 1];
	i++;
    } while (total < qual_tot && i < i_end);
    i--;

    if (win_len > 5) {
	int c = avg_clip(io, &conf[i - win_len/2],
			 win_len*2+1, win_len-2, avg_qual);
	if (c != -1)
	    i = i-win_len/2 + c;
    }
    
    return i;
}

/*
 * Clips the contig end back based on an average quality value.
 * We only look for the leftmost and rightmost reading. The furthest these
 * can be clipped back is to where the next reading is found (ie depth=2).
 */
void quality_clip_ends(GapIO *io, int contig, int avg_qual)
{
    int rnum, rnum1, rnum2;
    int rpos1, rpos2;
    GReadings r;
    int1 *conf = 0;
    int conf_len;
    int clip;
    int i;
    int win_len = 31;
    int printed_cname = 0;
    int clipped;

    /* --- Left end --- */
    /* Identify 2 leftmost sequences */
    rnum1 = io_clnbr(io, contig);
    rnum2 = io_rnbr(io, rnum1);

    /* Load confidence */
    gel_read(io, rnum1, r);
    if (NULL == (conf = (int1 *)xcalloc(r.length, sizeof(*conf))))
	return;
    conf_len = r.length;
    if (DataRead(io, r.confidence, conf, r.length * sizeof(*conf),
		 sizeof(*conf)))
	return;

    /* Clip sequence */
    clip = avg_clip(io, conf, conf_len, win_len, avg_qual)+2;
    if (clip-2 != -1 && clip > r.start && rnum2) {
	int apos = r.position - r.start + clip;
	if (apos > io_relpos(io, rnum2)) {
	    clip -= apos - io_relpos(io, rnum2);
	    apos = r.position - r.start + clip;
	}

	r.position = r.position - r.start + clip;
	r.start = clip;
	clipped = r.sequence_length - (r.end - r.start - 1);
	if (clipped) {
	    vmessage("Contig %s     ", io_rname(io, io_clnbr(io, contig)));
	    printed_cname = 1;
	    vmessage("clip %d from left     ", clipped);
	}
	r.sequence_length = r.end - r.start - 1;
	io_relpos(io, rnum1) = r.position;
	io_length(io, rnum1) = r.sequence_length * (r.sense ? -1 : 1);
	gel_write(io, rnum1, r);
    }
    xfree(conf);


    /* --- Right end --- */
    rnum1 = rnum = io_crnbr(io, contig);
    gel_read(io, rnum, r);
    rpos1 = r.position + r.sequence_length-1;
    rpos2 = 0;
    rnum2 = 0;

    /* Identify 2 rightmost sequences */
    while (rnum = io_lnbr(io, rnum)) {
	gel_read(io, rnum, r);
	/*
	 * Optimisation - REAL traces are not (yet) going to be longer than
	 * 2Kb, so stop looking at that point.
	 */
	if (io_clength(io, contig) - r.position >= 2000)
	    break;

	if (rpos1 <= r.position + r.sequence_length-1) {
	    rpos1  = r.position + r.sequence_length-1;
	    rnum2 = rnum1;
	    rpos2 = rpos1;
	    rnum1 = rnum;
	} else if (rpos2 < r.position + r.sequence_length-1) {
	    rpos2 = r.position + r.sequence_length-1;
	    rnum2 = rnum;
	}
    }
	
    /* Load and reverse confidence array */
    gel_read(io, rnum1, r);
    if (NULL == (conf = (int1 *)xcalloc(r.length, sizeof(*conf))))
	return;
    conf_len = r.length;
    if (DataRead(io, r.confidence, conf, r.length * sizeof(*conf),
		 sizeof(*conf)))
	return;
    for (i = 0; i < r.length/2; i++) {
	int tmp;
	tmp = conf[i];
	conf[i] = conf[r.length-1-i];
	conf[r.length-1-i] = tmp;
    }

    /* Compute clip point (in original orientation) */
    clip = avg_clip(io, conf, conf_len, win_len, avg_qual);
    if (clip != -1)
	clip = r.length - clip;

    /* Clip */
    if (clip != -1 && ++clip < r.end && rnum2) {
	int apos = r.position - r.start + clip;
	if (apos < io_relpos(io, rnum2) + ABS(io_length(io, rnum2))-1) {
	    clip-= apos - (io_relpos(io, rnum2) + ABS(io_length(io, rnum2))-1);
	    apos = r.position - r.start + clip;
	}

	r.end = clip+2;
	clipped = r.sequence_length - (r.end - r.start - 1);
	if (clipped) {
	    if (!printed_cname)
		vmessage("Contig %s     ", io_rname(io, io_clnbr(io, contig)));
	    vmessage("clip %d from right",
		     r.sequence_length - (r.end - r.start - 1));
	    printed_cname = 1;
	}
	r.sequence_length = r.end - r.start - 1;
	io_length(io, rnum1) = r.sequence_length * (r.sense ? -1 : 1);
	gel_write(io, rnum1, r);
    }

    if (printed_cname)
	vmessage("\n");

    xfree(conf);

    remove_contig_holes(io, contig);
}

/*

set io [open_db -name F35H8 -version 1 -access rw]
quality_clip -io $io -contigs {=6 =7} -quality 35
close_db -io $io
exit

*/
