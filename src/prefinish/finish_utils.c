#include <stdio.h>

#include "IO.h"
#include "misc.h"
#include "finish.h"
#include "finish_utils.h"
#include "finish_hash.h"
#include "tagUtils.h"

int finish_next_expt_id(int reset) {
    static int expt_num = 1;
    if (reset)
	expt_num = reset;

    return expt_num++;
}

int finish_next_group_id(int reset) {
    static int group_num = 1;
    if (reset)
	group_num = reset;

    return group_num++;
}

/*
 * Computes the average expected length for a sequence with position 'start'
 * to 'end' on a template with minimum and maximum size ranges from
 * 't_start1'..'t_end1' and 't_start2'..'t_end2'.
 * All ranges are inclusive.
 *
 * When the sequence lies entirely within the minimum template size then
 * the expected length is the length we've specified (start to end). Otherwise
 * there is a chance that the template may be short (> min, but < max) and so
 * the sequence may be terminated prematurely. We cannot possibly know when
 * this will happen, but we can guess at an average case.
 *
 * Returns:
 *	An average length as the return value. new_start and new_end are
 *	filled in with the new expected start and end values, which also
 *	takes into account template 'clipping'.
 */
int finish_avg_length(int start, int end, int dir,
		      int t_start1, int t_end1,
		      int t_start2, int t_end2,
		      int *new_start, int *new_end) {
    int i;
    double len = 0;

    /* Inefficient, but simple to understand */
    for (i = start; i <= end; i++) {
	if (i >= t_start1 && i <= t_end1) {
	    len++;
	} else if (i < t_start1 && i >= t_start2) {
	    len += (i-t_start2+1)/(double)(t_start1-t_start2+1);
	} else if (i > t_end1 && i <= t_end2) {
	    len += (t_end2-i+1)/(double)(t_end2-t_end1+1);
	} /* else off end of template */
    }
    
    /* Clip by template */
    if (start < t_start2) start = t_start2;
    if (start > t_end2)   start = t_end2;
    if (end   < t_start2) end   = t_start2;
    if (end   > t_end2)	  end   = t_end2;

    *new_start = start;
    *new_end   = end;

    /* Update new length */
    if (dir == 1) {
	*new_end   = start + (len-1);
    } else {
	*new_start = end - (len-1);
    }

    return (int)len;
}

/*
 * Computes the chance (0 to 1) of all the sequence at position start..end
 * existing on a template spanning between t_start1..t_end1 (minimum) and
 * t_start2..t_end2 (maximum).
 *
 * Returns chance (1 == always, 0 == never, plus between values)
 */
double template_exists_chance(int start, int end,
			      int t_start1, int t_end1,
			      int t_start2, int t_end2) {
    double err[2], emax;
    int i, j;

    /* Compute error at both ends */
    for (j = 0; j < 2; j++) {
	i = j ? end : start;

	if (i >= t_start1 && i <= t_end1) {
	    /* ok */
	} else if (i < t_start1 && i >= t_start2) {
	    err[j] = 1 - (i-t_start2+1)/(double)(t_start1-t_start2+1);
	} else if (i > t_end1 && i <= t_end2) {
	    err[j] = 1 - (t_end2-i+1)/(double)(t_end2-t_end1+1);
	} else {
	    err[j] = 1;
	}
    }

    /* Take max error */
    emax = MAX(err[0], err[1]);

    return emax > 1 ? 0 : 1 - emax;
}


/*
 * complement_seq_qual_mapping
 *
 * Complements a sequence, quality buffer, and integer array (which is
 * used externally for mapping unpadded to padded positions, but for the
 * purposes of this code it could contain anything).
 *
 * Arguments:
 *	len		Dimension of seq, qual and map arrays
 *	seq		The DNA sequence (does not need to be nul terminated)
 *	qual		The confidence values
 *	map		The integer mapping array.
 */
void complement_seq_qual_mapping(int len, char *seq, float *qual, int *map) {
    char tmp_seq;
    float tmp_qual;
    int tmp_map;
    int i, j;
    static char complementary_base[256];
    static int initialised = 0;
    
    if (!initialised) {
	for (i = 0; i < 256; i++)
	    complementary_base[i] = i;

	complementary_base['a'] = 't';
	complementary_base['A'] = 'T';
	complementary_base['c'] = 'g';
	complementary_base['C'] = 'G';
	complementary_base['g'] = 'c';
	complementary_base['G'] = 'C';
	complementary_base['t'] = 'a';
	complementary_base['T'] = 'A';
	complementary_base['u'] = 'a';
	complementary_base['U'] = 'A';

	initialised = 1;
    }

    for (i = 0, j = len-1; i <= j; i++, j--) {
	tmp_seq  = complementary_base[(unsigned)seq[i]];
	seq[i]   = complementary_base[(unsigned)seq[j]];
	seq[j]   = tmp_seq;

	tmp_qual = qual[i];
	qual[i]  = qual[j];
	qual[j]  = tmp_qual;

	tmp_map  = map[i];
	map[i]   = map[j];
	map[j]   = tmp_map;
    }
}

/*
 * Looks for cloning vector. It's useful to know if its there so that we
 * don't allow extending into vector. If svec_also is true then this also
 * treats SVEC as the same as CVEC. The reason for this is to treat sequencing
 * vector as the clone-end for EST projects.
 *
 * Returns:
 *	Updates contents of left and right to be 0 for no vector and 1 for
 * 	vector present.
 */
void find_cloning_vector(GapIO *io, int contig, int *left, int *right,
			 int svec_also) {
    int gel;
    int cvec_l = 0, cvec_r = 0;
    char *tag_types[] = {"CVEC", "SVEC"};
    GAnnotations *a;
    GReadings r;

    /* left end - within 5 bases */
    for (gel = io_clnbr(io, contig); gel; gel = io_rnbr(io, gel)) {
	gel_read(io, gel, r);
	if (r.position - r.start > 1)
	    break;

	a = vtagget(io, gel, 1 + (svec_also>0), tag_types);
	while (a && a != (GAnnotations *)-1) {
	    int start;

	    start = r.position - r.start +
		(r.sense
		 ? r.length - (a->length + a->position - 1)
		 : a->position - 1);

	    if (start <= 5) {
		cvec_l = 1;
		break;
	    }

	    a = vtagget(io, 0,  1 + (svec_also>0), tag_types);
	}
    }

    /* Right end - within 5 bases */
    for (gel = io_crnbr(io, contig); gel; gel = io_lnbr(io, gel)) {
	gel_read(io, gel, r);
	if (r.position < io_clength(io, contig) - max_gel_len(io))
	    break;

	a = vtagget(io, gel, 1 + (svec_also>0), tag_types);
	while (a && a != (GAnnotations *)-1) {
	    int end;

	    end = r.position - r.start +
		(r.sense
		 ? r.length - (a->length + a->position - 1)
		 : a->position - 1)
		+ a->length - 1;

	    if (end + 5 >= io_clength(io, contig))
		cvec_r = 1;

	    a = vtagget(io, 0, 1 + (svec_also>0), tag_types);
	}
    }

    if (cvec_l)
	puts("Cloning vector detected at left end of contig");

    if (cvec_r)
	puts("Cloning vector detected at right end of contig");

    *left = cvec_l;
    *right = cvec_r;
}

/*
 * Clip suggested start and end points for new sequences based on the known
 * vector found in an existing sequence. In the case where vector is found
 * in the middle of s_start to s_end, but not either side, s_start and s_end
 * are set to be just one end (the start<->vector).
 *
 * Arguments:
 *	io		Gap IO handle
 *	s_start		Current suggestion for new sequence start in contig
 *	s_end		Current suggestion for new sequence end in contig
 *	rnum		Reading number to compare against
 *
 * Returns:
 *	Void return, but updates s_start and s_end.
 */
void finish_clip_svec(GapIO *io, int *s_start, int *s_end, int rnum) {
    char *tag_types[] = {"SVEC"};
    GAnnotations *a;
    GReadings r;

    gel_read(io, rnum, r);

    a = vtagget(io, rnum, 1, tag_types);
    while (a && a != (GAnnotations *)-1) {
	int start, end;

	start = r.position - r.start +
	    (r.sense
	     ? r.length - (a->length + a->position - 1)
	     : a->position - 1);
	end = start + a->length - 1;

	if (start <= *s_start && end >= *s_end) {
	    /* ...VVVVVVVVVVVVVVVVVVVVVVVVV... */
	    /* --------start----end----------- */
	    *s_end = *s_start;

	} else if (start <= *s_start && end >= *s_start) {
	    /* ...VVVVVVVVVVV................. */
	    /* --------start----end----------- */
	    *s_start = end + 1;

	} else if (start <= *s_end && end >= *s_end) {
	    /* ..............VVVVVVVVVVV...... */
	    /* --------start----end----------- */
	    *s_end = start - 1;

	} else if (start >= *s_start && start <= *s_end &&
		   end >= *s_start && end <= *s_end) {
	    /* .............VVV............... */
	    /* --------start----end----------- */
	    *s_end = *s_start;
	}

	a = vtagget(io, 0, 1, tag_types);
    }
}

/*
 * Finds the readings covering a specific consensus base
 *
 * Arguments:
 *	io		Gap IO handle
 *	contig		Contig number
 *	pos		Contig position
 *
 * Returns:
 *	A malloc array of integer reading numbers. It is up to the caller
 *	to free this array with xfree.
 *	NULL for failure
 */
int *seqs_at_pos(GapIO *io, int contig, int pos) {
    int rnum;
    int *seqs = xmalloc(8 * sizeof(int)); /* 7 seqs plus 0 terminator */
    int seqs_dim = 8;
    int count = 0;

    if (!seqs)
	return NULL;

    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	if (io_relpos(io, rnum) + ABS(io_length(io, rnum)) - 1 < pos)
	    continue;
	if (io_relpos(io, rnum) > pos)
	    break;
	if (count >= seqs_dim-1 /* -term */) {
	    seqs_dim *= 2;
	    if (NULL == (seqs = xrealloc(seqs, seqs_dim * sizeof(int))))
		return NULL;
	}
	seqs[count++] = rnum;
    }
    seqs[count++] = 0;

    return seqs;
}

/*
 * From a template (zero for none), finds a reading appropriate for
 * tagging (positions start to end).
 *
 * Ideally we choose a reading using that template.
 * If none exists then we choose any reading.
 * If this doesn't work then the tag gets truncated.
 *
 * Arguments:
 *	io		Gap IO handle
 *	contig		Contig number
 *	template	Template number (0 for none)
 *	start		First base in tag (offset into contig)
 *	end		Last base in tag
 *
 * Returns:
 *	Success: Template number. Also may modify start and end values.
 *	Failure: 0
 */
int tag_template(GapIO *io, int contig, int template,
		 int *start, int *end) {
    int *seqs;
    int i;
    int best1 = 0; /* same template */
    int best2 = 0; /* tag fits */
    int best3 = 0; /* tag truncated */
    int new_end = *start;

    /* Find sequences covering start of tag */
    if (!(seqs = seqs_at_pos(io, contig, *start)))
	return 0;

    for (i = 0; seqs[i]; i++) {
	GReadings r;
	int in_seq;

	gel_read(io, seqs[i], r);

	in_seq = (r.position + r.sequence_length - 1 >= *end) ? 1 : 0;
	
	if (!best1 && in_seq && (r.template == template)) {
	    best1 = seqs[i];
	}

	if (!best2 && in_seq) {
	    best2 = seqs[i];
	}
	
	if (new_end < r.position + r.sequence_length - 1) {
	    best3 = seqs[i];
	    new_end = r.position + r.sequence_length - 1;
	}
    }

    xfree(seqs);

    if (best1)
	return best1;
    if (best2)
	return best2;
    *end = new_end;
    return best3;
}

/*
 * Finds where a template has a duplicate or not.
 *
 * Returns 1 for yes,
 *         0 for no.
 */
int template_is_dup(finish_t *fin, int *templates_picked,
		    int num_picked,
		    int template) {
    int dup = 0;
    int picked;
    int t;

    if (!fin->template_dup)
	return 0;

    for (picked = 0; picked < num_picked; picked++) {
	t = fin->template_dup[template];
	while (t != template) {
	    if (t == templates_picked[picked])
		dup = 1;
	    t = fin->template_dup[t];
	    if (!t) {
		fprintf(stderr, "Error: broken template_dup linked list\n");
		break;
	    }
	}
    }

    return dup;
}


/*
 * secondary_primer_match
 *
 * Identifies whether there is a high scoring match with a particular primer
 * elsewhere in a set of files, within the entire database sequence, or
 * within just a section of sequence.
 *
 * This uses fin->extern_seq buffer when check_external is set to true.
 * If check_contig is zero then do not compare against the consensus.
 * If check_contig >0 then compare against region of that contig (between
 * contig_start and contig_end).
 * If check_contig <0 then compare against all contigs in the database.
 *
 * If skip_self is set then the self_match (finding your own match in
 * the consensus) is ignored.
 *
 * Arguments:
 *	fin			'Finish' object
 *	check_contig		Contig to check against (0 => none, -1 => all)
 *	check_start		Start position in check_contig
 *	check_end		End position in check_contig
 *	self_match		Whether to ignore the self-match (true if so)
 *	self_strand		Which strand is self_match relevant to
 *	check_external		Should we check the external file list?
 *	prim			The primer to find matches to.
 *
 * Returns:
 * 	0	for no match (other than self)
 *	>=0	for match
 */
double secondary_primer_match(finish_t *fin,
			      int check_contig,
			      int check_start,
			      int check_end,
			      int self_match,
			      int self_strand,
			      int check_external,
			      char *prim) {
    size_t primer_len;
    char primer[100];
    double sc = 0;

    /* Take copy of primer */
    strncpy(primer, prim, 100);
    primer[99] = 0;
    primer_len = strlen(primer);

    if (check_contig < 0 && fin->all_cons_h) {
	/* All contigs */
	if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
	    printf("Check allcons self=%d strand %d\n",
		   self_match, self_strand);
	sc = hash_compare_primer(fin->all_cons_h, primer, primer_len,
				 fin->opts.pwalk_max_match,
				 self_match, self_strand);
    } else if (check_contig > 0) {
	/* Specific contig */
	if (check_contig != fin->contig) {
	    printf("Trying to check against the wrong 'specific contig'\n");
	    return 0;
	}

	if (check_start < 0)
	    check_start = 0;
	if (check_end >= io_clength(fin->io, check_contig)) {
	    check_end = io_clength(fin->io, check_contig)-1;
	}
	if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
	    printf("Check cons %d..%d self=%d strand %d\n",
		   check_start, check_end, self_match, self_strand);
	sc = compare_primer(fin->cons + check_start,
			    check_end - check_start + 1,
			    primer, primer_len,
			    fin->opts.pwalk_max_match,
			    self_match, self_strand);
    }
    /* else no contig check */

    /* Compare against external sequences */
    if (check_external && fin->external_seq) {
	double vsc;
	if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
	    printf("Check extern self=%d strand %d\n", 0, 0);
	vsc = hash_compare_primer(fin->external_seq_h, primer, primer_len,
				  fin->opts.pwalk_max_match, 0, 0);
	if (vsc > sc)
	    sc = vsc;
    }

#if !defined(USE_OSP)
    /* Use primer3 matching too */
    /* FIXME: to do. Cannot do this yet as it doesn't check its own sequence
     * and we cannot get it to skip 1 hit (self match).
     */
#endif

    return sc;
}
