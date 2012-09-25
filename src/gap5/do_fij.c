#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "fij.h"
#include "misc.h" 
#include "dna_utils.h"
#include "align_lib.h"
#include "newgap_structs.h"
#include "io_lib/hash_table.h"
/*#include "hash_lib.h"*/

/* 1/6/98 johnt - need to explicitly import globals from DLLs with Visual C++ */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif

extern DLL_IMPORT int unknown_char;

typedef struct {
    fij_arg *fij_args;
    HashTable *lib_hash;
    Contig_parms *contig_list1;
    Contig_parms *contig_list2;
    int number_of_contigs1;
    int number_of_contigs2;
    int *depad_to_pad1;
    int *depad_to_pad2;
#if 0
    int min_overlap;          /* fij_args->min_overlap */
    int max_percent_mismatch; /* fij_args->max_mis */
    int max_alignment;        /* fij_args->max_display */
#endif
    int seq2_len;
    int one_by_one;
} add_fij_t;

/*
 * Check a sequence overlap to see if enough read pairs support it.
 * If reverse is true, contig2 was complemented.  If this
 * is the case, s2l and s2r should have been adjusted by the caller to
 * the correct positions on the uncomplemented contig2.
 *
 * Returns 0 if enough pairs were found.
 *         1 if not enough pairs were found.
 *        -1 on failure.
 */

static int check_overlap_pairs(add_fij_t *cd, int reverse,
			       int contig1_num, int s1l, int s1r,
			       int contig2_num, int s2l, int s2r) {
    /* For reverse alignment, s2l,s2r have already been converted to
       positions on the original strand */
    GapIO    *io = cd->fij_args->io;
    tg_rec crec1 = cd->contig_list1[contig1_num].contig_number;
    tg_rec crec2 = cd->contig_list2[contig2_num].contig_number;
    contig_t        *c1;
    contig_t        *c2;
    contig_iterator *ci      = NULL;
    rangec_t        *r1, *r2;
    HashTable       *pairs   = NULL;
    pool_alloc_t    *rp_pool = NULL;
    int              good_pairs = 0;
    int              target = cd->fij_args->rp_min_freq;
    int es1l, es1r, es2l, es2r;

    if (target < 1) target = 1;

    /* Expand start/end a bit to find pairs that span out beyond the ends
       of the overlapping region */
    c1 = cache_search(io, GT_Contig, crec1);
    if (NULL == c1) return -1;
    es1l = c1->start < s1l - cd->fij_args->rp_end_size
	? s1l - cd->fij_args->rp_end_size : c1->start;
    es1r = c1->end > s1r + cd->fij_args->rp_end_size
	? s1r + cd->fij_args->rp_end_size : c1->end;
    c2 = cache_search(io, GT_Contig, crec2);
    if (NULL == c2) return -1;
    es2l = c2->start < s2l - cd->fij_args->rp_end_size
	? s2l - cd->fij_args->rp_end_size : c2->start;
    es2r = c2->end > s2r + cd->fij_args->rp_end_size
	? s2r + cd->fij_args->rp_end_size : c2->end;

    ci = contig_iter_new_by_type(io, crec1, 0, CITER_FIRST,
				 es1l, es1r, GRANGE_FLAG_ISSEQ);
    if (NULL == ci) return -1;
    pairs = HashTableCreate(1024, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);
    if (NULL == pairs) goto fail;
    rp_pool = pool_create(sizeof(rangec_t));
    if (NULL == rp_pool) goto fail;
    while (NULL != (r1 = contig_iter_next(io, ci))) {
	HashItem *hi;
	HashData hd;
	if (!r1->pair_rec) continue;  /* Unpaired */
	if (NULL != (hi = HashTableSearch(pairs, (char *)&r1->pair_rec,
					  sizeof(r1->pair_rec)))) {
	    /* internal read pair */
	    r2 = hi->data.p;
	    pool_free(rp_pool, r2);
	    HashTableDel(pairs, hi, 0);
	    continue;
	}
	r2 = pool_alloc(rp_pool);
	if (NULL == r2) goto fail;
	*r2 = *r1;
	hd.p = r2;
	if (!HashTableAdd(pairs, (char *)&r2->rec, sizeof(r2->rec), hd, 0))
	    goto fail;
    }
    contig_iter_del(ci);
    ci = contig_iter_new_by_type(io, crec2, 0, CITER_FIRST,
				 es2l, es2r, GRANGE_FLAG_ISSEQ);
    if (NULL == ci) goto fail;
    while (good_pairs < target && NULL != (r2 = contig_iter_next(io, ci))) {
	HashItem *hi;
	if (NULL != (hi = HashTableSearch(pairs, (char *)&r2->pair_rec,
					  sizeof(r2->pair_rec)))) {
	    seq_t *s;
	    int o1, o2, orient, dist, r1min, r1max, r2min, r2max;
	    tg_rec lib_rec;
	    library_t *lib;
	    int total_count;
	    
	    /* Found a possible link */
	    r1 = hi->data.p;
	    lib_rec = r1->library_rec;
	    /* Get read orientations */
	    o1 = r1->comp;
	    o2 = r2->comp ^ reverse;
	    s = cache_search(io, GT_Seq, r1->rec);
	    if (NULL == s) goto fail;
	    o1 ^= (s->len < 0);
	    if (!lib_rec && s->parent_type == GT_Library) {
		lib_rec = s->parent_rec;
	    }
	    s = cache_search(io, GT_Seq, r2->rec);
	    if (NULL == s) goto fail;
	    o2 ^= (s->len < 0);
	    if (!lib_rec && s->parent_type == GT_Library) {
		lib_rec = s->parent_rec;
	    }
	    if (!lib_rec) continue;
	    
	    /* Filter by library */
	    if (cd->lib_hash
		&& !HashTableSearch(cd->lib_hash,
				    (char *) &lib_rec, sizeof(lib_rec))) {
		continue;
	    }

	    /* Get approximate aligned position */
	    if (reverse) {
		r2->start = s1r - (r2->start - s2l);
		r2->end   = s1r - (r2->end   - s2l);
	    } else {
		r2->start += s1l - s2l;
		r2->end   += s1l - s2l;
	    }

	    /* Work out relative orientation and distance */
	    r1min = MIN(r1->start, r1->end);
	    r1max = MAX(r1->start, r1->end);
	    r2min = MIN(r2->start, r2->end);
	    r2max = MAX(r2->start, r2->end);
	    if (o1 == o2) {
		orient = LIB_T_SAME;
	    } else {
		if ((o1 == 0 && r1min < r2max) || (o1 == 1 && r2min < r1max)) {
		    orient = LIB_T_INWARD;
		} else {
		    orient = LIB_T_OUTWARD;
		}
	    }
	    dist = MAX(r1max, r2max) - MIN(r1min, r2min);

	    /* Get library stats */
	    lib = cache_search(io, GT_Library, lib_rec);
	    if (NULL == lib) goto fail;

	    if (lib->flags == 0) {
		if (0 != update_library_stats(io, lib_rec, 100,
					      NULL, NULL, NULL)) goto fail;
	    }

	    /* Work out if it's a good pair */
	    total_count = lib->counts[0] + lib->counts[1] + lib->counts[2];
	    if (lib->counts[orient] >= .05 * total_count &&
		ABS(dist - lib->insert_size[orient]) < 3 * lib->sd[orient]) {
		good_pairs++;
	    }
	}
    }

    HashTableDestroy(pairs, 0);
    pool_destroy(rp_pool);
    contig_iter_del(ci);

    return good_pairs >= target ? 0 : 1;

 fail:
    if (NULL != rp_pool) pool_destroy(rp_pool);
    if (NULL != pairs)   HashTableDestroy(pairs, 0);
    if (NULL != ci)      contig_iter_del(ci);
    return -1;
}

static void add_fij_overlap(OVERLAP *overlap, int contig1_num,
			    int contig2_num, void *clientdata) {
    double percent_mismatch;
    int seq1_start_f, seq2_start_f;
    static char buf[1024],name1[10],name2[10];
    add_fij_t *cd = (add_fij_t *)clientdata;
    int c1_s;

    c1_s = cd->one_by_one
	? 0 
	: cd->contig_list1[contig1_num].contig_start_offset;

    percent_mismatch = 100.0 - overlap->percent;

    if ((overlap->length < cd->fij_args->min_overlap) || 
	(percent_mismatch > cd->fij_args->max_mis)) return;

    /* note conversion depadded to padded coordinates */
    /* See comments in seq_utils/align_lib.h for how the OVERLAP struct works.
       overlap->left2 = no. pads to left of seq2 = start pos in seq1
       overlap->left1 = no. pads to left of seq1 = start pos in seq2 */
    seq1_start_f = cd->depad_to_pad1[overlap->left2 + c1_s]
	- cd->depad_to_pad1[c1_s]
	- cd->contig_list1[contig1_num].contig_left_extension + 1 ; 
    
    /* add 1 to get base number */
    seq2_start_f = cd->depad_to_pad2[overlap->left1]
	- cd->contig_list2[contig2_num].contig_left_extension + 1 ; 

    /* Check read pairs if screening on */
    if (cd->fij_args->rp_mode >= 0) {
	/* Overhang at right of seq1 = overlap->right1 - overlap->right
	   Overhang at right of seq2 = overlap->right2 - overlap->right
	   Hence end positions in seq1 and seq2 are... */
	int seq1e = overlap->seq1_len - 1 - overlap->right1 + overlap->right;
	int seq2e = overlap->seq2_len - 1 - overlap->right2 + overlap->right;
	int seq1_end_f = cd->depad_to_pad1[seq1e + c1_s]
	    - cd->depad_to_pad1[c1_s]
	    - cd->contig_list1[contig1_num].contig_left_extension + 1;
	int seq2_end_f = cd->depad_to_pad2[seq2e]
	    - cd->contig_list2[contig2_num].contig_left_extension + 1;
	
	if (check_overlap_pairs(cd, 0, contig1_num, seq1_start_f, seq1_end_f,
				contig2_num, seq2_start_f, seq2_end_f)) {
	    return;
	}
    }
    
    sprintf(name1,"%"PRIrec,
	    cd->contig_list1[contig1_num].contig_left_gel);
    sprintf(name2,"%"PRIrec,
	    cd->contig_list2[contig2_num].contig_left_gel);
    sprintf(buf,
	    " Possible join between contig =%"PRIrec
	    " in the + sense and contig =%"PRIrec"\n"
	    " Length %d",
	    cd->contig_list1[contig1_num].contig_number,
	    cd->contig_list2[contig2_num].contig_number,
	    overlap->length);
    
#if 0
    /* Oops.  The initial coordinates in list_alignment are
       padded, but then the alignment is depadded.  Hopefully
       no-one will notice ! */
    
    overlap->seq1_out[overlap->right+1] = '\0';
    overlap->seq2_out[overlap->right+1] = '\0';
    
    if (overlap->length <= cd->fij_args->max_display) {
	ret = list_alignment(&overlap->seq1_out[overlap->left],
			     &overlap->seq2_out[overlap->left],
			     name1,name2,seq1_start_f,seq2_start_f,buf);
    } else {
	vmessage("%s\n", buf);
	vmessage(" Percentage mismatch %5.1f\n\n",
		 percent_mismatch);
    }
#else
    vmessage("%s\n", buf);
    vmessage(" Percentage mismatch %5.1f\n\n",
	     percent_mismatch);
#endif
    
    buffij(cd->contig_list1[contig1_num].contig_left_gel,
	   seq1_start_f,
	   cd->contig_list2[contig2_num].contig_left_gel,
	   seq2_start_f,
	   overlap->length, (int)overlap->score,
	   percent_mismatch);
}

static void add_fij_overlap_r(OVERLAP *overlap, int contig1_num,
			      int contig2_num, void *clientdata) {
    double percent_mismatch;
    int seq1_start_r, seq1_end_r, seq2_start_r, seq2_end_r, seq1e, seq2e;
    static char buf[1024],name1[10],name2[10];
    add_fij_t *cd = (add_fij_t *)clientdata;
    int c1_s;

    c1_s = cd->one_by_one
	? 0 
	: cd->contig_list1[contig1_num].contig_start_offset;

    percent_mismatch = 100.0 - overlap->percent;

    if ((overlap->length < cd->fij_args->min_overlap) ||
	(percent_mismatch > cd->fij_args->max_mis)) return;

    /* note conversion depadded to padded coordinates */
    /* See comments in seq_utils/align_lib.h for how the OVERLAP struct works.
       overlap->left2 = no. pads to left of seq2 = start pos in seq1
       overlap->left1 = no. pads to left of seq1 = start pos in seq2 */
    seq1_start_r = cd->depad_to_pad1[overlap->left2 + c1_s]
	- cd->depad_to_pad1[c1_s]
	- cd->contig_list1[contig1_num].contig_left_extension + 1 ; 
    
    /* add 1 to get base number */
    seq2_start_r = cd->depad_to_pad2[overlap->left1]
	- cd->contig_list2[contig2_num].contig_right_extension + 1 ; 
#if 0
    {
	int p = overlap->left1 + overlap->length - 1;
	int diff = 0;
	if (p > cd->seq2_len-1) {
	    diff = p - (cd->seq2_len-1);
	    p = cd->seq2_len-1;
	}
	seq2_end_r = cd->depad_to_pad2[p] + diff;
    }
#endif

    /* Overhang at right of seq1 = overlap->right1 - overlap->right
       Overhang at right of seq2 = overlap->right2 - overlap->right
       Hence end positions in seq1 and seq2 are... */
    seq2e = overlap->seq2_len - 1 - overlap->right2 + overlap->right;
    seq2_end_r = cd->depad_to_pad2[seq2e]
	- cd->contig_list2[contig2_num].contig_right_extension + 1;

    /* Check read pairs if screening on */
    if (cd->fij_args->rp_mode >= 0) {
	int seq2len_pad = cd->depad_to_pad2[cd->seq2_len - 1];
	seq1e = overlap->seq1_len - 1 - overlap->right1 + overlap->right;
	seq1_end_r = cd->depad_to_pad1[seq1e + c1_s]
	    - cd->depad_to_pad1[c1_s]
	    - cd->contig_list1[contig1_num].contig_left_extension + 1;
	/* Note conversion of seq2_start, seq2_end to complementary strand */
	if (check_overlap_pairs(cd, 1, contig1_num, seq1_start_r, seq1_end_r,
				contig2_num, seq2len_pad - seq2_end_r,
				seq2len_pad - seq2_start_r)) {
	    return;
	}
    }

    sprintf(name1,"%"PRIrec,
	    cd->contig_list1[contig1_num].contig_left_gel);
    sprintf(name2,"%"PRIrec,
	    cd->contig_list2[contig2_num].contig_left_gel);
    sprintf(buf,
	    " Possible join between contig =%"PRIrec
	    " in the - sense and contig =%"PRIrec"\n"
	    " Length %d",
	    cd->contig_list1[contig1_num].contig_number,
	    cd->contig_list2[contig2_num].contig_number,
	    overlap->length);
    
    /* Oops.  The initial coordinates in list_alignment are
       padded, but then the alignment is depadded.  Hopefully
       no-one will notice ! */
    
#if 0
    overlap->seq1_out[overlap->right+1] = '\0';
    overlap->seq2_out[overlap->right+1] = '\0';
    
    if (overlap->length <= cd->fij_args->max_display) {
	ret = list_alignment(&overlap->seq1_out[overlap->left],
			     &overlap->seq2_out[overlap->left],
			     name1,name2,seq1_start_r,seq2_start_r,buf);
    } else {
	vmessage("%s\n", buf);
	vmessage(" Percentage mismatch %5.1f\n\n",
		 percent_mismatch);
    }
#else
    vmessage("%s\n", buf);
    vmessage(" Percentage mismatch %5.1f\n\n",
	     percent_mismatch);
#endif

    buffij(cd->contig_list1[contig1_num].contig_left_gel,
	   seq1_start_r,
	   -cd->contig_list2[contig2_num].contig_left_gel,
	   seq2_end_r,
	   overlap->length, (int)overlap->score,
	   percent_mismatch);
}

int do_it_fij(fij_arg *fij_args, char *seq, int seq_len,
	      int gap_open, int gap_extend,
	      int compare_mode, HashTable *lib_hash, HashTable *links,
	      Contig_parms *contig_list1, int number_of_contigs1,
	      Contig_parms *contig_list2, int number_of_contigs2,
	      int num_shared) {
    int ret, i;
    int max_contig1, max_contig2;
    int seq1_len, seq2_len, contig1_num, contig2_num;
    double comp[5];
    char *depad_seq1    = NULL, *depad_seq2    = NULL;
    int  *depad_to_pad1 = NULL, *depad_to_pad2 = NULL;
    int edge_mode, job, seq1_start, seq2_start;
    int compare_method;
    Hash *h = NULL;
    OVERLAP *overlap = NULL;
    ALIGN_PARAMS *params;
    Contig_parms *contig_list_depadded = NULL;
    add_fij_t add_fij_cd;
    int one_by_one;
    int c1_loop_size;
    int retval = -1;
    int last_contig1 = number_of_contigs1 - 1;
    int strand;
    int first_shared1 = number_of_contigs1 - num_shared;

    if (number_of_contigs1 == 0 || number_of_contigs2 == 0) return -2;

    /* 17 = fast mode; 31 = sensitive */
    compare_method = fij_args->min_match ? 17 : 31;

    /*
     * TODO: For now we have compare_b_bulk(), but no compare_a_bulk().
     * So if we're using the sensitive method then for now we need to revert
     * back to comparing individual contigs rather than a contig against
     * the entire combined consensus.
     */
    one_by_one = (compare_method != 17 || compare_mode == COMPARE_ONE_BY_ONE)
	? 1 : 0;

    if (NULL == (params = create_align_params())) return -1;

    edge_mode = 10;
    seq1_start = 0;
    seq2_start = 0;
    job = RETURN_NEW_PADS | RETURN_EDIT_BUFFERS | RETURN_SEQ;

    if (set_align_params (params, fij_args->band, gap_open, gap_extend,
			  edge_mode, job, seq1_start, seq2_start,0,0,0)) {
	goto out;
    }

    if (NULL == (overlap = create_overlap())) goto out;
    init_overlap (overlap, seq, seq, seq_len, seq_len);

    /* 
     * find longest contig 
     * initialise hashing
     * compare each contig in each orientation
     */

    if (one_by_one) {
	for ( i = 0, max_contig1 = 0; i < number_of_contigs1; i++) {
	    int l = (contig_list1[i].contig_end_offset
		     - contig_list1[i].contig_start_offset);
	    if (l > max_contig1) max_contig1 = l;
	}
	max_contig1++;
    } else {
	max_contig1 = contig_list1[last_contig1].contig_end_offset + 1;
    }

    for (i = 0, max_contig2 = 0; i < number_of_contigs2; i++) {
	int l = (contig_list2[i].contig_end_offset
		 - contig_list2[i].contig_start_offset);
	if (l > max_contig2) max_contig2 = l;
    }
    max_contig2++;

    if (NULL == (depad_seq1   = (char *)xmalloc(sizeof(char) * max_contig1))||
	NULL == (depad_seq2   = (char *)xmalloc(sizeof(char) * max_contig2))||
	NULL == (depad_to_pad1 = (int *)xmalloc(sizeof(int)  * max_contig1))||
	NULL == (depad_to_pad2 = (int *)xmalloc(sizeof(int)  * max_contig2))) {
	goto out;
    }

    if ( init_hash8n(MAX(max_contig1, max_contig2), max_contig2,
		     fij_args->word_len, 8, fij_args->min_match,
		     compare_method, &h )) {
	goto out;
    }

    h->fast_mode = fij_args->fast_mode;
    h->filter_words =  0;
    
    if ( HASH_JOB_EXPD & compare_method ) {

        p_comp(comp,seq,seq_len);

	if(poisson_diagonals(MINMAT, MIN(max_contig1, max_contig2),
			     h->word_length, fij_args->max_prob,
			     h->expected_scores, comp)) {
	    goto out;
	}
    }

    add_fij_cd.fij_args = fij_args;
    add_fij_cd.lib_hash = lib_hash;
    add_fij_cd.contig_list1 = contig_list1;
    add_fij_cd.contig_list2 = contig_list2;
    add_fij_cd.number_of_contigs1 = number_of_contigs1;
    add_fij_cd.number_of_contigs2 = number_of_contigs2;
    add_fij_cd.depad_to_pad1 = depad_to_pad1;
    add_fij_cd.depad_to_pad2 = depad_to_pad2;
#if 0
    add_fij_cd.min_overlap = fij_args->min_overlap;
    add_fij_cd.max_percent_mismatch = fij_args->max_mis;
    add_fij_cd.max_alignment = fij_args->max_display;
#endif
    add_fij_cd.one_by_one = one_by_one;

    c1_loop_size = one_by_one ? number_of_contigs1 : 1;

    for ( contig1_num = 0; contig1_num < c1_loop_size; contig1_num++ ) {
	seq1_len = (one_by_one
		    ? (contig_list1[contig1_num].contig_end_offset
		       - contig_list1[contig1_num].contig_start_offset + 1)
		    : (contig_list1[last_contig1].contig_end_offset + 1));
	
	if (seq1_len < fij_args->min_overlap)
	    continue;

	/* Depad seq1 */
	overlap->seq1 = (one_by_one
			 ? &seq[(contig_list1[contig1_num].contig_start_offset)]
			 : seq);
	copy_seq(depad_seq1, overlap->seq1, seq1_len);
	depad_seq(depad_seq1, &seq1_len, depad_to_pad1);
	overlap->seq1 = h->seq1 = depad_seq1;
	h->seq1_len = overlap->seq1_len = seq1_len;

	/* Convert padded to depadded coords in contig_list */
	if (!one_by_one) {
	    int p, cnum;
	    contig_list_depadded = malloc(number_of_contigs1 *
					  sizeof(Contig_parms));
	    if (NULL == contig_list_depadded) goto out;

	    memcpy(contig_list_depadded,
		   contig_list1,
		   number_of_contigs1 * sizeof(Contig_parms));


	    for (cnum = p = 0; p < seq1_len; p++) {
		if (depad_to_pad1[p] == contig_list1[cnum].contig_start_offset) 
		    contig_list_depadded[cnum].contig_start_offset = p;

		if (depad_to_pad1[p] >= contig_list1[cnum].contig_end_offset)
		    contig_list_depadded[cnum++].contig_end_offset = p;
	    }

	    add_fij_cd.contig_list1 = contig_list_depadded;
	}

	vmessage("Hashing sequences... ");
	//UpdateTextOutput();
	if (hash_seqn(h, 1)) {
	    verror(ERR_WARN, "find internal joins",
		   "hashing 1st sequence (#%"PRIrec")",
		   contig_list1[contig1_num].contig_left_gel); 
	    continue;
	}
	vmessage("done\n\n");
	//UpdateTextOutput();
    
	(void) store_hashn ( h );

	if (fij_args->filter_words > 0) {
	    /* Compute expected mean number of hits per hash key */
	    double m = (seq1_len - 20*number_of_contigs1)
		/ pow(4, h->word_length);
	    h->filter_words = fij_args->filter_words * m + 0.5;
	    if (h->filter_words < 2)
		h->filter_words = 2;
	} else if (fij_args->filter_words < 0) {
	    h->filter_words = -fij_args->filter_words;
	    if (h->filter_words < 2)
		h->filter_words = 2;
	} else {
	    h->filter_words = 0;
	}
	if (h->filter_words)
	    vmessage("Filtering words occuring more than %d times.\n",
		     h->filter_words);

	for (contig2_num = (contig1_num < first_shared1
			    ? 0 : contig1_num - first_shared1 + 1);
	     contig2_num < number_of_contigs2; contig2_num++) {

	    if (contig_list2[contig2_num].contig_number ==
		contig_list1[contig1_num].contig_number) {
		continue;
	    }
	    if (one_by_one && links) {
		HashItem    *hi;
		contig_pair cp;
		cp.c1 = MIN(contig_list1[contig1_num].contig_number,
			    contig_list2[contig2_num].contig_number);
		cp.c2 = MAX(contig_list1[contig1_num].contig_number,
			    contig_list2[contig2_num].contig_number);
		hi = HashTableSearch(links, (char *)&cp, sizeof(cp));
		if (NULL == hi) continue;
	    }
#if 0
	    if (one_by_one) {
		fprintf(stderr,
			"Aligning %d of %d =%"PRIrec
			" against %d of %d =%"PRIrec"\n",
			contig1_num + 1, number_of_contigs1,
			contig_list1[contig1_num].contig_number,
			contig2_num + 1, number_of_contigs2,
			contig_list2[contig2_num].contig_number);
	    } else {
		fprintf(stderr,
			"Aligning all against %d of %d =%"PRIrec"\n",
			contig2_num + 1, number_of_contigs2,
			contig_list2[contig2_num].contig_number);
	    }
#endif
	    //UpdateTextOutput();

	    for (strand = 0; strand < 2; strand++) {
		seq2_len = contig_list2[contig2_num].contig_end_offset 
		    - contig_list2[contig2_num].contig_start_offset + 1;

		if (seq2_len < fij_args->min_overlap)
		    continue;

		/* Depad seq2 */
		overlap->seq2 =
		    &seq[(contig_list2[contig2_num].contig_start_offset)];
		copy_seq(depad_seq2, overlap->seq2, seq2_len);
		if (1 == strand) { /* Reverse strand */
		    complement_seq(depad_seq2, seq2_len);
		}
		depad_seq(depad_seq2, &seq2_len, depad_to_pad2);
		overlap->seq2 = h->seq2 = depad_seq2;
		h->seq2_len = overlap->seq2_len = add_fij_cd.seq2_len = seq2_len;

		/* fflush(stdout); */
		if (hash_seqn(h, 2)) {
		    verror(ERR_WARN, "find internal joins",
			   "hashing 2nd sequence (#%"PRIrec")",
			   contig_list2[contig2_num].contig_left_gel); 
		    continue;
		}

		ret = 0;
		if ( compare_method != 17 ) {
		    /*
		     * Accept that some matches may just be too slow
		     * or memory hungry for sensitive mode, so we fall
		     * back to fast methods when this fails.
		     */
		    ret = compare_a ( h, params, overlap );
		    if (ret < 0) {
			verror(ERR_WARN, "find internal joins",
			       "alignment too large for sensitive mode;"
			       " falling back to quick alignment"); 
		    }
		    if (ret > 0) {
			if (0 == strand) { /* forward */
			    add_fij_overlap(overlap, contig1_num, contig2_num,
					    &add_fij_cd);
			} else { /* reverse */
			    add_fij_overlap_r(overlap, contig1_num, contig2_num,
					      &add_fij_cd);
			}
		    }

		}
		if ( ret < 0 || compare_method == 17 ) {
		    if (one_by_one) {
			ret = compare_b ( h, params, overlap );

			if ( ret < 0 ) {
			    verror(ERR_WARN, "find internal joins", "hashing" );
			    continue;
			}
		    
			if ( ret ) {
			    if (0 == strand) {
				add_fij_overlap(overlap, contig1_num,
						contig2_num, &add_fij_cd);
			    } else {
				add_fij_overlap_r(overlap, contig1_num,
						  contig2_num, &add_fij_cd);
			    }
			}
		    } else {
			int lim = (contig2_num < num_shared
				   ? contig_list_depadded[number_of_contigs1
							  - num_shared
							  + contig2_num - 1].contig_end_offset
				   : (contig_list_depadded[last_contig1].contig_end_offset + 1));
			ret = compare_b_bulk ( h, params, overlap,
					       contig2_num,
					       contig_list2[contig2_num].contig_number,
					       contig_list_depadded,
					       number_of_contigs1,
					       lim, links,
					       (0 == strand
						? add_fij_overlap
						: add_fij_overlap_r),
					       &add_fij_cd);
		    }
		}
		free_overlap(overlap);
	    }
	}
    }

    retval = 0;

 out:
    if (NULL != depad_seq1)    xfree(depad_seq1);
    if (NULL != depad_seq2)    xfree(depad_seq2);
    if (NULL != depad_to_pad1) xfree(depad_to_pad1);
    if (NULL != depad_to_pad2) xfree(depad_to_pad2);
    if (NULL != h) free_hash8n(h);
    destroy_alignment_params (params);
    if (NULL != overlap) destroy_overlap (overlap);
    if (NULL != contig_list_depadded) xfree(contig_list_depadded);

    return retval;
}
