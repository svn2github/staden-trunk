#include <string.h>
#include <stdio.h>

#include "tg_gio.h"
#include "tg_index_common.h"


typedef struct {
    int rec;
    int bin;
    int idx;
    int crec;
    int pos;
} pair_loc_t;

/* debugging functions */
static void print_pair(pair_loc_t *p) {
    fprintf(stderr, "rec:%d bin:%d idx:%d crec:%d pos:%d\n", p->rec, p->bin, p->idx, p->crec, p->pos);
} 

static void print_range(range_t *r) {
    fprintf(stderr, "start:%d end:%d rec:%d mqual:%d pair_rec:%d flags:%d\n", r->start, r->end, r->rec, r->mqual, r->pair_rec, r->flags);
}


/* save sequence, returns recno */
int save_sequence(GapIO *io, seq_t *seq, bin_index_t *bin, range_t *r_out) {

    seq->bin = bin->rec;
    seq->bin_index = r_out - ArrayBase(range_t, bin->rng);
    
    return sequence_new_from(io, seq);
}


void find_pair(GapIO *io, HacheTable *pair, int recno, char *tname, bin_index_t *bin, contig_t *c,
                  seq_t *seq, tg_args *a, range_t *r_out, library_t *lib){		
    int new = 0;
    HacheData hd;
    pair_loc_t *pl = NULL;
    HacheItem *hi  = NULL;

    /* Add data for this end */
    pl = (pair_loc_t *)malloc(sizeof(*pl));
    pl->rec  = recno;
    pl->bin  = bin->rec;
    pl->crec = c->rec;
    pl->pos  = seq->len >= 0 ? seq->pos : seq->pos - seq->len - 1;
    pl->idx  = seq->bin_index;
    hd.p = pl;

    hi = HacheTableAdd(pair, tname, strlen(tname), hd, &new);
    
    /* Pair existed already */
    if (!new) {
	pair_loc_t *po = (pair_loc_t *)hi->data.p;
	
	bin_index_t *bo;
	range_t *ro;
	
	/* We found one so update r_out now, before flush */
	r_out->flags &= ~GRANGE_FLAG_TYPE_MASK;
	r_out->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	r_out->pair_rec = po->rec;

	if (!a->fast_mode) {
	    /* fprintf(stderr, "Get other side\n"); */
	    /* Link other end to 'us' too */
	    bo = (bin_index_t *)cache_search(io, GT_Bin, po->bin);
	    bo = cache_rw(io, bo);
	    bo->flags |= BIN_RANGE_UPDATED;
	    ro = arrp(range_t, bo->rng, po->idx);
	    ro->flags &= ~GRANGE_FLAG_TYPE_MASK;
	    ro->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	    ro->pair_rec = pl->rec;
	}
	
	if (lib) {
	    /* Increment insert size in library */
	    if (po->crec == pl->crec) {
		int isize = pl->pos - po->pos;
		/*
		 * We can get +ve isize via:
		 * |------->     <-------|
		 *
		 * and -ve isize via:
		 * <-------|
		 *    |------->
		 *
		 * We know that 's' is the right-most sequence so
		 * when this_pos as the input is sorted.
		 * Therefore we can tell which case it is by the orientation
		 * of this sequence, and negate isize for the 2nd case.
		 */
		if (!(r_out->flags & GRANGE_FLAG_COMP1))
		    isize = -isize;

		lib = cache_rw(io, lib);
		
		if (a->fast_mode || (r_out->flags & GRANGE_FLAG_COMP1) != 
		    (ro->flags & GRANGE_FLAG_COMP1)) {
		    accumulate_library(io, lib,
				       isize >= 0
				       ? LIB_T_INWARD
				       : LIB_T_OUTWARD,
				       ABS(isize));
		} else {
		    accumulate_library(io, lib, LIB_T_SAME, ABS(isize));
		}
	    }
	}	    

	/* And, making an assumption, remove from hache */
	HacheTableDel(pair, hi, 1);
	free(pl);
    }
}


int save_range_sequence(GapIO *io, seq_t *seq, uint8_t mapping_qual, HacheTable *pair, int is_pair, char *tname, contig_t *c,
                  	 tg_args *a, int flags, library_t *lib) {
    range_t r, *r_out;
    int recno;
    bin_index_t *bin;

    /* Create range */
    r.start = seq->pos;
    r.end   = seq->pos + ABS(seq->len) - 1;
    r.rec   = 0;
    r.mqual = mapping_qual;
    r.pair_rec = 0;
    r.flags = flags;

    /* Add the range to a bin, and see which bin it was */
    bin = bin_add_range(io, &c, &r, &r_out, NULL);

    /* Save sequence */
    recno = save_sequence(io, seq, bin, r_out);


    if (is_pair) {
	find_pair(io, pair, recno, tname, bin, c, seq, a, r_out, lib);
    }

    if (!a->no_tree)
	sequence_index_update(io, seq->name, seq->name_len, recno);
	
    free(seq->data);

    /* Link bin back to sequence too before it gets flushed */
    r_out->rec = recno;
    
    return recno;
}

void create_new_contig(GapIO *io, contig_t **c, char *cname, int merge) {

    if (*c) {
	cache_decr(io, *c);
    }	    

    if (!merge || (NULL == (*c = find_contig_by_name(io, cname)))) {
    	*c = contig_new(io, cname);
	contig_index_update(io, cname, strlen(cname), (*c)->rec);
    }
    
    cache_incr(io, *c);
}    






