/*
 * tg_index_common.c - common functions for tg_index
 *
 * Andrew Whitwham, October 2009
 * Wellcome Trust Sanger Institute
 *
 */

#include <string.h>
#include <stdio.h>
#include <fcntl.h>

#include "tg_gio.h"
#include "tg_index_common.h"


typedef struct {
    int rec;
    int bin;
    int idx;
    int crec;
    int pos;
    int orient;
} pair_loc_t;

/* --------------------------------------------------------------------------
 * Temporary file handling for storing name + record.
 * This is used in the name B+Tree generation code as it's more efficient
 * to delay generation of the B+Tree until after adding all the sequence
 * records.
 */

bttmp_t *bttmp_file_open(void) {
    bttmp_t *tmp = malloc(sizeof(*tmp));
    int fd;

    if (!tmp)
	return NULL;

    /*
     * This emits a warning from gcc:
     *     the use of `tmpnam' is dangerous, better use `mkstemp'
     *
     * The problem is, mkstemp isn't standard while tmpnam is. Instead
     * we use tmpnam in a safe manner via open with O_CREAT|O_EXCL and then
     * convert this to a FILE pointer. This is basically what tmpfile does, 
     * but in our case we need to know the file name too so we can sort
     * it later on.
     */
    if (NULL == tmpnam(tmp->name)) {
	free(tmp);
	return NULL;
    }
    
    if (-1 == (fd = open(tmp->name, O_RDWR|O_CREAT|O_EXCL, 0666))) {
	free(tmp);
	return NULL;
    }

    tmp->fp = fdopen(fd, "wb+");

    return tmp;
}

void bttmp_file_close(bttmp_t *tmp) {
    if (tmp->fp) {
	fclose(tmp->fp);
	tmp->fp = NULL;
    }

    unlink(tmp->name);
}

/*
 * Stores a name and record in a temporary file suitable for sorting and
 * adding to the name index at a later stage.
 */
void bttmp_file_store(bttmp_t *tmp,  size_t name_len, char *name, int rec) {
    fprintf(tmp->fp, "%.*s %d\n", name_len, name, rec);
}

/* Sort the temporary file, and rewind to start */
void bttmp_file_sort(bttmp_t *tmp) {
    char new_tmp[L_tmpnam];
    char buf[100+2*L_tmpnam];

    tmpnam(new_tmp);
    sprintf(buf, "sort < %s > %s", tmp->name, new_tmp);
    fclose(tmp->fp);

    /* Use unix sort for now */
    printf("buf=%s\n", buf);
    system(buf);
    printf("done\n");

    unlink(tmp->name);
    strcpy(tmp->name, new_tmp);
    tmp->fp = fopen(tmp->name, "rb+");
}

/*
 * Repeatedly fetch lines from the temp file.
 * NB: non-reentrant. Value is valid only until the next call to this
 * function.
 *
 * Return name on success and fills out rec.
 *       NULL on EOF (*rec==0) or failure (*rec==1)
 */
char *bttmp_file_get(bttmp_t *tmp, int *rec) {
    static char line[8192];
    static int recno;

    if (fscanf(tmp->fp, "%s %d\n", line, &recno) == 2) {
	*rec = recno;
	return line;
    }

    *rec = feof(tmp->fp) ? 0 : 1;
	
    return NULL;
}

/* debugging functions */
static void print_pair(pair_loc_t *p) {
    fprintf(stderr, "rec:%d bin:%d idx:%d crec:%d pos:%d\n", p->rec, p->bin, p->idx, p->crec, p->pos);
} 

static void print_range(range_t *r) {
    fprintf(stderr, "start:%d end:%d rec:%d mqual:%d pair_rec:%d flags:%d\n", r->start, r->end, r->rec, r->mqual, r->pair_rec, r->flags);
}



/* temp file handling */
static FILE *fp = NULL;
static int max_bin = 0;

int open_tmp_file(void) {
    
    fp = tmpfile();
    
    return fp ? 1 : 0;
}

void close_tmp_file(void) {
    fclose(fp);
}

/* --------------------------------------------------------------------------
 * Read-pair and sequence storing functions. Common to all file format
 * parsers.
 */

/* save sequence, returns recno */
int save_sequence(GapIO *io, seq_t *seq, bin_index_t *bin, range_t *r_out) {

    seq->bin = bin->rec;
    seq->bin_index = r_out - ArrayBase(range_t, bin->rng);
    
    return sequence_new_from(io, seq);
}


void find_pair(GapIO *io, HacheTable *pair, int recno, char *tname,
	       bin_index_t *bin, contig_t *c, seq_t *seq, tg_args *a,
	       range_t *r_out, library_t *lib){		
    int new = 0;
    HacheData hd;
    pair_loc_t *pl  = NULL;
    HacheItem *hi   = NULL;
 
    /* Add data for this end */
    pl = (pair_loc_t *)malloc(sizeof(*pl));
    pl->rec    = recno;
    pl->bin    = bin->rec;
    pl->crec   = c->rec;
    pl->pos    = seq->len >= 0 ? seq->pos : seq->pos - seq->len - 1;
    pl->idx    = seq->bin_index;
    pl->orient = seq->len < 0;
    hd.p = pl;

    hi = HacheTableAdd(pair, tname, strlen(tname), hd, &new);
    
    /* Pair existed already */
    if (!new) {
	pair_loc_t *po = (pair_loc_t *)hi->data.p;
	
	bin_index_t *bo;
	
	/* We found one so update r_out now, before flush */
	r_out->flags &= ~GRANGE_FLAG_TYPE_MASK;
	r_out->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	r_out->pair_rec = po->rec;
	
	if (!a->fast_mode) {
	    /* TEMP - move later*/
	    fprintf(fp, "%d %d %d\n", po->bin, po->idx, pl->rec);
	
	    if (po->bin > max_bin) max_bin = po->bin;
	    
	    /* fprintf(stderr, "Get other side\n"); */
	    /* Link other end to 'us' too */
	    /*
	    bo = (bin_index_t *)cache_search(io, GT_Bin, po->bin);
	    bo = cache_rw(io, bo);
	    bo->flags |= BIN_RANGE_UPDATED;
	    ro = arrp(range_t, bo->rng, po->idx);
	    ro->flags &= ~GRANGE_FLAG_TYPE_MASK;
	    ro->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	    ro->pair_rec = pl->rec;
	    */
	}
	
	if (lib) {
	    /* Increment insert size in library */
	    if (po->crec == pl->crec) {
		int isize = pl->pos - po->pos;
		int ltype;

		/*
		 * We know that 'seq' is the right-most sequence so
		 * this position minus previous position gives us the
		 * insert size and comparing this vs previous orient
		 * we can work out the type of the library. Note although
		 * this read is further right than the previous one, when
		 * overlapping it's possible for the 5' of this minus the
		 * 5' of previous to appear as a negative size.
		 *
		 * Types of pair orientations:
		 *
		 *                           isize   ltype
		 * |------->     <-------|   +ve     IN
		 *
		 * <-------|                 -ve     IN
		 *    |------->
		 *
		 * <-------|     |------->   +ve     OUT
		 *
		 * |------->                 +ve     OUT
		 *    <-------|
		 *
		 * <-------|     <-------|   +ve     SAME
		 *
		 * |------->     |------->   +ve     SAME
		 *
		 * <-------|                 +ve     SAME
		 *    <-------|
		 *
		 * |------->                 +ve     SAME
		 *    |------->
		 */
		if (pl->orient == po->orient) {
		    ltype = LIB_T_SAME;
		} else {
		    if ((isize >= 0 && pl->orient == 1 /* <----| */) ||
			(isize <  0 && pl->orient == 0 /* |----> */))
			ltype = LIB_T_INWARD;
		    else
			ltype = LIB_T_OUTWARD;
		}
		
		lib = cache_rw(io, lib);
		
		accumulate_library(io, lib, ltype, ABS(isize));
	    }
	}	    

	/* And, making an assumption, remove from hache */
	HacheTableDel(pair, hi, 1);
	free(pl);
    }
}


int save_range_sequence(GapIO *io, seq_t *seq, uint8_t mapping_qual,
			HacheTable *pair, int is_pair, char *tname,
			contig_t *c, tg_args *a, int flags, library_t *lib) {
    range_t r, *r_out;
    int recno;
    bin_index_t *bin;
    HacheItem *hi;
    static int fake_recno = 1;

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
    if (a->data_type == DATA_BLANK)
	recno = fake_recno++;
    else
	recno = save_sequence(io, seq, bin, r_out);

    if (is_pair) {
	find_pair(io, pair, recno, tname, bin, c, seq, a, r_out, lib);
    }

    if (a->tmp && (a->data_type & DATA_NAME))
	bttmp_file_store(a->tmp, seq->name_len, seq->name, recno);

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

static void sort_file (FILE *old_files[], int div) {
    FILE *new_files[10];
    char line[100];
    int i;

    memset(new_files, 0, sizeof(FILE *) * 10);
    
    for (i = 0; i < 10; i++) {
	new_files[i] = tmpfile(); /* should auto-delete */
    }
    
    for (i = 0; i < 10; i++) {
    	if (old_files[i]) {
    	    rewind(old_files[i]);

	    while (fgets(line, 100, old_files[i])) {
    		int bin;
		int mod;

		sscanf(line, "%d", &bin);

		bin = bin / div;

		if (bin) {
	    	    mod = bin % 10;
		} else {
	    	    mod = 0;
		}

		fputs(line, new_files[mod]);
	    }
	    
	    fclose(old_files[i]);
	}
	
	old_files[i] = new_files[i];
    }
}


int sort_pair_file (void) {
    FILE  *old_files[11];
    int div = 1;
    int max_div = 10; /* temp, needs to be variable */
    int i = 0;
    FILE *final;
    
    memset(old_files, 0, sizeof(FILE *) * 11);

    old_files[0] = fp;
    
    
    while ((max_bin % max_div) != max_bin) {
    	max_div *= 10;
    }
    
    while (div < max_div) {
    	sort_file(old_files, div);
	div = div * 10;
    }
    
    /* gather files together here */
    
    final = tmpfile();
    
    while (old_files[i]) {
    	char line[100];
     	rewind(old_files[i]);
	
	while (fgets(line, 100, old_files[i])) {
  	    fputs(line, final);
	}
 	
     	fclose(old_files[i++]);
    }
    
    fp = final;
    
    return 1;
}


void complete_pairs(GapIO *io) {
    bin_index_t *bo;
    range_t *ro;
    int current_bin = -1;
    char line[100];
    int rec_count = 0;
    int flush = 1;
    
    rewind(fp);
    
    while (fgets(line, 100, fp)) {
    	int bin, idx, rec;
	
        sscanf(line, "%d %d %d", &bin, &idx, &rec);
	
	if (bin != current_bin) {

	    if (rec_count > 50000) {
	    	cache_flush(io);
		rec_count = 0;
	    }

	    bo = (bin_index_t *)cache_search(io, GT_Bin, bin);
	    bo = cache_rw(io, bo);
	    bo->flags |= BIN_RANGE_UPDATED;
    	    current_bin = bin;
	}

	ro = arrp(range_t, bo->rng, idx);
	ro->flags &= ~GRANGE_FLAG_TYPE_MASK;
	ro->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	ro->pair_rec = rec;
	
	rec_count++;

    }
    
    cache_flush(io);
    
}
