/*
 * tg_index_common.c - common functions for tg_index
 *
 * Andrew Whitwham, October 2009
 * Wellcome Trust Sanger Institute
 *
 */

#include <staden_config.h>

#include <string.h>
#include <stdio.h>
#include <fcntl.h>

#include "tg_gio.h"
#include "tg_index_common.h"
#include "zfio.h"

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
    free(tmp);
}

/*
 * Stores a name and record in a temporary file suitable for sorting and
 * adding to the name index at a later stage.
 */
void bttmp_file_store(bttmp_t *tmp,  size_t name_len, char *name, tg_rec rec) {
    fprintf(tmp->fp, "%.*s %"PRIrec"\n", (int)name_len, name, rec);
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
char *bttmp_file_get(bttmp_t *tmp, tg_rec *rec) {
    static char line[8192];
    int64_t recno;

    if (fscanf(tmp->fp, "%s %"PRId64"\n", line, &recno) == 2) {
	*rec = recno;
	return line;
    }

    *rec = feof(tmp->fp) ? 0 : 1;
	
    return NULL;
}

#if 0
/* debugging functions */
static void fprint_pair(FILE *tp, char *name, pair_loc_t *p) {
    fprintf(tp, "%s %"PRIrec" %"PRIrec" %d %"PRIrec" %d\n",
	    name, p->rec, p->bin, p->idx, p->crec, p->pos);
} 

static void print_pair(pair_loc_t *p) {
    fprintf(stderr, "rec:%"PRIrec" bin:%"PRIrec" idx:%d crec:%"PRIrec
	    " pos:%d\n",
	    p->rec, p->bin, p->idx, p->crec, p->pos);
} 

static void print_range(range_t *r) {
    fprintf(stderr, "start:%d end:%d rec:%"PRIrec" mqual:%d pair_rec:%"PRIrec
	    " flags:%d\n",
	    r->start, r->end, r->rec, r->mqual, r->pair_rec, r->flags);
}
#endif

/* --------------------------------------------------------------------------
 * Read-pair and sequence storing functions. Common to all file format
 * parsers.
 */

/* save sequence, returns recno */
tg_rec save_sequence(GapIO *io, seq_t *seq, bin_index_t *bin, range_t *r_out) {

    seq->bin = bin->rec;
    seq->bin_index = r_out - ArrayBase(range_t, bin->rng);
    
    return sequence_new_from(io, seq);
}


typedef struct {
    char *tname;
    char *data;
} pair_data;


static int cmp_pair(const void *p1, const void *p2) {
    pair_data *first = (pair_data *)p1; 
    pair_data *second = (pair_data *)p2;
    
    return strcmp(first->tname, second->tname);
}


static char *pair_to_pooled_string(string_alloc_t *s_pool, pair_loc_t *p) {
    char holder[255];
    
    sprintf(holder, " %"PRIrec" %"PRIrec" %d %"PRIrec" %d %d %d", 
    	    p->rec, p->bin, p->idx, p->crec, p->pos, p->orient, p->flags);
	    
    return string_dup(s_pool, holder);
}


static pair_queue_t *add_pair_queue(tg_pair_t *pair) {
    pair_queue_t *tmp;
    int cur = pair->que_size;
    
    pair->que_size++;
    
    tmp = (pair_queue_t *)realloc(pair->que, sizeof(pair_queue_t) * pair->que_size);
    
    if (tmp) {
    	pair->que = tmp;
    } else {
    	return NULL;
    }
    
    if (NULL == (pair->que[cur].fp = tmpfile())) {
    	fprintf(stderr, "Cannot open tmp file in add_pair_queue\n");
    	return NULL;
    }
    
    pair->que[cur].pair      = NULL;
    pair->que[cur].name_pool = NULL;
    pair->que[cur].index     = 0;
    pair->que[cur].pair_size = 0;
    
    fprintf(stderr, "New queue added, no. %d\n", pair->que_size);
    
    return &pair->que[cur];
}


static void save_pair_data(tg_pair_t *pair) {
    HacheIter *iter;
    HacheItem *hi;
    pair_data *save_pair;
    string_alloc_t *str_pool;
    int i = 0, i_max;
    pair_queue_t *que;
    
    if (NULL == (save_pair = (pair_data *)malloc(sizeof(pair_data) * pair->write_size))) {
    	fprintf(stderr, "Can't allocate memory in save_pair_data\n");
	return;
    }
    
    str_pool = string_pool_create(1024 * 1024);
    
    if (NULL == str_pool) {
    	fprintf(stderr, "Can't allocate string pool memory in save_pair_data\n");
	return;
    }
    
    iter = HacheTableIterCreate();

    while ((hi = HacheTableIterNext(pair->phache, iter))) {
	save_pair[i].tname = string_dup(str_pool, hi->key);
	save_pair[i].data  = pair_to_pooled_string(str_pool, (pair_loc_t *)hi->data.p);
    	i++;	
    }
    
    i_max = i;
    
    HacheTableIterDestroy(iter);
    
    qsort(save_pair, pair->count, sizeof(pair_data), cmp_pair);
    
    que = add_pair_queue(pair);
    
    if (NULL == que) {
    	fprintf(stderr, "Can't create new pair queue\n");
	return;
    }

    for (i = 0; i < i_max; i++) {
    	fprintf(que->fp, "%s %s\n", save_pair[i].tname, save_pair[i].data);
    }

    if (HacheTableEmpty(pair->phache, 1)) {
    	// TEST TEST TEST need to put in proper fail
    	fprintf(stderr, "save_pair_data failed on HacheTableEmpty\n");
    }
    
    string_pool_destroy(str_pool);
    free(save_pair);
    fflush(que->fp);
}
    	

static void find_pair(GapIO *io, tg_pair_t *pair, tg_rec recno, char *tname,
	       bin_index_t *bin, contig_t *c, seq_t *seq, tg_args *a,
	       range_t *r_out, library_t *lib) {		
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
    pl->flags  = seq->flags;
    pl->mq     = seq->mapping_qual;
    hd.p = pl;
    
    hi = HacheTableAdd(pair->phache, tname, strlen(tname), hd, &new);
        
    if (new) pair->count++;
    
    /* Pair existed already */
    if (!new) {
	pair_loc_t *po = (pair_loc_t *)hi->data.p;
	//bin_index_t *bo;
	
	/* We found one so update r_out now, before flush */
	r_out->flags &= ~GRANGE_FLAG_TYPE_MASK;
	r_out->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	r_out->pair_rec = po->rec;
	if ((po->flags & SEQ_END_MASK) == SEQ_END_REV)
	    r_out->flags |= GRANGE_FLAG_PEND_REV;
	if (po->flags & SEQ_COMPLEMENTED)
	    r_out->flags |= GRANGE_FLAG_COMP2;
	
	if (!a->fast_mode) {
	    /* TEMP - move later*/
	    fprintf(pair->finish, "%"PRIrec" %d %"PRIrec" %d\n",
		    po->bin, po->idx, pl->rec, pl->flags);
	
	    if (po->bin > pair->max_bin) pair->max_bin = po->bin;
	    
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

		//if (pl->mq > 10 && po->mq > 10) /* good only */
		accumulate_library(io, lib, ltype, ABS(isize));
	    }
	}	    

	/* And, making an assumption, remove from hache */
	HacheTableDel(pair->phache, hi, 1);
	pair->count--;
	free(pl);
    }
    
    if (pair->write_size && pair->count >= pair->write_size) {
    	// too many open pairs, store some away till later
    	fprintf(stderr, "Stored pairs %d\n", pair->count);    
	save_pair_data(pair);
	pair->count = 0;
    }
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


/*
 * Turns a comma separated list of data types into a bit-field.
 */
int parse_data_type(char *type) {
    char *cp;
    int data_type = 0;

    do {
	cp = strchr(type, ',');

	if (0 == strncmp(type, "seq", 3))
	    data_type |= DATA_SEQ;
	else if (0 == strncmp(type, "qual", 4))
	    data_type |= DATA_QUAL;
	else if (0 == strncmp(type, "name", 4))
	    data_type |= DATA_NAME;
	else if (0 == strncmp(type, "anno", 4))
	    data_type |= DATA_ANNO;
	else if (0 == strncmp(type, "all",  3))
	    data_type = DATA_ALL;
	else if (0 == strncmp(type, "none", 4))
	    data_type = 0;
	else if (0 == strncmp(type, "blank", 4))
	    data_type = DATA_BLANK;
	else
	    fprintf(stderr, "Ignoring unknown data_type '%.*s'\n",
		    (int)(cp ? cp-type : strlen(type)), type);

	type = cp ? cp+1 : NULL;
    } while (type);

    return data_type;
}

/* ------------------------------------------------------------------------ */
/* Auto file type detection */
int tg_index_file_type (char *fn) {
    char *suffix = strrchr(fn, '.');
    char data[11];
    zfp *fp;

    /* By standard suffix */
    if (suffix) {
	if (0 == strcmp(suffix, ".gz")) {
	    char *suffix2, tmp;
	    tmp = *suffix;
	    *suffix = 0;
	    suffix2 = strrchr(fn, '.');
	    *suffix = tmp;
	    if (suffix2)
		suffix = suffix2;
	}

	suffix++;
	if (0 == strcmp(suffix, "bam") ||
	    0 == strcmp(suffix, "BAM"))
	    return 'b';

	if (0 == strcmp(suffix, "sam") ||
	    0 == strcmp(suffix, "sam.gz") ||
	    0 == strcmp(suffix, "SAM"))
	    return 's';
	
	if (0 == strcmp(suffix, "ace") ||
	    0 == strcmp(suffix, "ace.gz") ||
	    0 == strcmp(suffix, "ACE"))
	    return 'A';
	
	if (0 == strcmp(suffix, "baf") ||
	    0 == strcmp(suffix, "baf.gz") ||
	    0 == strcmp(suffix, "BAF"))
	    return 'B';
	
	if (0 == strcmp(suffix, "map") ||
	    0 == strcmp(suffix, "MAP") ||
	    0 == strcmp(suffix, "maq"))
	    return 'm';

	if (0 == strcmp(suffix, "caf") ||
	    0 == strcmp(suffix, "CAF"))
	    return 'C';

	if (0 == strcmp(suffix, "fna") ||
	    0 == strcmp(suffix, "FNA") || 
	    0 == strcmp(suffix, "fasta") || 
	    0 == strcmp(suffix, "FASTA")) {
	    return 'F';
	}

	if (0 == strcmp(suffix, "fnq") ||
	    0 == strcmp(suffix, "FNQ") || 
	    0 == strcmp(suffix, "fastq") || 
	    0 == strcmp(suffix, "FASTQ")) {
	    return 'Q';
	}
    }

    /* By contents */
    if (NULL == (fp = zfopen(fn, "rb"))) {
	perror(fn);
	return '?';
    }

    if (NULL == zfgets(data, 10, fp)) {
	zfclose(fp);
	return '?';
    }
    zfclose(fp);

    if (0 == strncmp(data, "BAM\001", 4))
	return 'b'; /* bam */

    if (0 == strncmp(data, "AS ", 3))
	return 'A'; /* ace */

    /* Gets trickier to detect from here on */
    if (0 == strncmp(data, "CO=", 3))
	return 'B'; /* baf */

    if (0 == strncmp(data, "@HD\t", 3) ||
	0 == strncmp(data, "@SQ\t", 3) ||
	0 == strncmp(data, "@RG\t", 3) ||
	0 == strncmp(data, "@PG\t", 3))
	return 's'; /* sam */


    /* These are now pretty tenuous */
    if (*data == '>') {
	return 'F'; /* fasta */
    }

    if (*data == '@') {
	return 'Q'; /* fastq */
    }

    /*
     * And if still not found, well it maybe maq, but is just as likely
     * a differently formatted sam or baf. Give up at this point.
     */

    return '?';
}


/* ------------------------------------------------------------------------ */

/*
 * Relaces \n with newline and \\ with \.
 * Modifies the line in-situ as it can never grow.
 */
void unescape_line(char *txt) {
    char *cp;
    for (cp = txt; *txt; txt++) {
	if (*txt != '\\') {
	    *cp++ = *txt;
	} else {
	    if (*++txt == 'n')
		*cp++ = '\n';
	    else
		*cp++ = *txt;
	    if (!*txt)
		break;
	}
    }
    *cp++ = 0;
}


tg_pair_t *create_pair(int queue) {
    tg_pair_t *pair;
    
    if (NULL == (pair = (tg_pair_t *)malloc(sizeof(tg_pair_t)))) {
    	return NULL;
    }
    
    pair->que          = NULL;
    pair->que_size     = 0;
    pair->working_size = 1000;  
    pair->write_size   = queue; 
    pair->count        = 0;
    pair->phache       = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
    pair->phache->name = "pair";
    
    if (NULL == (pair->finish = tmpfile())) {
    	free(pair);
	return NULL;
    }
    
    pair->max_bin = 0;
    
    return pair;
}

    
tg_rec save_range_sequence(GapIO *io, seq_t *seq, uint8_t mapping_qual,
			   tg_pair_t *pair, int is_pair, char *tname,
			   contig_t *c, tg_args *a, int flags, library_t *lib) {
    range_t r, *r_out;
    tg_rec recno;
    bin_index_t *bin;
    static tg_rec fake_recno = 1;
    int comp;

    /* Update sequence library, aka read-group */
    if (lib && !seq->parent_type) {
	seq->parent_type = GT_Library;
	seq->parent_rec = lib->rec;
    }

    /* Create range */
    r.start = seq->pos;
    r.end   = seq->pos + ABS(seq->len)-1;
    r.rec   = 0;
    r.mqual = mapping_qual;
    r.pair_rec = 0;
    r.flags = flags;

    /* Add the range to a bin, and see which bin it was */
    bin = bin_add_range(io, &c, &r, &r_out, &comp, 1);

    /* Save sequence */
    if (a->data_type == DATA_BLANK) {
	recno = fake_recno++;
    } else {
	if (comp) {
	    complement_seq_t(seq);
	    seq->len = -seq->len;
	}

	recno = save_sequence(io, seq, bin, r_out);
    }

    if (is_pair) {
	find_pair(io, pair, recno, tname, bin, c, seq, a, r_out, lib);
    }

    if (a->tmp && (a->data_type & DATA_NAME))
	bttmp_file_store(a->tmp, seq->name_len, seq->name, recno);
    
    if (seq->name)
	free(seq->name);

    /* Link bin back to sequence too before it gets flushed */
    r_out->rec = recno;
    
    return recno;
}
    
    

   
static int sort_pair_file(tg_pair_t *pair) {
    FILE  *old_files[11];
    int div = 1;
    int max_div = 10; /* temp, needs to be variable */
    int i = 0;
    FILE *final;
    tg_rec max_bin = pair->max_bin; 
    
    memset(old_files, 0, sizeof(FILE *) * 11);

    old_files[0] = pair->finish;
    
    
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
    
    pair->finish = final;
    
    return 1;
}


static void complete_pairs(GapIO *io, tg_pair_t *pair) {
    bin_index_t *bo;
    range_t *ro;
    tg_rec current_bin = -1;
    char line[100];
    int rec_count = 0;
    int total_count = 0;
    
    rewind(pair->finish);
    
    while (fgets(line, 100, pair->finish)) {
	int idx, flags;
	tg_rec bin, rec;
	
        sscanf(line, "%"PRIrec" %d %"PRIrec" %d", &bin, &idx, &rec, &flags);
	
	if (bin != current_bin) {
	    if (rec_count > 50000) {
	    	total_count += rec_count;
	    	cache_flush(io);
		
		fprintf(stderr, "%d pairs finished so far\n", total_count);
	    
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
	if ((flags & SEQ_END_MASK) == SEQ_END_REV)
	    ro->flags |= GRANGE_FLAG_PEND_REV;
	if (flags & SEQ_COMPLEMENTED)
	    ro->flags |= GRANGE_FLAG_COMP2;
	
	rec_count++;

    }
    
    total_count += rec_count;
    
    fprintf(stderr, "%d pairs finished in total.\n", total_count);
    
    cache_flush(io);
    
}


/* an alternative to the Gnu specific getline */
long tg_get_line(char **line, long *length, FILE *fp) {
    char *in_line;
    long len;
    long offset = 0;
    long in_chars;
    
    if (line == NULL || fp == NULL || length == NULL) {
    	return -1;
    }
    
    if (*line == NULL || *length <= 0) {
    	if ((*line = malloc(256 * sizeof(char))) == NULL) {
	    return -1;
	}
	
	*length = 256;
    }
    
    in_line = *line;
    len     = *length;
    
    while ((fgets((in_line + offset), (len - offset), fp))) {
    	char *tmp = NULL;
	
    	in_chars = strlen(in_line);
	
    	offset = in_chars;
	
	// see if we have our full line
	if (*(in_line + offset - 1) == '\n') {
	    break;
	}
    	
    	len += len;
	tmp = realloc(in_line, len * sizeof(char));
	    
    	if (tmp) {
	    in_line = tmp;
	} else {
	    fprintf(stderr, "Memory error in get_line\n");
	    return -1;
	}
    }
    
    *line = in_line;
    *length = len;
    
    return offset;
}


static int load_data(pair_queue_t *pq) {
    int i;
    char *line_in = NULL;
    long line_size = 0;
    
    if (pq->name_pool) {
    	string_pool_destroy(pq->name_pool);
    }
    
    pq->name_pool = string_pool_create(1024);
    
    for (i = 0; i < pq->pair_size; i++) {
    	char name[1024];
    	int found;
    	int ret = tg_get_line(&line_in, &line_size, pq->fp);
	
	if (ret <= 0) break;
	
	found = sscanf(line_in, "%s %"PRId64" %"PRId64" %d %"PRId64" %d %d %d",
	    name, &pq->pair[i].rec, &pq->pair[i].bin, &pq->pair[i].idx,
	    &pq->pair[i].crec, &pq->pair[i].pos, &pq->pair[i].orient,
	    &pq->pair[i].flags);
	    
	if (found != 8) {
	    fprintf(stderr, "Error found in line: %s\n", line_in);
	    break;
	}
	
	pq->pair[i].tname = string_dup(pq->name_pool, name);
    }
    
    pq->pair_size = i;
    pq->index = 0;
    
    if (line_in) free(line_in);
    
    return pq->pair_size;
}


static int initialise_queues(tg_pair_t *pair) {
    int i;
    
    for (i = 0; i < pair->que_size; i++) {
    	rewind(pair->que[i].fp);
	
	if (NULL == (pair->que[i].pair = (pair_loc_t *)malloc(sizeof(pair_loc_t) * pair->working_size))) {
	   fprintf(stderr, "Out of memory allocating pairs in initialise_queues\n");
	   return -1;
	}
	
	pair->que[i].name_pool = NULL;
	pair->que[i].index = 0;
	pair->que[i].pair_size = pair->working_size;
	
	if (0 == load_data(&pair->que[i])) {
	    fprintf(stderr, "Initial data load failed on file %d\n", i);
	    return -1;
	}
    }
    
    return 0;
}


static void get_next(pair_queue_t *que) {
    que->index++;
    
    if (que->index == que->pair_size) {
    	load_data(que);
    }
}


static void save_match_pair(tg_pair_t *pair, pair_loc_t *p1, pair_loc_t *p2) {
    
    fprintf(pair->finish, "%"PRIrec" %d %"PRIrec" %d\n",
    	p1->bin, p1->idx, p2->rec, p2->flags);
	
    fprintf(pair->finish, "%"PRIrec" %d %"PRIrec" %d\n",
    	p2->bin, p2->idx, p1->rec, p1->flags);
}
	


static int find_saved_pairs(GapIO *io, tg_pair_t *pair) {
    int num_found = 0;
    int i;
    int still_looking = 1;

    // scan the files until there is nothing left to see
    
    while (still_looking) {
    	char *compare = NULL;
	int file = 0;
	int done = 0, match = 0;
	
	for (i = 0; i < pair->que_size; i++) {
	    pair_queue_t *que = &pair->que[i];
	
	    if (que->pair_size) {
	    	done++;
	    
	    	if (compare == NULL) {
		    compare = que->pair[que->index].tname;
		    file = i;
		} else {
		    int cmp = strcmp(compare, que->pair[que->index].tname);
		    
		    if (cmp > 0) {
		    	compare = que->pair[que->index].tname;
		    	file = i;
		    } else if (cmp == 0) {
		    	match = i;
			break;
		    }
		}
	    }
	}
	
	if (done) {
	    if (match) {
		save_match_pair(pair, &pair->que[match].pair[pair->que[match].index], &pair->que[file].pair[pair->que[file].index]);
		get_next(&pair->que[match]);
		num_found++;
	    }

	    get_next(&pair->que[file]);
	} else {
	    still_looking = 0;
	}
    }
    
    return num_found;
}


void finish_pairs(GapIO *io, tg_pair_t *pair) {
    int ret = 1;
    
    if (pair->que_size) {
    	// deal with left over pairs
    	save_pair_data(pair);

    	fprintf(stderr, "Resolving pair queues. Total is %d\n", pair->que_size);
    
	initialise_queues(pair);
	ret = find_saved_pairs(io, pair);
    
    	fprintf(stderr, "%d pairs found\n", ret);
    } else {
    	fprintf(stderr, "No pair queue found\n");
    }
    
    ret = sort_pair_file(pair);
    
    if (ret) {
    	fprintf(stderr, "Run complete pairs.\n");
    	complete_pairs(io, pair);
    } else {
    	fprintf(stderr, "sort_pair_file failed");
    }
    
    fprintf(stderr, "Pairs complete\n");
}


void delete_pair(tg_pair_t *pair) {
    int i;
    
    for (i = 0; i < pair->que_size; i++) {
    	if (pair->que[i].fp) fclose(pair->que[i].fp);
	
	if (pair->que[i].pair) free(pair->que[i].pair);
	
	if (pair->que[i].name_pool) string_pool_destroy(pair->que[i].name_pool);
    }
    
    if (pair->que) free(pair->que);
    
    if (pair->phache) HacheTableDestroy(pair->phache, 1);
    
    if (pair->finish) fclose(pair->finish);
    
    free(pair);
}
    


