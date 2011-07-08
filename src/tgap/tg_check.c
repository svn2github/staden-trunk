/*
 * Checks the validity of data in the database.
 *
 * These can be split into actual corruptions and credibility warnings.
 * Eg a sequence length not matching the start/end range_t is a fatal error,
 * but a sequence being clipped to zero bases is valid, just rather
 * unbelievable.
 */
#include <staden_config.h>
#include <string.h>
#include <stdlib.h>

#include "tg_gio.h"
#include "tg_check.h"
#include "hache_table.h"

/* Enable debugging, which may help track down the location of some errors */
//#define DEBUG_CHECK

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

typedef struct {
    /* Contig start/end */
    int cstart;
    int cend;
    
    /* Item counts */
    int nseq;
    int nanno;
    int nref;
} bin_stats;

/*
 * Checks a sequence struct is consistent, both internally and with the
 * range that points to it.
 */
static int check_seq(GapIO *io, int fix, bin_index_t *bin, range_t *r,
		     HacheTable *lib_hash, int is_cons) {
    int err = 0;
    int i, len;
    seq_t *s = cache_search(io, GT_Seq, r->rec);

    static char isbase[256], init_done = 0;

    if (!init_done) {
	memset(isbase, 0, 256);
	isbase['A'] = isbase['a'] = 1;
	isbase['C'] = isbase['c'] = 1;
	isbase['G'] = isbase['g'] = 1;
	isbase['T'] = isbase['t'] = 1;
	isbase['N'] = isbase['n'] = 1;
	isbase['*'] = isbase['-'] = 1;
	isbase['0'] = isbase['1'] = 1;
	isbase['2'] = isbase['3'] = 1;
	init_done = 1;
    }

    if (!s) {
	vmessage("Seq %"PRIrec": failed to load\n", r->rec);
	if (fix && !io->base) {
	    /* Returns bin if io->base is NULL */
	    cache_rw(io, bin);
	    r->rec      = 0;
	    r->flags    = GRANGE_FLAG_UNUSED;
	    bin->flags |= BIN_RANGE_UPDATED;
	}
	return 1;
    }

    if (!bin->rng ||
	s->bin_index != r - ArrayBase(range_t, bin->rng)) {
	vmessage("Seq %"PRIrec": bin_index does not match range index\n",
		 s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->bin_index = r - ArrayBase(range_t, bin->rng);
	}
    }

    if (s->bin != bin->rec) {
	vmessage("Seq %"PRIrec": bin does not match observed bin\n", s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->bin = bin->rec;
	}
    }
    if (r->end - r->start + 1 != ABS(s->len)) {
	vmessage("Seq %"PRIrec": length does not match bin-range\n", s->rec);
	err++;
	if (fix && !io->base) {
	    cache_rw(io, bin);
	    r->end = ABS(s->len)-1 + r->start;
	    bin->flags |= BIN_RANGE_UPDATED;
	}
    }

    if (s->left < 1 || s->right > ABS(s->len)) {
	vmessage("Seq %"PRIrec": left/right clips outside of sequence "
		 "bounds.\n", s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    if (s->left < 1)
		s->left = 1;
	    if (s->right > ABS(s->len))
		s->right = s->len;
	}
    }

    if (s->right < s->left) {
	vmessage("Seq %"PRIrec": right clip starts before left clip\n",
		 s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->left = s->right;
	}
    }

    if (s->mapping_qual != r->mqual) {
	vmessage("Seq %"PRIrec": mapping_qual disagrees with range\n",
		 s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->mapping_qual = r->mqual;
	}
    }

    /* TODO: Check r->pair_rec? */

    /* TODO: Check r->flags match s->flags */

    if (s->seq_tech > STECH_454) {
	vmessage("Seq %"PRIrec": Unknown seq_tech value\n", s->rec);
	err++;
    }

    if (s->format > 3) {
	vmessage("Seq %"PRIrec": Unknown format value\n", s->rec);
	err++;
    }

    if (s->format != SEQ_FORMAT_CNF1) {
	vmessage("Seq %"PRIrec": Warning - outdated seq format\n", s->rec);
    }

    /* Parent rec / type */
    if (s->parent_rec && !cache_exists(io, s->parent_type, s->parent_rec)) {
	vmessage("Seq %"PRIrec": parent_rec/type do not match\n", s->rec);
	err++;
    }

    if (lib_hash && s->parent_type == GT_Library) {
	HacheItem *hi = HacheTableQuery(lib_hash, (char *)&s->parent_rec,
					sizeof(s->parent_rec));
	if (!hi) {
	    HacheData hd;
	    vmessage("Seq %"PRIrec": parent_rec %"PRIrec" of type "
		     "GT_Library is not in db->library array\n",
		     s->rec, s->parent_rec);
	    hd.i = 0;
	    HacheTableAdd(lib_hash, (char *)&s->parent_rec,
			  sizeof(s->parent_rec), hd, NULL);
	    err++;
	}
    }
    

    /* Check valid sequence characters */
    len = ABS(s->len);
    for (i = 0; i < len; i++) {
	if (!isbase[(unsigned char)s->seq[i]] &&
	    !(s->seq[i] == ' ' && is_cons)) {
	    vmessage("Seq %"PRIrec": seq has unexpected char '%c' (%d)\n",
		     s->rec, isprint(s->seq[i]) ? s->seq[i] : '?', s->seq[i]);
	    err++;
	    break;
	}
    }

    if ((s->name && s->name_len != strlen(s->name)) ||
	(!s->name && s->name_len != 0)) {
	vmessage("Seq %"PRIrec": name and name_len do not match\n", s->rec);
	err++;
    }

    return err;
}

/*
 * Checks a anno_ele struct is consistent, both internally and with the
 * range that points to it.
 */
static int check_anno(GapIO *io, bin_index_t *bin, range_t *r) {
    int err = 0;
    anno_ele_t *a = cache_search(io, GT_AnnoEle, r->rec);

    cache_incr(io, a);

    /* Bin records match */
    if (a->bin != bin->rec) {
	vmessage("Anno %"PRIrec": bin does not match observed bin\n", a->rec);
	err++;
    }

    /* Types match */
    switch (a->obj_type) {
    case GT_Contig:
	if (r->flags & GRANGE_FLAG_TAG_SEQ) {
	    vmessage("Anno %"PRIrec": range flags indicate sequence, but anno "
		     "obj_type claims contig\n", a->rec);
	    err++;
	}
	break;

    case GT_Seq:
	if (!(r->flags & GRANGE_FLAG_TAG_SEQ)) {
	    vmessage("Anno %"PRIrec": range flags indicate contig, but anno "
		     "obj_type claims sequence\n", a->rec);
	    err++;
	}
	break;

    default:
	vmessage("Anno %"PRIrec": unrecognised value (%d) for obj_type\n",
		 a->rec, a->obj_type);
	err++;
    }

    /* That anno links to the correct data type */
    if (a->obj_rec && !cache_exists(io, a->obj_type, a->obj_rec)) {
	vmessage("Anno %"PRIrec": obj_rec/type do not match\n", a->rec);
	err++;
    }

    /* Consistency of cached portions in range */
    if (a->tag_type != r->mqual) {
	vmessage("Anno %"PRIrec": type does not match copy in "
		"range_t struct.\n", a->rec);
	err++;
    }
    if (a->obj_type != GT_Contig && a->obj_rec != r->pair_rec) {
	vmessage("Anno %"PRIrec": obj_rec does not match copy "
		 "in range_t struct.\n", a->rec);
	err++;
    }

    /*
     * Also check position of item being annotated - it should overlap
     * the tag and be in the same contig. This is a time consuming check
     * though.
     */
    {
	int astart, aend, ostart, oend;
	tg_rec acontig, ocontig;

	if (-1 == bin_get_item_position(io, GT_AnnoEle, a->rec,
					&acontig, &astart, &aend, NULL,
					NULL, NULL, NULL)) {
	    cache_decr(io, a);
	    return err+1;
	}

	if (a->obj_type == GT_Contig) {
	    contig_t *c;

	    ocontig = a->obj_rec;
	    if (ocontig && ocontig != acontig) {
		/* Previously an error, but removed assumptions on this */
		vmessage("Anno %"PRIrec": non-zero obj_rec with obj_type "
			 "GT_Contig is deprecated\n", a->rec);
	    }
	    ocontig = acontig;
	    c = cache_search(io, GT_Contig, ocontig);
	    if (!c) {
		cache_decr(io, a);
		return err+1;
	    }
	    ostart = c->start;
	    oend = c->end;
	} else {
	    if (-1 == bin_get_item_position(io, a->obj_type, a->obj_rec,
					    &ocontig, &ostart, &oend, NULL,
					    NULL, NULL, NULL)) {
		cache_decr(io, a);
		return err+1;
	    }
	}

	if (acontig != ocontig || astart < ostart || aend > oend) {
	    vmessage("Anno %"PRIrec": does not overlap annotated object:\n",
		     a->rec);
	    vmessage("\tTag Ctg %"PRIrec" at %d..%d\n", acontig, astart, aend);
	    vmessage("\tObj Ctg %"PRIrec" at %d..%d\n", ocontig, ostart, oend);
	    
	    err++;
	}
    }

    cache_decr(io, a);
    return err;
}

/*
 * Checks a range of type GRANGE_FLAG_ISREFPOS
 */
static int check_refpos(GapIO *io, range_t *r) {
    int err = 0;

    if (r->start != r->end) {
	vmessage("RefPos %"PRIrec": start and end positions differ\n", r->rec);
	err++;
    }

    return err;
}

/*
 * Walks a contig bin tree, executing callbacks per bin.
 *
 * Returns 0 on success
 *         number of errors on failure
 */
static int bin_walk(GapIO *io, int fix, tg_rec rec, int offset, int complement,
		    int level, HacheTable *lib_hash, bin_stats *bs) {
    bin_index_t *bin;
    int i, f_a, f_b, err = 0;
    bin_stats child_stats;
    int start, end, cstart, cend, nthis_seq = 0;
    int valid_range;
    int db_vers = io->base ? io->base->db->version : io->db->version;

    if (!rec)
	return 1;

    if (NULL == (bin = cache_search(io, GT_Bin, rec)))
	return 1;

    /* Set f_a & f_b for NMIN/NMAX macros & handle complementing */
    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }
    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

#ifdef DEBUG_CHECK
    printf("Bin %"PRIrec"  %d..%d  nseq=%d\n",
	   bin->rec,
	   NMIN(bin->pos, bin->pos + bin->size-1),
	   NMAX(bin->pos, bin->pos + bin->size-1),
	   bin->nseqs);
#endif

    cache_incr(io, bin);


    /* Recurse */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;
	ch = get_bin(io, bin->child[i]);
	if (!ch) {
	    vmessage("bin %"PRIrec" failed to load.\n",
		     bin->child[i]);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->child[i] = 0;
	    }
	    continue;
	}

	child_stats.cstart = INT_MAX;
	child_stats.cend   = INT_MIN;
	child_stats.nseq   = 0;
	child_stats.nanno  = 0;
	child_stats.nref   = 0;

	err += bin_walk(io, fix, bin->child[i],
			NMIN(ch->pos, ch->pos + ch->size-1) /* offset */,
			complement, level, lib_hash, &child_stats);

	bs->nseq  += child_stats.nseq;
	bs->nanno += child_stats.nanno;
	bs->nref  += child_stats.nref;

	if (child_stats.nseq) {
	    if (bs->cstart > child_stats.cstart)
		bs->cstart = child_stats.cstart;
	    if (bs->cend   < child_stats.cend)
		bs->cend   = child_stats.cend;
	}
    }


    /* Check range items */
    start  = INT_MAX;
    end    = INT_MIN;
    cstart = INT_MAX;
    cend   = INT_MIN;
    valid_range = 0;
    if (bin->rng) {
	char *loop;
	int last_range = -1;

	/* Check rng_free list */
	if (bin->rng_free < -1 ||
	    (bin->rng_free != -1 && bin->rng_free >= ArrayMax(bin->rng))) {
	    vmessage("bin %"PRIrec": rng_free has invalid value %d\n",
		     bin->rec, bin->rng_free);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->rng_free = -1;
		bin->flags |= BIN_BIN_UPDATED;
	    } else {
		goto skip_checks;
	    }
	}

	loop = (char *)calloc(ArrayMax(bin->rng), sizeof(char));
	i = bin->rng_free;
	while (i != -1) { 
	    range_t *r = arrp(range_t, bin->rng, i);

	    if (!(r->flags & GRANGE_FLAG_UNUSED)) {
		vmessage("bin %"PRIrec" free range item %d: "
			 "flagged as in use\n", bin->rec, i);
		err++;
		if (fix) {
		    bin = cache_rw(io, bin);
		    /* Terminate rng_free list and just accept the
		     * chance of leaking.
		     */
		    if (last_range != -1) {
			r = arrp(range_t, bin->rng, last_range);
			r->rec = -1;
			bin->flags |= BIN_RANGE_UPDATED;
		    } else {
			bin->rng_free = -1;
			bin->flags |= BIN_BIN_UPDATED;
		    }
		}
		break;
	    }

	    if ((int)r->rec < -1 || 
		((int)r->rec != -1 && (int)r->rec >= ArrayMax(bin->rng))) {
		vmessage("bin %"PRIrec" free range item %d: "
			 "next rec (r->rec) is invalid\n", bin->rec, i);
		err++;
		if (fix) {
		    bin = cache_rw(io, bin);
		    r->rec = -1;
		    bin->flags |= BIN_RANGE_UPDATED;
		}
		break;
	    }

	    /* Loop detection */
	    if (loop[i]) {
		vmessage("bin %"PRIrec": loop detected in free list\n",
			 bin->rec);
		err++;
		if (fix) {
		    bin = cache_rw(io, bin);
		    if (last_range != -1) {
			r = arrp(range_t, bin->rng, last_range);
			r->rec = -1;
			bin->flags |= BIN_RANGE_UPDATED;
		    } else {
			bin->rng_free = -1;
			bin->flags |= BIN_BIN_UPDATED;
		    }
		}
		break;
	    }
	    loop[i] = 1;

	    last_range = i;
	    i = (int)r->rec;
	}
	free(loop);

    skip_checks:

	/* Iterate through USED bin items and check them */
	for (i = 0; i < ArrayMax(bin->rng); i++) {
	    range_t *r = arrp(range_t, bin->rng, i);

	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    valid_range = 1;

#ifdef DEBUG_CHECK
	    printf("Range item %d (%"PRIrec" flag %d): %d..%d (abs %d..%d)\n",
		   i, r->rec, r->flags, r->start, r->end,
		   NMIN(r->start, r->end), NMAX(r->start, r->end));
#endif
	    if (start > r->start)
		start = r->start;
	    if (end   < r->end)
		end   = r->end;

	    switch (r->flags & GRANGE_FLAG_ISMASK) {
	    case GRANGE_FLAG_ISSEQ: {
		bs->nseq++;
		nthis_seq++;
		if (level > 1)
		    err += check_seq(io, fix, bin, r, lib_hash, 0);

		if (cstart > r->start)
		    cstart = r->start;
		if (cend   < r->end)
		    cend   = r->end;
		break;
	    }

	    case GRANGE_FLAG_ISANNO:
		bs->nanno++;
		if (level > 1)
		    err += check_anno(io, bin, r);
		break;

	    case GRANGE_FLAG_ISREFPOS:
		bs->nref++;
		if (level > 1)
		    err += check_refpos(io, r);
		break;

	    case GRANGE_FLAG_ISCONS:
		if (level > 1)
		    err += check_seq(io, fix, bin, r, NULL, 1);
		break;

	    case GRANGE_FLAG_ISUMSEQ:
	    case GRANGE_FLAG_ISREF:
		/* FIXME: check */
		break;

	    default:
		/* Unknown ISMASK flag */
		vmessage("bin %"PRIrec" range item %d: "
			 "Unknown GRANGE_FLAG_IS? flag: %d\n",
			 bin->rec, i, r->flags & GRANGE_FLAG_ISMASK);
		err++;
	    }
	}
    }


    /* Check count validity to ensure this + children are correct */
    if (bin->nseqs != bs->nseq) {
	vmessage("bin %"PRIrec": nseqs does not match observed counts\n",
		 bin->rec);
	if (fix) {
	    bin = cache_rw(io, bin);
	    bin->flags |= BIN_BIN_UPDATED;
	    bin->nseqs = bs->nseq;
	}
	err++;
    }
    if (db_vers > 1 && bin->nanno != bs->nanno) {
	vmessage("bin %"PRIrec": nanno does not match observed counts\n",
		 bin->rec);
	if (fix)
	    bin->nanno = bs->nanno;
	err++;
    }
    if (db_vers > 1 && bin->nrefpos != bs->nref) {
	vmessage("bin %"PRIrec": nrefpos does not match observed counts\n",
		 bin->rec);
	if (fix) {
	    bin = cache_rw(io, bin);
	    bin->flags |= BIN_BIN_UPDATED;
	    bin->nrefpos = bs->nref;
	}
	err++;
    }


    /*
     * Check used start/end range, and accumulate absolute positions so we
     * can check the contig later.
     */
    if (valid_range) {
	if (start != bin->start_used ||
	    end   != bin->end_used) {
	    vmessage("bin %"PRIrec": used start/end range are incorrect\n",
		     bin->rec);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->start_used = start;
		bin->end_used = end;
	    }
	}

	if (nthis_seq) {
	    if (bs->cstart > NMIN(cstart, cend))
		bs->cstart = NMIN(cstart, cend);
	    if (bs->cend   < NMAX(cstart, cend))
		bs->cend   = NMAX(cstart, cend);
	}
    } else {
	if (bin->start_used != 0 || bin->end_used != 0) {
	    vmessage("bin %"PRIrec": used start/end are non-zero "
		     "in an empty bin\n", bin->rec);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->start_used = 0;
		bin->end_used = 0;
	    }
	}
    }

    cache_decr(io, bin);

    return err;
}


/*
 * Checks a single contig. lib_hash may be NULL, but if not then it is a
 * HacheTable keyed on library record numbers. We then validate all sequences
 * against this to check that no unknown libraries are present.
 *
 * Returns the number of errors found
 *         0 on success.
 */
int check_contig(GapIO *io, tg_rec crec, int fix, int level,
		 HacheTable *lib_hash) {
    contig_t *c;
    bin_stats bs;
    int err;

    if (!cache_exists(io, GT_Contig, crec)) {
	vmessage("Record %"PRIrec" is not a contig, but in contig order\n",
		 crec);
	return 1;
    }

    c = cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    bs.cstart = INT_MAX;
    bs.cend   = INT_MIN;
    bs.nseq   = 0;
    bs.nanno  = 0;
    bs.nref   = 0;

    err = bin_walk(io, fix, c->bin, contig_offset(io, &c), 0, level,
		   lib_hash, &bs);
    if (bs.cstart != c->start ||
	bs.cend   != c->end) {
	vmessage("Contig %"PRIrec": used start/end range are incorrect\n",
		 crec);
	err++;
	if (fix) {
	    c = cache_rw(io, c);
	    c->start = bs.cstart;
	    c->end   = bs.cend;
	}
    }
    
    cache_decr(io, c);

    return err;
}


/*
 * Performs a thorough internal consistency check of all on disk data
 * structures. It's therefore quite slow, but can highlight algorithm
 * problems during development.
 *
 * Returns the number of errors found
 *         or 0 on success.
 */
int check_database(GapIO *io, int fix, int level) {
    database_t *db;
    ArrayStruct *contig_order, *library;
    int i;
    int err = 0;
    HacheTable *hash = NULL;

    vfuncheader("Check Database");

    /* Load low level db structs; already in GapIO, but do full reload */
    db = cache_search(io, GT_Database, 0);
    if (!db) {
	vmessage("Failed to read GT_Database record 0\n");
	return ++err;
    }
    cache_incr(io, db);


    /* Check contig_order: 1 per contig and no more than 1 */
    contig_order = cache_search(io, GT_RecArray, db->contig_order);
    if (!contig_order) {
	vmessage("Failed to read contig order array\n");
	cache_decr(io, db);
	return ++err;
    }
    cache_incr(io, contig_order);

    hash = HacheTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    if (db->Ncontigs != ArrayMax(contig_order)) {
	vmessage("Contig order array is not the same size as db->Ncontigs\n");
	err++;
	if (fix) {
	    cache_rw(io, io->contig_order);
	    ArrayMax(io->contig_order) = io->db->Ncontigs;
	    ArrayMax(contig_order) = io->db->Ncontigs;
	}
    }

    for (i = 0; i < ArrayMax(contig_order); i++) {
	tg_rec crec = arr(tg_rec, contig_order, i);
	HacheData hd;
	int new;

	hd.i = 0;
	HacheTableAdd(hash, (char *)&crec, sizeof(crec), hd, &new);
	if (!new) {
	    vmessage("Contig %"PRIrec" occurs more than once in the "
		     "contig_order array\n", crec);
	    err++;
	}
    }
    HacheTableDestroy(hash, 0);
    cache_decr(io, contig_order);


    /* Also check contig btree - every contig should have name in index */
    /* TODO */


    /* Check libraries, keep hash for later checking too */
    library = cache_search(io, GT_RecArray, db->library);
    if (!library) {
	vmessage("Failed to read library array\n");
	cache_decr(io, db);
	return ++err;
    }
    cache_incr(io, library);


    hash = HacheTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    if (db->Nlibraries != ArrayMax(library)) {
	vmessage("library array is not the same size as db->Nlibraries\n");
	err++;
    }
    for (i = 0; i < ArrayMax(library); i++) {
	tg_rec lrec = arr(tg_rec, library, i);
	HacheData hd;
	int new;

	hd.i = 0;
	HacheTableAdd(hash, (char *)&lrec, sizeof(lrec), hd, &new);
	if (!new) {
	    vmessage("Library %"PRIrec" occurs more than once in the "
		     "library array\n", lrec);
	    err++;
	}
    }
    cache_decr(io, db);
    cache_decr(io, library);


    /* Check each contig in turn */
    for (i = 0; i < ArrayMax(contig_order); i++) {
	tg_rec crec = arr(tg_rec, contig_order, i);
	vmessage("--Checking contig #%"PRIrec" (%d of %d)\n",
		 crec, i+1, (int)ArrayMax(contig_order));
	UpdateTextOutput();
	err += check_contig(io, crec, fix, level, hash);
    }

    HacheTableDestroy(hash, 0);

    vmessage("\n*** Total number of errors: %d ***\n", err);

    return err;
}
