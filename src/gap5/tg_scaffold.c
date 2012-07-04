/*

An implementation of scaffolds in Gap5.

Scaffolds are collections of contigs, just as contigs are collections
of sequences. Traditionally the scaffold defines the ordering,
orientation and some minor placement (gap size) of contigs. In this
implementation we ignore orientation as we will complement the contig
objects themselves to represent this information.

Contigs have a parent object now - the scaffold. They also maintain a
list of links to other contigs. These links form a (directed?)  graph.
Generally this will be linear if we've loaded an AGP file, but it
could be a true graph if we've read data from a FASTG format file.

Scaffolds have lists of member objects (contigs) with gap sizes. The
order of items is largely irrelevant. The list is simply to permit
fast extraction and identification of siblings. The actual flattened
contig order is the old database contig_order array. However a
function update_scaffold_order() exists to coerce the scaffold contig
into to be the same as the contig_order. 

Users can manually manipulate the main io->contig_order array and
therefore implicitly the order of members within a scaffold. It is
expected that there will also be ways of drawing a scaffold on top of
the contig_link structure, but note this latter is the evidence used
to produce a scaffold and not the actual scaffold itself.

The main database object now holds a scaffold record array in addition
to the old contig_order array. The latter still needs updating as it
is the canonical order for iterating through flattened contigs, and
the user may even have elected to ignore scaffolds and sort by contig
size or contig name.  The order of this array however is irrelevant at
present. (db->contig_order is the master.)

Like sequences and more recently contigs, scaffolds are internally
stored in blocks of up to 1024 scaffolds to a GT_ScaffoldBlock.

Scaffold names must be unique and distinct from contig names. This
will permit specifying scaffold names to dialogues that can take lists
of contigs, as an indication to mean "all the contigs in this
scaffold". For simplicity we hold these names in the same B+Tree we
use for contigs, but prefix them with "/scaffold/" to force the
scaffolds to end up in their own branch of the tree and to explicitly
avoid clashing just incase we're forced down this road at some stage.

 */

#include <staden_config.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <string.h>
#include <errno.h>

#include "tg_gio.h"
#include "list_proc.h"

/*
 * Sets a scaffold name.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_set_name(GapIO *io, scaffold_t **f, char *name) {
    scaffold_t *n;
    GapIO *iob = io->base ? io->base : io;
    char cname[8192];

    sprintf(cname, "/scaffold/%s", name);

    if (!(n = cache_rw(io, *f)))
	return -1;

    /* Delete old name */
    if (n->name) {
	tg_rec r = iob->iface->contig.index_del(iob->dbh, n->name, n->rec);
	if (r != -1 && r != io->db->contig_name_index) {
	    io->db = cache_rw(io, io->db);
	    io->db->contig_name_index = r;
	}
    }

    if (NULL == (n = cache_item_resize(n, sizeof(*n) + strlen(name)+1)))
	return -1;

    *f = n;

    /* Add new name */
    n->name   = (char *)(&n->data);
    strcpy(n->name, name);

    if (*name) {
	tg_rec r = iob->iface->contig.index_add(iob->dbh, name, n->rec);
	if (r != -1 && r != io->db->contig_name_index) {
	    io->db = cache_rw(io, io->db);
	    io->db->contig_name_index = r;
	}
    }

    return 0;
}


/*
 * Creates a new named scaffold.
 *
 * Returns scaffold pointer on success.
 *         NULL on failure
 */
scaffold_t *scaffold_new(GapIO *io, char *name) {
    tg_rec rec;
    scaffold_t *f, init_f;

    if (!io->db->scaffold)
	return NULL;

    memset(&init_f, 0, sizeof(scaffold_t));
    init_f.name = name;

    /* Allocate our contig */
    rec = cache_item_create(io, GT_Scaffold, &init_f);

    /* Initialise it */
    f = (scaffold_t *)cache_search(io, GT_Scaffold, rec);
    f = cache_rw(io, f);

    if (name)
        scaffold_set_name(io, &f, name);
    else
        f->name = NULL;

    /* Add it to the scaffold order too */
    io->scaffold = cache_rw(io, io->scaffold);
    io->db = cache_rw(io, io->db);
    ARR(tg_rec, io->scaffold, io->db->Nscaffolds++) = rec;

    /* Add to the new contigs list */
    if (name)
	add_to_list("new_scaffolds", name);

    return f;
}


/*
 * Adds a contig to a scaffold array.
 * Gap size, type and evidence refer to the gap between this and the
 * "previous" contig - ie the last in the scaffold. More complex
 * scaffold manipulations will be handled elsewhere.
 *
 * Set these fields to 0 if you do not know them.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_add(GapIO *io, tg_rec scaffold, tg_rec contig,
		 int gap_size, int gap_type, int evidence) {
    scaffold_t *f;
    contig_t *c;
    scaffold_member_t *m;
    int i;
    
    if (!(f = cache_search(io, GT_Scaffold, scaffold)))
	return -1;

    /* Check if it already exists */
    for (i = 0; i < ArrayMax(f->contig); i++) {
	m = arrp(scaffold_member_t, f->contig, ArrayMax(f->contig)-1);
	if (m->rec == contig)
	    return 0;
    }

    /* Append */
    f = cache_rw(io, f);
    m = ArrayRef(f->contig, ArrayMax(f->contig)); // extend
    m->rec = contig;
    m->gap_size = ArrayMax(f->contig) > 1 ? gap_size : 0;
    m->gap_type = gap_type;
    m->evidence = evidence;

    /* Update the contig record too */
    c = cache_search(io, GT_Contig, contig);
    c = cache_rw(io, c);
    c->scaffold = scaffold;

#if 0
    /* Add a scaffold link to the contig graph too */
    if (ArrayMax(f->contig) >= 2) {
	m = arrp(scaffold_member_t, f->contig, ArrayMax(f->contig)-2);
	contig_link_t lnk;

	lnk.rec1 = contig;
	lnk.rec2 = m->rec;
	/* Best guess */
	lnk.pos1 = 0; lnk.end1 = 1;
	lnk.pos2 = 0; lnk.end2 = 0;
	lnk.orientation = 0;
	lnk.size = 100;
	lnk.type = CLINK_TYPE_SCAFFOLD;
	lnk.score = 0;

	contig_add_link(io, &lnk);
    }
#endif

    return 0;
}

tg_rec scaffold_index_query(GapIO *io, char *name) {
    char cname[8192];

    sprintf(cname, "/scaffold/%s", name); 
    return io->iface->contig.index_query(io->dbh, name, 0);
}

/*
 * Adds a contig named ctg_name to a scaffold named scaf_name. The names are
 * looked up in the B+Tree index.
 */
int scaffold_add_by_name(GapIO *io, char *scaf_name, char *ctg_name,
			 int gap_size, int gap_type, int evidence) {
    tg_rec srec, crec;

    if ((crec = contig_index_query(io, ctg_name)) <= 0)
	return -1;

    if ((srec = scaffold_index_query(io, scaf_name)) <= 0) {
	scaffold_t *f = scaffold_new(io, scaf_name);
	srec = f->rec;
    }

    return scaffold_add(io, srec, crec, gap_size, gap_type, evidence);
}

/*
 * Removes a contig from a scaffold.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_remove(GapIO *io, tg_rec scaffold, tg_rec contig) {
    scaffold_t *f;
    scaffold_member_t *m, *m2;
    contig_t *c;
    int i;

    c = cache_search(io, GT_Contig, contig);
    f = cache_search(io, GT_Scaffold, scaffold);

    if (!c || !f)
	return -1;

    if (c->scaffold != scaffold) {
	verror(ERR_WARN, "scaffold_remove", "Attempted to remove contig #%"
	       PRIrec" from a scaffold #%"PRIrec" it is not a member of",
	       contig, scaffold);
	return -1;
    }

    c = cache_rw(io, c);
    c->scaffold = 0;

    f = cache_rw(io, f);
    for (i = 0; i < ArrayMax(f->contig); i++) {
	m = arrp(scaffold_member_t, f->contig, i);
	if (m->rec == contig) {
	    /* Shuffle array down */
	    for (i++; i < ArrayMax(f->contig); i++) {
		m2 = arrp(scaffold_member_t, f->contig, i);
		*m = *m2;
		m = m2;
	    }
	    ArrayMax(f->contig)--;
	}
    }

    return 0;
}

/*
 * Loads a new scaffold from an AGP file.
 * Contigs that are not listed in this will keep their old scaffold data.
 * Contigs which are listed and have an existing scaffold will be amended.
 *
 * AGP format is http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
 *
 * It has multiple columns with the major ones for our purpose being
 * 1: Scaffold name
 * 5: N or U for a gap, anything else for a sequence object
 *
 * If a gap:
 * 6: gap size
 * 7: gap type
 * 8: linkage y/n
 * 9: linkage evidence
 *
 * If a sequence:
 * 6: contig name
 * 9: orientation
 *
 * There are also a host of start/end coordinates in AGP, but these are
 * ignored as they should match the contig sizes. (And if not it's an error
 * and we would want to take the contig sizes as the definitive versions.)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int scaffold_from_agp(GapIO *io, char *fn) {
    FILE *fp;
    char line[8192];
    char *scaf_name, *ctg_name, *cp;
    int lno = 0;
    int linkage = 0, gap_size = 0, gap_type = 0, evidence = 0, orientation = 0;

    if (NULL == (fp = fopen(fn, "r"))) {
	verror(ERR_WARN, "scaffold_from_agp", "%s: %s",
	       fn, strerror(errno));
	return -1;
    }

    while (fgets(line, 8192, fp)) {
	lno++;
	if (!(scaf_name = strtok(line, "\t"))) {
	    verror(ERR_WARN, "scaffold_from_agp", "Invalid data on line %d",
		   lno);
	    return -1;
	}

	strtok(NULL, "\t"); /* start */
	strtok(NULL, "\t"); /* end */
	strtok(NULL, "\t"); /* element number */
	cp = strtok(NULL, "\t"); /* type */

	if (*cp == 'N' || *cp == 'U' || *cp == 'n' || *cp == 'u') {
	    /* Gap */
	    cp = strtok(NULL, "\t"); gap_size = atoi(cp);
	    cp = strtok(NULL, "\t"); gap_type = (cp[0]<<8) + cp[1];
	    cp = strtok(NULL, "\t"); linkage  = (*cp=='y'||*cp=='Y');
	    cp = strtok(NULL, "\t"); evidence = 0; /* FIXME */
	} else {
	    ctg_name = strtok(NULL, "\t");
	    strtok(NULL, "\t"); /* component begin */
	    strtok(NULL, "\t"); /* component end */
	    cp = strtok(NULL, "\t"); orientation = cp && *cp=='-' ? 1 : 0;

	    scaffold_add_by_name(io, scaf_name, ctg_name,
				 gap_size, gap_type, evidence);
	}
    }

    fclose(fp);
    return 0;
}

/* ------------------------------------------------------------------------ */
typedef struct {
    tg_rec scaffold;
    int ctg_idx;
} scaf_ctg_t;

static int scaf_ctg_sort(const void *p1, const void *p2) {
    const scaf_ctg_t *s1 = (scaf_ctg_t *)p1;
    const scaf_ctg_t *s2 = (scaf_ctg_t *)p2;

    if (s1->scaffold < s2->scaffold)
	return -1;
    else if (s1->scaffold > s2->scaffold)
	return 1;
    else
	return s1->ctg_idx - s2->ctg_idx;
}

/*
 * Given a contig order and a set of current scaffolds, this updates the
 * order of entries within each scaffold to match the contig order.
 *
 * For example if we have contigs in order 1 3 5 2 6 8 4 7 9 and
 * scaffolds {1 2 3 4} {5 6 7 8 9} we would shuffle the scaffold members
 * to        {1 3 2 4} {5 6 8 7 9}
 *
 * The purpose is for integration with contig shuffling in the Contig List
 * or Contig Selector. The master contig order array is what gets shuffled
 * manually by the user and it is also the definitive order to use when
 * outputting data (so it is completely under users control whether they
 * sort by name, size or scaffold).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int update_scaffold_order(GapIO *io) {
    int i, j, ret = -1;
    int nc;
    int ns;
    tg_rec *crecs;

    if (!io->scaffold)
	return 0; /* Not supported, but considered success */

    nc = ArrayMax(io->contig_order);
    ns = ArrayMax(io->scaffold);

    scaf_ctg_t *a = (scaf_ctg_t *)malloc(nc * sizeof(*a));
    if (!a)
	return -1;

    /*
     * Produce an array of scaffold and contig recs, so we can sort on
     * both fields.
     */
    crecs = ArrayBase(tg_rec, io->contig_order);
    for (i = 0; i < nc; i++) {
	contig_t *c = cache_search(io, GT_Contig, crecs[i]);
	if (!c)
	    goto err;

	a[i].ctg_idx = i;
	a[i].scaffold = c->scaffold;
    }

    qsort(a, nc, sizeof(*a), scaf_ctg_sort);

    /*
     * Now recreate scaffold orders from the sorted contig list.
     */
    for (i = 0; i < nc; i++) {
	scaffold_t *f;
	int k;

	if (!a[i].scaffold)
	    continue;

	j = i;
	while (i < nc && a[i].scaffold == a[j].scaffold)
	    i++;

	/* j .. i-1 share the same scaffold */
	f = cache_search(io, GT_Scaffold, a[j].scaffold);
	if (!f)
	    goto err;

	if (!f->contig || ArrayMax(f->contig) != i-j) {
	    verror(ERR_WARN, "update_scaffold_order", "Scaffold %"PRIrec
		   "has different number of entries than contigs claim.",
		   f->rec);
	    goto err;
	}

	/* Only mark r/w and update if they differ */
	for (k = 0; k < ArrayMax(f->contig); k++) {
	    if ((arrp(scaffold_member_t, f->contig, k))->rec
		!= crecs[a[j+k].ctg_idx])
		break;
	}
	
	if (k != ArrayMax(f->contig)) {
	    f = cache_rw(io, f);
	    for (k = 0; k < ArrayMax(f->contig); k++)
		(arrp(scaffold_member_t, f->contig, k))->rec
		    = crecs[a[j+k].ctg_idx];
	}
	
	i--;
    }

    ret = 0;
 err:
    free(a);
    return ret;
}