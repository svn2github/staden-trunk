#define USE_HASH_TABLES

#include <tcl.h>

#include "xalloc.h"
#include "IO.h"
#include "vseqs.h"
#include "text_output.h"
#include "misc.h"

/*
 * new_vcontig
 *
 * Creates a virtual contig. This creates the linked list of vrseqs,
 * which initially just contain references to the real sequences in this
 * contig.
 * The vcontig_t should destroyed by calling del_vcontig.
 *
 * Arguments:
 *	io		Gap4 IO handle
 *	db		A vcontig_t database structure
 *	contig		The contig number to copy into db
 *
 * Returns:
 *	The virtual contig (vcontig_t *) on success.
 *	NULL on failure.
 */
vcontig_t *new_vcontig(GapIO *io, int contig) {
    int rnum;
    vcontig_t *db;
    vrseq_t *vrseq, *last;

    if (NULL == (db = (vcontig_t *)xmalloc(sizeof(*db))))
	return NULL;

#ifdef USE_HASH_TABLES
    /* Hash table to map reading numbers to vrseq_t pointers */
    Tcl_InitHashTable(&db->num_hash, TCL_ONE_WORD_KEYS);
#endif

    /* Allocate vrseqs (one per real read) and link them together */
    last = NULL;
    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	if (NULL == (vrseq = (vrseq_t *)xmalloc(sizeof(*vrseq))))
	    return NULL;

	vrseq->vseq = NULL;
	vrseq->rnum = rnum;
	vrseq->position = io_relpos(io, rnum);
	vrseq->left = last;
	if (last)
	    last->right = vrseq;
	else
	    db->left = vrseq;

#ifdef USE_HASH_TABLES
	{
	    Tcl_HashEntry *hash;
	    int new_entry;
	    hash = Tcl_CreateHashEntry(&db->num_hash, (char *)rnum,
				       &new_entry);
	    Tcl_SetHashValue(hash, (ClientData)vrseq);
	}
#endif
	last = vrseq;
    }
    vrseq->right = NULL;
    db->right = vrseq;

    /* Initialise other db records */
    db->io = io;
    db->contig = contig;
    db->next_rnum = NumReadings(io)+1;
    db->cons = NULL;

    return db;
}

/*
 * del_vcontig
 *
 * This frees the memory allocated within (and including) the vcontig, as
 * allocated by new_vcontig and subsequent functions.
 *
 * Arguments:
 *	vc		The virtual contig pointer to free
 *
 * Returns:
 * 	None
 */
void del_vcontig(vcontig_t *vc) {
    vrseq_t *v, *v2;

    v = vc->left;
    while (v) {
	v2 = v->right;
	if (v->vseq) {
	    if (v->vseq->seq)
		xfree(v->vseq->seq);
	    if (v->vseq->conf)
		xfree(v->vseq->conf);
	    xfree(v->vseq);
	}
	xfree(v);
	v = v2;
    }

#ifdef USE_HASH_TABLES
    Tcl_DeleteHashTable(&vc->num_hash);
#endif

    xfree(vc);
}

/*
 * new_vrseq
 *
 * Allocates and blanks a new virtual sequence within a virtual contig.
 * Although we return a vrseq_t (virtual or real) this contains a reference
 * to the new virtual sequence.
 *
 * Arguments:
 *	vc		The virtual contig to place this virtual sequence in
 *
 * Returns:
 *	On success; A vrseq_t struct referencing the new virtual sequence.
 *	On failure; NULL
 */
vrseq_t *new_vrseq(vcontig_t *vc) {
    vrseq_t *vrseq;
    vseq_t *vseq;

    /* Allocate a new vrseq */
    if (NULL == (vrseq = (vrseq_t *)xmalloc(sizeof(*vrseq))))
	return NULL;

    vrseq->left = NULL;
    vrseq->right = NULL;

    /* Allocate a new vseq */
    if (NULL == (vseq = (vseq_t *)xcalloc(1, sizeof(*vseq))))
	return NULL;

    vseq->seq = NULL;
    vseq->conf = NULL;

    /* Link the vseq to vrseq */
    vrseq->rnum = vc->next_rnum++;
    vrseq->vseq = vseq;
    vrseq->position = 0;

#ifdef USE_HASH_TABLES
    /* Add hash table entry */
    {
	Tcl_HashEntry *hash;
	int new;

	hash = Tcl_CreateHashEntry(&vc->num_hash, (char *)(vrseq->rnum), &new);
	Tcl_SetHashValue(hash, (ClientData)vrseq);
    }
#endif

    /* printf("Add virtual sequence %d\n", vrseq->rnum); */
    
    return vrseq;
}

/*
 * del_vrseq
 *
 * Deletes a virtual sequence from the database. Although we pass in a vrseq,
 * we cannot delete this if it references a real sequence.
 * This also relinks any neighbours of this vrseq as appropriate.
 *
 * Arguments:
 *	vc		The virtual contig containing the sequence.
 *	vrseq		The real/virtual sequence referencing the virtual seq.
 *
 * Returns:
 *	None
 */
void del_vrseq(vcontig_t *vc, vrseq_t *vrseq) {
    if (!vrseq || !vc)
	return;

#if 0
    if (!vrseq->vseq) {
	verror(ERR_WARN, "del_vrseq",
	       "virtual attempt to delete a real sequence\n");
	return;
    }
#endif

    /* printf("Del virtual sequence %d\n", vrseq->rnum); */

    /* Check contig ends */
    if (vc->left == vrseq) {
	vc->left = vrseq->right;
    }
    if (vc->right == vrseq) {
	vc->right = vrseq->left;
    }

    /* Relink neighbours */
    if (vrseq->left)
	vrseq->left->right = vrseq->right;
    if (vrseq->right)
	vrseq->right->left = vrseq->left;

#ifdef USE_HASH_TABLES
    /* Delete hash table entry */
    {
	Tcl_HashEntry *hash;

	hash = Tcl_FindHashEntry(&vc->num_hash, (char *)(vrseq->rnum));
	if (hash)
	    Tcl_DeleteHashEntry(hash);
    }
#endif


    /* Deallocate */
    if (vrseq->vseq) {
	if (vrseq->vseq->seq)
	    xfree(vrseq->vseq->seq);
	if (vrseq->vseq->conf)
	    xfree(vrseq->vseq->conf);
	xfree(vrseq->vseq);
    }
    xfree(vrseq);
}

/*
 * link_vrseq
 *
 * This links a vrseq into a virtual contig. In order to do this it needs
 * to find the correct offset within the contig and insert vrseq into the
 * doubly linked list.
 * (FIXME: For speed, we may wish to keep an array of positions so that we
 * can binary search this to identify the insertion point quickly)
 *
 * If the sequence being linked in doesn't have a sequence or confidence
 * then we create them based on the consensus and expected confidence
 * values.
 *
 * Arguments:
 *	vc		The virtual contig
 *	vrseq		The real/virtual sequence
 *	position	Where to place the real/virtual sequence
 *
 * Returns:
 *	None
 */
void link_vrseq(vcontig_t *vc, vrseq_t *vrseq, int position) {
    vrseq_t *vr;
    int pos;
    int len;

    /* printf("Link seq %d\n", vrseq->rnum); */

    vrseq->position = position;

    /*
     * Find insertion point.
     * After this we can insert in between vr and vr->left.
     */
    for (vr = vc->left; vr; vr = vr->right) {
	if (vr->position >= position)
	    break;
    }

   
    /* Link in sequence */
    if (!vr) {
	/* At end of contig */
	vc->right->right = vrseq;
	vrseq->left = vc->right;
	vrseq->right = NULL;
	vc->right = vrseq;
    } else if (!vr->left) {
	/* At start of contig */
	vrseq->left = NULL;
	vrseq->right = vr;
	vr->left = vrseq;
	vc->left = vrseq;
    } else {
	/* Middle of contig */
	vrseq->left = vr->left;
	vrseq->right = vr;
	vr->left->right = vrseq;
	vr->left = vrseq;
    }
    
    if (!vrseq->vseq)
	return;

    pos = vrseq->position;
    len = vrseq->vseq->r.end - vrseq->vseq->r.start - 1;

    /* Create a virtual sequence from the consensus */
    if (!vrseq->vseq->seq) {
	int i;

	if (!vc->cons) {
	    fprintf(stderr, "No consensus - hence no virtual sequence");
	    return;
	}

	vrseq->vseq->seq = (char *)xmalloc(len+1);
	if (pos > 0 && pos + len <= io_clength(vc->io, vc->contig)) {
	    for (i = 0; i < len; i++) {
		if (vc->cons[pos-1+i] != '-' &&
		    vc->cons[pos-1+i] != 'N')
		    vrseq->vseq->seq[i] = vc->cons[pos-1+i];
		else
		    vrseq->vseq->seq[i] = 'A';
	    }
	} else {
	    for (i = 0; i < len; i++) {
		if (pos+i <= 0 ||
		    pos+i > io_clength(vc->io, vc->contig)) {
		    vrseq->vseq->seq[i] = 'A'; /* Make something up! */
		} else {
		    vrseq->vseq->seq[i] = vc->cons[pos+i-1];
		    if (vrseq->vseq->seq[i] == '-' ||
			vrseq->vseq->seq[i] == 'N')
			vrseq->vseq->seq[i] = 'A'; /* made up! */
		}
	    }
	}
    }

    /* Generate expected confidence values */
    if (!vrseq->vseq->conf) {
	int i, j;
	vseq_qdist_t qdist;

	/* FIXME: Hard code an expected quality distribution */
	qdist.q_start[0] = 15;
	qdist.q_end  [0] = 40;
	qdist.p_start[0] = 0;
	qdist.p_end  [0] = 10;

	qdist.q_start[1] = 40;
	qdist.q_end  [1] = 40;
	qdist.p_start[1] = 10;
	qdist.p_end  [1] = 50;

	qdist.q_start[2] = 40;
	qdist.q_end  [2] = 35;
	qdist.p_start[2] = 50;
	qdist.p_end  [2] = 70;

	qdist.q_start[3] = 35;
	qdist.q_end  [3] = 15;
	qdist.p_start[3] = 70;
	qdist.p_end  [3] = 100;

	qdist.n_items    = 4;
	
	vrseq->vseq->conf = (int1 *)xmalloc(sizeof(int1) * (len+1));
	for (j = 0; j < qdist.n_items; j++) {
	    int i_start, i_end;
	    double q, q_inc;
	    /*
	     * FIXME: A "Total Hack"(TM) incoming...
	     *
	     * For fragments of sequence (eg where we know it goes into
	     * vector) we do not want to apply the quality distribution
	     * squashed down into len bases. Hence we have a minimum length
	     * of 400 to apply it to.
	     */
	    i_start = MAX(400,len) * qdist.p_start[j]/100.0;
	    i_end   = MAX(400,len) * qdist.p_end[j]  /100.0;

	    if (i_start >= i_end)
		continue;

	    q       = qdist.q_start[j];
	    q_inc   = (qdist.q_end[j] - qdist.q_start[j]) /
		(double)(i_end - i_start);
	    
	    for (i = i_start; i < i_end && i < len; i++) {
		vrseq->vseq->conf[i] = q;
		q += q_inc;
	    }
	}
	/* Just in case of rounding errors */
	for (; i < len; i++)
	    vrseq->vseq->conf[i] = 0;

	if (vrseq->vseq->r.sense) {
	    /* Complemented; reverse confidence */
	    for (i = 0, j = len-1; i < j; i++, j--) {
		int1 tmp;
		tmp = vrseq->vseq->conf[i];
		vrseq->vseq->conf[i] = vrseq->vseq->conf[j];
		vrseq->vseq->conf[j] = tmp;
	    }
	}
    }
}

/* FIXME: This is horribly slow! */
vrseq_t *vrseq_index2ptr(vcontig_t *vc, int num) {
#ifdef USE_HASH_TABLES
    Tcl_HashEntry *hash;

    hash = Tcl_FindHashEntry(&vc->num_hash, (char *)num);
    if (hash)
	return (vrseq_t *)Tcl_GetHashValue(hash);
    else
	return NULL;

#else

    vrseq_t *vr;

    vr = vc->left;
    for (vr = vc->left; vr; vr = vr->right)
	if (vr->rnum == num)
	    break;

    return vr;
#endif
}

/*
 * virtual_info_func
 *
 * This is a query function suitable to pass as an argument to the consensus
 * calculate (calc_consensus) and fragment finder (find_fragments) functions.
 *
 * This is based on the database_info function found in qualIO.c, but it
 * extends this by allowing the addition of sequences that are not in the
 * database; "fabricated sequences" if you wish. The purpose of this is to
 * allow the consensus calculation and find_fragments functions to operate
 * on our expected experimental results without actually doing them.
 *
 * Arguments:
 *	job		The type of information we are requesting.
 *	mydata		Client data used by this function (through parent) 
 *	theirdata	Both input and output: where to put results of job
 *
 * Returns:
 *	No return value, but the clientdata (depth_t) is modified.
 */
int virtual_info_func(int job, void *mydata, info_arg_t *theirdata) {
    vcontig_t *vc = (vcontig_t *)mydata;
    GapIO *io;

    if (vc == NULL || (io = vc->io) == NULL)
	return -1;

    switch (job) {
    case GET_SEQ:
	{
	    gel_seq_t *gel_seq = &theirdata->gel_seq;
	    vrseq_t *vr;

	    vr = vrseq_index2ptr(vc, gel_seq->gel);
	    if (!vr)
		return -1;

	    if (vr->vseq) {
		char *g_seq;
		int1 *g_conf;
		int pos, len;

		pos = vr->position;
		len = vr->vseq->r.end - vr->vseq->r.start - 1;
		g_seq = (char *)xmalloc(len+1);
		g_conf = (int1 *)xmalloc(sizeof(int1) * (len+1));

		memcpy(g_seq, vr->vseq->seq, len);
		memcpy(g_conf, vr->vseq->conf, len);
		

		gel_seq->gel_start  = 0;
		gel_seq->gel_end    = len+1;
		gel_seq->gel_seq    = g_seq;
		gel_seq->gel_conf   = g_conf;
		gel_seq->gel_opos   = NULL;
		gel_seq->gel_length = len;

	    } else {
		char *g_seq  = NULL;
		int1 *g_conf = NULL;
		int g_start, g_end, g_len;
	    
		/*printf("GET_INFO real seq %d/%d\n", vr->rnum, gel_seq->gel);*/

		if (0 != io_aread_seq(io, gel_seq->gel, &g_len, &g_start,
				     &g_end, &g_seq, &g_conf,
				      /*g_opos*/ NULL, 0)) {
		    if (g_seq)  xfree(g_seq);
		    if (g_conf) xfree(g_conf);
		    return -1;
		}

		gel_seq->gel_start  = g_start;
		gel_seq->gel_end    = g_end;
		gel_seq->gel_seq    = g_seq;
		gel_seq->gel_conf   = g_conf;
		gel_seq->gel_opos   = /*g_opos*/NULL;
		gel_seq->gel_length = g_len;
	    }
	    
	    return 0;
	}
    case DEL_SEQ:
	{
	    gel_seq_t *gel_seq = &theirdata->gel_seq;

	    if (gel_seq->gel_seq) xfree(gel_seq->gel_seq);
	    if (gel_seq->gel_conf) xfree(gel_seq->gel_conf);
	    /*if (gel_seq->gel_opos) xfree(gel_seq->gel_opos);*/
	    
	    return 0;
	}
    case GET_CONTIG_INFO:
	{
	    contig_info_t *contig_info = &theirdata->contig_info;

	    contig_info->length = io_clength(io, vc->contig);
	    contig_info->leftgel = vc->left->rnum;
	    /*
	    GContigs c;

	    GT_Read(io, arr(GCardinal, io->contigs, contig_info->contig-1),
		    &c, sizeof(c), GT_Contigs);

	    contig_info->length  = c.length;
	    contig_info->leftgel = c.left;
	    */
	    
	    return 0;
	}
    case DEL_CONTIG_INFO:
	{
	    return 0;
	}
    case GET_GEL_INFO:
	{
	    gel_info_t *gel_info = &theirdata->gel_info;
	    vrseq_t *vr;

	    vr = vrseq_index2ptr(vc, gel_info->gel);

	    if (!vr) {
		printf("GET_INFO: No seq %d\n", gel_info->gel);
		return -1;
	    }

	    gel_info->next_right   = vr->right ? vr->right->rnum : 0;

	    if (vr->vseq) {
		/*printf("GET_INFO virtual seq %d/%d\n", vr->rnum,
		  gel_info->gel);*/
		gel_info->unclipped_len= vr->vseq->r.length;
		gel_info->length       = vr->vseq->r.end - vr->vseq->r.start - 1;
		gel_info->complemented = vr->vseq->r.sense;
		gel_info->position     = vr->position;
		gel_info->start	       = 1;
		gel_info->as_double    =
		    vr->vseq->r.chemistry & GAP_CHEM_TERMINATOR;
		gel_info->template     = vr->vseq->r.template;
	    } else {
		GReadings r;
		gel_read(io, gel_info->gel, r);

		/*printf("GET_INFO real seq %d/%d\n", vr->rnum, gel_info->gel);*/

		gel_info->length       = r.end - r.start - 1;
		gel_info->unclipped_len= r.length;
		gel_info->complemented = r.sense;
		gel_info->position     = r.position;
		gel_info->as_double    = r.chemistry & GAP_CHEM_TERMINATOR;
		gel_info->start	       = r.start;
		gel_info->template     = r.template;
		/*
		  gel_info->as_double	   =
		  ((r.chemistry & GAP_CHEM_TERMINATOR) == GAP_CHEM_TERMINATOR) &&
		  ((r.chemistry & GAP_CHEM_TYPE_MASK) == GAP_CHEM_TYPE_BIGDYE);
		*/
	    }
	    
	    return 0;
	}
    case DEL_GEL_INFO:
	{
	    return 0;
	}
    case GET_GEL_LEN:
	{
	    return max_gel_len(io);
	}
    default:
	verror(ERR_FATAL, "database_info", "Unknown job number (%d)", job);
	return -1;
    }
}


