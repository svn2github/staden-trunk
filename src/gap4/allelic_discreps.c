/*
 * FIXME
 *
 * o  Reversing a node does not necessarily reverse the cost as there are
 *    more than two possible base calls (ie it may still not match)
 *
 * o  new_node and new_edge are horribly inefficient. Asking for N nodes
 *    takes O(N^2) time. We need to deal with free nodes/edges better.
 *    Maybe by linking them to each other in a list and having a pointer
 *    to the first free item (or NULL if none).
 *
 * o  In seqs at position, filter out matches that aren't in the top two.
 *    Ie if AACAGCCA then remove the G. This is so that 'reverse' node works
 *    properly.
 *
 * o  To solve the previous problem in a better way, handle reverses by
 *    moving to another graph. That way we can potentially have more than
 *    2 graphs, indicating more than 2 populations.
 */

/*
 * IDEAS:
 *
 * Instead of "reversing" nodes, try taking nodes out of the graph into
 * a new graph. New graphs containing just that node can be created.
 *
 * (Implementing this is simply a matter of adding a graph
 * number to a node. Scoring a node only scores edges when both ends have
 * the same graph number.)
 */

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>

#include "IO.h"
#include "qual.h"
#include "gap_globals.h"
#include "array.h"
#include "misc.h"
#include "template.h"
#include "allelic_discreps.h"

typedef struct {
    int seqnum;
    int tmplate;
    char base;
    int1 conf;
    int cost;
} seq_base_t;

typedef struct {
    int pos;
    double score;
    seq_base_t *seqs;
    int nseqs;
} snp_t;


/* ------------------------------------------------------------------------ */
/* BASIC GRAPH PRIMITIVES */

typedef struct {
    /* Generic */
    Array edge; /* ad_edge* */
    double edge_sum;
    int index;

    /* Specific */
    int checked;
    int orient;
    int tnum;
} ad_node;

typedef struct {
    /* Generic */
    ad_node *n1;
    ad_node *n2;
    float weight;
    int index;

    /* Specific */
    snp_t *snp;
    int pos;
} ad_edge;

typedef struct {
    /* Generic */
    Array nodes; /* ad_node* */
    Array edges; /* ad_edge* */
    int free_node; /* 0 => index, -1 == unknown, -2 == none */
    int free_edge; /* 0 => index, -1 == unknown, -2 == none */
} ad_graph;


/* ------------------------------------------------------------------------ */

/*
 * Creates a new ad_edge and initialises all the fields to zero / NULL.
 * The edge should be freed using the del_edge() function.
 *
 * Returns ad_edge* on success
 *         NULL on failure.
 */
static ad_edge *new_edge(ad_graph *g) {
    ad_edge **e, **base;
    int sz, i;

    if (!g)
	return NULL;

    sz = ArrayMax(g->edges);
    base = ArrayBase(ad_edge *, g->edges);
    for (i = 0; i < sz; i++) {
	if (!base[i]) {
	    break;
	}
    }
    e = (ad_edge **)ArrayRef(g->edges, i);

    if (!e)
	return NULL;

    if (i == sz || !*e) {
	if (NULL == (*e = (ad_edge *)xmalloc(sizeof(ad_edge))))
	    return NULL;
    }

    memset(*e, 0, sizeof(**e));
    (*e)->n1 = NULL;
    (*e)->n2 = NULL;
    (*e)->weight = 0;
    (*e)->index = i;
    (*e)->snp = NULL;

    return *e;
}


/*
 * Deallocates an ad_edge 'e' from graph 'g'.
 * No attempt is made to check whether the nodes that this edge has
 * references to still maintain a reference to this edge.
 */
static void del_edge(ad_graph *g, ad_edge *e) {
    if (!e)
	return;

    arr(ad_edge *, g->edges, e->index) = NULL;
    xfree(e);
}


/* ------------------------------------------------------------------------ */

/*
 * Creates a new ad_node and initialises all the fields to zero / NULL.
 * The node should be freed using the del_node() function.
 *
 * Returns ad_node* on success
 *         NULL on failure.
 */
static ad_node *new_node(ad_graph *g) {
    ad_node **n, **base;
    int sz, i;

    if (!g)
	return NULL;

    sz = ArrayMax(g->nodes);
    base = ArrayBase(ad_node *, g->nodes);
    for (i = 0; i < sz; i++) {
	if (!base[i]) {
	    break;
	}
    }
    n = (ad_node **)ArrayRef(g->nodes, i);

    if (!n)
	return NULL;

    if (i == sz) {
	if (NULL == (*n = (ad_node *)xmalloc(sizeof(ad_node))))
	    return NULL;
    }

    memset(*n, 0, sizeof(**n));
    (*n)->edge = NULL;
    (*n)->index = i;

    return *n;
}

/*
 * Removes all the edges for a specific node, and if del_edge is true
 * the edge gets deallocated too if it no longer references a node.
 */
static void del_node_edges(ad_graph *g, ad_node *n, int del_edges) {
    int i, sz;

    if (!n || !n || !n->edge)
	return;

    sz = ArrayMax(n->edge);
    for (i = 0; i < sz; i++) {
	ad_edge *e = arr(ad_edge *, n->edge, i);

	/* Unlink edge from this node */
	if (e->n1 == n)
	    e->n1 = NULL;

	if (e->n2 == n)
	    e->n2 = NULL;

	/* Deallocate edge if it no longer references any nodes */
	if (!e->n1 && !e->n2 && del_edges) {
	    del_edge(g, e);
	}
    }

    ArrayDestroy(n->edge);
    n->edge = 0;
}

/* 
 * Deallocates an ad_node created via new_node().
 * If del_edges is true then any edges referenced by this node have their
 * link back to this node removed. If that means the edge has no link to a
 * valid node then the edge is also deallocated.
 */
static void del_node(ad_graph *g, ad_node *n, int del_edges) {
    if (!n || !g)
	return;

    if (n->edge) {
	del_node_edges(g, n, del_edges);
    }

    arr(ad_node *, g->nodes, n->index) = NULL;
    xfree(n);

    return;
}


/*
 * Returns a node pointer with a specific index.
 */
ad_node *get_node(ad_graph *g, int index) {
    return arr(ad_node *, g->nodes, index);
}


#if 0
/* DEBUG */
static void dump_node(ad_node *n) {
    float score;
    int j, l;
    ad_edge **e;

    printf("TEMPLATE %d\n", n->tnum);
    if (!n->edge)
	return;

    l = ArrayMax(n->edge);
    e = ArrayBase(ad_edge *, n->edge);
    score = 0;
    for (j = 0; j < l; j++) {
	printf("  Edge %d(%p) weight %f pos %d template #%d\n",
	       j, e[j], e[j]->weight, e[j]->pos,
	       e[j]->n1->tnum != n->tnum ? e[j]->n1->tnum : e[j]->n2->tnum);
	score += e[j]->weight;
    }
    printf("  Total weight %f\n", score);
}
#endif

/* ------------------------------------------------------------------------ */

/*
 * Creates and initialises an ad_graph structure.
 * This should be deallocated using the del_graph() function.
 *
 * Returns ad_graph pointer on success
 *         NULL on failure.
 */
static ad_graph *new_graph(void) {
    ad_graph *g = (ad_graph *)xcalloc(1, sizeof(*g));
    if (!g)
	return NULL;

    g->nodes = ArrayCreate(sizeof(ad_node *), 0);
    g->edges = ArrayCreate(sizeof(ad_node *), 0);
    g->free_node = -2;
    g->free_edge = -2;

    return g;
}


/*
 * Deallocates an ad_graph.
 * If del_nodes is true then the nodes and edges within the graph are also
 * deallocated.
 */
static void del_graph(ad_graph *g, int del_nodes) {
    if (!g)
	return;

    if (g->nodes) {
	if (del_nodes) {
	    int i, sz;

	    sz = ArrayMax(g->nodes);
	    for (i = 0; i < sz; i++) {
		del_node(g, arr(ad_node *, g->nodes, i), 1);
	    }
	}
	ArrayDestroy(g->nodes);
    }

    if (g->edges) {
	int i, sz;

	sz = ArrayMax(g->edges);
	for (i = 0; i < sz; i++) {
	    del_edge(g, arr(ad_edge *, g->edges, i));
	}
	ArrayDestroy(g->edges);
    }

    xfree(g);
}

/* ------------------------------------------------------------------------ */
/* GRAPH MANIPULATION FUNCTIONS */

/*
 * Links an edge between two nodes.
 */
void link_edge(ad_edge *e, ad_node *n1, ad_node *n2) {
    int n;

    e->n1 = n1;
    if (!n1->edge)
	n1->edge = ArrayCreate(sizeof(ad_edge *), 0);
    n = ArrayMax(n1->edge);
    ArrayRef(n1->edge, n);
    arr(ad_edge *, n1->edge, n) = e;
    /* ARR(ad_edge *, n1->edge, n) = e; */

    e->n2 = n2;
    if (!n2->edge)
	n2->edge = ArrayCreate(sizeof(ad_edge *), 0);
    n = ArrayMax(n2->edge);
    ArrayRef(n2->edge, n);
    arr(ad_edge *, n2->edge, n) = e;
    /* ARR(ad_edge *, n2->edge, n) = e; */
}


/*
 * Finds and returns an edge between two nodes.
 * If it doesn't exist it creates one first.
 *
 * Returns: ad_edge pointer on success
 *          NULL on failure (only if out of memory).
 */
ad_edge *get_edge(ad_graph *g, ad_node *n1, ad_node *n2,
		  int pos, snp_t *snp) {
    int i;
    ad_edge *e;

    /* Search through template with fewest edges */
    if (n1->edge && n2->edge) {
	Array edges = (ArrayMax(n1->edge) < ArrayMax(n2->edge))
	    ? n1->edge
	    : n2->edge;
	int sz = ArrayMax(edges);

	for (i = 0; i < sz; i++) {
	    e = arr(ad_edge *, edges, i);
	    if ((e->n1 == n1 && e->n2 == n2 && e->snp == snp) ||
		(e->n1 == n2 && e->n2 == n1 && e->snp == snp)) {
		return e;
	    }
	}
    }

    /* No edge between n1 and n2 at pos, so add one */
    e = new_edge(g);
    e->weight = 0;
    link_edge(e, n1, n2);
    e->pos = pos;
    e->snp = snp;

    return e;
}


/* ------------------------------------------------------------------------ */
/* SNP/HAPLOTYPE code from here on */

/*
 * Use calc_discrepancies() to identify candidate SNP locations within
 * a specific contig.
 *
 * Returns a malloced array of snp_t structures holding the SNPs.
 *         NULL on failure.
 */
snp_t *candidate_snps(GapIO *io, int contig, int *nsnps_p) {
    float *q1, *q2;
    int clen;
    int i;
    snp_t *snps = NULL;
    int nsnps = 0;
    int *costs;

    /* Calculate discrepancies matrix */
    clen = io_clength(io, contig);
    if (NULL == (q1 = (float *)xmalloc(clen * sizeof(float))))
	return NULL;

    if (NULL == (q2 = (float *)xmalloc(clen * sizeof(float)))) {
	xfree(q1);
	return NULL;
    }

    calc_discrepancies(contig, 1, clen, q1, q2,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)io);


    /* Search for high scoring SNPs and record position/score */
    for (i = 0; i < clen; i++) {
	float combined;

	/* combined = (q1[i] + log(MAX(.00001, 1-q2[i]/3))*-500)/2; */
	combined = MIN(300, q1[i] * q2[i]);
	if (combined > 40) {
	    /* Grow snps array in blocks of 32 at a time */
	    if (nsnps % 32 == 0) {
		snps = (snp_t *)xrealloc(snps, (nsnps + 32) * sizeof(snp_t));
		if (!snps)
		    return NULL;
	    }
	    snps[nsnps].pos   = i+1;
	    snps[nsnps].score = combined;
	    snps[nsnps].seqs  = NULL;
	    snps[nsnps].nseqs = 0;
	    nsnps++;
	}
    }

    xfree(q1);
    xfree(q2);

    /*
     * Reduce score of neighbouring snps. Each SNP is a triangular impulse
     * of 0, 1, 2, 3, 2, 1, 0 centred on the SNP location. We then divide
     * the SNP score by 3/costs[?] so that neighbouring SNPs get reduced.
     */
    if (NULL == (costs = (int *)xcalloc(clen+1, sizeof(int))))
	return NULL;

    for (i = 0; i < nsnps; i++) {
	int j;
	for (j = 0; j < 3; j++) {
	    if (snps[i].pos - j >= 1 && snps[i].pos - j <= clen)
		costs[snps[i].pos - j] += 3-j;
	    if (j != 0 && (snps[i].pos + j >= 1 && snps[i].pos + j <= clen))
		costs[snps[i].pos + j] += 3-j;
	}
    }

    for (i = 0; i < nsnps; i++) {
	snps[i].score /= costs[snps[i].pos] / 3.0;
    }
    xfree(costs);


    *nsnps_p = nsnps;
    return snps;
}


/*
 * Produce a list of all templates overlapping base 'pos' in 'contig'.
 * If a sequence at the specific position has a quality lower than the global
 * quality_cutoff then it is ignored. If a template is represented more than
 * once (ie by more than one sequence) then the sequences are checked for
 * matches and mismatches. Mismatching templates will be rejected if the
 * match is high. Matching sequences on a templates get coalesced into one
 * record.
 *
 * Returns an array of seq_base_t structs on success, of length *nseqs_p.
 *            NULL on failure
 */
seq_base_t *seqs_at_region(GapIO *io, int contig, template_c **tarr,
			   int pos, int *nseqs_p) {
    int rnum;
    seq_base_t *seqs = NULL;
    int nseqs = 0;
    GReadings r;
    int i, j;

    /*
     * Loop through all sequences in the contig. Could be optimised as
     * this gets called multiple times, but for now take the easy
     * approach.
     */
    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	int cost;
	int1 conf;
	char *seq;

	/* Read pos are sorted - optimise */
	if (io_relpos(io, rnum) > pos)
	    break;

	if (io_relpos(io, rnum) + ABS(io_length(io, rnum)) < pos)
	    continue;

	/* Read reading and sequence details */
	gel_read(io, rnum, r);
	if (r.confidence) {
	    int1 *c = (int1 *)xmalloc(r.length);
	    if (c) {
		DataRead(io, r.confidence, c, r.length, 1);
		conf = c[pos - r.position + r.start];
		xfree(c);
	    } else {
		conf = 0;
	    }
	} else {
	    conf = 0;
	}
	seq = SeqReadStatic(io, r.sequence, r.length);

	/* Filter out sequences of poor quality */
	if (conf < quality_cutoff)
	    continue;

#if 0
	/* FIXME: Filter out known bad templates. A test for whether a
	 * cyclic approach will work.
	 */
	if (/* Definitely chimeric */
	    r.template == 37 ||
	    r.template == 46 ||
	    r.template == 48 ||
	    r.template == 188||

	    /* Possibly chimeric, or poorly place PCR errors */
	    r.template == 68 ||
	    r.template == 263||
	    r.template == 279)
	    continue;
#endif

	/* Penalise for poor quality templates */
	cost = (r.template && tarr[r.template]->consistency) ? 4 : 1;

	/* Grow in blocks of 8 */
	if (nseqs % 8 == 0) {
	    seqs = (seq_base_t *)xrealloc(seqs,
					  (nseqs + 8) * sizeof(seq_base_t));
	    if (!seqs)
		return NULL;
	}

	seqs[nseqs].seqnum  = rnum;
	seqs[nseqs].tmplate = r.template;
	seqs[nseqs].cost    = cost;
	seqs[nseqs].conf    = conf;
	seqs[nseqs].base    = toupper(seq[pos - r.position + r.start]);
	nseqs++;
    }

    /*
     * We now have a list of sequences, but we want to remove any duplicate
     * templates. We mark sequences for removal by setting cost to 0.
     */
    for (i = 0; i < nseqs; i++) {
	int t = seqs[i].tmplate;

	for (j = i+1; j < nseqs; j++) {
	    if (seqs[j].tmplate == t) {
		if (seqs[i].base == seqs[j].base) {
		    if (seqs[i].conf > seqs[j].conf)
			seqs[j].cost = 0;
		    else
			seqs[i].cost = 0;
		} else {
		    /* Disagreement, both with conf > quality_cutoff. */
		    /* => reject this template */
		    seqs[i].cost = 0;
		    seqs[j].cost = 0;
		}
	    }
	}
    }

    /* Finally remove the templates with cost == 0 */
    for (i = j = 0; i < nseqs; i++) {
	if (seqs[i].cost) {
	    seqs[j++] = seqs[i];
	}
    }

    *nseqs_p = j;
    return seqs;
}


/*
 * Given a set of snps within a contig, this adds edges between all
 * pair-wise combinations of templates that have sequence cover each
 * SNP.
 */
void add_snp_edges(GapIO *io, int contig, template_c **tarr,
		   ad_graph *g, snp_t *snp, int nsnps) {
    int i, j, k;

    for (i = 0; i < nsnps; i++) {
	/* Identify sequences overlapping this SNP position */
	seq_base_t *seqs;
	int nseqs;
	double weight;

	/* Skip previously filtered SNPs */
	if (snp[i].pos <= 0)
	    continue;

	seqs = seqs_at_region(io, contig, tarr, snp[i].pos, &nseqs);
	snp[i].seqs  = seqs;
	snp[i].nseqs = nseqs;
	weight = snp[i].score;

	printf("SNP %d at %d qual_score %f. Seqs: ",
	       i, snp[i].pos, snp[i].score);

	/* Iterate through all pairs */
	for (j = 0; j < nseqs; j++) {
	    ad_node *n1, *n2;

	    printf("%c", seqs[j].base);
	    n1 = get_node(g, seqs[j].tmplate);

	    for (k = j+1; k < nseqs; k++) {
		ad_edge *e;

		/* Should not happen now, but double check */
		if (seqs[j].tmplate == seqs[k].tmplate) {
		    printf("Warning, duplicate template in graph\n");
		    continue;
		}

		n2 = get_node(g, seqs[k].tmplate);
		e = get_edge(g, n1, n2, snp[i].pos, &snp[i]);

		if (seqs[j].base == seqs[k].base) {
		    e->weight += weight / seqs[j].cost;
		} else {
		    e->weight -= weight / seqs[j].cost;
		}
	    }
	}
	puts("");
    }
}


/*
 * For all nodes, sum the weights on the edges to compute a node weight.
 * "node" is a node number to recalculate, or -1 for all.
 *
 * Returns the total score for all nodes or the specific node requested.
 */
int calculate_node_weights(ad_graph *g, int node) {
    int i, i_start, i_end;
    double total_score = 0;

    /* Single node or all nodes */
    i_start = (node == -1) ? 0 : node;
    i_end   = (node == -1) ? ArrayMax(g->nodes) : node+1;

    for (i = i_start; i < i_end; i++) {
	ad_node *n = arr(ad_node *, g->nodes, i);
	ad_edge **e;
	int j, sz;
	double score;

	if (!n->edge)
	    continue;

	e = ArrayBase(ad_edge *, n->edge);
	sz = ArrayMax(n->edge);

	/* printf("Node %d: Edges ", i); */
	for(score = j = 0; j < sz; j++) {
	    /* printf("%f ", e[j]->weight); */
	    score += e[j]->weight;
	}

	/* printf(" => sum %f\n", score); */
	n->edge_sum = score;
	total_score += score;
    }
    
    return total_score;
}


/*
 * Reverse a node by toggling its orient and reversing the sums on the edge.
 * FIXME: There may be more than two base types, so perhaps negating the
 * edge is not correct.
 */
void reverse_node(ad_node *node) {
    int i, n;
    ad_edge **e;
    int last_pos = -1;

    node->orient ^= 1;
    if (!node->edge)
	return;

    printf("  => Reversing template %d\n", node->tnum);

    n = ArrayMax(node->edge);
    e = ArrayBase(ad_edge *, node->edge);
    for (i = 0; i < n; i++) {
	e[i]->weight = -e[i]->weight;
	if (e[i]->pos != last_pos) {
	    last_pos = e[i]->pos;
	}
    }

    /* FIXME: we only need to recompute node weights linked to the edges
     * we change, but for now we recompute all!
     */
}


/* Qsort sort func for node scores */
int sort_nodes(const void *n1p, const void *n2p) {
    ad_node **n1 = (ad_node **)n1p;
    ad_node **n2 = (ad_node **)n2p;
    return (*n1)->edge_sum - (*n2)->edge_sum;
}


/*
 * Optimises a graph by 'reversing' nodes until we can no longer obtain an
 * improvement by doing so.
 * We avoid loops by disallowing any node to be reversed more than a fixed
 * number of times.
 */
void optimise_graph(ad_graph *g) {
    ad_node **sorted_nodes = NULL;
    int i;

    puts("=== Optimise graph ===");

    /* Build a table of node pointers so we can sort the pointers */
    sorted_nodes = (ad_node **)xcalloc(ArrayMax(g->nodes), sizeof(ad_node *));
    for (i = 0; i < ArrayMax(g->nodes); i++)
	sorted_nodes[i] = arr(ad_node *, g->nodes, i);

    /* Loop until we can loop no more! (ie no beneficial change) */
    for(;;) {
	/*
	 * At present we only attempt to reverse the first node.
	 * This is because we can guarantee that the first node has the best
	 * improvement in score by reversing it as the differencer is just
	 * score - -score.
	 *
	 * When reverse_node is improved to recalculate properly, rather than
	 * assuming that all mismatches turn into matches (which is not true
	 * when the alphabet is 4 instead of 2), then we can set rev_count
	 * to a higher figure.
	 */
	int rev_count = 1;
	int best_node = -1;
	int best_score = 0;
	int i;

	printf("> New round\n");

	/* Compute edge sum for each node and sort by this field */
	calculate_node_weights(g, -1);
	qsort(sorted_nodes, ArrayMax(g->nodes), sizeof(ad_node *),
	      sort_nodes);

	/*
	 * Check the impact on the score of reversing nodes. We try
	 * to find the first 'rev_count' nodes that give a +ve score
	 * difference, and then pick the best of those.
	 */
	for (i = 0; i < ArrayMax(g->nodes) && rev_count; i++) {
	    if (sorted_nodes[i]->checked) {
		int before, after;
		
		before = sorted_nodes[i]->edge_sum;
		if (before >= 0)
		    continue;

		reverse_node(sorted_nodes[i]);
		calculate_node_weights(g, sorted_nodes[i]->tnum);
		reverse_node(sorted_nodes[i]);
		after = sorted_nodes[i]->edge_sum;

		printf("    Index %d (tnum %d). %d -> %d\n",
		       i, sorted_nodes[i]->tnum, before, after);
		    
		if (after > before)
		    rev_count--;
		
		if (best_score < after-before) {
		    best_score = after-before;
		    best_node = i;
		}
	    }
	}

	/* If best_node == -1 then we've run out of improvements */
	if (-1 == best_node)
	    break;

	/* Reverse the best node */
	reverse_node(sorted_nodes[best_node]);
	sorted_nodes[best_node]->checked--;
    }

    xfree(sorted_nodes);
}


/*
 * Set the edge weight for any edges at position 'pos' to zero, effectively
 * removing their impact from the graph.
 */
static void remove_edges_at_pos(ad_graph *g, int pos) {
    int i, sz;
    ad_edge **base;

    sz = ArrayMax(g->edges);
    base = ArrayBase(ad_edge *, g->edges);
    for (i = 0; i < sz; i++) {
	if (base[i] && base[i]->pos == pos)
	    base[i]->weight = 0;
    }
}


/*
 * Filters SNPs based on the layout of the current graph. If a SNP
 * appears to contribute to too many negative edges then it is probably
 * a false SNP (assuming that the graph is layed out correctly).
 *
 * Filtering is made simply by setting the SNP position to zero.
 */
void filter_snps(GapIO *io, int contig, template_c **tarr,
		 ad_graph *g, snp_t *snp, int nsnps) {
    int i, j, k;
    int base_lookup[256];

    memset(base_lookup, 5, 256*sizeof(*base_lookup));
    base_lookup['A'] = 0;
    base_lookup['C'] = 1;
    base_lookup['G'] = 2;
    base_lookup['T'] = 3;
    base_lookup['*'] = 4;

    for (i = 0; i < nsnps; i++) {
	seq_base_t *seqs;
	int nseqs;
	int match, mismatch;

	if (!snp[i].pos)
	    continue;

	seqs = snp[i].seqs;
	nseqs = snp[i].nseqs;
	match = mismatch = 0;

	/* FIXME
	 * Alternative method is to loop through all edges in the graph
	 * building up + and - totals for each SNP position.
	 * SNPs with large sums of both (or maybe just large -ve sum) are
	 * most likely inaccurate.
	 *
	 * Similarly we can find nodes (templates) which seem to be suspect,
	 * or perhaps chimeric. Should we filter these out first?
	 */

	
	/* Count the number of matches/mismatches between all pairs */
	for (j = 0; j < nseqs; j++) {
	    ad_node *n1 = get_node(g, seqs[j].tmplate);
	    for (k = j+1; k < nseqs; k++) {
		ad_node *n2 = get_node(g, seqs[k].tmplate);

		/* FIXME: compute %age identity for all + and all -
		 * orient sequences, and see if both are high?
		 */

		/* same base & orient or opposite base & orient => match */
		if ((seqs[j].base == seqs[k].base) ==
		    (n1->orient   == n2->orient))
		    match++;
		else
		    mismatch++;
	    }
	}

	printf("SNP %d at %d, match/mis = %d/%d\n",
	       i, snp[i].pos, match, mismatch);

	if ((double)match/(mismatch+match) < 0.6 /* || base[0] == base[1] */) {
	    printf("   Filtering SNP\n");
	    remove_edges_at_pos(g, snp[i].pos);
	    snp[i].pos = 0;
	}
    }

    calculate_node_weights(g, -1);
}


/*
 * Displays a list of templates(nodes) in the graph and their relevant
 * scores.
 */
static void dump_templates(ad_graph *g) {
    int i, j, sz;
    ad_node **base;

    sz = ArrayMax(g->nodes);
    base = ArrayBase(ad_node *, g->nodes);
    for (i = 0; i < sz; i++) {
	int last_pos = -1;
	double score = 0;
	double total_score = 0;
	int ne;
	ad_edge **e;

	printf("Node %d\n", i);
	if (!base[i] || !base[i]->edge)
	    continue;

	e = ArrayBase(ad_edge *, base[i]->edge);
	ne = ArrayMax(base[i]->edge);
	for (j = 0; j < ne; j++) {
	    if (last_pos != e[j]->pos) {
		if (last_pos != -1 && score != 0)
		    printf("    SNP at %5d score %f\n", last_pos, score);
		total_score += score;
		last_pos = e[j]->pos;
		score = 0;
	    }
	    score += e[j]->weight;
	}

	if (last_pos != -1 && score != 0) {
	    printf("    SNP at %5d score %f\n", last_pos, score);
	    total_score += score;
	}

	printf("    Total %f\n", total_score);
    }
}


static void dump_snps(ad_graph *g, snp_t *snp, int nsnp) {
    int i, j, sz;
    ad_edge **base;

    sz = ArrayMax(g->edges);
    base = ArrayBase(ad_edge *, g->edges);
    for (j = 0; j < nsnp; j++) {
	int pos = snp[j].pos;
	double score;
	int count;

	score = 0;
	count = 0;
	for (i = 0; i < sz; i++) {
	    if (!base[i] || base[i]->pos != pos)
		continue;
	    
	    score += base[i]->weight;
	    count++;
	}
	printf("SNP %d pos %5d count %4d score %f avg %f\n",
	       j, snp[j].pos, count, score, score/count);

	if (score/count < 10) {
	    printf("SNP %d: REMOVED\n", j);
	    remove_edges_at_pos(g, snp[j].pos);
	    snp[j].pos = -snp[j].pos;
	}
    }

    calculate_node_weights(g, -1);
}


/*
 * Assign templates, and hence sequences to each of the separate alleles
 * based on node->orient.
 *
 * It returns a tcl-style list of two lists (one for each list), with
 * each sub list being a list of readings in the that allele. Not all
 * readings will make their way into an allele as some have "unknown"
 * placements.
 *
 * Returns dstring_t for success
 *         NULL for failure
 */
dstring_t *assign_alleles(GapIO *io, int conitg, ad_graph *g) {
    int rnum;
    dstring_t *ds0, *ds1, *ds;

    ds0 = dstring_create(NULL);
    ds1 = dstring_create(NULL);

    for (rnum = io_clnbr(io, conitg); rnum; rnum = io_rnbr(io, rnum)) {
	GReadings r;
	ad_node *n;

	gel_read(io, rnum, r);
	n = get_node(g, r.template);

	if (!n || n->edge_sum <= 0) {
	    printf("READ %s UNKNOWN allele\n",
		   io_rname(io, rnum));
	    continue;
	}

	if (n->orient == 0)
	    dstring_appendf(ds0, "%s ", io_rname(io, rnum));
	if (n->orient == 1)
	    dstring_appendf(ds1, "%s ", io_rname(io, rnum));

	printf("READ %s allele %d (score %f)\n",
	       io_rname(io, rnum), n->orient, n->edge_sum);
    }

    ds = dstring_create(NULL);
    dstring_insertf(ds, 0, "{%s} {%s}", dstring_str(ds0), dstring_str(ds1));
    dstring_destroy(ds0);
    dstring_destroy(ds1);

    return ds;
}


/*
 * Generates text output in dot format for use with graphviz/dotty.
 */
static void dump_graph(ad_graph *g) {
    int i, sz;
    ad_edge **base;

    fprintf(stderr, "graph ad_graph {\n");
    sz = ArrayMax(g->edges);
    base = ArrayBase(ad_edge *, g->edges);
    for (i = 0; i < sz; i++) {
	if (!base[i] || base[i]->pos <= 0)
	    continue;

	fprintf(stderr, "    p%p -- p%p\n",
		base[i]->n1, base[i]->n2);
    }
    fprintf(stderr, "}\n");
}

/* ------------------------------------------------------------------------ */
/*
 * MAIN ENTRY POINT.
 *
 * Given a contig containing a mixed assembly from two alleles, this
 * function attempts to split templates into two sets corresponding
 * to each of the two alleles.
 */
dstring_t *allelic_discreps(GapIO *io, int contig) {
    ad_graph *g = NULL;
    int ntemplates, i;
    template_c **tarr;
    snp_t *snps;
    int nsnps;
    int nchecks = 2;
    int cycle;
    dstring_t *ds = NULL;

    puts("Attach debugger now");
    sleep(1);


    /* Compute the discrepancies and obtain possible SNP locations & scores */
    snps = candidate_snps(io, contig, &nsnps);
    if (!nsnps)
	return NULL;

    /* Compute template consistency status; used for scoring SNPs */
    if (NULL == (tarr = init_template_checks(io, 1, &contig, 1)))
	goto error;

    check_all_templates(io, tarr);

    /*
     * Create a graph and add one node per template. We use node 0 too so
     * that node indicies match template indicies.
     */
    if (NULL == (g = new_graph()))
	goto error;

    ntemplates = Ntemplates(io);
    for (i = 0; i <= ntemplates; i++) {
	ad_node *n;
	n = new_node(g);
	n->tnum = i;
	n->checked = nchecks; /* number of remaining checks allowed */
    }

    /*
     * 2 cycles:
     * Add edges, optimise, find false looking SNPs and chimeric templates,
     * rebuild edges, and optimise again.
     */
    for (cycle = 0; cycle < 2; cycle++) {
	printf("=== CYCLE %d ===\n", cycle);
	
	/* Add edges between all nodes that reside at this snp. */
	add_snp_edges(io, contig, tarr, g, snps, nsnps);

	/* Optimise graph by reversing nodes */
	optimise_graph(g);

	/* Loop throgh SNPs once more determining probably false ones */
	dump_snps(g, snps, nsnps); /* Also does filtering */
	filter_snps(io, contig, tarr, g, snps, nsnps);

	dump_templates(g);

	/* 1 cycle only for now... FIXME :) */
	break;

	/* Clear edges from all nodes and reset score and check count */
	for (i = 0; i <= ntemplates; i++) {
	    ad_node *n = get_node(g, i);
	    del_node_edges(g, n, 1);
	    n->checked = nchecks;
	    n->orient = 0;
	}
    }

    /*
     * Spot potential chimeric clones.
     * If we've sequenced through the chimeric point then it'll be obvious
     * in the assembly and invariably it will have been clipped.
     * More likely though is that the chimera has joined in the middle
     * somewhere and we have each end sequenced.
     * We can search through each template looking for the SNPs that are on
     * it. If all the start SNPs match and the end ones mismatch, or vice
     * versa, and the split from matching to mismatching spans a section
     * of template for which we have no sequence data for, then label
     * this clone as chimeric. For extra checking, we may want a minimum of
     * 2 SNPs at each end to avoid simple PCR errors.
     */

    /* List SNPs */

    /* List templates/reads on each allele */
    ds = assign_alleles(io, contig, g);

    /* dump_graph(g); */

 error:
    if (snps)
	xfree(snps);

    if (g)
	del_graph(g, 1);

    if (tarr)
	uninit_template_checks(io, tarr);

    return ds;
}



/*
 * Tests on gap4 assembled CMT1A.2.
 * CONFLICTS are:
 *
 * PASS: CMT1APCR-01g02.q1k wrong allele (1). CMT1APCR-01g02.p1k is
 *       correct allele (2). => correct location (2 to 1 vote).
 *
 * PASS: CMT1APCR-04f05.q1k correct (1) wrong (1).
 *       CMT1APCR-04f05.p1k correct (2)
 *
 * FAIL: CMT1APCR-01d06.p1k wrong (3). CMT1APCR-01d06.q1k wrong (3).
 *
 * FAIL: CMT1APCR-01e05.p1k wrong (2). CMT1APCR-01e05.q1k correct(?) (2
 *       from a polyA).
 *
 * ????: CMT1APCR-01d05.q1k wrong (3). CMT1APCR-01d05.p1k correct (3).
 *
 * ???: CMT1APCR-03b09.q1k wrong (3). CMT1APCR-03b09.p1k correct (2).
 *
 * ???: CMT1APCR-01e03.p1k correct allele (2). CMT1APCR-01e03.q1k
 *      wrong (1 4-base deletion)
 *
 * ???: CMT1APCR-04d10.q1k. One match on each strand => PCR mutation and not
 *      possible to see where it belongs.
 *
 *
 * A test with CMT1A.P (phrap assembly) gives much the same ??? reads,
 * but only has one FAIL.
 *
 * Suspected chimeric clones: 01e10, 01d05, 03b09, 01g02?
 */


/* e05(t48), d06(t38) */
