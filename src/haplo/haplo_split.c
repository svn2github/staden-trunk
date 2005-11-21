/*
 * I've concluded that link_score is suspect at best. The idea that node A
 * and node B should be merged because A and B have a shared similarlity with 
 * node set X and a shared dissimilarlity with node set Y does not necessarily
 * mean that A and B are in the same set.
 *
 * If the edge score between A and B is bad it may mean that A or B are a
 * chimera or it simply be that there are 3 total sets to derive and not 2.
 *
 * For the haploid splitting case the link score idea may still be valid,
 * but for a repeat splitting problem we do not know the copy number.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "IO.h"
#include "haplo.h"
#include "xalloc.h"

struct _node;
struct _edge;

static int verbosity = 0;

/* For uninitialised edge and linkage scores */
#define UNDEF_SCORE -9999999
#define NO_LINK     -9999998

/* Array of Node pointers */
typedef struct {
    struct _node **node;
    int nnodes;
    int allocated;
} node_array;

/* Array of Edge pointers */
typedef struct {
    struct _edge **edge;
    int nedges;
    int allocated;
} edge_array;

/* A graph Node */
typedef struct _node {
    int number; /* -1 for unused */
    edge_array *edges;
    char *tname;
    double tscore;
    int (*matrix)[6];
    node_array *merged; /* order of merging */
    double chimeric_score;
} node;

/* A graph Edge */
typedef struct _edge {
    /* struct _edge **n1e */
    /* struct _edge **n2e */
    node *n1;
    node *n2;
    double edge_score;
    double linkage_score;
} edge;

/* The total graph */
typedef struct _graph {
    node_array *nodes;
    edge_array *edges;
    int (*matrix)[6];
    double *snp_scores;
    int nsnps;
    int ntemplates;
    double correlation_offset;
} graph;

#ifndef ABS
#    define ABS(a) ((a) > 0 ? (a) : -(a))
#endif

/* NODE_ARRAYS ------------------------------------------------------------- */

/**
 * Allocates and initialises a new node_array.
 */
node_array *node_array_create(void) {
    node_array *n = (node_array *)malloc(sizeof(*n));
    if (!n)
	return NULL;

    n->node = NULL;
    n->nnodes = 0;
    n->allocated = 0;

    return n;
}


/**
 * Deallocates memory used by a node_array.
 */
void node_array_destroy(node_array *n) {
    if (n) {
	if (n->node)
	    free(n->node);
	free(n);
    }
}

/**
 * Adds (and returns) a new node to a node_array, allocating memory as
 * appropriate. The node 'np' is added to the node array.
 *
 * Returns a pointer to the node pointer on success
 *         NULL on failure.
 */
node **node_array_add(node_array *n, node *np) {
    node **npp;

    if (n->nnodes >= n->allocated) {
	n->allocated = n->allocated ? n->allocated * 2 : 4;
	npp = (node **)realloc(n->node, n->allocated * sizeof(node **));
	if (NULL == npp)
	    return NULL;
	n->node = npp;
    }
    n->nnodes++;
    n->node[n->nnodes-1] = np;

    return &n->node[n->nnodes-1];
}


/**
 * Computes the intersection of a pair of node_arrays.
 * The returned value should be free using node_array_destroy
 *
 * Assumption: the nodes are sorted incrementally by node number.
 *
 * Returns an allocated node_array struct on success.
 *         NULL on failure.
 */
node_array *node_array_intersection(node_array *n1, node_array *n2) {
    int i, j;
    node_array *inter;

    if (NULL == (inter = node_array_create()))
	return NULL;

    for (i = j = 0; i < n1->nnodes; i++) {
	int nnum = n1->node[i]->number;
	while (j < n2->nnodes && n2->node[j]->number < nnum)
	    j++;
	if (j < n2->nnodes && n2->node[j]->number == nnum) {
	    if (NULL == node_array_add(inter, n1->node[i]))
		return NULL;
	}
    }

    return inter;
}


/**
 * Computes the union of a pair of node_arrays.
 * The returned value should be free using node_array_destroy
 *
 * Assumption: the nodes are sorted incrementally by node number.
 *
 * Returns an allocated node_array struct on success.
 *         NULL on failure.
 */
node_array *node_array_union(node_array *n1, node_array *n2) {
    int i, j;
    node_array *onion;

    if (NULL == (onion = node_array_create()))
	return NULL;

    i = j = 0;
    while (i < n1->nnodes && j < n2->nnodes) {
	if (n1->node[i]->number < n2->node[j]->number) {
	    while (i < n1->nnodes &&
		   n1->node[i]->number < n2->node[j]->number) {
		if (NULL == node_array_add(onion, n1->node[i]))
		    return NULL;
		i++;
	    }
	} else if (n1->node[i]->number > n2->node[j]->number) {
	    while (j < n2->nnodes &&
		   n2->node[j]->number < n1->node[i]->number) {
		if (NULL == node_array_add(onion, n2->node[j]))
		    return NULL;
		j++;
	    }
	} else {
	    if (NULL == node_array_add(onion, n1->node[i]))
		return NULL;
	    i++;
	    j++;
	}
    }

    while (i < n1->nnodes) {
	if (NULL == node_array_add(onion, n1->node[i]))
	    return NULL;
	i++;
    }

    while (j < n2->nnodes) {
	if (NULL == node_array_add(onion, n2->node[j]))
	    return NULL;
	j++;
    }

    return onion;
}

void node_array_print(char *msg, node_array *na) {
    int i;
    if (msg)
	printf("%s", msg);
    for (i = 0; i < na->nnodes; i++) {
	printf(" %d", na->node[i]->number);
    }
    printf("\n");
}


/* EDGE_ARRAYS ------------------------------------------------------------- */

/**
 * Allocates and initialises a new edge_array.
 */
edge_array *edge_array_create(void) {
    edge_array *e = (edge_array *)malloc(sizeof(*e));
    if (!e)
	return NULL;

    e->edge = NULL;
    e->nedges = 0;
    e->allocated = 0;

    return e;
}


/**
 * Deallocates memory used by a edge_array.
 */
void edge_array_destroy(edge_array *e) {
    if (e) {
	if (e->edge)
	    free(e->edge);
	free(e);
    }
}

/**
 * Adds (and returns) a new edge to a edge_array, allocating memory as
 * appropriate. The edge 'ep' is added to the edge array.
 *
 * Returns a pointer to the edge pointer on success
 *         NULL on failure.
 */
edge **edge_array_add(edge_array *e, edge *ep) {
    edge **epp;

    if (e->nedges >= e->allocated) {
	e->allocated = e->allocated ? e->allocated * 2 : 4;
	epp = (edge **)realloc(e->edge, e->allocated * sizeof(edge **));
	if (NULL == epp)
	    return NULL;
	e->edge = epp;
    }
    e->nedges++;
    e->edge[e->nedges-1] = ep;

    return &e->edge[e->nedges-1];
}


/* NODE -------------------------------------------------------------------- */

/**
 * Allocates an initialises a new node.
 */
node *node_create(void) {
    node *n = (node *)malloc(sizeof(*n));
    if (!n)
	return NULL;

    n->number = 0;
    n->edges = edge_array_create();
    n->merged = NULL;

    return n;
}

/**
 * Deallocates memory used by a node.
 */
void node_destroy(node *n) {
    if (n) {
	if (n->edges)
	    edge_array_destroy(n->edges);
	if (n->tname)
	    free(n->tname);
	if (n->merged)
	    node_array_destroy(n->merged);
	free(n);
    }
}


/**
 * qsort callback function for node_sort_edges()
 *
 * Hack - there's no "clientdata" to pass into qsort, so we use a
 * global variable to achieve the same result.
 */
static int qsort_node_num = -1;
static int qsort_edge(const void *v1, const void *v2) {
    edge **e1 = (edge **)v1;
    edge **e2 = (edge **)v2;
    int n1 = (*e1)->n1->number == qsort_node_num ?
	(*e1)->n2->number : (*e1)->n1->number;
    int n2 = (*e2)->n1->number == qsort_node_num ?
	(*e2)->n2->number : (*e2)->n1->number;

    return n1-n2;
}

/**
 * Sorts an edge_array by the node numbers linked to, in ascending order.
 */
void node_sort_edges(node *n) {
    qsort_node_num = n->number;
    qsort(n->edges->edge, n->edges->nedges, sizeof(edge *), qsort_edge);
}


/**
 * Returns a node_array of nodes connected to this node.
 * The node_array should be deallocated by the caller using
 * node_array_destroy().
 */
node_array *node_neighbours(node *n) {
    int i;
    node_array *na = node_array_create();

    for (i = 0; i < n->edges->nedges; i++) {
	edge *e = n->edges->edge[i];
	if (!e)
	    continue;
	node_array_add(na, e->n1 == n ? e->n2 : e->n1);
    }

    return na;
}

/* EDGE -------------------------------------------------------------------- */

/**
 * Allocates an initialises a new edge.
 */
edge *edge_create(void) {
    edge *e = (edge *)malloc(sizeof(*e));
    if (!e)
	return NULL;

    e->n1 = NULL;
    e->n2 = NULL;
    e->edge_score = UNDEF_SCORE;
    e->linkage_score = UNDEF_SCORE;

    return e;
}

/**
 * Deallocates memory used by a edge.
 */
void edge_destroy(edge *e) {
    if (e) {
	free(e);
    }
}

/* GRAPH ------------------------------------------------------------------- */

/**
 * Allocates and initialises a blank graph.
 */
graph *graph_create(void) {
    graph *g = (graph *)malloc(sizeof(*g));
    if (!g)
	return NULL;
    
    g->nodes = node_array_create();
    g->edges = edge_array_create();
    g->matrix = NULL;
    g->nsnps = 0;
    g->ntemplates = 0;
    g->correlation_offset = 0.9;

    return g;
}

/**
 * Recursively deallocates a node and all the nodes in its merged array.
 */
void node_recursive_destroy(node *n) {
    int i;
    for (i = 0; n->merged && i < n->merged->nnodes; i++) {
	node_recursive_destroy(n->merged->node[i]);
    }
    node_destroy(n);
}


/**
 * Deallocates a graph, including deallocating all of the nodes and edges
 * associated with the graph.
 */
void graph_destroy(graph *g) {
    if (!g)
	return;

    if (g->nodes) {
	int i;
	for (i = 0; i < g->nodes->nnodes; i++) {
	    if (g->nodes->node[i])
		node_recursive_destroy(g->nodes->node[i]);
	}
	node_array_destroy(g->nodes);
    }

    if (g->edges) {
	int i;
	for (i = 0; i < g->edges->nedges; i++)
	    edge_destroy(g->edges->edge[i]);
	edge_array_destroy(g->edges);
    }

    if (g->matrix)
	free(g->matrix);

    free(g);
}

/**
 * Creates a new node in the graph and returns the node pointer.
 */
node *graph_add_node(graph *g) {
    node *n = node_create();

    if (!n || NULL == node_array_add(g->nodes, n))
	return NULL;

    return n;
}

edge *graph_add_edge(graph *g, node *n1, node *n2, double edge_score) {
    edge *e = edge_create();

    if (!e || NULL == edge_array_add(g->edges, e))
	return NULL;

    e->n1 = n1;
    e->n2 = n2;
    e->edge_score = edge_score;
    e->linkage_score = UNDEF_SCORE;

    edge_array_add(n1->edges, e);
    edge_array_add(n2->edges, e);

    return e;
}


/**
 * Doesn't really work as a good measure as it has lots of false
 * positives, but at least the real chimeric clones also get a low
 * score so we will still do most of the good clones first when
 * merging.
 *
 * Computes the total score for the edge divided by the total of the
 * absolute scores for this edge. The idea being that mostly +ve or
 * -ve scoring edges are fine, but a mixture of the two implies maybe
 * one node is chimeric or simply "poor".
 */
double chimeric_score(graph *g, edge *e) {
    int (*M1)[6] = e->n1->matrix;
    int (*M2)[6] = e->n2->matrix;
    int i, j, k, score = 0, total = 0;
    int count = 0;
    double dscore;

    for (i = 0; i < g->nsnps; i++) {
	for (j = 1; j < 6; j++) {
	    for (k = 1; k < 6; k++) {
		if (M1[i][j] && M2[i][k]) {
		    if (j == k)
			score += g->snp_scores[i];
		    else
			score -= g->snp_scores[i];
		    total += g->snp_scores[i];
		    count++;
		}
	    }
	}
    }

    if (score < 0) score = -score;
    dscore = (score + 500) / (double)(total + 500);
    return dscore * dscore;
}

void graph_calc_chimeric_scores(graph *g) {
    double *score;
    double *minscore;
    int *count;
    int i;

    if (verbosity)
	puts("Chimera checking; low scores *may* indicate chimeras");

    score = (double *)malloc(g->nodes->nnodes * sizeof(*score));
    minscore = (double *)malloc(g->nodes->nnodes * sizeof(*score));
    count = (int *)malloc(g->nodes->nnodes * sizeof(*score));
    for (i = 0; i < g->nodes->nnodes; i++) {
	minscore[i] = 1;
	score[i] = 0;
	count[i] = 0;
    }

    for (i = 0; i < g->edges->nedges; i++) {
	edge *e = g->edges->edge[i];
	double sc;

	if (NULL == e) 
	    continue;

	sc = chimeric_score(g, e);
	if (minscore[e->n1->number] > sc)
	    minscore[e->n1->number] = sc;
	if (minscore[e->n2->number] > sc)
	    minscore[e->n2->number] = sc;
	score[e->n1->number] += sc;
	score[e->n2->number] += sc;
	count[e->n1->number]++;
	count[e->n2->number]++;
    }

    for (i = 0; i < g->nodes->nnodes; i++) {
	g->nodes->node[i]->chimeric_score =
	    minscore[i] * (score[i]+5)/(count[i]+5);
	if (verbosity >= 2)
	    printf("CHIMERIC %f %s\n",
		   g->nodes->node[i]->chimeric_score,
		   g->nodes->node[i]->tname);
    }
    
    free(score);
    free(minscore);
    free(count);
}

void graph_print(graph *g, int verbose) {
    int i, j;
    for (i = 0; i < g->nodes->nnodes; i++) {
	node *n1 = g->nodes->node[i];
	node *n2;

	if (!n1)
	    continue;

	printf("Node %d :", n1->number);
	for (j = 0; j < n1->edges->nedges; j++) {
	    if (NULL == n1->edges->edge[j])
		continue;

	    n2 = n1->edges->edge[j]->n1 == n1
		? n1->edges->edge[j]->n2
		: n1->edges->edge[j]->n1;
	    if (verbose)
		printf(" (%d=%+3f,%+3f)",
		       n2->number,
		       n1->edges->edge[j]->edge_score,
		       n1->edges->edge[j]->linkage_score);
	    else
		printf(" %d/%d", n2->number, (int)(n1->edges->edge[j]->edge_score/100));
	}
	puts("");
    }
}

/* File I/O ---------------------------------------------------------------- */

/**
 * Find an edge between node n1 and n2.
 * Returns NULL if no edge is found.
 */
edge *edge_find(node *n1, node *n2) {
    int i, nedges;
    edge **edges;

    /* Search through the node with smallest number of edges */
    if (n1->edges->nedges > n2->edges->nedges) {
	nedges = n2->edges->nedges;
	edges = n2->edges->edge;
    } else {
	nedges = n1->edges->nedges;
	edges = n1->edges->edge;
    }

    for (i = 0; i < nedges; i++) {
	if (!edges[i])
	    continue;

	if ((edges[i]->n1 == n1 && edges[i]->n2 == n2) ||
	    (edges[i]->n1 == n2 && edges[i]->n2 == n1)) {
	    return edges[i];
	}
    }

    return NULL;
}


/**
 * Unlinks an edge from the nodes that use it.
 */
void edge_unlink(edge *e) {
    int i, j;

    for (j = 0; j < 2; j++) {
	node *n = j ? e->n1 : e->n2;

	for (i = 0; i < n->edges->nedges; i++) {
	    if (n->edges->edge[i] == e) {
		n->edges->edge[i] = NULL;
		break;
	    }
	}
    }

    e->n1 = NULL;
    e->n2 = NULL;
    e->edge_score = NO_LINK;
    e->linkage_score = NO_LINK;
}


/**
 * Computes a "link" score between nodes n1 and n2. If recalc is zero
 * we attempt to use the cached copy, if known. Otherwise we force a
 * recalculation of this figure.
 *
 * Let E(n1,n2) = edge score between n1 and n2.
 * If S(n) is the set of edges connecting node n.
 * Let I(n1,n2) = Intersection of S(n1) and S(n2).
 * "|A|" = Absolute value of A.
 *                          __
 * Link(n1,n2) = E(n1,n2) + \   ( |(E(n1,X)+E(n1,X)| - |(E(n1,X)-E(n2,X)| )
 *                          /_ 
 *                       X member of I(n1,n2)
 *
 * Returns the link_score if n1 and n2 are directly connected,
 *         or NO_LINK if n1 and n2 are not connected.
 */
int link_score(node *n1, node *n2, int recalc) {
    double score = 0;
    edge *e, *e1, *e2;
    node_array *na1, *na2, *na;
    int i;

    if (NULL == (e = edge_find(n1, n2)))
	return NO_LINK;

    if (!recalc && e->linkage_score != UNDEF_SCORE)
	return e->linkage_score;

    /* Obtain the array of nodes linked to by both n1 and n2 */
    na1 = node_neighbours(n1);
    na2 = node_neighbours(n2);
    na = node_array_intersection(na1, na2);

    /* Compute score */
    score = e->edge_score;

#if 1
    if (score >= 0) {
	for (i = 0; i < na->nnodes; i++) {
	    e1 = edge_find(n1, na->node[i]);
	    e2 = edge_find(n2, na->node[i]);
	    score += ABS(e1->edge_score + e2->edge_score) / 100;
	    score -= ABS(e1->edge_score - e2->edge_score) / 100;
	}
    }
#endif

    node_array_destroy(na);
    node_array_destroy(na1);
    node_array_destroy(na2);

    e->linkage_score = score * n1->chimeric_score * n2->chimeric_score *
	n1->tscore * n2->tscore;
    return score;
}

void graph_calc_link_scores(graph *g, int recalc) {
    int i, j;
    node *n1, *n2;
    node_array *na;

    for (i = 0; i < g->nodes->nnodes; i++) {
	n1 = g->nodes->node[i];

	if (!n1)
	    continue;

	na = node_neighbours(n1);

	for (j = 0; j < na->nnodes; j++) {
	    n2 = na->node[j];

	    if (n2->number < n1->number)
		continue; /* will do this from the other node */

	    link_score(n1, n2, recalc);
	}

	node_array_destroy(na);
    }
}

/**
 * Finds the edge with the highest linkage score and returns it.
 */
edge *best_edge(graph *g) {
    int i, best_score = -1e6;
    edge *best_edge = NULL;

    for (i = 0; i < g->edges->nedges; i++) {
	edge *e = g->edges->edge[i];
	if (NULL == e) 
	    continue;

	if (e->linkage_score == UNDEF_SCORE)
	    link_score(e->n1, e->n2, 0 /* cached */);

	if (best_score < e->linkage_score) {
	    best_score = e->linkage_score;
	    best_edge = e;
	}
    }

    return best_edge;
}

/**
 * Computes the edge score between matrix M1 and M2. These matrices
 * have "not present", A, C, G, T, pad as the elements 0-5
 * inclusive. We therefore loop from 1-5 inclusive and skip the
 * not-present case when computing the correlation between the
 * matrices.
 *
 * The correlation scores for each SNP positions are them summed
 * together to produce the overall edge score.
 */
double calc_edge_score(int (*M1)[6], int (*M2)[6], double *scores,
		    int nsnps, int *countp, double offset) {
    int i, j, count = 0;
    double score = 0;

    for (i = 0; i < nsnps; i++) {
	double avg1, avg2;

	/* Compute correlation between M1[i] to M2[i] */
	avg1 = (M1[i][1] + M1[i][2] + M1[i][3] + M1[i][4] + M1[i][5]) / 5.0;
	avg2 = (M2[i][1] + M2[i][2] + M2[i][3] + M2[i][4] + M2[i][5]) / 5.0;
	double numerator = 0;
	double denom1 = 0, denom2 = 0;
	
	for (j = 1; j < 6; j++) {
	    double diff1 = M1[i][j] - avg1;
	    double diff2 = M2[i][j] - avg2;
	    numerator += diff1 * diff2;
	    denom1 += diff1 * diff1;
	    denom2 += diff2 * diff2;
	}
	if (denom1 * denom2) {
	    count++;
	    double corr = numerator / sqrt(denom1 * denom2) - offset;
	    /*
	    printf("R({%d,%d,%d,%d,%d}, {%d,%d,%d,%d,%d}) = %f\n",
		   M1[i][1], M1[i][2], M1[i][3], M1[i][4], M1[i][5], 
		   M2[i][1], M2[i][2], M2[i][3], M2[i][4], M2[i][5], 
		   numerator / sqrt(denom1 * denom2));
	    */
	    score += 100 * corr * scores[i];
	}
    }

    if (countp)
	*countp = count;

    return score;
}



/**
 * Merges the nodes linked to by edge 'e' and sets any neighbouring
 * edge / link scores to be invalid.
 */
void merge_node(graph *g, edge *e) {
    node *n1, *n2;
    node_array *na1, *na2, *na;
    int i, j;

    if (verbosity >= 2)
	printf("Merging %d / %d (score %8.2f, link %8.2f)   %s / %s\n",
	       e->n1->number, e->n2->number,
	       e->edge_score,
	       e->linkage_score,
	       e->n1->tname, e->n2->tname);

    /*
      print_matrix_node(g, e->n1);
      print_matrix_node(g, e->n2);
    */

    /* Find all the nodes linked to either n1 or n2 where n1 and n2
     * are the nodes in this edge.
     */
    n1 = e->n1;
    n2 = e->n2;
    na1 = node_neighbours(n1);
    na2 = node_neighbours(n2);
    na = node_array_union(na1, na2);
    node_array_destroy(na1);
    node_array_destroy(na2);

    /* Attach n2 to the node_array in n1 - allows traceback */
    if (!n1->merged) {
	n1->merged = node_array_create();
    }
    node_array_add(n1->merged, n2);

    /*
     * Merge the matrix rows.
     */
    for (i = 0; i < g->nsnps; i++) {
	for (j = 0; j < 6; j++)
	    n1->matrix[i][j] += n2->matrix[i][j];
    }

    /*
     * Forall nodes in our set, find the edges between that and n1
     * and/or n2. If it links with both then set the edge score to be
     * the average of the two values, otherwise use the single score.
     * Set the linkage score to be undefined.
     * Reset the edge to be between this node and n1 always (as n2
     * will then be disconnected and considered to be merged with n1).
     */
    for (i = 0; i < na->nnodes; i++) {
	node *n = na->node[i];
	edge *e1, *e2;

	if (n == n1 || n == n2)
	    continue;

	e1 = edge_find(n, n1);
	e2 = edge_find(n, n2);

	if (!e1 && !e2)
	    continue;

	if (e1 && e2) {
	    /* links to both, so remove edge to n2 */
	    e1->edge_score = (e1->edge_score + e2->edge_score) / 2;
	    edge_unlink(e2);
	} else if (e2) {
	    /* links only to n2, so relink edge to n1 */
	    if (e2->n1 == n)
		e2->n2 = n1;
	    else
		e2->n1 = n1;

	    edge_array_add(n1->edges, e2);

	    e1 = e2;
	}

	e1->linkage_score = UNDEF_SCORE;
	e1->edge_score = UNDEF_SCORE;
    }
    node_array_destroy(na);

    edge_unlink(e);

    for (i = 0; i < g->nodes->nnodes; i++) {
	if (g->nodes->node[i] == n2) {
	    g->nodes->node[i] = NULL;
	    break;
	}
    }

    /* Recompute all undefined edge scores */
    for (i = 0; i < g->edges->nedges; i++) {
	edge *e;

	if (!(e = g->edges->edge[i]))
	    continue;

	if (!e->n1 || !e->n2)
	    continue;

	e->edge_score =
	    calc_edge_score(e->n1->matrix, e->n2->matrix, g->snp_scores,
			    g->nsnps, NULL, g->correlation_offset);
    }
}

static void print_matrix_node(graph *g, node *n) {
    int j, k;
    printf("%s :\n", n->tname);
    for (k = 1; k < 6; k++) {
	int line = 0;
	for (j = 0; j < g->nsnps; j++) {
	    if (n->matrix[j][k]) {
		line = 1;
		break;
	    }
	}
	if (line || 1) {
	    printf("Seq %d:%c ", n->number, "-ACGT*"[k]);
	    for (j = 0; j < g->nsnps; j++) {
		printf("%c", n->matrix[j][k] + '0');
	    }
	    puts("");
	}
    }
}

static void print_matrix(graph *g) {
    int i;
    node *n;

    puts("===Matrix===");
    for (i = 0; i < g->nodes->nnodes; i++) {
	if (!(n = g->nodes->node[i]))
	    continue;

	printf("%d ", i);
	print_matrix_node(g, n);
    }
}

static void print_group_recurse(node *n, int depth) {
    int i;

    if (!n->merged)
	return;

    for (i = 0; i < n->merged->nnodes; i++) {
	node *sub = n->merged->node[i];
	printf("%.*s%d %s\n",
	       depth,
	       "                                                            ",
	       sub->number, sub->tname);
	print_group_recurse(sub, depth+1);
    }
}

void print_groups(graph *g) {
    int i, ngroups = 0;
    node *n;

    puts("++groups");
    for (i = 0; i < g->nodes->nnodes; i++) {
	if (!(n = g->nodes->node[i]))
	    continue;

	printf("Group %d\n", ngroups++);
	printf(">%d %s\n", n->number, n->tname);
	print_group_recurse(n, 2);
    }
    puts("--groups");
}


/**
 * Adds edges between every remaining node and gives these edges a
 * zero-value edge score.
 *
 * The purpose of this is to compute the link score between nodes with
 * no direct edge. This allows for inference of which nodes go
 * together based on their correlation of disagreement with other
 * nodes rather than agreement, which in turn means the 2-haplotype
 * problem can be solved.
 *
 * This is probably not a desirable function to use when dealing with
 * repeat splitting though as in that case the number of repeat copies
 * may well not be exactly two.
 */
void add_zero_edges(graph *g) {
    int i, j;
    node *n1, *n2;
    edge *e;

    for (i = 0; i < g->nodes->nnodes; i++) {
	if (!(n1 = g->nodes->node[i]))
	    continue;

	for (j = i+1; j < g->nodes->nnodes; j++) {
	    if (!(n2 = g->nodes->node[j]))
		continue;

	    if ((e = edge_find(n1, n2)))
		continue;

	    /* No edge between i and j, so make one */
	    graph_add_edge(g, n1, n2, 0);
	}
    }
}

/* Qsort comparison function */
static int int_compare(const void *vp1, const void *vp2) {
    return (*(const int *)vp1) - (*(const int *)vp2);
}

/**
 * Builds a snp graph from the snps array
 */
graph *graph_from_snps(GapIO *io, snp_t *snp, int nsnps, double c_offset) {
    graph *g = NULL;
    int i, j;
    int ntemplates, *templates = NULL;
    node **tmpl_map = NULL;
    int lookup[256];

    if (verbosity)
	puts("Building graph");

    /* lookup table for fast base to index conversion */
    for (i = 0; i < 256; i++)
	lookup[i] = 0;
    lookup['A'] = lookup['a'] = 1;
    lookup['C'] = lookup['c'] = 2;
    lookup['G'] = lookup['g'] = 3;
    lookup['T'] = lookup['t'] = 4;
    lookup['*'] = 5;

    g = graph_create();
    g->correlation_offset = c_offset;

    /*
     * Obtain a list of unique template numbers, by collecting and sorting.
     * For each unique template we initialise the graph node.
     *
     * At the end of this we have ntemplates unique templates and
     * tmpl_map[] indexed by Gap4 template number maps us to the
     * node in the graph.
     */
    for (ntemplates = i = 0; i < nsnps; i++) {
	ntemplates += snp[i].nseqs;
    }
    templates = (int *)malloc(ntemplates * sizeof(int));
    for (ntemplates = i = 0; i < nsnps; i++) {
	for (j = 0; j < snp[i].nseqs; j++) {
	    templates[ntemplates++] = snp[i].seqs[j].tmplate;
	}
    }
    tmpl_map = (node **)xcalloc(Ntemplates(io)+1, sizeof(*tmpl_map));
    qsort(templates, ntemplates, sizeof(int), int_compare);
    for (i = j = 0; i < ntemplates; j++) {
	tmpl_map[templates[i]] = graph_add_node(g);
	tmpl_map[templates[i]]->number = j;
	tmpl_map[templates[i]]->tname  =
	    strdup(get_template_name(io, templates[i]));

	do {
	    i++;
	} while (i < ntemplates && templates[i] == templates[i-1]);
    }
    xfree(templates);
    ntemplates = j;

    g->nsnps = nsnps;
    g->ntemplates = ntemplates;


    /*
     * The matrix is of size nsnps by ntemplates. It initially will consist
     * of entirely empty vectors.
     */
    g->matrix = (int (*)[6])malloc(ntemplates * nsnps * sizeof(*g->matrix));
    memset(&g->matrix[0][0], 0, ntemplates * nsnps * 6 * sizeof(int));
    for (i = j = 0; j < ntemplates; i++) {
	if (!tmpl_map[i])
	    continue;

	tmpl_map[i]->matrix = &g->matrix[j * nsnps];
	j++;
    }


    /* Create the snp_score array */
    g->snp_scores = (double *)malloc(nsnps * sizeof(*g->snp_scores));
    for (i = 0; i < nsnps; i++) {
	g->snp_scores[i] = snp[i].score;
    }

    /* Now add basecalls to the matrix from snp info */
    for (i = 0; i < nsnps; i++) {
	for (j = 0; j < snp[i].nseqs; j++) {
	    node *n = tmpl_map[snp[i].seqs[j].tmplate];
	    n->matrix[i][lookup[snp[i].seqs[j].base]]++; 
	    n->tscore = snp[i].seqs[j].tscore; 
	    /* FIXME: confidence not yet used */
	}

	for (j = 0; j < ntemplates; j++) {
	    if (!g->nodes->node[j]->matrix[i][1] &&
		!g->nodes->node[j]->matrix[i][2] &&
		!g->nodes->node[j]->matrix[i][3] &&
		!g->nodes->node[j]->matrix[i][4] &&
		!g->nodes->node[j]->matrix[i][5])
 		g->nodes->node[j]->matrix[i][0]++;
	}
    }

    return g;
}

void graph_add_edges(graph *g) {
    node *n1, *n2;
    int i, j;

    for (i = 0; i < g->ntemplates; i++) {
	n1 = g->nodes->node[i];
	for (j = i+1; j < g->ntemplates; j++) {
	    int count;
	    double score;
	    n2 = g->nodes->node[j];

	    score = calc_edge_score(n1->matrix, n2->matrix,
				    g->snp_scores, g->nsnps, &count,
				    g->correlation_offset);
	    if (count) {
		graph_add_edge(g, g->nodes->node[i], g->nodes->node[j], score);
	    }
	}
    }

    /* Sort edges-out for each node by node number */
    for (i = 0; i < g->ntemplates; i++) {
	node_sort_edges(g->nodes->node[i]);
    }
}

static void list_group_recurse(dstring_t *ds, node *n) {
    int i;

    if (!n->merged)
	return;

    for (i = 0; i < n->merged->nnodes; i++) {
	node *sub = n->merged->node[i];
	dstring_appendf(ds, " %s", sub->tname);
	list_group_recurse(ds, sub);
    }
}

static dstring_t *list_groups(graph *g) {
    int i;
    node *n;
    dstring_t *ds;

    ds = dstring_create(NULL);

    /* We follow this by each group of nodes */
    for (i = 0; i < g->nodes->nnodes; i++) {
	if (!(n = g->nodes->node[i]))
	    continue;

	dstring_appendf(ds, "{%s", n->tname);
	list_group_recurse(ds, n);
	dstring_appendf(ds, "} ");
    }

    return ds;
}

int count_groups(graph *g) {
    int i;
    int ngroups = 0;

    /* We follow this by each group of nodes */
    for (i = 0; i < g->nodes->nnodes; i++) {
	if (!g->nodes->node[i])
	    continue;

	ngroups++;
    }

    return ngroups;
}

dstring_t *haplo_split(GapIO *io, snp_t *snp, int nsnps, int verbose,
		       double min_score, int two_pass, int fast_mode,
		       double c_offset, int max_sets) {
    graph *g;
    edge *e;
    dstring_t *ds;

    verbosity = verbose;
    g = graph_from_snps(io, snp, nsnps, c_offset);
    if (verbosity >= 3)
	print_matrix(g);

    graph_add_edges(g);
    graph_calc_chimeric_scores(g);
    graph_calc_link_scores(g, 1);
    if (verbosity >= 3)
	graph_print(g, 0);

    if (verbosity)
	puts("Merging graph nodes");

    while ((e = best_edge(g)) && (e->linkage_score > min_score)) {
	if (verbosity >= 1) {
	    putchar('.');
	    fflush(stdout);
	}
	merge_node(g, e);
	graph_calc_link_scores(g, fast_mode ? 0 : 1);
	if (verbosity >= 4) {
	    print_matrix(g);
	    graph_print(g, 1);
	}
    }
    if (verbosity >= 1)
	puts("");

    /* graph_print(g, 1); */

    if (two_pass) {
	/* Add fake zero-score edges if we want just 2-haplotypes */
	add_zero_edges(g);
	graph_calc_link_scores(g, 1);
	if (verbosity >= 4)
	    graph_print(g, 1);

	puts("===pass 2===");
	while ((e = best_edge(g)) && (e->linkage_score > min_score)) {
	    merge_node(g, e);
	    graph_calc_link_scores(g, fast_mode ? 0 : 1);
	    /* graph_print(g, 1); */
	}
	/* graph_print(g, 1); */
    }

    /* Force number of groups to be X? */
    if (max_sets) {
	int ngroups = count_groups(g);
	add_zero_edges(g);
	for (; ngroups > max_sets; ngroups--) {
	    e = best_edge(g);
	    if (!e) {
		printf("Bailed out as no edge connecting groups\n");
		break;
	    }
	    merge_node(g, e);
	    graph_calc_link_scores(g, fast_mode ? 0 : 1);
	}
    }

    /* print_groups(g); */

    ds = list_groups(g);
    graph_destroy(g);

    return ds;
}
