#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>
#include "QueueInterface.h"
#include "IO.h"
#include "template.h"
#include "io-reg.h"
#include "misc.h"

/* contig ordering program 

   initially for "read-pairs" so the input is a list of contigs
   and the id's of the contigs that they are linked to by read-pairs.
   For each contig, the contigs it is linked to are termed "mates".

   1. the input data is converted to the "graph" structure. The graph
      is described by adjacency records but within the same structure
      we later build trees (called Spanning Paths, or SP's) using linked
      lists. An SP is a list of vertices that cover the graph from one
      end to the other, but does not include all the vertices. Having
      obtained an SP we try to include other vertices that its members
      overlap.
   2. fill in the addresses of each vertex's mates
   3. from a randomly selected vertex ([0]) create a tree and return its
      final vertex
   4. from the final vertex of the previous tree create another tree. This
      is our SP.
   5. Our SP (defined as a singly linked list) only has left links. Add
      the right links to make it doubly linked
   6. Now we try to add in as many of the unused vertices (contigs) as
      possible. First find the vertices that are not part of the spanning 
      path and mark them unvisited. Next we set the weights for those that 
      are visited (their weight is the same as their position, ie the first
      has weight 1, the next 2, etc. Then we calculate weights for the 
      unvisited vertices dependent on the visited vertices they overlap.
      eg if an unvisited vertex overlaps the vertices with weights 3 and 4,
      then it will be given weight 3.5. We want to intercalate these vertices
      so we sort them into left to right order and add them to the SP. From
      here we return one end of the SP.

   7. The contigs in the db should now be ordered using the SP. (Not done).

   NOTES

     I have not yet written the freeing routines for these data structures.
     What can go wrong or has not been programmed: In general the routine
     will be used when there is not enough read-pair data to order all the
     contigs, so really it should return several SP's. Some read-pair data
     may be false (eg because data have been misnamed) and there is no
     checking for this. The algorithm cannot be certain about the relative
     order of the two ends of the SP - ie the last two contigs could be in
     the wrong order, ditto the first two. The actual lengths of the templates
     (and contigs) has not been used, so if there are a few very long templates
     mixed in with much smaller ones I am not sure of the outcome.
     The standard gap4 error routines are needed in a few places marked FIXME.

*/


void InitializeQueue(Queue *Q)
{                  
    Q->Count = 0;           /* Count == number of items in the queue */
    Q->Front = 0;           /* Front == location of item to remove next */
    Q->Rear = 0;            /* Rear == place to insert next item */
    Q->Items = ArrayCreate(sizeof(ItemType), 0); /* Queued items */
}

void DestroyQueue(Queue *Q) {
    ArrayDestroy(Q->Items);
}

int Empty(Queue *Q)
{
    return ( Q->Count == 0 );
}


int Full(Queue *Q)
{
    return 0;
}


void Insert(ItemType R, Queue *Q)

{
    ArrayRef(Q->Items, Q->Rear);
    arr(ItemType, Q->Items, Q->Rear) = R;
    Q->Rear++;
    ++(Q->Count);
}

void Remove(Queue *Q, ItemType *F)
{
    if (Q->Count == 0) {
	/* FIXME */
	printf("attempt to remove item from empty Queue\n");
    } else {
	*F = arr(ItemType, Q->Items, Q->Front);
	Q->Front++;
	--(Q->Count);
    }
}


typedef enum {false, true} Boolean;

/* Search graph G beginning at vertex v */
void 
GraphSearch(Graph *g, 
	    int cur_sp, 
	    AdjacencyRec *vin, 
	    AdjacencyRec **vout)  
{
    Queue Q, *q;
    int i;
    AdjacencyRec *a = NULL,*x,*w, *p;

#ifdef DEBUG
    printf("GraphSearch\n");
#endif

    q = &Q;

    (void) InitializeQueue(&Q);

/*
    printf("front %d rear %d count %d\n",q->Front,q->Rear,q->Count);
    printf("degree %d\n",vin->degree);
*/

	/* mark each vertex x in V as being unvisited */

    for (i = 0; i < g->number_of_verts; i++) {
	/* second call to GraphSearch */
	if (cur_sp == g->recs[i]->visited) {
	    g->recs[i]->visited = false;
	    g->recs[i]->left = NULL;
	    g->recs[i]->right = NULL;
	} else {
	    /* first call to GraphSearch of new SP */
	    if (!g->recs[i]->visited) {
		g->recs[i]->left = NULL;
		g->recs[i]->right = NULL;
	    }    
	}
    }

   /* Use vertex v in V as a starting point, and put v in container C */
    (void) Insert(vin, &Q );

    p = NULL;

    while (! Empty(&Q) ) {
	
	/*   (Remove a vertex x from container C); */
	(void) Remove( &Q, &x);
	
	if (!(x->visited)) {     /* if vertex hasn't been visited already */
	    
	    x->visited = cur_sp;            /* mark x as having been visited */
	    if ( ! ( x->left ) ) x->left = p;
	    p = x;
	    a = x;	                  /* remember last vertex visited */
	    
	    /*         
	     * for (each vertex w in Vx) {    
	     * Enter all unvisited vertices   
	     * if (!(w.Visited))  (Put w into C);         
	     */
	    /* printf("x %d degree %d \n", x->id, x->degree); */
	    for ( i = 0; i < x->degree; i++ ) {
		/* printf("i %d x->id %d mates[i] %d\n", i,x->id,x->mates[i]); */
		w = x->mates_ad[i]; 
		if ( ! ( w->visited ) ) {
		    (void) Insert( w, &Q );
		    if ( ! ( w->left ) ) w->left = p;
		    /* printf("queueing id %d parent %d\n", w->id,w->left->id); */
		}
	    }
	}
	else {
	    /* printf("been here %d parent %d\n", x->id, x->left->id); */
	}
     }
    *vout = a;

    DestroyQueue(&Q);
}


#ifdef REMOVE
AdjacencyRec *adjacency_record ( IdType items[], int *start ) {
  int i, num_elts;
  AdjacencyRec *a, **ai;
  IdType *mates, t;

  /* when we come in *start is  id element */

  t = items[(*start)++];
  num_elts = *start;
  while ( items[num_elts] ) { 
    num_elts++;
  }

  num_elts = num_elts - *start;
  a = (AdjacencyRec *) malloc ( sizeof ( AdjacencyRec ));
  mates = (IdType *) malloc ( sizeof ( IdType ) * num_elts);
  ai = (AdjacencyRec **) malloc ( sizeof ( AdjacencyRec * ) * num_elts );
  for ( i = 0; i < num_elts; i++, (*start)++ ) {
    mates[i] = items[*start];
  }
  (*start)++;
  a->id = t;
  a->degree = num_elts;
  a->mates = mates;
  a->mates_ad = ai;
  a->left = NULL;
  a->right = NULL;
  a->weight = 0.0;
  return a;
}
#endif

void AddMateAddresses ( Graph *g ) {

    int i, j, k;
    AdjacencyRec *ai;

    /* it is convenient to work with the addresses so here we add them */

    /* for each AdjacencyRec i take its address ai
       for each AdjacencyRec k access each mate j
       if the mate[j] == id[i] then mates_ad[j] = ai
    */
  for ( i=0;i<g->number_of_verts; i++ ) {
      ai = g->recs[i];
#ifdef DEBUG
      printf("array num %d vertex %d degree %d\n",i,g->recs[i]->id,(g->recs[i])->degree); 
#endif
      for ( k=0;k<g->number_of_verts; k++ ) {
	  for (j=0;j<(g->recs[k])->degree; j++ ) {
	      if (abs(g->recs[k]->mates[j].m) == g->recs[i]->id) {
		  g->recs[k]->mates_ad[j] = ai;
#ifdef DEBUG
		 printf("i %d id %d mates %d\n", i, g->recs[i]->id, g->recs[k]->mates[j]);
#endif
	      }
	  }
      }
  }
}

void 
AddRightLinks (Graph *g, 
	       AdjacencyRec *v) 
{
    int i;
    AdjacencyRec *aj, *ai;

    /* starting from vertex v add right links to a left linked graph */

  i = 0;
  aj = NULL;
  ai = v;

  while ( (ai->left) && (i<g->number_of_verts) ) {
      ai->right = aj;
      aj = ai;
      ai = ai->left;
      i++;
  }
    ai->right = aj;
}

int 
compare_double (AdRecSort *r1, 
		AdRecSort *r2) {
    if ( (*r1).weight < (*r2).weight ) {
	return (-1);
    }
    else if ( (*r1).weight == (*r2).weight ) {
	return (0);
    }
    else {
	return (1);
    }
}

/* 
 * find the vertices that are not part of the spanning path 
 * and mark them unvisited 
 * set the weights for those that are visited
 * then add in the others
 * SP starts at vin. When done set vout to start of SP
 */
int
Augment_SP(Graph *g, 
	   int cur_sp,
	   AdjacencyRec *vin, 
	   AdjacencyRec **vout)  
{
    int i, j, num_wt, num_to_add;
    double temp_wt;
    AdjacencyRec *ai, *atl, *atr;
    AdRecSort *sort_array, *arp;

    /* remember right end of chain */
    atr = vin;

    if (NULL == (sort_array = (AdRecSort *) xmalloc (sizeof (AdRecSort) * 
						     g->number_of_verts)))
	return -1;

    arp = sort_array;

    /* 
     * NB the strategy of setting all unvisited 
     * will affect ALL SP's we may have found so beware
     */
    for ( i = 0; i < g->number_of_verts; i++ ) {
	/* only set to false those verts that I have visited in this SP */
	if (g->recs[i]->visited == cur_sp) {
	    g->recs[i]->visited = false;
	}
    }

    /* set weights for members of SP */

    ai = atr;
    i = 1;
    while ( ai->left ) {
	ai->visited = cur_sp;
	ai->weight = i++;
	ai = ai->left;
    }
    ai->visited = cur_sp;
    ai->weight = i;
/*
    printf("connected vertices\n");
    for ( i=0;i<g->number_of_verts; i++ ) {
	if (g->recs[i]->visited) {
	    printf("id %d weight %f",g->recs[i]->id,g->recs[i]->weight);
	    if(g->recs[i]->left) {
		printf(" %d ",g->recs[i]->left->id);
	    }
	    else {
		printf(" 0 ");
	    }
	    if(g->recs[i]->right) {
		printf(" %d ",g->recs[i]->right->id);
	    }
	    else {
		printf(" 0 ");
	    }
	    printf("\n");
	}
    }
*/
    /* printf("unconnected vertices\n"); */

    for ( i = 0, num_to_add = 0; i < g->number_of_verts; i++ ) {
	if (!(g->recs[i]->visited)) {
 	    /* printf("%d\n",g->recs[i]->id); */
	    temp_wt = 0.0;
	    num_wt = 0;
	    for (j = 0; j < (g->recs[i])->degree; j++ ) {
		/* find those mates which have been visited for this SP */
		if ( g->recs[i]->mates_ad[j]->visited == cur_sp) {
		    /*
		     * temp_wt += g->recs[i]->mates_ad[j]->weight;
		     * num_wt++;
		     */
		    temp_wt += (g->recs[i]->mates_ad[j]->weight *
				g->recs[i]->mates[j].weight);
		    num_wt += g->recs[i]->mates[j].weight;

/*		    
		       printf("i %d id %d mates %d\n", i, g->recs[i]->id, g->recs[i]->mates[j].m);
		       printf("weight %f \n", g->recs[i]->mates[j].weight);
*/		      
		}
	    }
	    if ( num_wt ) {
		temp_wt /= num_wt;
#ifdef DEBUG
		printf("temp_wt %f\n", temp_wt);
#endif
		num_to_add++;
		arp->weight = temp_wt;
		arp->adrec = g->recs[i];
		arp++;
	    }  
	}
    }

    if ( num_to_add ) {

	/* sort_array is an array of Adjacency record pointers. 
	   We need to copy the addresses of the vertices to add into this
	   array, then sort them on their weight, then intercalate them
	   into the spanning path
	*/

	qsort ( (void *) sort_array, num_to_add, sizeof(AdRecSort),
		(int (*)(const void *, const void *))compare_double );
/*
	for ( i=0; i<num_to_add; i++ ) {
	    printf("i %d sort_array[i].weight %f\n",i,sort_array[i].weight);
	}
*/
	/* move right to left looking for the spot to drop each new item.
	   needs to work when it is off either end as well.
	*/

	/* do right end */

	ai = atr;
	i = 0;
/*	sort_array[0].weight = 0.5; */
	if ( sort_array[0].weight < ai->weight ) {
	    sort_array[0].adrec->right = NULL;
	    sort_array[0].adrec->left = ai;
	    sort_array[0].adrec->visited = cur_sp; 
	    ai->right = sort_array[0].adrec;
	    
	    i++;
	    /*printf("id %d  %d %d %d\n",ai->id,ai->left->id,g->recs[2]->left->id,ai->right->id); */

	    /* reset right end of SP atr */

	    atr = ai->right;
	}

	/* do the left end */

	/* need to find left end from SP */

	ai = atr;
	while ( ai->left ) {
	    ai = ai->left;
	}
	atl = ai;
/*	sort_array[2].weight = 10.5; */
	if ( sort_array[num_to_add-1].weight > ai->weight ) {
	    sort_array[num_to_add-1].adrec->left = NULL;
	    sort_array[num_to_add-1].adrec->right = ai;
	    sort_array[num_to_add-1].adrec->visited = cur_sp; 
	    ai->left = sort_array[num_to_add-1].adrec;
	    atl = ai->left;
	    num_to_add--;
	}

/*
    for ( j=0;j<g->number_of_verts; j++ ) {
	if ( g->recs[j]->left && g->recs[j]->right ) {
	printf("id %d weight %f %d %d\n",g->recs[j]->id,g->recs[j]->weight,g->recs[j]->left->id,g->recs[j]->right->id);
    }
    }
*/
	/* do those between the SP elements */

	ai = atr;
	while ( (ai->left) && (i < num_to_add) ) {
	    /* printf("ai %d i %d aiw %f saw %f\n",ai->id,i,ai->weight, sort_array[i].weight ); */
	    if ( ai->weight > sort_array[i].weight ) {
		/*
		   printf("found home for %d",sort_array[i].adrec->id);
		   printf(" right neighbour %d ",ai->right->id);
		   printf("left neighbour %d\n",ai->id);
		   */
		sort_array[i].adrec->visited = cur_sp; 
		sort_array[i].adrec->right = ai->right;
		sort_array[i].adrec->left = ai;
		(ai->right)->left = sort_array[i].adrec;
		ai->right = sort_array[i].adrec;
		i++;
 	    } else {
		ai = ai->left;
	    }
	}
/*
    printf("connected vertices\n");
    for ( j=0;j<g->number_of_verts; j++ ) {
	if ( g->recs[j]->left && g->recs[j]->right ) {
	printf("id %d weight %f %d %d\n",g->recs[j]->id,g->recs[j]->weight,g->recs[j]->left->id,g->recs[j]->right->id);
    }
    }
*/
	/* list the final order in both directions */
/*
	ai = atl;
	while ( ai->right ) {
	    if ( ai->right && ai->left ) {
		printf("id %d l %d r %d\n",ai->id,ai->left->id,ai->right->id);
	    }
	    ai = ai->right;
	}
*/
	ai = atr;
	while ( ai->left ) {
	    if ( ai->left && ai->right ) {
/*
		printf("id %d l %d r %d\n",ai->id,ai->left->id,ai->right->id);
*/
	    }
	    ai = ai->left;
	}
    }
    *vout = atr;

    xfree(sort_array);
    return 0;
}

/* I have added all this new stuff to try to get the ordered contigs
 * emerging from the process in the correct relative orientations.
 * i.e. to work out which ones need to be complemented. Now each contig
 * has a direction value coded as 1 or -1 - the -1 contigs need to be
 * complemented. The algorithm assumes that the adjacency records
 * now store the mate data so that +ve mates are those to the right
 * and -ve mates are those to the left of the contig. 
 */

#define SIGN(i) ( ( i < 0 ) ? -1 : (( i == 0 ) ? 0 : 1 ) ) /* return 1,-1,0 */

int sign_mates_array_elt( Mates *array, int array_size, int value ) {

    /* find the sign of value value in array array 
     * return -1,+1 or 0 for error (also 0 for value=0)
     */

    int i, abs_value;

    abs_value = abs ( value );

    for ( i=0; i<array_size; i++ ) {
	if ( abs_value  == abs ( array[i].m ) ) {
	    return SIGN ( array[i].m );
	}
    }
    return 0;
}

void FindContigDirections ( Graph *g, AdjacencyRec *vin ) {
    int aid, ajd, aisign;
    AdjacencyRec *aj, *ai, *ak = NULL;

    /*  come in with left-right ordered contigs and mates coded
     *  +ve for right link, -ve for left link.
     *  Leave with the direction value of each contig signifying
     *  if it needs complementing: 1 = no change, -1 = complement.
     */

    ai = vin;

    while ( aj = ai->right ) {

	aid = ai->id;
	ajd = aj->id;

	/* the contigs should be in the correct left to right order
	 * so if a contig thinks the link to its right neighbour is
	 * to the left we must complement it
	 */

	aisign = sign_mates_array_elt( ai->mates, ai->degree, ajd );

	if ( aisign == -1) ai->direction = -1;

	ak = ai; 	/* save ai for last contig check */
	ai = aj;

    }

    if (ak) {
	/* when we get to the right end the last link should have
	 * been to the left. If not complement.
	 */
	aisign = sign_mates_array_elt( ai->mates, ai->degree, ak->id );
	if ( aisign == 1) ai->direction = -1;
    } else {
	ai->direction = 1; /* one contig => no change */
    }
}

AdjacencyRec *add_adjacency_record(id_wt *contig,
				   int num_elts)
{
  int i;
  AdjacencyRec *a, **ai;
  Mates *mates;
    
  a = (AdjacencyRec *) xmalloc ( sizeof ( AdjacencyRec ));

  mates = (Mates *) xmalloc ( sizeof ( Mates ) * num_elts);

  ai = (AdjacencyRec **) xmalloc ( sizeof ( AdjacencyRec * ) * num_elts );

  for (i = 0; i < num_elts; i++) {
    mates[i].m = contig[i].id;
    mates[i].weight = contig[i].weight;
#ifdef DEBUG
    printf("id %d weight %f num_elts %d\n", mates[i].m, mates[i].weight,
	   num_elts);
#endif
  }
  a->id = contig[0].id;
  a->direction = 1;
  a->degree = num_elts;
  a->mates = mates;
  a->mates_ad = ai;
  a->left = NULL;
  a->right = NULL;
  a->weight = 0.0;
  return a;
}

void
print_adjacency_record(AdjacencyRec *a)
{
    int i;

    printf("id %d degree %d weight %f \n", a->id, a->degree, a->weight);
    for (i =0; i < a->degree; i++) {
	printf("i %d mate %d \n", i, a->mates[i].m);
    }
}

/*
 * return 1 if the reading is within "dist" from the start or end of the contig
 * otherwise return 0
 */
int
TemplateDistance(GapIO *io,
		 gel_cont_t *gc,
		 int dist)
{
    GReadings r;
    int r_start;
    int r_end;

    gel_read(io, gc->read, r);
    r_start = r.position;
    r_end = r.position + (r.end - r.start - 2);

/*
    printf("r %d st %d end %d \n", gc->read, r_start, r_end);
    printf("c %d end %d \n", gc->contig, io_clength(io, gc->contig));
*/
    if ((r_start > dist) && (r_end < (io_clength(io, gc->contig) - dist))) {
	return 0;
    } else {
	return 1;
    }
}

/*
 * determine the direction of the template
 * return 0 if the template points to the right
 * return 1 if the template points to the left
 */
int
template_direction(template_c *t)
{
    int direction = 0; /* Default to zero in the unknown case */
    
    /* forward reading, guessed end */
    /*   ---->                      */
    /* st------------------ end     */

    /* reverse reading, guessed start */
    /*    ---->                       */
    /* end------------------start     */
/*    
    printf ("st %d end %d flags guessed st %d end %d\n", t->start, t->end,
	    t->flags & TEMP_FLAG_GUESSED_START, 
	    t->flags & TEMP_FLAG_GUESSED_END);
*/
    if (((t->start <= t->end) && !(t->flags & TEMP_FLAG_GUESSED_START))||
	((t->start >= t->end) && !(t->flags & TEMP_FLAG_GUESSED_END))) {
	direction = 0;
    }

    /* forward reading, guessed start */
    /*               <-----           */
    /* st------------------ end       */
    
    /* reverse reading, guessed end   */
    /*               <------          */
    /* end------------------start     */
    
    if (((t->start <= t->end) && !(t->flags & TEMP_FLAG_GUESSED_END)) ||
	    ((t->start >= t->end) && !(t->flags & TEMP_FLAG_GUESSED_START))) {
	direction = 1;
    }

    /* printf("direction %d \n", direction); */
    return direction;
}

/*
 * determine which end of the contig a template starts
 * return 0 for left 
 * return 1 for right
 */
int 
TemplateEnd(GapIO *io,
	    template_c *t,
	    int read,
	    int contig)
{
    int len;
    GReadings r;
    int start;
    
    len = io_clength(io, contig)/2;
    gel_read(io, read, r);

#ifdef DEBUG
    printf("len %d start %d end %d guessed start %d end %d\n", 
	   len, t->start, t->end, t->flags & TEMP_FLAG_GUESSED_START,
	   t->flags & TEMP_FLAG_GUESSED_END);
#endif
    start = r.position;
#ifdef DEBUG
    printf("len %d len/2 %d start %d \n", io_clength(io, contig), len,
	   start);
#endif

    if (start < len) {
	return 0;
    } else {
	return 1;
    }
}

/*
 * check whether the template is pointing outwards from the contig
 * return 1 if it points outwards (therefore correct)
 * return 0 if it points inwards (therefore incorrect)
 */
int
TemplateDirection(GapIO *io,
		  template_c *t,
		  int contig,
		  int read) 
{
    int t_start, t_end;

    get_template_positions(io, t, contig);
    t_start = MIN(MIN(t->start, t->end), t->min);
    t_end = MAX(MAX(t->start, t->end), t->max);

    /* for very small contigs, return true */
    if (io_clength(io, contig) < t_end - t_start) {
#ifdef DEBUG
	printf("clen %d tlen %d\n", io_clength(io, contig), 
	       t_end - t_start);
	printf("****************************HERE *****************\n");
#endif
	return 1;
    }

    /* left end of contig */
    if (TemplateEnd(io, t, read, contig) == 0) {
	/* pointing leftwards */
	if (template_direction(t) == 1) {
	    /* ok */
	    return 1;
	} else {
	    return 0;
	}
    } else {
	/* right end of contig */
	/* pointing rightwards */
	if (template_direction(t) == 0) {
	    /* ok */
	    return 1;
	} else {
	    return 0;
	}
    }
}

/*
 * find those contigs with readpairs and create adjacency structures
 * remove readpairs which are further than "dist" from the ends of the contig
 * and those which are pointing in an inconsistent direction (ie inwards)
 */
int
init_contig_order(GapIO *io,
		  AdjacencyRec ***a,
		  int *n_vert)
{
    int i, j;
    template_c **tarr;
    template_c *t;
    item_t *item;
    gel_cont_t *gc;
    gel_cont_t *first;
    item_t *top;
    AdjacencyRec **adj;
    int vert; 
    cnt_comp **contig;
    id_wt *contigs;
    int num_contigs;
    int cnt;
    GReadings r1, r2;

    int dist = 1000;

#define DO_DIST
#define DO_DIR

    if (Ntemplates(io) == 0)
	return -1;
    
    num_contigs = NumContigs(io);
    if (NULL == (adj = (AdjacencyRec **)xmalloc((NumContigs(io) + 1) *
					       sizeof(AdjacencyRec*)))){
	return -1;
    }
    if (NULL == (contig = (cnt_comp **)xmalloc((NumContigs(io) + 1) * 
					       sizeof(cnt_comp*)))){
	return -1;
    }
    for (i = 1; i <= NumContigs(io); i++) {
	if (NULL == (contig[i] = (cnt_comp *)xcalloc((NumContigs(io) + 1), 
					      sizeof(cnt_comp)))){
	    return -1;
	}
    }
    if (NULL == (contigs = (id_wt *)xmalloc((NumContigs(io) + 1) * 
					  sizeof(id_wt)))){
	return -1;
    }
	
    
    /* Initialise templates */
    if (NULL == (tarr = init_template_checks(io, 0, NULL, 1)))
	return -1;

    check_all_templates(io, tarr); 
    contig_spanning_templates(io, tarr);

#ifdef DEBUG
    for (i = 1; i <= Ntemplates(io); i++){
	if (tarr[i]) {
	    t = tarr[i];
	    printf("template %d \n", t->num);
	    for (item = head(t->gel_cont); item; item=item->next) {
		gc = (gel_cont_t *)(item->data);
		printf("contig %d read %d \n", gc->contig, gc->read);
	    }
	}
    }
    return -1;
#endif

    for (i = 1; i <= Ntemplates(io); i++){
	if (tarr[i]) {
	    t = tarr[i];
	    top = head(t->gel_cont);
	    first = (gel_cont_t *)(top->data);
#ifdef DEBUG
	    printf("READING1 %d contig %d \n", first->read, first->contig);
#endif
#ifdef DO_DIST
	    if (!TemplateDistance(io, first, dist)) {
		/* printf("REMOVED on distance \n"); */
		continue;
	    }
#endif
#ifdef DO_DIR
	    if (!TemplateDirection(io, t, first->contig, first->read)) {
		/* printf("REMOVED on direction \n"); */
		continue;
	    }
#endif
	    for (item = top->next; item; item=item->next) {
		gc = (gel_cont_t *)(item->data);

		/* may have multiple readings on the same end of template */
		if (first->contig == gc->contig)
		    continue;
#ifdef DEBUG
		printf("    READING2 %d contig %d \n", gc->read, gc->contig);
#endif		

#ifdef DO_DIST
		if (TemplateDistance(io, gc, dist)) {
#ifdef DO_DIR
		    if (TemplateDirection(io, t, gc->contig, gc->read)) {
			
			/*
			 * find the sense of each reading on contig
			 */
			gel_read(io, first->read, r1);
			gel_read(io, gc->read, r2);

#ifdef DEBUG
			printf("r1: c1 %d c2 %d sense %d\n",
			       first->contig, gc->contig, r1.sense);

			printf("r2: c1 %d c2 %d sense %d\n",
			       gc->contig, first->contig, r2.sense);
#endif
			if (r1.sense == 1) {
			    contig[first->contig][gc->contig].comp++;
			} else {
			    contig[first->contig][gc->contig].orig++;
			}

			if (r2.sense == 1) {
			    contig[gc->contig][first->contig].comp++;
			} else {
			    contig[gc->contig][first->contig].orig++;
			}
			contig[first->contig][gc->contig].cnt++;

		    } else {
			/* printf("REMOVED on direction \n"); */
		    }
#endif
		} else {
			/* printf("REMOVED on distance \n"); */
		}
#else
		contig[first->contig][gc->contig]++;
#endif

	    }
	    /* printf("\n"); */
	}
    }

    if (tarr) uninit_template_checks(io, tarr);
    
    for (i = 1; i <= num_contigs; i++) {
	for (j = 1; j <= num_contigs; j++) {
	    if (contig[i][j].cnt) {
		contig[j][i].cnt = contig[i][j].cnt;
	    }
	}
    }

    vert = 0;
    for (i = 1; i <= num_contigs; i++) {
	cnt = 1;
	for (j = 1; j <= num_contigs; j++) {
	    if (contig[i][j].cnt) {
		/* printf("i %d j %d rp %d \n", i, j, contig[i][j]); */
		contigs[0].id = i;
		contigs[0].weight = 0;

		/* 
		 * hopefully either comp or orig is 0 at this point - if
		 * neither are, the data is inconsistent.
		 * Use the largest value in direction calculation
		 */
		if (contig[i][j].comp > contig[i][j].orig) {
		    contigs[cnt].id = j * -1;
		} else {
		    contigs[cnt].id = j;
		}
		contigs[cnt++].weight = (double)contig[i][j].cnt; 
		/* printf("cnt %d weight %f \n", cnt-1, contigs[cnt-1].weight); */
	    }
	    
	}
	if (cnt > 1) {
	    adj[vert++] = add_adjacency_record(contigs, cnt);
	}
    }
    *n_vert = vert;
    *a = adj;

    xfree(contigs);
    for (i = 1; i <= num_contigs; i++) {
	xfree(contig[i]);
    }
    xfree(contig);

    return 0;
}

void free_contig_order(AdjacencyRec **a, int nvert) {
    int i;

    for (i = 0; i < nvert; i++) {
	xfree(a[i]->mates);
	xfree(a[i]->mates_ad);
	xfree(a[i]);
    }

    xfree(a);
}

/*
 * add list of contigs to list box
 */
void
UpdateAutomaticContigOrder(Tcl_Interp *interp,
			   GapIO *io,
			   id_dir *contig,
			   int num_contigs)
{
    Tcl_DString c_list;
    Tcl_DString d_list;
    Tcl_DString d_cmd;
    char buf[30];
    int i;

    Tcl_DStringInit(&c_list);
    Tcl_DStringInit(&d_list);
    Tcl_DStringInit(&d_cmd);
    
    for (i = 0; i < num_contigs; i++) {
	Tcl_DStringAppendElement(&c_list, get_contig_name(io, 
							  ABS(contig[i].id)));
	sprintf(buf, "%d", contig[i].dir);
	Tcl_DStringAppendElement(&d_list, buf);
    }

    Tcl_DStringAppendElement(&d_cmd, "create_contig_order_list");
    sprintf(buf, "%d", *handle_io(io));
    Tcl_DStringAppendElement(&d_cmd, buf);
    Tcl_DStringAppendElement(&d_cmd, Tcl_DStringValue(&c_list));
    Tcl_DStringAppendElement(&d_cmd, Tcl_DStringValue(&d_list));
    if (TCL_ERROR == Tcl_Eval(interp, Tcl_DStringValue(&d_cmd)))
	printf("UpdateAutomaticContigOrder %s\n", Tcl_GetStringResult(interp)); 

    Tcl_DStringFree(&c_list);
    Tcl_DStringFree(&d_list);
    Tcl_DStringFree(&d_cmd);
}

int
FindSpanningPath(Graph *g, 
		 int cur_sp,
		 id_dir *contig,
		 int *num_visited)
{
    AdjacencyRec *ai, *vout, *vout1, *vout2;
    int cnt;
    int num_contigs;
    int i;
    int start = -1;

    /* select first non-visited vertex */
    for (i = 0; i < g->number_of_verts; i++) {
	if (!g->recs[i]->visited) {
	    start = i;
	    break;
	}
    }
    if (start == -1)
	return -1;

    /* from a random vertex create the first tree */
    (void) GraphSearch(g,  cur_sp, g->recs[start], &vout);
 
    /* from the end of the first tree find a spanning path */
    (void) GraphSearch(g, cur_sp, vout, &vout1);

    /*  printf("parent %d\n",vout1->left->id); */
    
    /* 
     * when we get here vout1 is the end of the chain and 
     * we only have links in one direction so put in the reverse links
     */

    (void) AddRightLinks(g,  vout1 );
    
    /* list out our SP */
/*
   ai = vout1;
   i = 0;
   
    while ( (ai->left) && (i<num_verts) ) {
	if ((ai->id)&&(ai->left)&&(ai->right)) {
	    printf("id %d left %d right %d \n",ai->id,ai->left->id,ai->right->id);
	}
	ai = ai->left;
	i++;
    }
*/
/*
    while ( ai->right ) {
	if ( ai->right && ai->left ) {
	    printf("%d ", ai->id);
	}
	ai = ai->right;
    }
*/
    /* intercalate the vertices that are not part of the spanning path */
    if (-1 == Augment_SP(g, cur_sp, vout1, &vout))
	return -1;

    vout = vout1;
  
    /* list the final order in both directions */
#ifdef DEBUG
    printf("final path\n");
#endif
    ai = vout;

    cnt = 0;

#ifdef DEBUG
    printf("id %d \n", ai->id);
#endif
    while ( ai->left ) {
	if ( ai->left && ai->right ) {
#ifdef DEBUG
	    printf("id %d l %d r %d\n",ai->id,ai->left->id,ai->right->id);
#endif
	    cnt++;
	}
	ai = ai->left;
    }

#ifdef DEBUG
    printf("%d ", ai->id);
#endif

    vout2 = ai;
    FindContigDirections(g, vout2);

    contig[0].id = ai->id;
    contig[0].dir = ai->direction;
    num_contigs = 1;
    while ( ai->right ) {
	ai = ai->right;
#ifdef DEBUG
	printf("final %d %d\n", ai->id, ai->direction);
#endif
	contig[num_contigs].id = ai->id;
	contig[num_contigs++].dir = ai->direction;
    }
#ifdef DEBUG
    printf("\n");
#endif

    *num_visited = num_contigs;

    return 0;
}

int
find_contig_order(Tcl_Interp *interp,
		  GapIO *io) 
{
    char cmd[1024];
    IdType data[]= {17, 31,33,36,0,
			19, 21,27,28,0,
			21, 19,28,0,
			23, 24,30,0,
			24, 23,0,
			27, 19,30,32,0,
			28, 19,21,36,0,
			30, 23,27,32,0,
			31, 17,33,0,
			32, 27,30,0,
			33, 17,31,0,
			36, 17,28,0};

    /* AdjacencyRec *a[MAX_VERTICES], *ai, *aj, *vout, *vout1; */
    AdjacencyRec **a;
    Graph *g;
    int num_verts, num_elts, vert, i;
    int start, end;
    int non_visited, num_visited;
    int cur_sp;
    id_dir *contig;

    vert = 0;
    num_elts = sizeof ( data ) / sizeof (IdType);
    
    /* find read pairs and enter into structure */
    init_contig_order(io, &a, &vert); 
    
    /* copy the overlap data to the graph structure */

    start = 0;
    end = 0;

#ifdef HACK
    if (NULL == (a = (AdjacencyRec **)xmalloc(MAX_VERTICES *
					      sizeof(AdjacencyRec*)))){
	return;
    }
    vert = 0;
    while ( start < num_elts ) {
	a[vert] = adjacency_record ( data, &start );
	vert++;
    }
#endif
    num_verts = vert;
    
#ifdef DEBUG
    for (i = 0; i < num_verts; i++) {
	print_adjacency_record(a[i]);
    }
#endif    

    if (NULL == (g = (Graph *) xmalloc ( sizeof ( Graph ))))
	return -1;

    g->number_of_verts = num_verts;
    g->recs = a;

/*
  for ( i=0;i<g->number_of_verts; i++ ) {
      int j;
      printf("array num %d vertex %d degree %d\n",i,g->recs[i]->id,(g->recs[i])->degree);
      for (j=0;j<(g->recs[i])->degree; j++ ) {
	  printf("mates %d\n",g->recs[i]->mates[j]);
      }
  }
*/
    /* fill in the addresses of each vertex's mates */
    
    (void) AddMateAddresses( g );

    if (TCL_ERROR == Tcl_VarEval(interp, "init_contig_order_list", NULL))
	verror(ERR_WARN, "init_c_order_list",  "%s \n", Tcl_GetStringResult(interp));

    non_visited = g->number_of_verts;
    cur_sp = 1;

    /* initialise all vertices */
    for (i = 0; i < g->number_of_verts; i++) {
	g->recs[i]->visited = false;
	g->recs[i]->left = NULL;
	g->recs[i]->right = NULL;
    }

    if (NULL == (contig = (id_dir *)xmalloc((NumContigs(io) + 1) * 
					    sizeof(id_dir)))){
	return -1;
    }
	
    /*
     * find all spanning paths 
     */
    while (non_visited) {
	if (-1 == FindSpanningPath(g, cur_sp++, contig, &num_visited)) 
	    return -1;
	non_visited -= num_visited;
	if (num_visited > 1) {
	    UpdateAutomaticContigOrder(interp, io, contig, num_visited);
	}
    }

    sprintf(cmd, "contig_order_listbox %d ", *handle_io(io));
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "find_contig_order", " %s\n", Tcl_GetStringResult(interp));

    free_contig_order(g->recs, num_verts);
    xfree(g);

    xfree(contig);
    return 0;
}




