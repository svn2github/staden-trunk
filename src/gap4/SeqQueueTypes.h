/*    "SeqQueueTypes.h"    */
#include <stdio.h>
#include <stdlib.h>
#include "array.h"

typedef int IdType;

#define MAX_VERTICES 100     /* ie the number of contigs */
 
/* most of the time we need the address of the AdjacencyRec for each mate
   so store them in mates_ad
*/
 
typedef struct _Mates {
    IdType m;
    double weight;
} Mates;

typedef struct _id_wt {
    int id;
    double weight;
} id_wt;

typedef struct _cnt_comp {
  int cnt;
  int comp;
  int orig;
} cnt_comp;

typedef struct _id_dir {
  int id;
  int dir;
} id_dir;

typedef struct  _AdjacencyRec {
    IdType  id;                      /* contig id */
    int     direction;               /* the final direction +/- 1 */
    int     degree;                  /* number of mates */
    Mates   *mates;                  /* array of mate ids */
    struct  _AdjacencyRec **mates_ad;/* array of pointers to mates */
    int     visited;                 /* flag to say we have been here */
    double  weight;                  /* weight to define position in SP */
    struct  _AdjacencyRec *left;     /* pointer to left neighbour in tree */
    struct  _AdjacencyRec *right;    /* pointer to right neighbour in tree */
} AdjacencyRec;
 
typedef struct {                     /* our graph structure */
    int     number_of_verts;         /* number of vertices (contigs) */
    AdjacencyRec **recs;             
} Graph;
 
typedef struct {                     /* a temporary structure for sorting */
    AdjacencyRec *adrec;             /*   the non-SP contigs into left-right */
    double       weight;             /*   order, prior to intercalation */
} AdRecSort;
 
typedef  AdjacencyRec *ItemType;
 

typedef  struct {  
  /* our queue for the breadth */
  /* first search */
  int       Count;                  /* number of queue items */
  int       Front;                  /* front of queue */
  int       Rear;                   /* rear of queue */
  Array     Items;		    /* array of 'ItemType', queued items */
} Queue;



