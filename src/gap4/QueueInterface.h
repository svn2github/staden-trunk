#include "SeqQueueTypes.h"     /* imports the data type definitions of */
                               /* ItemType and Queue */

   /* defined operations */

extern void InitializeQueue(Queue *Q);
   /* Initialize the queue Q to be the empty queue */

extern int Empty(Queue *Q);
   /* Returns TRUE == 1 if and only if the queue Q is empty */

extern int Full(Queue *Q);
   /* Returns TRUE == 1 if and only if the queue Q is full */

extern void Insert(ItemType R, Queue *Q);
   /* If Q is not full, insert a new item R onto the rear of Q */

extern void Remove(Queue *Q, ItemType *F);
   /* If Q is non-empty, remove the frontmost item of Q and put it in F */

