#include "os.h"

/*
 *-----------------------------------------------------------------------------
 * Debugging for the BUBBL3 routine
 *-----------------------------------------------------------------------------
 */
#if 0
f_proc_ret bubbl3x_(int *val, int *list, int *listel, int *listal, int *idim) {
    int i;

    printf("---- %d ----\n", *val);
    for (i=0; i<*idim; i++) {
	printf("%d: %d %d %d\n", i, list[i], listel[i], listal[i]);
    }
}
#endif


/*
 *-----------------------------------------------------------------------------
 * Our own quick sort implementation to sort three arrays simultaneously, using
 * the first as the key.
 *-----------------------------------------------------------------------------
 */
#if 0
static void sort3(int *a1, int *a2, int *a3, int left, int right) {
    int i, last, tmp, ind;

    if (left >= right)
        return;

    /* Swap(arrays[left], arrays[(left+right)/2]) */
    ind = (left + right)/2;
    tmp = a1[left]; a1[left] = a1[ind]; a1[ind] = tmp;
    tmp = a2[left]; a2[left] = a2[ind]; a2[ind] = tmp;
    tmp = a3[left]; a3[left] = a3[ind]; a3[ind] = tmp;

    last = left;

    for (i= left+1; i <= right; i++){
        if (a1[i] >= a1[left]) {
	    ++last;
	    /* Swap(arrays[last], arrays[i]) */
	    tmp = a1[last]; a1[last] = a1[i]; a1[i] = tmp;
	    tmp = a2[last]; a2[last] = a2[i]; a2[i] = tmp;
	    tmp = a3[last]; a3[last] = a3[i]; a3[i] = tmp;
	}
    }
    /* Swap(arrays[left], arrays[last]) */
    tmp = a1[last]; a1[last] = a1[left]; a1[left] = tmp;
    tmp = a2[last]; a2[last] = a2[left]; a2[left] = tmp;
    tmp = a3[last]; a3[last] = a3[left]; a3[left] = tmp;

    /* Recurse */
    sort3(a1, a2, a3, left, last-1);
    sort3(a1, a2, a3, last+1, right);
}

f_proc_ret bubbl3_(int *list, int *listel, int *listal, int *idim) {
    sort3(list, listel, listal, 0, *idim-1);
    f_proc_return();
}
#endif


/*
 *-----------------------------------------------------------------------------
 * A BUBBL3 implementation using the system qsort function. This sorts a single
 * array which is an index into the list. All three lists are then reordered
 * using this index array.
 *-----------------------------------------------------------------------------
 */
#if 0
static int *lpointer;

int three_sort(const int *i1, const int *i2) {
    return (lpointer[*i2] - lpointer[*i1]);
}


f_proc_ret bubbl3_(int *list, int *listel, int *listal, int *idim) {
    int *ind, *listn;
    int i;

    ind = (int *)malloc(sizeof(int) * *idim);
    listn = (int *)malloc(sizeof(int) * *idim);
    for (i=0; i<*idim; i++) {
	ind[i] = i;
    }
    lpointer = list;
    qsort(ind, *idim, sizeof(int), three_sort);

    for (i=0; i<*idim; i++)
	listn[i] = list[ind[i]];
    memcpy(list, listn, sizeof(int) * *idim);
    for (i=0; i<*idim; i++)
	listn[i] = listel[ind[i]];
    memcpy(listel, listn, sizeof(int) * *idim);

    for (i=0; i<*idim; i++)
	listn[i] = listal[ind[i]];
    memcpy(listal, listn, sizeof(int) * *idim);

    xfree(ind);
    xfree(listn);

    f_proc_return();
}
#endif


/*
 *-----------------------------------------------------------------------------
 * Finally a direct copy of the original fortran bubble sort routine. Despite
 * the fact that bubble sort is a terrible sort (except for nearly sorted
 * data), it seems to work as well as either of the above two for our data.
 * Also note that as this is the original algorithm assemblies will be the
 * same. With another sort routine the order of items may differ (when
 * multiple items have the same sort key) which gives marginally differing
 * alignments.
 *
 * The reason that this has been recoded in C is simply that when unoptimised
 * it is around 6-7 times faster! Naturally when optimised both Fortran and
 * C give similar speeds.
 *
 * The original routine follows:
 *
 * C     BUBBL3
 * C   SUBROUTINE TO SORT INTEGER ARRAY (LIST) INTO ASCENDING  ORDER
 * C
 *       SUBROUTINE BUBBL3(LIST,LISTEL,LISTAL,IDIM)
 * C   AUTHOR: RODGER STADEN
 *       INTEGER LIST(IDIM),LISTEL(IDIM),LISTAL(IDIM)
 * C      CALL BUBBL3X(0, LIST, LISTEL, LISTAL, IDIM)
 * C
 * C   SET POINTERS TO ZERO
 *       I=0
 *       J=0
 * C
 * 10    CONTINUE
 * C
 * C   SET I=J IF WE HAVE JUST CORRECTLY POSITIONED AN ELEMENT
 *       IF(J.GT.I)I=J
 * C
 * C   INCREMENT POINTER TO NEXT ELEMENT
 *       I=I+1
 * C   TEST FOR END OF ARRAY
 *       IF(I.EQ.IDIM) THEN
 * C         CALL BUBBL3X(1, LIST, LISTEL, LISTAL, IDIM)
 *          RETURN
 *       ENDIF
 * C
 * 20    CONTINUE
 * C
 * C   COMPARE ADJACENT ELEMENTS
 *       IF(LIST(I).GE.LIST(I+1))GO TO 10
 * C
 * C   FIRST MOVE THIS ELEMENT? IF SO SET POINTER TO ITS INITIAL POSITION
 *       IF(J.LT.I)J=I
 * C
 * C   EXCHANGE ADJACENT ELEMENTS
 *       ITEMP=LIST(I)
 *       LIST(I)=LIST(I+1)
 *       LIST(I+1)=ITEMP
 * C
 *       ITEMP=LISTEL(I)
 *       LISTEL(I)=LISTEL(I+1)
 *       LISTEL(I+1)=ITEMP
 *       ITEMP=LISTAL(I)
 *       LISTAL(I)=LISTAL(I+1)
 *       LISTAL(I+1)=ITEMP
 * C
 * C
 * C   DECREMENT BACK THRU LIST WITH THIS ELEMENT
 *       IF(I.GT.1)I=I-1
 * C
 *       GO TO 20
 *       END
 *
 *-----------------------------------------------------------------------------
 */
f_proc_ret bubbl3_(int *list, int *listel, int *listal, int *idim) {
    int i = 0, j = 0, temp;

    /* Not ANSI, but good for a quick hack */
    list--;
    listel--;
    listal--;

 label10:
    if (j > i)
	i = j;

    i++;

    if (i == *idim) {
	return;
    }

 label20:
    if (list[i] >= list[i+1])
	goto label10;
	
    if (j < i)
	j = i;

    temp=list[i];
    list[i]=list[i+1];
    list[i+1]=temp;

    temp=listel[i];
    listel[i]=listel[i+1];
    listel[i+1]=temp;

    temp=listal[i];
    listal[i]=listal[i+1];
    listal[i+1]=temp;

    if (i > 1)
	i--;

    goto label20;
}
