#include <stdio.h>
#include <string.h>

#include "IO.h"
#include "io-reg.h"
#include "alter_rel.h"
#include "misc.h"
#include "array.h"
#include "gap-error.h"
#include "fort.h"
#include "tagUtils.h"
#include "dis_readings.h"


/* ------------------------------------------------------------------------- */
/* Annotation list handling */

/*
 * Search for an annotation number 'search' within a list starting at 'anno'.
 * Assume memory of correct size in 'loop' to clear and use for checking
 * looping.
 * Returns N for found, where N is 1 for found once, 2 for loop.
 *         0 for not found.
 */
static int in_list(GapIO *io, int anno, int search, char *loop) {
    GAnnotations a;
    int count = 0;
    
    memset(loop, 0, Nannotations(io)+1);

    while (anno && !loop[anno]) {
        if (anno == search)
	    count++;

        if (0 == GT_Read(io, arr(GCardinal, io->annotations, anno-1),
                         &a, sizeof(a), GT_Annotations)) {
            loop[anno]++;
            anno = a.next;
        } else
            anno = 0;
    }

    return anno ? 2 : count;
}

/*
 * Looks for the referee of an annotation. The type (0=reading, 1=contig,
 * 2=freelist) and element number are stored.
 * If init==1 then we initialise the search (and do one search). If
 * init==0 then we continue from the last found. If init==2 then we quit
 * the search. 
 *
 * We return N for found (where N is 1 for found, 2 for found in loop),
 *           0 for not found, and -1 for error.
 */
int annotation_address(GapIO *io, int init, int search, int *type, int *num) {
    int i, count;
    GContigs c;
    GReadings r;
    static char *loop = NULL;
    static int current_mode;    /* 0==reading, 1==contig, 2==freelist */
    static int current_point;   /* contig or reading number */

    if (2 == init) {
	if (loop)
	    xfree(loop);

	return 0;
    }

    if (1 == init) {
	current_mode = 0;
	current_point = 0;

	if (NULL == (loop = (char *)xcalloc(Nannotations(io)+1, 1))) {
	    return -1;
	}
    }


    /* Contigs */
    if (0 == current_mode) {
	for (i = current_point; i < NumContigs(io); i++) {
	    GT_Read(io, arr(GCardinal, io->contigs, i),
		    &c, sizeof(c), GT_Contigs);

	    if (count = in_list(io, c.annotations, search, loop)) {
		*type = 1;
		*num = i+1;
		current_point = i+1;

		return count;
	    }
	}

	current_mode++;
	current_point = 0;
    }

    /* Readings */
    if (1 == current_mode) {
	for (i = current_point; i < NumReadings(io); i++) {
	    gel_read(io, i+1, r);

	    if (count = in_list(io, r.annotations, search, loop)) {
		*type = 0;
		*num = i+1;
		current_point = i+1;
		
		return count;
	    }
	}

	current_mode++;
	current_point = 0;
    }

    /* Freelist */
    if (2 == current_mode) {
	if (count = in_list(io, io->db.free_annotations, search, loop)) {
	    *type = 2;
	    *num = 0;

	    current_mode++;
	    return count;
	}

	current_mode++;
    }

    return 0;
}

/* ------------------------------------------------------------------------- */
/* Reinitialising of the contig order */
int reset_contig_order(GapIO *io) {
    int i;

    if (io->db.contig_order == 0) {
	if (-1 == (io->db.contig_order = allocate(io, GT_Array))) {
	    GAP_ERROR_FATAL("Initialising contig order array");
	    return -1;
	}
	
	io->contig_order = ArrayCreate(sizeof(GCardinal), io->db.Ncontigs);
	ArrayDelay(io, io->db.contig_order, io->db.Ncontigs, io->contig_order);
    }
    
    ArrayRef(io->contig_order, io->db.Ncontigs-1);
    
    for (i = 0; i < io->db.Ncontigs; i++)
	arr(GCardinal, io->contig_order, i) = i+1;
    
    ArrayDelay(io, io->db.contig_order, io->db.Ncontigs, io->contig_order);
    DBDelayWrite(io);

    flush2t(io);
    return 0;
}
