/*
 * Routines to read the blast matrix file formats.
 *
 * This is a lash-up designed to work specifically with the DNA matrices.
 */

#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "read_matrix.h"
#include "misc.h"

/*
 * Creates a score matrix by loading the blast format matrix from file 'fn'.
 * The axis of the matrix are the characters in base_order, with the matrix
 * size being strlen(base_order) in each direction.
 * Elements listed in base_order that do not appear in the file are entered
 * in the matrix as zero. (Should compute average and use this?)
 *
 * The matrix is allocated as an array of strlen(base_order) pointers to
 * arrays of strlen(base_order) integers.
 *
 * Returns the matrix for success, or NULL for failure.
 */
int **create_matrix(char *fn, char *base_order) {
    int **matrix;
    int len = strlen(base_order);
    int i, j, ncols = 0;
    FILE *fp;
    signed char lookup[256];
    int first_line = 1;
    char line[1024], *linep;
    signed char cols[256];

    /* Open the file */
    if (NULL == (fp = fopen(fn, "r")))
	return NULL;

    /* Allocate matrix */
    if (NULL == (matrix = (int **)xmalloc(len * sizeof(int *))))
	return NULL;

    for (i = 0; i < len; i++) {
	if (NULL == (matrix[i] = (int *)xcalloc(len, sizeof(int))))
	    return NULL;
    }

    /* Initialise our character value to base_order index translation. */
    memset(lookup, -1, 256);
    for (i = 0; i < len; i++) {
	lookup[toupper(base_order[i])] = i;
	lookup[tolower(base_order[i])] = i;
    }

    /*
     * Read the file itself.
     * Lines starting with '#' are comments.
     * After that, the first line indicates the order of elements in the cols.
     * Then we have the rows with data values.
     */

    while(fgets(line, 1024, fp)) {
	/* Skip comments */
	if (line[0] == '#')
	    continue;
	    
	if (first_line) {
	    /*
	     * For the first line, we initialise a translation from
	     * column characters in file to indices into 'matrix'.
	     */
	    first_line = 0;
	    for (j = i = 0, linep = line; *linep; linep++) {
		if (isspace(*linep))
		    continue;
		cols[j++] = lookup[(unsigned int)*linep];
	    }
	    ncols = j;
	} else {
	    /*
	     * Otherwise, we read the first character, and then fill in
	     * 'matrix' column by column.
	     */
	    int row_index, element;

	    for (linep = line; *linep && isspace(*linep); linep++)
		;
	    row_index = lookup[(unsigned int)*linep++];
	    if (row_index != -1) {
		for (j = 0; j < ncols; j++) {
		    element = strtol(linep, &linep, 10);
		    if (cols[j] != -1)
			matrix[row_index][cols[j]] = element;
		}
	    }
	}
    }

    /* Tidy up */
    fclose(fp);

#if 0
    /* Debug - print out the matrix */
    putchar(' ');
    for (j = 0; j < len; j++) {
	printf("   %c", base_order[j]);
    }
    putchar('\n');
    for (j = 0; j < len; j++) {
	printf("%c ", base_order[j]);
	for (i = 0; i < len; i++) {
	    printf("%3d ", matrix[i][j]);
	}
	putchar('\n');
    }
#endif
    
    return matrix;
}


/*
 * Free a matrix which has been allocated by create_matrix().
 * We still need to pass in base_order as this is how we indicate the size of
 * the matrix.
 */
void free_matrix(int **matrix, char *base_order) {
    int i, len = strlen(base_order);

    if (!matrix)
	return;

    for (i = 0; i < len; i++) {
	if (matrix[i]) xfree(matrix[i]);
    }
    xfree(matrix);

    return;
}
