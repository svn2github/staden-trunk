#ifndef _READ_MATRIX_H_
#define _READ_MATRIX_H_

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
int **create_matrix(char *file_name, char *base_order);

/*
 * Free a matrix which has been allocated by create_matrix().
 * We still need to pass in base_order as this is how we indicate the size of
 * the matrix.
 */
void free_matrix(int **matrix, char *base_order);

#endif

