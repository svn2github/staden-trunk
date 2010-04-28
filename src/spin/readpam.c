/* routines to read the following matrix file formats:

   blast

*/
#include <staden_config.h>

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "sequence_formats.h"
#include "readpam.h"
#include "misc.h" /* MIN MAX */
#include "dna_utils.h"

/* int score_matrix[MAX_SCORE_MATRIX][MAX_SCORE_MATRIX]; */

int **score_matrix = {0};

#define DUMMY_SCORE 99999
#define ROUND(x)	((x) < 0.0 ? ceil((x) - 0.5) : floor((x) + 0.5))

/************************************************************/
/* 
 * default dna matrix 
 */
int dna_matrix[5][5] = {
    { 1,0,0,0,0 },
    { 0,1,0,0,0 },
    { 0,0,1,0,0 },
    { 0,0,0,1,0 },
    { 0,0,0,0,0 }
};
void
identity_dna_matrix(int ***matrix)
{
    int i, j;

    for (i = 0; i < 5; i++) {
	for (j = 0; j < 5; j++) {
	    (*matrix)[i][j] = dna_matrix[i][j];
	}
    }
}

void 
identity_prot_matrix (int ***matrix)
{
    int i, j;
    int *p = get_protein_lookup();

    for (i = 0; i < char_set_size; i++) {
	for (j = 0; j < char_set_size; j++) {
	    if ((i == j) && (i != p['X']) && (i != p['*']) && (i != p['-']))
		(*matrix)[i][j] = 1;
	    else
		(*matrix)[i][j] = 0;
	/* printf("%d ", (*matrix)[i][j]);  */
	}
	/* printf("\n"); */
    }

}

/* read in blast and staden format matrix files */

/* Deal with 2 special line types: comments that have "#" in column 0
 * and the first data record has " ". 
 * Return codes: 0 OK, else is an offset off the end of vector 
*/

int get_matrix (int *vector, 
		int max_vector, 
		int *dimension1, 
		int *dimension2, 
		FILE *fp)
{

#define MAX_LINE 256
    
  /*    char **matrix_pt;*/
    char line[MAX_LINE];
    int i, j, k, row, rows, column, columns, looking_for_start = 1;
    int lookup[100],offset;
    column = columns = rows = 0;

    
    /*   if(NULL == (matrix_pt = (char **)malloc(MAX_SCORE_MATRIX *sizeof(char *))))
   printf("erroe in memory allocation\n");
        for (i = 0; i < MAX_SCORE_MATRIX; i++) {
   if(NULL == (matrix_pt[i] = (char *) malloc(MAX_SCORE_MATRIX *
   sizeof(char))))
    printf("erroe in memory allocation\n");
    }*/

 
    for (i = 0; i < 100; i++) {
	lookup[i] = char_lookup[char_set_size - 1];
    }

    for (i = 0; i < max_vector; i++) {
	vector[i] = DUMMY_SCORE;
    }
    while ( fgets( line,sizeof(line),fp ) ) {

	/* Check for special lines of type "#"*/

	if ( '#' != line[0] ) {

	    if ( looking_for_start ) {

		if ( ' ' == line[0] ) {
		    looking_for_start = 0;
		    for (j = 0; j < MAX_LINE && line[j]; j++) {
			if ( isgraph ( (int) line[j]) ) {
			    lookup[column] = char_lookup[line[j]];
			    column++;
			}
		    }
		}
	    }
	    else {

		/* we are into the data */

		/*
		 * need to check rows < columns - this will ignore anything at
		 * the bottom of the matrix
		 */
		if (rows >= column) {
		    break;
		}
		rows++;
		/*  printf("\n");*/
		row = char_lookup[line[0]];
		i = 1;

		for(k = 0; k < column; k++) {
		
		    for(;line[i]&&(line[i]==' ');i++)
			;
		         j = atoi(&line[i]);
			 /*     j = strtol(&line[i],matrix_pt,10);
                     if (isgraph((int)matrix_pt[0][0] && (int)matrix_pt[0][0]<48||(int)matrix_pt[0][0]>57))
		     return -1;*/
	
			 /*    printf("%3d", j);*/
		         offset = (row * column) + lookup[k];

		      if (offset >= max_vector) return offset;
		      vector[offset] = j;


		      for(;line[i]&&(line[i]!=' ');i++){
		      if(isgraph((int)line[i]) && ((int)line[i] <48 &&(int)line[i]!=45 || (int)line[i]==45 &&(int)line[i-1]!= ' ' || (int)line[i] > 57))
		      return -1;
		    }

		    
		}
	    }
	}
    }
    if (column >=MAX_SCORE_MATRIX || rows >=MAX_SCORE_MATRIX) 
	return -1;
    
    *dimension1 = column;
    *dimension2 = rows;
     return looking_for_start;
}

int
find_matrix_average(int **matrix)
{
    int i, j;
    int avg = 0;

    for (i = 0; i < 22; i++) {
	for (j = 0; j < 22; j++) {
	    avg += matrix[i][j];
	}
    }
    /* avg = (int) (0.5 + (double) avg / 484.0 ); */
    avg = (int)ROUND((double) avg / 484.0);
    return avg;
}

void
set_score_matrix(int **matrix)
{
    score_matrix = matrix;
}

int 
create_pam_matrix(char *file_name_ptr,
		  int ***matrix)
{
#define MAX_VECTOR 1000
    int i, j, max_vector;
    FILE *fp;
    int *vector, dimension1, dimension2;
    int avg;

    #ifdef PAM_DEBUG
    int b[]={
     0,15,12, 3, 2,14, 4, 6, 7, 8,10, 9,11, 5,13,16,17,19,20,18, 1,21,22,23,24};

/*   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * -
*/
    #endif
#ifdef SIP_DEBUG
    int b[]={
     2,16,17,13, 0, 6,12, 3, 4,14, 1,21, 7,15, 9,11, 8,10,18, 5,20,19,24,22,23};
/*
     C  S  T  P  A  G  N  D  E  Q  B  Z  H  R  K  M  I  L  V  F  Y  W  -  X  *
*/
#endif

    if (NULL == (fp = fopen(file_name_ptr,"r"))) {
	verror(ERR_WARN, "file open", "Unable to open file %s", file_name_ptr);
	return -1;
    }

    if ( ! ( vector = (int *) malloc(sizeof(int)*MAX_VECTOR))) { 
        return -1;
    }
    /* normally the matrix array will be global.
       here we get the matrix data in an array vector
       and copy it to the matrix array. get_matrix
       also returns the row and column sizes in
       dimension1 and dimension2. The vector array
       should be freed. I guess max_vector should
       be set to the size of the max matrix. decide
       if the caller or get_matrix should handle
       file opening?
    */
    max_vector = MAX_VECTOR;
    i = get_matrix ( vector, max_vector, &dimension1, &dimension2, fp );
    if (i) {
	free(vector);
	return -1;
    }

    /* can't deal with non-square matrices */
    if (dimension1 != dimension2) {
	free(vector);
	return -1;
    }

    /* copy 1 dimensional vector array into 2 dimensional matrix array */
    for(i = 0; i < dimension1; i++) {
	for(j = 0; j < dimension2; j++) {
	      (*matrix)[i][j] = vector[(dimension1*i)+j];
	     
		}
    }
    /*   printf("\ndim1 %d dim2 %d \n", dimension1, dimension2);
        for(i = 0; i < dimension1; i++) {
	for(j = 0; j < dimension2; j++) {
	    printf("%3d", (*matrix)[i][j]);
	}
	printf("\n");
	}

	printf("\n\n\n");*/

    avg = find_matrix_average(*matrix);

    /* printf("DIM1 %d DIM2 %d \n",dimension1, dimension2); */
    /* replace 'missing' values with the matrix average */
    for(i = 0; i < dimension1; i++) {
	for(j = 0; j < dimension2; j++) {
	    if ((*matrix)[i][j] == DUMMY_SCORE) {
		(*matrix)[i][j] = avg;
	    }
	}
    }
    

#ifdef DEBUG
    printf("\ndim1 %d dim2 %d \n", dimension1, dimension2);
    for(i = 0; i < dimension1; i++) {
	for(j = 0; j < dimension2; j++) {
	    printf("%3d", (*matrix)[i][j]);
	}
	printf("\n");
    }

    printf("\n\n\n");
    /*
    for(i=0;i<dimension1;i++) {
	for(j=0;j<dimension2;j++) {
	    printf("%3d",(*matrix)[b[i]][b[j]]);
	}
	printf("\n");
    }
    */
#endif
    fclose ( fp );
    free ( vector );
    return 0;
}

void print_matrix(int **matrix)
{
    int i, j;
    
    for(i = 0; i < MAX_SCORE_MATRIX; i++) {
	for(j = 0; j < MAX_SCORE_MATRIX; j++) {
	    printf("%3d", matrix[i][j]);
	}
	printf("\n");
    }

    printf("\n\n\n\n");
}




