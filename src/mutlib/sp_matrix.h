#ifndef _SP_MATRIX_H_
#define _SP_MATRIX_H_


#include <cstdio>      /* For FILE* */


namespace sp {



/* Data */
typedef struct
{
    int** data;
    int   rows;
    int   cols;

}matrix_t;



/* Services */
int   matrix_create( matrix_t* m, int rows, int cols );
int*  matrix_row( matrix_t* m, int row );
int   matrix_rows( matrix_t* m );
int   matrix_cols( matrix_t* m );
int** matrix_data( matrix_t* m );
void  matrix_fill( matrix_t* m, int value );
void  matrix_destroy( matrix_t* m );
void  matrix_print( matrix_t* m, std::FILE* s );


}


#endif
