#include <cassert>
#include <sp_matrix.h>
extern "C" {
#include <Misc/misc.h>      /* For xmalloc() */
}

namespace sp {


int matrix_create( matrix_t* m, int rows, int cols )
{
    int r;
    assert(m != NULL);
    assert(rows>0);
    assert(cols>0);


    /* Allocate array for rows */
    m->data = (int**) xmalloc( rows * sizeof(int*) );
    if( !m->data )
        return -1;

    
    /* Initialise it */
    for( r=0; r<rows; r++ )
        m->data[r] = 0;
    m->rows = rows;
    m->cols = 0;


    /* Allocate each row */
    for( r=0; r<rows; r++ )
    {
        m->data[r] = (int*) xmalloc( cols * sizeof(int) );
        if( !m->data[r] )
        {
            matrix_destroy(m);
            return -1;
        }
    }
    m->cols = cols;
    return 0;
}



void matrix_destroy( matrix_t* m )
{
    int r;
    assert(m != NULL);


    /* If nothing to do, exit */
    if( !m->data )
        return;


    /* Free rows */
    for( r=0; r<m->rows; r++ )
    {
        if( m->data[r] )
            xfree( m->data[r] );
    }


    /* Free array for rows */
    xfree( m->data );
    m->data = 0;
    m->rows = 0;
    m->cols = 0;
}


int* matrix_row( matrix_t* m, int row )
{
    assert(m != NULL);
    assert(row<m->rows);
    if( !m || (row<0) || (row>=m->rows) )
        return 0;
    else
        return m->data[row];
}



int matrix_rows( matrix_t* m )
{
    assert(m != NULL);
    return m->rows;
}



int matrix_cols( matrix_t* m )
{
    assert(m != NULL);
    return m->cols;
}


int** matrix_data( matrix_t* m )
{
    assert(m != NULL);
    return m->data;
}


void matrix_fill( matrix_t* m, int value )
{
    int r;
    int c;
    assert(m != NULL);
    assert(m->data != NULL);
    

    /* Valid input? */
    if( !m || !m->data )
        return;


    /* Fill the matrix */
    for( r=0; r<m->rows; r++ )
        for( c=0; c<m->cols; c++ )
            m->data[r][c] = value;
}


void matrix_print( matrix_t* m, std::FILE* s )
{
    int r;
    int c;
    assert(m != NULL);
    assert(s != NULL);
    assert(m->data != NULL);
    

    /* Valid input? */
    if( !m || !s || !m->data )
        return;


    /* Print the matrix to a stream */
    for( r=0; r<m->rows; r++ )
    {
        for( c=0; c<m->cols; c++ )
            std::fprintf( s, "%6d ", m->data[r][c] );
        std::fprintf( s, "\n" );
    }
}



#ifdef TEST_SP_MATRIX


#include <stdio.h>
int main( void )
{
    matrix_t m;
    std::printf("Creating matrix...\n" );
    matrix_create(&m, 10, 10);
    std::printf("Rows = %d\nCols = %d\n", matrix_rows(&m), matrix_cols(&m) );
    matrix_fill( &m, 123456 );
    matrix_print( &m, stdout );
    std::printf("Destroying matrix...\n" );
    matrix_destroy(&m);
    return 0;
}



#endif



}
