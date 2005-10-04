#include <stdio.h>
#include <math.h>

#include "matrix.h"

void make_coord(m_coords *point, 
		double x, 
		double y) 
{
    point->x = x;
    point->y = y;
    point->w = 1.0;
}

void print_coords(m_coords *coords[], 
		  int num) {
    int i;
    for (i = 0; i < num; i++) {
	printf("%d %f %f\n", i, coords[i]->x, coords[i]->y);
    }
}

void make_zoom_matrix(double M[3][3], 
		      ZoomValues *z) 
{  
    M[0][0] = z->scale_x;
    M[0][1] = 0.0;
    M[0][2] = z->origin_x * (1 - z->scale_x);
    M[1][0] = 0.0;
    M[1][1] = z->scale_y;
    M[1][2] = z->origin_y * (1 - z->scale_y);
    M[2][0] = M[2][1] = 0.0;
    M[2][2] = 1.0;
}

void make_identity_matrix(double I[3][3]) 
{
  int i,j;

  for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
	  I[i][j] = 0.0;
      }
  }
  I[0][0] = I[1][1] = I[2][2] = 1.0;
}

void make_rotation_matrix(double M[3][3],
			  double angle)
{

    M[0][0] = cos(angle);
    M[0][1] = -1 * sin(angle);
    M[0][2] = 0.0;
    M[1][0] = sin(angle);
    M[1][1] = cos(angle);
    M[1][2] = 0.0;
    M[2][0] = 0.0;
    M[2][1] = 0.0;
    M[2][2] = 1.0;
}

void make_translation_matrix(double M[3][3],
			     double x,
			     double y)
{

    M[0][0] = 1.0;
    M[0][1] = 0.0;
    M[0][2] = x;
    M[1][0] = 0.0;
    M[1][1] = 1.0;
    M[1][2] = y;
    M[2][0] = 0.0;
    M[2][1] = 0.0;
    M[2][2] = 1.0;
}

void copy_matrix(double to[3][3], 
		 double from[3][3]) 
{
    int i, j;

    for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) {
	    to[i][j] = from[i][j];
	}
    }
}

/*
 * multiply B x A
 */
void multiply_matrices(double out[3][3], 
		       double b[3][3], 
		       double a[3][3]) 
{
    int i,j,k;
    
    for (i = 0; i < 2; i++) {
	for (j = 0; j < 3; j++) {
	    out[i][j] = 0.0;
	    for (k = 0; k < 3; k++) {
		out[i][j] = out[i][j] + b[i][k] * a[k][j];
	    }
	}
    }
    out[2][0] = out[2][1] = 0.0;
    out[2][2] = 1.0;
}

/*
 * multiply COORDS x Matrix
 */
void multiply_coords_by_matrix(m_coords *out, 
			       double m[3][3], 
			       m_coords *in) 
{
    
    out->x = m[0][0] * in->x + m[0][1] * in->y + m[0][2] * in->w;
    out->y = m[1][0] * in->x + m[1][1] * in->y + m[1][2] * in->w;
}

