#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct ZoomValues_ {
    double origin_x;
    double origin_y;
    double scale_x;
    double scale_y;
} ZoomValues;

typedef struct m_coords_ {
  double x;
  double y;
  double w;
} m_coords;


void make_coord(m_coords *point, double x, double y);

void make_identity_matrix(double I[3][3]) ;

void make_zoom_matrix(double M[3][3], ZoomValues *z);
void make_translation_matrix(double M[3][3], double x, double y);

void copy_matrix(double to[3][3], 
		 double from[3][3]); 

void multiply_matrices(double out[3][3], 
		       double b[3][3], 
		       double a[3][3]); 

void multiply_coords_by_matrix(m_coords *out, double m[3][3], m_coords *in);

void print_matrix(double S[3][3]);

#endif
