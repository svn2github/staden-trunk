#include <limits.h>
#include <float.h>
#include <math.h>

#include "misc.h"

/* general array stuff */

void fill_double_array ( double array[], int num_elements, double filler ) {
    int i;

    for (i=0;i<num_elements;i++) {
	array[i] = filler;
    }
}

void fill_int_array ( int array[], int num_elements, int filler ) {
    int i;

    for (i=0;i<num_elements;i++) {
	array[i] = filler;
    }
}

int max_int_array ( int array[], int num_elements ) {
    int i, max = INT_MIN;

    for (i=0;i<num_elements;i++) {
	max = MAX( max, array[i] );
    }
    return max;
}

int min_int_array ( int array[], int num_elements ) {
    int i, min = INT_MAX;

    for (i=0;i<num_elements;i++) {
	min = MIN( min, array[i] );
    }
    return min;
}

double max_double_array ( double array[], int num_elements ) {
    int i;
    double max = -DBL_MAX;

    for (i=0;i<num_elements;i++) {
	max = MAX( max, array[i] );
    }
    return max;
}

double min_double_array ( double array[], int num_elements ) {
    int i;
    double min = DBL_MAX;
    for (i=0;i<num_elements;i++) {
	min = MIN( min, array[i] );
    }
    return min;
}

void reset_zeroes ( double table[], int num_elements, double correction ) {
    int i;
    double small = DBL_EPSILON;

    for ( i=0; i<num_elements; i++ ) {
	if ( table[i] <= small ) table[i] = correction;
    }
}

double sum_double_array ( double array[], int num_elements ) {
    int i;
    double t=0.0;
    for ( i=0; i<num_elements; i++ ) {
	t += array[i];
    }
    return t;
}

void div_double_array ( double array[], int num_elements, double divisor ) {
    int i;
    double small = DBL_EPSILON;

/* sanity check */
    if (( divisor <= small ) && ( divisor >= 0.0 )) return;
    for ( i=0; i<num_elements; i++ ) {
	array[i] /= divisor;
    }
}

void mult_double_array ( double array[], int num_elements, double mult ) {
    int i;

    for ( i=0; i<num_elements; i++ ) {
	array[i] *= mult;

    }
}

void ratio_double_arrays ( double denominator[], double numerator[],
			  int num_elements ) {
  /* return the ratio of the elements in the arrays in denominator */

  int i;

  for(i=0;i<num_elements;i++) {
    if ( numerator[i] > DBL_EPSILON ) denominator[i] /= numerator[i];
  }
}

void scale_double_array ( double array[], int num_elements, double total ) {
    int i;
    double small = DBL_EPSILON;
    double count;

    /* scaling by division */

/* sanity check */
    if ( total <= small ) return;
    for ( i=0,count=0.0; i<num_elements; i++ ) {
	count += array[i];
    }

    if ( count < small ) return;
    count = total/count;
    for ( i=0; i<num_elements; i++ ) {
	array[i] *= count;
    }
}

void scale_double_array1 ( double array[], int num_elements, double total ) {

  /* want sum of array elements to be total, so add them and then subtract
     sum/num_elements from each
     */
  int i;
  double sum, mean;

  for ( i=0, sum=0.0; i<num_elements; i++ ) {
    sum += array[i];
  }
  mean = sum / num_elements;
  for ( i=0, sum=0.0; i<num_elements; i++ ) {
    array[i] -= mean;
  }
}


void log_double_array ( double array[], int num_elements ) {
    int i;

    for ( i=0; i<num_elements; i++ ) {
	if ( array[i] > 0.0 ) array[i] = log ( array[i] );
    }
}

void log10_double_array ( double array[], int num_elements ) {
    int i;

    for ( i=0; i<num_elements; i++ ) {
	if ( array[i] > 0.0 ) array[i] = log10 ( array[i] );
    }
}

void exp_double_array ( double array[], int num_elements ) {
    int i;

    for ( i=0; i<num_elements; i++ ) {
	array[i] = exp ( array[i] );
    }
}

