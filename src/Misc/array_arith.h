
/* general array stuff */

void fill_int_array ( int array[], int num_elements, int filler );

void fill_double_array ( double array[], int num_elements, double filler );

int max_int_array ( int array[], int num_elements );

int min_int_array ( int array[], int num_elements );

double max_double_array ( double array[], int num_elements );

double min_double_array ( double array[], int num_elements );

void reset_zeroes ( double table[], int num_elements, double correction );

double sum_double_array ( double array[], int num_elements );

void div_double_array ( double array[], int num_elements, double divisor );

void scale_double_array ( double array[], int num_elements, double total );

void scale_double_array1 ( double array[], int num_elements, double total );

void ratio_double_arrays ( double denominator[], double numerator[], int num_elements );

void mult_double_array ( double array[], int num_elements, double mult );

void log_double_array ( double array[], int num_elements );
