int Plot_Base_Comp(int win_len, 
		   int *score,
		   char *seq,
		   int seq_len,
		   int start,
		   int end,
		   int *match,
		   int *max_match);

int 
get_base_comp_res(char seq[], int seq_length, 
		  int window_length, /* sum over all windows this length */
		  int user_start,    /* seq start numbering from 1 */
		  int user_end,      /* seq end numbering from 1 */
		  double score[],    /* score for each base type */
		  double result[],   /* put results here */
		  double *min,       /* min result */
		  double *max);      /* max result */

double get_base_comp_mass(int a, int c, int g, int t);

void get_aa_comp (char *seq, int seq_length, double aa_comp[25]);
void get_aa_comp_mass(double *aa_comp, double *mass);

