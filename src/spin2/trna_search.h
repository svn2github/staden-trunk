#ifndef _TRNA_SEARCH_
#define _TRNA_SEARCH_

#define MAX_TRNA 100

typedef struct _TrnaRes {
    char *seq;
    int seq_length;
    int aa_right;
    int aa_left;
    int ac_left;
    int ac_right;
    int tu_right;
    int tu_left;
    int intron_length;
    int aa_score;
    int ac_score;
    int tu_score;
    int d_score;
    int total_bp_score;
    int total_cb_score;
} TrnaRes;

typedef struct _TrnaSpec {
    int max_trna_length;
    int min_trna_length;
    int max_intron_length;
    int min_intron_length;
    int max_tu_loop_length;
    int min_tu_loop_length;
    int min_aa_to_ac_length;
    int max_aa_to_ac_length;
    int min_acs_to_ace_length;
    int max_var_loop_length;
    int min_aa_score;
    int min_ac_score;
    int min_tu_score;
    int min_d_score;
    int min_total_bp_score;
    int min_total_cb_score;
    int con_bases1[18];
    int con_bases2[18];
    int con_base_offsets[18];
    int con_base_scores[18];
} TrnaSpec;


int trna_search ( char seq[], int seq_length, int user_start, int user_end,
		  TrnaRes ***results, int *nmatch, int *max_total_bp_score, 
		  TrnaSpec **t);

void draw_trna ( TrnaRes *r );

#endif
