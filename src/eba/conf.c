#include "Read.h"
#include "misc.h"
#include "conf.h"

#define WINDOW_SIZE 15
/*#define WINDOW_SIZE 9*/

static int eba2phred[] = {
    /*  0 */ 4,
    /*  1 */ 5,
    /*  2 */ 5,
    /*  3 */ 5,
    /*  4 */ 5,
    /*  5 */ 6,
    /*  6 */ 6,
    /*  7 */ 6,
    /*  8 */ 6,
    /*  9 */ 7,
    /* 10 */ 7,
    /* 11 */ 7,
    /* 12 */ 7,
    /* 13 */ 7,
    /* 14 */ 7,
    /* 15 */ 7,
    /* 16 */ 7,
    /* 17 */ 8,
    /* 18 */ 8,
    /* 19 */ 8,
    /* 20 */ 8,
    /* 21 */ 8,
    /* 22 */ 8,
    /* 23 */ 8,
    /* 24 */ 8,
    /* 25 */ 8,
    /* 26 */ 8,
    /* 27 */ 8,
    /* 28 */ 8,
    /* 29 */ 8,
    /* 30 */ 8,
    /* 31 */ 8,
    /* 32 */ 8,
    /* 33 */ 8,
    /* 34 */ 8,
    /* 35 */ 8,
    /* 36 */ 9,
    /* 37 */ 9,
    /* 38 */ 9,
    /* 39 */ 9,
    /* 40 */ 9,
    /* 41 */ 9,
    /* 42 */ 9,
    /* 43 */ 9,
    /* 44 */ 10,
    /* 45 */ 10,
    /* 46 */ 10,
    /* 47 */ 10,
    /* 48 */ 10,
    /* 49 */ 11,
    /* 50 */ 11,
    /* 51 */ 11,
    /* 52 */ 11,
    /* 53 */ 12,
    /* 54 */ 12,
    /* 55 */ 12,
    /* 56 */ 13,
    /* 57 */ 13,
    /* 58 */ 14,
    /* 59 */ 14,
    /* 60 */ 15,
    /* 61 */ 15,
    /* 62 */ 16,
    /* 63 */ 16,
    /* 64 */ 17,
    /* 65 */ 17,
    /* 66 */ 18,
    /* 67 */ 19,
    /* 68 */ 19,
    /* 69 */ 20,
    /* 70 */ 21,
    /* 71 */ 22,
    /* 72 */ 22,
    /* 73 */ 23,
    /* 74 */ 24,
    /* 75 */ 25,
    /* 76 */ 26,
    /* 77 */ 27,
    /* 78 */ 29,
    /* 79 */ 30,
    /* 80 */ 31,
    /* 81 */ 32,
    /* 82 */ 33,
    /* 83 */ 34,
    /* 84 */ 36,
    /* 85 */ 37,
    /* 86 */ 38,
    /* 87 */ 39,
    /* 88 */ 40,
    /* 89 */ 41,
    /* 90 */ 42,
    /* 91 */ 43,
    /* 92 */ 44,
    /* 93 */ 45,
    /* 94 */ 45,
    /* 95 */ 46,
    /* 96 */ 48,
    /* 97 */ 50,
    /* 98 */ 52,
    /* 99 */ 43,
    /*100 */ 43
};

/*
 * Calibration from -filtered 1 -non_filtered 0 -average 0 -offset 20.
 */
static int ebanew2phred[] = {
    /*   0 */ 5, 
    /*   1 */ 5, 
    /*   2 */ 5, 
    /*   3 */ 6, 
    /*   4 */ 6, 
    /*   5 */ 6, 
    /*   6 */ 6, 
    /*   7 */ 6, 
    /*   8 */ 6, 
    /*   9 */ 6, 
    /*  10 */ 6, 
    /*  11 */ 6, 
    /*  12 */ 6, 
    /*  13 */ 6, 
    /*  14 */ 7, 
    /*  15 */ 7, 
    /*  16 */ 7, 
    /*  17 */ 7, 
    /*  18 */ 7, 
    /*  19 */ 7, 
    /*  20 */ 7, 
    /*  21 */ 7, 
    /*  22 */ 7, 
    /*  23 */ 7, 
    /*  24 */ 7, 
    /*  25 */ 7, 
    /*  26 */ 7, 
    /*  27 */ 7, 
    /*  28 */ 8, 
    /*  29 */ 8, 
    /*  30 */ 8, 
    /*  31 */ 8, 
    /*  32 */ 8, 
    /*  33 */ 8, 
    /*  34 */ 8, 
    /*  35 */ 8, 
    /*  36 */ 8, 
    /*  37 */ 8, 
    /*  38 */ 8, 
    /*  39 */ 8, 
    /*  40 */ 8, 
    /*  41 */ 8, 
    /*  42 */ 8, 
    /*  43 */ 9, 
    /*  44 */ 9, 
    /*  45 */ 9, 
    /*  46 */ 9, 
    /*  47 */ 9, 
    /*  48 */ 9, 
    /*  49 */ 9, 
    /*  50 */ 9, 
    /*  51 */ 9, 
    /*  52 */ 10,
    /*  53 */ 10,
    /*  54 */ 10,
    /*  55 */ 10,
    /*  56 */ 10,
    /*  57 */ 10,
    /*  58 */ 10,
    /*  59 */ 10,
    /*  60 */ 10,
    /*  61 */ 11,
    /*  62 */ 11,
    /*  63 */ 11,
    /*  64 */ 11,
    /*  65 */ 11,
    /*  66 */ 11,
    /*  67 */ 11,
    /*  68 */ 12,
    /*  69 */ 12,
    /*  70 */ 12,
    /*  71 */ 13,
    /*  72 */ 13,
    /*  73 */ 13,
    /*  74 */ 13,
    /*  75 */ 14,
    /*  76 */ 14,
    /*  77 */ 14,
    /*  78 */ 15,
    /*  79 */ 15,
    /*  80 */ 15,
    /*  81 */ 16,
    /*  82 */ 17,
    /*  83 */ 17,
    /*  84 */ 17,
    /*  85 */ 17,
    /*  86 */ 18,
    /*  87 */ 19,
    /*  88 */ 19,
    /*  89 */ 20,
    /*  90 */ 20,
    /*  91 */ 21,
    /*  92 */ 21,
    /*  93 */ 22,
    /*  94 */ 23,
    /*  95 */ 24,
    /*  96 */ 25,
    /*  97 */ 28,
    /*  98 */ 31,
    /*  99 */ 34,
    /* 100 */ 34
};

static int get_conf(Read *r, int i) { 
    int conf; 
 
    switch (r->base[i]) { 
    case 'A': 
    case 'a': 
        conf = r->prob_A[i]; 
        break; 
         
    case 'C': 
    case 'c': 
        conf = r->prob_C[i]; 
        break; 
         
    case 'G': 
    case 'g': 
        conf = r->prob_G[i]; 
        break; 
 
    case 'T': 
    case 't': 
        conf = r->prob_T[i]; 
        break; 
 
    default: 
        conf = (r->prob_A[i] + r->prob_C[i] + r->prob_G[i] + r->prob_T[i]) / 4; 
    } 
 
    return conf; 
} 

static void set_conf(Read *r, int i, int val) { 
    switch (r->base[i]) { 
    case 'A': 
    case 'a': 
        r->prob_A[i] = val; 
        break; 
         
    case 'C': 
    case 'c': 
        r->prob_C[i] = val; 
        break; 
         
    case 'G': 
    case 'g': 
        r->prob_G[i] = val; 
        break; 
 
    case 'T': 
    case 't': 
        r->prob_T[i] = val; 
        break; 
 
    default: 
	r->prob_A[i] = val;
	r->prob_C[i] = val;
	r->prob_G[i] = val;
	r->prob_T[i] = val;
    } 
} 

/*
 * Rescales eba scores to phred-style scores.
 */
void rescale_scores(Read *r, int table) {
    int i;
    for (i = 0; i < r->NBases; i++) {
	int conf = get_conf(r, i);
	if (conf < 0) conf = 0;
	if (conf > 100) conf = 100;
	set_conf(r, i, table == 0 ? eba2phred[conf] : ebanew2phred[conf]);
    }
}

/*
 * MODULE    SeqQual12
 *
 * area of the called divided by the area of the max non-called
 * for range from 1/2 between this base and the lastbase to 1/2
 * between this base and the next base
 * 
 * So find the area for each base....divide area of the called
 * by other max area
 */

float get_area(TRACE *trace, int startp, int endp, int offset)
{
    int i;
    float sum=1.0e-10;
    
    for (i=startp; i<endp; i++)
	sum += trace[i];
  
    return(sum + offset);
}

float max_area(TRACE *tx, TRACE *ty, TRACE *tz , int stp, int endp, int offset)
{
    float x,y,z;
    x = get_area(tx,stp,endp,offset);
    y = get_area(ty,stp,endp,offset);
    z = get_area(tz,stp,endp,offset);
    
    if (x > y) {
	if (z > x) {
	    return (z);
	} else {
	    return (x);
	}
    } else {
	if (z > y) {
	    return (z);
	} else {
	    return (y);
	}
    }
}

/*
 * Like get_area, but we use just a 5-wide peak and apply a raised
 * cosine profile to work our which samples have the highest significance.
 */
float get_cosa(TRACE *trace, int pos, int offset)
{
    int i;
    float sum=1.0e-10;
    /* raised cosine of N*pi/5 where N<-1..4 */
    double rcos[5] = {1.000, .9045, .6545, .3455, .0955};
    
    sum += trace[pos];
    for (i = 1; i < 5; i++) {
	sum += rcos[i]*trace[pos+i];
	sum += rcos[i]*trace[pos-i];
    }

    return(sum + offset);
}

float max_cosa(TRACE *tx, TRACE *ty, TRACE *tz, int pos, int offset)
{
    float x,y,z;
    x = get_cosa(tx,pos,offset);
    y = get_cosa(ty,pos,offset);
    z = get_cosa(tz,pos,offset);
    
    if (x > y) {
	if (z > x) {
	    return (z);
	} else {
	    return (x);
	}
    } else {
	if (z > y) {
	    return (z);
	} else {
	    return (y);
	}
    }
}

char probFromQual(float qual)
{
    return (char)(100.0*(1.0 - min(qual,1.0)));
}


/*
 * Applies a simple filter to the trace data, producing new traces.
 * Equivalent to 1st-order differentiation followed by fitting a +-1 square
 * wave of width 11.
 */
void filter_trace(Read *r) {
    int *tmp, i, j;
    
    tmp = (int *)calloc(r->NPoints, sizeof(int));
    for (j = 0; j < 4; j++) {
	TRACE *t;
	switch (j) {
	case 0:
	    t = r->traceA;
	    break;
	case 1:
	    t = r->traceC;
	    break;
	case 2:
	    t = r->traceG;
	    break;
	case 3:
	    t = r->traceT;
	    break;
	}

	for (i = 5; i < r->NPoints-5; i++)
	    tmp[i] = -t[i-5] + 2*t[i] - t[i+5];
	for (i = 5; i < r->NPoints-5; i++)
	    t[i] = tmp[i] > 0 ? tmp[i] : 0;
    }
    xfree(tmp);
}


/*
 * Averages the confidence arrays in r over a window of length
 * WINDOW_SIZE (#define at top).
 */
int average_conf(Read *r) {
    char val, *conf_buf;
    double total;
    int i;
    
    /* Create a temp. copy of the confidence values */
    if (NULL == (conf_buf = (char *)xcalloc(r->NBases, 1))) {
	return -1;
    }

    for (i = 0; i < r->NBases; i++) {
	switch((r->base)[i]) {
	case 'A': case 'a':
	    conf_buf[i] = r->prob_A[i];
	    break;

	case 'C': case 'c':
	    conf_buf[i] = r->prob_C[i];
	    break;

	case 'G': case 'g':
	    conf_buf[i] = r->prob_G[i];
	    break;

	case 'T': case 't':
	    conf_buf[i] = r->prob_T[i];
	    break;

	default:
	    conf_buf[i] = 0;
	}
    }

    /*
     * Calculate average over a specified window size.
     * For the ends of the array we cheat and copy the nearest average.
     * It's easy to do and this is only temporary anyway.
     */
	
    /* Initialise total */
    total = 0;
    for (i = 0; i < WINDOW_SIZE && i < r->NBases; i++) {
	total += conf_buf[i];
    }
    
    /* Loop around for inner section of quality buffer */
    for (i = WINDOW_SIZE/2; i < r->NBases - WINDOW_SIZE/2 - 1; i++) {
	switch((r->base)[i]) {
	case 'A':
	case 'a':
	    (r->prob_A)[i] = total / WINDOW_SIZE;
	    break;
		
	case 'C':
	case 'c':
	    (r->prob_C)[i] = total / WINDOW_SIZE;
	    break;
		
	case 'G':
	case 'g':
	    (r->prob_G)[i] = total / WINDOW_SIZE;
	    break;

	case 'T':
	case 't':
	    (r->prob_T)[i] = total / WINDOW_SIZE;
	    break;
	}
	total += conf_buf[i + WINDOW_SIZE/2+1]
	    - conf_buf[i - WINDOW_SIZE/2];
    }

    /* do left end - extend from first done */
    i = WINDOW_SIZE/2+1 < r->NBases
	? WINDOW_SIZE/2+1 : (r->NBases > 0 ? r->NBases-1 : 0);
    val = conf_buf[i];
    for (i = 0; i < WINDOW_SIZE/2 && i < r->NBases; i++) {
	switch ((r->base)[i]) {
	case 'A':
	case 'a':
	    (r->prob_A)[i] = val;
	    break;
	case 'C':
	case 'c':
	    (r->prob_C)[i] = val;
	    break;
	case 'G':
	case 'g':
	    (r->prob_G)[i] = val;
	    break;
	case 'T':
	case 't':
	    (r->prob_T)[i] = val;
	    break;
	}
    }
    

    /* do right end - extend from last done */
    i = r->NBases - WINDOW_SIZE/2 - 2 >= 0 ? r->NBases - WINDOW_SIZE/2 - 2: 0;
    val = conf_buf[i];
    for (; i < r->NBases; i++) {
	switch ((r->base)[i]) {
	case 'A':
	case 'a':
	    (r->prob_A)[i] = val;
	    break;
	case 'C':
	case 'c':
	    (r->prob_C)[i] = val;
	    break;
	case 'G':
	case 'g':
	    (r->prob_G)[i] = val;
	    break;
	case 'T':
	case 't':
	    (r->prob_T)[i] = val;
	    break;
	}
    }
    
    xfree(conf_buf);

    return 0;
}


void calc_conf_values(Read *r, int phred_scale, int cosa, int offset) {
    int i;
    int pos,start_pos,end_pos;

    /*
     * Set confidence values for first and last bases to zero to simplify
     * the loop. We reset them later on so it's not too vital.
     */
    
    (r->prob_A)[0] = (r->prob_C)[0] = (r->prob_G)[0] = (r->prob_T)[0] = 0;
    if (r->NBases > 1)  {
	i = r->NBases - 1;
	(r->prob_A)[i] = (r->prob_C)[i] = (r->prob_G)[i] = (r->prob_T)[i] = 0;
    }
    
    
    for (i = 1; i < r->NBases-1; i++) {

	/* Clear the probability arrays */
	(r->prob_A)[i] = (r->prob_C)[i] =
	    (r->prob_G)[i] = (r->prob_T)[i] = 0;
	pos = (r->basePos)[i];

	if (pos >= r->NPoints-5 || pos < 5)
	    continue;

	start_pos = (r->basePos)[i]-((r->basePos)[i] - (r->basePos[i-1])) /2;
	end_pos   = (r->basePos)[i]+((r->basePos)[i+1] - (r->basePos[i])) /2;

	switch ((r->base)[i]) {
	case 'A':
	case 'a':
	    (r->prob_A)[i] = cosa
		? probFromQual((max_cosa(r->traceC, r->traceG, r->traceT, pos,
					 offset)/
				get_cosa(r->traceA, pos, offset)))
		: probFromQual((max_area(r->traceC, r->traceG, r->traceT,
					 start_pos, end_pos, offset) /
				get_area(r->traceA, start_pos, end_pos,
					 offset)));
	    break;

	case 'C':
	case 'c':
	    (r->prob_C)[i] = cosa
		? probFromQual((max_cosa(r->traceA, r->traceG, r->traceT, pos,
					 offset)/
				get_cosa(r->traceC, pos, offset)))
		: probFromQual((max_area(r->traceA, r->traceG, r->traceT,
					 start_pos, end_pos, offset) /
				get_area(r->traceC, start_pos, end_pos,
					 offset)));
	    break;
	    
	case 'G':
	case 'g':
	    (r->prob_G)[i] = cosa
		? probFromQual((max_cosa(r->traceC, r->traceA, r->traceT, pos,
					 offset)/
				get_cosa(r->traceG, pos, offset)))
		: probFromQual((max_area(r->traceC, r->traceA, r->traceT,
					 start_pos, end_pos, offset) /
				get_area(r->traceG, start_pos, end_pos,
					 offset)));
	    break;

	case 'T':
	case 't':
	    (r->prob_T)[i] = cosa
		? probFromQual((max_cosa(r->traceC, r->traceG, r->traceA, pos,
					 offset)/
				get_cosa(r->traceT, pos, offset)))
		: probFromQual((max_area(r->traceC, r->traceG, r->traceA,
					 start_pos, end_pos, offset) /
				get_area(r->traceT, start_pos, end_pos,
					 offset)));
	    break;

	default:
	    break;
	}
    }

    if (r->NBases > 1) {
	(r->prob_A)[0] = (r->prob_A)[1];
	(r->prob_C)[0] = (r->prob_C)[1];
	(r->prob_G)[0] = (r->prob_G)[1];
	(r->prob_T)[0] = (r->prob_T)[1];

	(r->prob_A)[r->NBases-1] = (r->prob_A)[r->NBases-2];
	(r->prob_C)[r->NBases-1] = (r->prob_C)[r->NBases-2];
	(r->prob_G)[r->NBases-1] = (r->prob_G)[r->NBases-2];
	(r->prob_T)[r->NBases-1] = (r->prob_T)[r->NBases-2];
    }

    return;
}

