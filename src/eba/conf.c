#include "Read.h"
#include "misc.h"

#define AVERAGE_QUAL
#define WINDOW_SIZE 15

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
static void rescale_scores(Read *r) {
    int i;
    for (i = 0; i < r->NBases; i++) {
	int conf = get_conf(r, i);
	if (conf < 0) conf = 0;
	if (conf > 100) conf = 100;
	set_conf(r, i, eba2phred[conf]);
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

float get_area(TRACE *trace, int startp, int endp)
{
    int i;
    float sum=1.0e-10;
    
    for (i=startp; i<endp; i++)
	sum += trace[i];
  
    return(sum);
}

float max_area(TRACE *tx, TRACE *ty, TRACE *tz , int stp, int endp)
{
    float x,y,z;
    x = get_area(tx,stp,endp);
    y = get_area(ty,stp,endp);
    z = get_area(tz,stp,endp);
    
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

void calc_conf_values(Read *r, int phred_scale) {
    int i;
    int pos,start_pos,end_pos;

#ifdef AVERAGE_QUAL
    char *conf_buf, val;
    double total;
    
    if (NULL == (conf_buf = (char *)xcalloc(r->NBases, 1))) {
	return;
    }
#endif

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
	start_pos = (r->basePos)[i]-((r->basePos)[i] - (r->basePos[i-1])) /2;
	end_pos   = (r->basePos)[i]+((r->basePos)[i+1] - (r->basePos[i])) /2;

	switch ((r->base)[i]) {
	case 'A':
	case 'a':
#ifdef AVERAGE_QUAL
	    conf_buf[i] =
#else
	    (r->prob_A)[i] =
#endif
		probFromQual((max_area(r->traceC, r->traceG, r->traceT,
				       start_pos, end_pos) /
			      get_area(r->traceA, start_pos, end_pos)));
	    break;

	case 'C':
	case 'c':
#ifdef AVERAGE_QUAL
	    conf_buf[i] =
#else
	    (r->prob_C)[i] =
#endif
		probFromQual((max_area(r->traceA, r->traceG, r->traceT,
				       start_pos, end_pos) /
			      get_area(r->traceC, start_pos, end_pos)));
	    break;
	    
	case 'G':
	case 'g':
#ifdef AVERAGE_QUAL
	    conf_buf[i] =
#else
	    (r->prob_G)[i] =
#endif
		probFromQual((max_area(r->traceC, r->traceA, r->traceT,
				       start_pos, end_pos) /
			      get_area(r->traceG, start_pos, end_pos)));
	    break;

	case 'T':
	case 't':
#ifdef AVERAGE_QUAL
	    conf_buf[i] =
#else
	    (r->prob_T)[i] =
#endif
		probFromQual((max_area(r->traceC, r->traceG, r->traceA,
				       start_pos, end_pos) /
			      get_area(r->traceT, start_pos, end_pos)));
	    break;

	default:
	    break;
	}
    }

    if (r->NBases > 1) {
#ifdef AVERAGE_QUAL
	conf_buf[0] = conf_buf[1];
	conf_buf[r->NBases-1] = conf_buf[r->NBases-2];
#else
	(r->prob_A)[0] = (r->prob_A)[1];
	(r->prob_C)[0] = (r->prob_C)[1];
	(r->prob_G)[0] = (r->prob_G)[1];
	(r->prob_T)[0] = (r->prob_T)[1];

	(r->prob_A)[r->NBases-1] = (r->prob_A)[r->NBases-2];
	(r->prob_C)[r->NBases-1] = (r->prob_C)[r->NBases-2];
	(r->prob_G)[r->NBases-1] = (r->prob_G)[r->NBases-2];
	(r->prob_T)[r->NBases-1] = (r->prob_T)[r->NBases-2];
#endif
    }

    /* Average confidence values */
#ifdef AVERAGE_QUAL
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
	total += conf_buf[i + WINDOW_SIZE/2+1] - conf_buf[i - WINDOW_SIZE/2];
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
    i = r->NBases - WINDOW_SIZE/2 >= 0 ? r->NBases - WINDOW_SIZE/2 : 0;
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
#endif

    if (phred_scale) {
	rescale_scores(r);
    }

    return;
}

