/*
 * Copyright (c) 1996, Whitehead Institute for Biomedical Research. All rights
 * reserved.  Please see full software use agreement in primer3_main.c or by
 * executing primer3 with -h.
 */

#include <limits.h>
#include <math.h>
#include <string.h>
#include "oligotm.h"
#include "primer3_release.h"


#define SANTALUCIA_1998
/*
 * Define this is you wish to use the newer NN parameters listed in the
 * SantaLucia paper: http://dx.doi.org/10.1073/pnas.95.4.1460
 *
 * "A unified view of polymer, dumbbell and oligonucleotide DNA
 * nearest-neighbor thermodynamics", SantaLucia JR (1998),
 * Proc Natl Acad Sci 95:1460-65.
 */


#ifndef SANTALUCIA_1998

/* 
 * Tables of nearest-neighbor thermodynamics for DNA bases.  See Breslauer,
 * Frank, Bloecker, and Markey, Proc. Natl. Acad. Sci. USA, vol 83, page 3748,
 * table 2.
 */
#define S_A_A 240
#define S_A_C 173
#define S_A_G 208
#define S_A_T 239
#define S_A_N 215
  
#define S_C_A 129
#define S_C_C 266
#define S_C_G 278
#define S_C_T 208
#define S_C_N 220  
  
#define S_G_A 135
#define S_G_C 267
#define S_G_G 266
#define S_G_T 173
#define S_G_N 210
  
#define S_T_A 169
#define S_T_C 135
#define S_T_G 129
#define S_T_T 240
#define S_T_N 168
  
#define S_N_A 168
#define S_N_C 210
#define S_N_G 220
#define S_N_T 215
#define S_N_N 203


#define H_A_A  91
#define H_A_C  65
#define H_A_G  78
#define H_A_T  86
#define H_A_N  80

#define H_C_A  58
#define H_C_C 110
#define H_C_G 119
#define H_C_T  78
#define H_C_N  91

#define H_G_A  56
#define H_G_C 111
#define H_G_G 110
#define H_G_T  65
#define H_G_N  85

#define H_T_A  60
#define H_T_C  56
#define H_T_G  58
#define H_T_T  91
#define H_T_N  66

#define H_N_A  66
#define H_N_C  85
#define H_N_G  91
#define H_N_T  80
#define H_N_N  80

/* Delta G's of disruption * 1000. */
#define G_A_A  1900
#define G_A_C  1300
#define G_A_G  1600
#define G_A_T  1500
#define G_A_N  1575

#define G_C_A  1900 
#define G_C_C  3100
#define G_C_G  3600
#define G_C_T  1600
#define G_C_N  2550

#define G_G_A  1600
#define G_G_C  3100
#define G_G_G  3100
#define G_G_T  1300
#define G_G_N  2275

#define G_T_A   900
#define G_T_C  1600
#define G_T_G  1900
#define G_T_T  1900
#define G_T_N  1575

#define G_N_A  1575
#define G_N_C  2275
#define G_N_G  2550
#define G_N_T  1575
#define G_N_N  1994

#else /* SANTALUCIA_1998 */

/*-----------------------------------------------------------------------------
 * Update using http://dx.doi.org/10.1073/pnas.95.4.1460
 * SantaLucia (1998) paper.
 */

/* Delta H and delta S from table 2. */
#define S_A_A 222
#define S_A_C 224
#define S_A_G 210
#define S_A_T 204
#define S_A_N 215
  
#define S_C_A 227
#define S_C_C 199
#define S_C_G 272
#define S_C_T 210
#define S_C_N 227
  
#define S_G_A 222
#define S_G_C 244
#define S_G_G 199
#define S_G_T 224
#define S_G_N 222
  
#define S_T_A 213
#define S_T_C 222
#define S_T_G 227
#define S_T_T 222
#define S_T_N 221
  
#define S_N_A 221
#define S_N_C 222
#define S_N_G 227
#define S_N_T 215
#define S_N_N 221

#define S_TERM_GC  28
#define S_TERM_AT -41
#define S_TERM_N   -6

#define S_SYM      14

#define H_A_A  79
#define H_A_C  84
#define H_A_G  78
#define H_A_T  72
#define H_A_N  78

#define H_C_A  85
#define H_C_C  80
#define H_C_G 106
#define H_C_T  78
#define H_C_N  87

#define H_G_A  82
#define H_G_C  98
#define H_G_G  80
#define H_G_T  84
#define H_G_N  86

#define H_T_A  72
#define H_T_C  82
#define H_T_G  85
#define H_T_T  79
#define H_T_N  79

#define H_N_A  79
#define H_N_C  86
#define H_N_G  87
#define H_N_T  78
#define H_N_N  83

#define H_TERM_GC  -1
#define H_TERM_AT -23
#define H_TERM_N  -12

#define H_SYM      0

/*
 * Unified (ref 22) column of table 1.
 */
/* Delta G's of disruption * 1000. */
#define G_A_A  1000
#define G_A_C  1440
#define G_A_G  1280
#define G_A_T   880
#define G_A_N  1150 /* Differs with p3-2.2.1 - they use min value */

#define G_C_A  1450
#define G_C_C  1840
#define G_C_G  2170
#define G_C_T  1280
#define G_C_N  1685 /* differs */

#define G_G_A  1300
#define G_G_C  2240
#define G_G_G  1840
#define G_G_T  1440
#define G_G_N  1705 /* differs */

#define G_T_A   580
#define G_T_C  1300
#define G_T_G  1450
#define G_T_T  1000
#define G_T_N  1082 /* differs */

#define G_N_A  1082 /* these 4 differ too */
#define G_N_C  1705
#define G_N_G  1685
#define G_N_T  1150
#define G_N_N  1405

#define G_TERM_GC  -980
#define G_TERM_AT -1003
#define G_TERM_N   -991

#endif /* SANTALUCIA_1998 */

#define A_CHAR 'A'
#define G_CHAR 'G'
#define T_CHAR 'T'
#define C_CHAR 'C'
#define N_CHAR 'N'

#define CATID5(A,B,C,D,E) A##B##C##D##E
#define CATID2(A,B) A##B
#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
     dh += CATID5(H,_,LAST,_,THIS); \
     ds += CATID5(S,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE);      \
  }

#define STATE(LAST)     \
   CATID2(LAST,_STATE): \
   c = *s; s++;         \
   DO_PAIR(LAST,A)      \
   else DO_PAIR(LAST,T) \
   else DO_PAIR(LAST,G) \
   else DO_PAIR(LAST,C) \
   else DO_PAIR(LAST,N) \
   else if ('\0' == c)  \
             goto DONE; \
   else goto ERROR \

/*
 * Query - should AAAACTTTT be considered symmetric? If not then
 * it implies all odd length strings are not. (We treat this as symmetric.)
 *
 * Returns 1 if 's' is a self complement.
 *         0 if not.
 */
static int is_sym(const char *s) {
    const char *e = s + strlen(s)-1;
    int comp[256];
    comp['A'] = 'T';
    comp['C'] = 'G';
    comp['G'] = 'C';
    comp['T'] = 'A';

    while (s < e) {
	if (comp[(const unsigned char)*s] != *e)
	    return 0;
	s++;
	e--;
    }

    return 1;
}

double 
oligotm(s, DNA_nM, K_mM, Mg_mM, dNTP_mM)
     const  char *s;
     double DNA_nM;
     double K_mM;
     double Mg_mM;
     double dNTP_mM;
{
    register int dh = 0, ds = 0;
    register char c;
    double delta_H, delta_S, Ct, salt;
    size_t len = strlen(s);
    int symmetric;

    /* const char *orig=s; */

#ifndef SANTALUCIA_1998
    ds += 108;
#else
    /* SantaLucia method */

    /* Terminal AT/GC scoring */
    if (*s == 'A' || *s == 'T') {
	ds += S_TERM_AT;
	dh += H_TERM_AT;
    } else if (*s == 'G' || *s == 'C') {
	ds += S_TERM_GC;
	dh += H_TERM_GC;
    } else {
	ds += S_TERM_N;
	dh += H_TERM_N;
    }

    if (s[len-1] == 'A' || s[len-1] == 'T') {
	ds += S_TERM_AT;
	dh += H_TERM_AT;
    } else if (s[len-1] == 'G' || s[len-1] == 'C') {
	ds += S_TERM_GC;
	dh += H_TERM_GC;
    } else {
	/* guess, pick avg */
	ds += S_TERM_N;
	dh += H_TERM_N;
    }

    /* Symmetry adjustment */
    if ((symmetric = is_sym(s))) {
	ds += S_SYM;
	dh += H_SYM;
    }
#endif /* SANTALUCIA_1998 */

    /* Use a finite-state machine (DFA) to calucluate dh and ds for s. */
    c = *s; s++;
    if (c == 'A') goto A_STATE;
    else if (c == 'G') goto G_STATE;
    else if (c == 'T') goto T_STATE;
    else if (c == 'C') goto C_STATE;
    else if (c == 'N') goto N_STATE;
    else goto ERROR;
    STATE(A);
    STATE(T);
    STATE(G);
    STATE(C);
    STATE(N);

 DONE:  /* dh and ds are now computed for the given sequence. */
    delta_H = dh * -100.0;  /* 
			     * Nearest-neighbor thermodynamic values for dh
			     * are given in 100 cal/mol of interaction.
			     */
    delta_S = ds * -0.1;     /*
			      * Nearest-neighbor thermodynamic values for ds
			      * are in in .1 cal/K per mol of interaction.
			      */

#ifndef SANTALUCIA_1998
    /* 
     * See Rychlik, Spencer, Rhoads, Nucleic Acids Research, vol 18, no 21,
     * page 6410, eqn (ii).
     */
    return delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0))
    	- 273.15 + 16.6 * log10(K_mM/1000.0);
    
#else
    /* SantaLucia salt concentrations: see equation 8:
     * dS(oligomer) = dS + 0.368 * N * ln(Na+)
     * dG(oligomer) = dG - 0.114 * N * ln(Na+)
     * N is number of NN pairs (ie length-1).
     */

    /* von Ahsen et al, "Oligonucleotide Melting Temperatures under
     * PCR Conditions: Nearest-Neighbor Corrections for Mg2+,
     * Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide
     * Concentrations with Comparison to Alternative Empirical Formulas"
     * Clinical Chemistry 47: 1956-1961, 2001; 
     *
     * Suggests non-linear equiv of 120*sqrt(Mg_mM)
     *
     * eg (units correct?)
     *
     *   salt = K_mM + 120 * sqrt(Mg_mM - dNTP_mM);
     *   delta_S += 0.368 * (len-1) * log(salt/1000.0); 
     *
     * Paper also mentions older work demonstrating that Mg_mM may
     * have a 140 fold effect compared to K.ie:
     *
     *   delta_S += 0.368 * (len-1) * log(K_mM/1000.0 + 140*Mg_mM/1000.0);
     */

    if (Mg_mM - dNTP_mM >= 0)
	salt = K_mM + 120 * sqrt(Mg_mM - dNTP_mM);
        //salt = K_mM + 140 * Mg_mM;
    else
	salt = K_mM;

    delta_S += 0.368 * (len-1) * log(salt/1000.0);

    /* Equation 3 */
    Ct = log(DNA_nM / (symmetric ? 1000000000. : 4000000000.));
    /* printf("%s dh=%f ds=%f tm=%f\n", orig, delta_H, delta_S,
     *        delta_H / (delta_S + 1.987 * Ct) - 273.15);
     */
    return delta_H / (delta_S + 1.987 * Ct) - 273.15;
#endif

 ERROR:  /* 
	  * length of s was less than 2 or there was an illegal character in
	  * s.
	  */
    return OLIGOTM_ERROR;
}
#undef DO_PAIR

#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
     dg += CATID5(G,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE);      \
  }

double 
oligodg(s)
    const char *s;       /* The sequence. */
{
    register int dg = 0;
    register char c;

    /* Use a finite-state machine (DFA) to calucluate dg s. */
    c = *s; s++;

#ifdef SANTALUCIA_1998
    /* Terminal AT/GC scoring */
    if (c == 'A' || c == 'T')
	dg += G_TERM_AT;
    else if (c == 'G' || c == 'C')
	dg += G_TERM_GC;
    else
	dg += G_TERM_N;
#endif

    if (c == 'A') goto A_STATE;
    else if (c == 'G') goto G_STATE;
    else if (c == 'T') goto T_STATE;
    else if (c == 'C') goto C_STATE;
    else if (c == 'N') goto N_STATE;
    else goto ERROR;
    STATE(A);
    STATE(T);
    STATE(G);
    STATE(C);
    STATE(N);

#ifdef SANTALUCIA_1998
    /* Terminal AT/GC scoring */
    c = s[-2]; /* last  base */
    if (c == 'A' || c == 'T')
	dg += G_TERM_AT;
    else if (c == 'G' || c == 'C')
	dg += G_TERM_GC;
    else
	dg += G_TERM_N;
#endif

 DONE:  /* dg is now computed for the given sequence. */
    return dg / 1000.0;

 ERROR:  /* 
	  * length of s was less than 2 or there was an illegal character in
	  * s.
	  */
    return OLIGOTM_ERROR;
}

double end_oligodg(s, len)
  const char *s;
  int len; /* The number of characters to return. */
{
  int x = strlen(s);
  return x < len ? oligodg(s) : oligodg(s + (x - len));
}

double seqtm(seq, dna_conc, salt_conc, Mg_conc, dNTP_conc, nn_max_len)
  const  char *seq;
  double dna_conc;
  double salt_conc;
  double Mg_conc;
  double dNTP_conc;
  int    nn_max_len;
{
  int len = strlen(seq);
  return (len > nn_max_len)
    ? long_seq_tm(seq, 0, len, salt_conc) : oligotm(seq, dna_conc, salt_conc,
						    Mg_conc, dNTP_conc);
}

/* See oligotm.h for documentation on this function and the formula it
   uses. */
double
long_seq_tm(s, start, len, salt_conc)
  const char *s;
  int start, len;
  double salt_conc;
{
  int GC_count = 0;
  const char *p, *end = &s[len];
  if (len <= 0) return OLIGOTM_ERROR;
  /* Length <= 0 is nonsensical. */
  for (p = &s[0]; p < end; p++) {
    if ('G' == *p || 'g' == *p || 'C' == *p || 'c' == *p)
      GC_count++;
  }

  return
    81.5 
    + (16.6 * log10(salt_conc / 1000.0))
    + (41.0 * (((double) GC_count) / len))
    - (600.0 / len);

}
