/*
 * This is a modified algorithm, from Huang, originally from Myers and Miller.
 *
 * Our changes are simply to massage into a format suitable for external use.
 * We also changed the C types used to FastInt - a typedef defined elsewhere
 * to be an fast integer type appropriate to the local system. This change was
 * suggested by DEC.
 *
 * The original comments follow, although not all of these will be relevant
 * now.
 */

/*  A GLOBAL ALIGNMENT PROGRAM (GAP):

    copyright (c) 1992 Xiaoqiu Huang
    The distribution of the program is granted provided no charge is made
    and the copyright notice is included.
    E-mail: huang@cs.mtu.edu

     Proper attribution of the author as the source of the software would
     be appreciated: "On global sequence alignment" (to appear in CABIOS).
	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931

    The GAP program computes a global alignment of two sequences
    without penalizing terminal gaps. It delivers the alignment in
    linear space, so long sequences can be aligned. 

    Users supply scoring parameters. In the simplest form, users just
    provide 3 integers: ms, q and r, where ms is the score of a mismatch
    and the score of an i-symbol indel is -(q + r * i). Each match
    automatically receives score 10. This simple scoring scheme may be
    used for DNA sequences. NOTE: all scores are integers.

    In general, users can define an alphabet of characters appearing
    in the sequences and a matrix that gives the substitution score
    for each pair of symbols in the alphabet. The 127 ASCII characters
    are eligible. The alphabet and matrix are given in a file, where
    the first line lists the characters in the alphabet and the lower
    triangle of the matrix comes next. An example file looks as follows:

    ARNDC	       
     13
    -15  19
    -10 -22  11
    -20 -10 -20  18
    -10 -20 -10 -20  12

    Here the -22 at position (3,2) is the score of replacing N by R.
    This general scoring scheme is useful for protein sequences where the
    set of protein characters and Dayhoff matrix are specified in the file.

    The GAP program is written in C and runs under Unix systems on
    Sun workstations and under DOS systems on PCs.
    We think that the program is portable to many machines.

    Sequences to be analyzed are stored in separate files.
    An input file contains all characters of a sequence, separated by
    newline characters, in linear order. No other characters are allowed.
    Since upper case and lower case characters are different, use the same
    case consistently. A sample sequence file of 4 lines is shown below.

GAATTCTAATCTCCCTCTCAACCCTACAGTCACCCATTTGGTATATTAAA
GATGTGTTGTCTACTGTCTAGTATCCCTCAAGTAGTGTCAGGAATTAGTC
ATTTAAATAGTCTGCAAGCCAGGAGTGGTGGCTCATGTCTGTAATTCCAG
CACTGGAGAGGTAGAAGTG

    To find the best alignment of two sequences in files A and B,
    use a command of form

	   gap  A  B  gs  ms  q  r > result

    where gap is the name of the object code, gs is the minimum length
    of any gap in the short sequence receiving a constant gap penalty,
    ms is a negative integer specifying mismatch weight, q and r are
    non-negative integers specifying gap-open and gap-extend penalties,
    respectively. Output alignment is saved in the file "result".

    For using a scoring matrix defined in file S, use a command of form

	   gap  A  B  gs  S  q  r > result

    Note that ms is replaced by the file S.

    Acknowledgments
    The functions diff2() and display() were originally written by Gene Myers.
    We made the following modifications: similarity weights (integer), instead of
    distance weights (float), are used, terminal gaps are not penalized, and
    any gap of length at least gs in the short sequence is given a constant
    penalty.
*/

#include <stdio.h>
#include "align.h"

#define NMAX 6000
static FastInt *CC, *DD;	/* saving matrix scores */
static FastInt *RR, *SS;	/* saving start-points */


static FastInt (*v)[128];				/* v = W */
static FastInt q, r;       /* gap penalties */
static FastInt qr;         /* qr = q + r */
static FastInt gaplen;     /* minimum length for constant-cost insertion */
static FastInt pay;	/* constant-cost for long insertion */

static int  zero = 0;				/* int type zero        */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

#define gap2(k)  ((k) <= 0 ? 0 : ((k) <= gaplen ? q+r*(k) : pay))
/* k-symbol insertion score */

static FastInt *sapp;				/* Current script append ptr */
static FastInt last;				/* Last script op appended */

static FastInt no_mat; 				/* number of matches */ 
static FastInt no_mis; 				/* number of mismatches */ 
static FastInt al_len; 				/* length of alignment */
						/* Append "Delete k" op */
#define DEL(k)				\
{ al_len += k;				\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ al_len += k;				\
  if (last > 0)				\
    last = sapp[-1] += (k);		\
  else					\
    last = *sapp++ = (k);		\
}
						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
  al_len += 1;				\
}

static int align(A,B,M,N,tb,te,sc,sr,ec,er)
unsigned char *A, *B; FastInt M, N, tb, te, sc, sr, ec, er;

{ FastInt midi, midj, type;	/* Midpoint, type, and cost */
  FastInt midc;
  FastInt ss,cc;

{ register FastInt i, j;
  register FastInt c, e, d, s;
           FastInt t;
           FastInt *va;
	   FastInt g, temp;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      if ( !sc || !ec )
	return 0;
      else
        return - gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
          if ( !sr || !er )
    	    return 0;
          else
            return - gap2(N);
        }
      midc = - (sc * (tb + r) + er * gap2(N) );
      midj = -1;
      if ( midc < ( c =  - (ec * (te + r) + sr * gap2(N) ) ) )
	{ midc = c;
	  midj = 0;
	}
      va = v[A[1]];
      for (j = 1; j <= N; j++)
	{ c = va[B[j]] - ( sr * gap2(j-1) + er * gap2(N-j) );
          if (c > midc)
           { midc = c;
             midj = j;
           }
	}
      if (midj == -1)
        { DEL(1) INS(N) }
      else
      if (midj == 0)
        { INS(N) DEL(1) }
      else
        { if (midj > 1) INS(midj-1)
          REP
	  if ( A[1] == B[midj] )
	     no_mat += 1;
	  else
	     no_mis += 1;
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = - q * sr;
  if ( N <= gaplen )
    for (j = 1; j <= N; j++)
      { CC[j] = t = (t-r) * sr;
        DD[j] = t-q;
      }
  else
   { for (j = 1; j <= gaplen; j++)
      { CC[j] = t = (t-r) * sr;
        DD[j] = t-q;
      }
     for (j = gaplen+1; j <= N; j++)
      { CC[j] = t = -pay * sr;
        DD[j] = t - q;
      }
   }
  if ( !ec ) DD[N] += q;
  t = -tb * sc;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = (t-r) * sc;
      e = t-q;
      g = t - pay;
      va = v[A[i]];
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
	  if ( j == N && !ec )
            { if ((c = CC[j] ) > (d = DD[j] )) d = c;}
	  else
            if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  c = s+va[B[j]];
          if (c < d) c = d;
          if (c < e) c = e;
	  if ( j - gaplen > 0 )
	    { if ( g < ( temp = CC[j-gaplen-1] - pay ) )
		g = temp;
	      if ( c < g ) c = g;
	    }
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  t = -q * er;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  if ( N <= gaplen )
    for (j = N-1; j >= 0; j--)
      { RR[j] = t = (t-r) * er;
        SS[j] = t-q;
      }
  else
   { temp = N - gaplen;
     for (j = N-1; j >= temp; j--)
      { RR[j] = t = (t-r) * er;
        SS[j] = t-q;
      }
     for (j = temp-1; j >= 0; j--)
      { RR[j] = t = -pay * er;
        SS[j] = t - q;
      }
   }
  if ( !sc ) SS[0] += q;
  t = -te * ec;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = (t-r) * ec;
      g = t - pay;
      e = t-q;
      va = v[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
	  if ( !j && !sc )
            { if ((c = RR[j] ) > (d = SS[j] )) d = c;}
	  else
            if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  c =  s+va[B[j+1]];
          if (c < d) c = d;
          if (c < e) c = e;
	  if ( j + gaplen < N )
	    { if ( g < ( temp = RR[j+gaplen+1] - pay ) )
		g = temp;
	      if ( c < g ) c = g;
	    }
          s = RR[j];
          RR[j] = c;
          SS[j] = d;
        }
    }
  SS[N] = RR[N];

  midc = CC[0]+RR[0];		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = CC[j] + RR[j]) >= midc)
      if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
        { midc = c;
          midj = j;
        }
  for (j = N; j >= 0; j--)
   { if ( j == N )
       d = q * ec;
     else
       if ( j == 0 )
         d = q * sc;
       else
	 d = q;
     if ((c = DD[j] + SS[j] + d) > midc)
       { midc = c;
         midj = j;
         type = 2;
       }
   }
}

/* Conquer: recursively around midpoint */

  cc = midj == N ? ec : 1;
  ss = midj == 0 ? sc : 1;
  if (type == 1)
    { (void) align(A,B,midi,midj,tb,q,sc,sr,cc,1);
      (void) align(A+midi,B+midj,M-midi,N-midj,q,te,ss,1,ec,er);
    }
  else
    { (void) align(A,B,midi-1,midj,tb,zero,sc,sr,cc,1);
      DEL(2);
      (void) align(A+midi+1,B+midj,M-midi-1,N-midj,zero,te,ss,1,ec,er);
    }
  return midc;
}

/* Interface and top level of comparator */
FastInt align_ss2(A,B,M,N,low,up,W,G,H,S,s1,s2,e1,e2)
char A[],B[]; FastInt M,N; FastInt W[][128],G,H; FastInt S[];
FastInt low,up;
FastInt s1,s2,e1,e2;
{ 
  FastInt c;
  int j;

  A--;B--;

  v = W;			/* Setup global parameters */
  q = G;
  r = H;
  qr = G+H;
  gaplen = 200; /* minimum length for constant cost insertion */
  pay = q + r * gaplen;
  sapp = S;
  last = 0;
  al_len = 0;
  no_mat = 0;
  no_mis = 0;
  j = sizeof(FastInt) * (MAX(N,M)+1);
  CC = ( FastInt * ) xmalloc(j);
  DD = ( FastInt * ) xmalloc(j);
  RR = ( FastInt * ) xmalloc(j);
  SS = ( FastInt * ) xmalloc(j);

  c = align(A,B,M,N,q,q,s1,s2,e1,e2);   /* OK, do it */
  
  free(CC);
  free(DD);
  free(RR);
  free(SS);

  return c;
}

/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];

void display_ss2(A,B,M,N,S,AP,BP)
    char A[], B[]; FastInt M, N; FastInt S[], AP, BP;
{ register char *a, *b, *c;
  register FastInt   i,  j, op;
           FastInt   lines, ap, bp;

  i = j = op = lines = 0;
  A--;
  B--;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  while (i < M || j < N)
    { if (op == 0 && *S == 0)
        { op = *S++;
          *a = A[++i];
          *b = B[++j];
          *c++ = (*a++ == *b++) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = B[++j];
              op--;
            }
          else
            { *a++ = A[++i];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N)
        { *a = *b = *c = '\0';
          vmessage("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
	      vmessage("    .    :");
          if (b <= a+5)
	      vmessage("    .");
          vmessage("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}

#ifdef never_used
/* CHECK_SCORE - return the score of the alignment stored in S */

static FastInt CHECK_SCORE(A,B,M,N,S,EG)
char A[], B[]; FastInt M, N; FastInt S[]; char EG;
{ 
  register FastInt   i,  j, op;
  FastInt score;

  score = i = j = op = 0;
  while (i < M || j < N) {
	op = *S++;
	if (EG == 1 && i == 0 && j == 0 && op != 0) {
		if (op > 0) j = j+op;
		else i = i-op;
	} else if (EG == 1 && (i == M || j == N)) {
		i = M;
		j = N;
	} else if (op == 0) 
		score = w[A[++i]][B[++j]] + score;
	else if (op > 0) {
		score = score - (g+op*h);
		j = j+op;
	} else {
		score = score - (g-op*h);
		i = i-op;
	}
  }
  return(score);
}
#endif
