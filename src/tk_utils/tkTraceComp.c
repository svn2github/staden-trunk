/*
 * Title:  opp.c
 *
 * File:   opp.c
 * Purpose: code for complementing sequences
 *
 * Last update: Mon Oct 31, 1994
 *
 * 15.01.90 SD  Taken from seqIOEdit.c
 * 31.10.94 JKB Convert to use Read structure
*/


#include "Read.h"
#include "tkTrace.h"
#include "tkTraceIO.h"

static char opp[256]; /* complement of any given base */

static void oppInitialize(void)
{
    int i;

    for (i = 0; i<256; i++) opp[i]='-';

    /* RMD 31/12/90 'N' -> '-' above.
     * removed 'N' and 'n' entries below and added reciprocal
     * 'K' and 'N' entries as for full Staden table
     */

    opp['A'] = 'T';
    opp['G'] = 'C';
    opp['T'] = 'A';
    opp['C'] = 'G';
    opp['a'] = 't';
    opp['g'] = 'c';
    opp['t'] = 'a';
    opp['c'] = 'g';
    opp['D'] = 'H';
    opp['H'] = 'D';
    opp['V'] = 'B';
    opp['B'] = 'V';
    opp['K'] = 'N';
    opp['N'] = 'K';
    opp['L'] = 'M';
    opp['M'] = 'L';
    opp['5'] = '6';
    opp['6'] = '5';
    opp['R'] = 'Y';
    opp['Y'] = 'R';
    opp['7'] = '7';
    opp['8'] = '8';
}


/* ---- Exports ---- */

/*
 * Complement and reverse bases and traces
 */
void complement_read(Read *read, int len)
{
    static int initialised = 0;
    int_2 temp_int2;
    TRACE *temp_TRACEptr;
    char temp_char;
    int i;

    if (!initialised) {
	oppInitialize();
	initialised = 1;
    }

    /* swap */
#define swap(A,B,I) ( (I)=(A), (A)=(B), (B)=(I) )

    /* complement and reverse traces */
    /* swap traces A<->T and C<->G */
    swap(read->traceA,read->traceT,temp_TRACEptr);
    swap(read->traceC,read->traceG,temp_TRACEptr);

    /* reverse points in traces */
    if (read->traceA) {
	for (i=0;i<read->NPoints/2;i++) {
	    swap(read->traceA[i],read->traceA[read->NPoints-i-1],temp_int2);
	    swap(read->traceC[i],read->traceC[read->NPoints-i-1],temp_int2);
	    swap(read->traceG[i],read->traceG[read->NPoints-i-1],temp_int2);
	    swap(read->traceT[i],read->traceT[read->NPoints-i-1],temp_int2);
	}
    }

    /* complement the sequence */
    for (i=0;i<read->NBases;i++) {
	signed int tpos;

	read->base[i] = opp[(unsigned)read->base[i]];
	tpos = read->NPoints - read->basePos[i] - 1;
	read->basePos[i] = tpos > 0 ? tpos : 0;

	swap(read->prob_A[i], read->prob_T[i], temp_char);
	swap(read->prob_C[i], read->prob_G[i], temp_char);
    }

    /* reverse sequence */
    for (i=0;i<read->NBases/2;i++) {
	swap(read->base[i],read->base[read->NBases-i-1],temp_char);
	swap(read->basePos[i],read->basePos[read->NBases-i-1],temp_int2);
	swap(read->prob_A[i], read->prob_A[read->NBases-i-1], temp_char);
	swap(read->prob_C[i], read->prob_C[read->NBases-i-1], temp_char);
	swap(read->prob_G[i], read->prob_G[read->NBases-i-1], temp_char);
	swap(read->prob_T[i], read->prob_T[read->NBases-i-1], temp_char);
    }

    /* swap cutoffs */
    i = read->leftCutoff;
    if (read->rightCutoff)
        read->leftCutoff = len - read->rightCutoff + 1;
    else
        read->leftCutoff = 0;

    if (i)
        read->rightCutoff = len - i + 1;
    else
        read->rightCutoff = 0;
}

/*
 * Complement and reverse a ted trace.
 */
void complement_trace(DNATrace *t) {
    int i;
    char temp_char;
    int_2 temp_int2;

    if (!t->read)
	return;

    complement_read(t->read, t->Ned);

    i = t->leftVector;
    if (t->rightVector != -1)
        t->leftVector = t->Ned - t->rightVector + 1;
    else
        t->leftVector = -1;

    if (i != -1)
        t->rightVector = t->Ned - i + 1;
    else
        t->rightVector = -1;

    /* complement the edited sequence */
    for (i=0;i<t->Ned;i++) {
	t->edBases[i] = opp[(unsigned)t->edBases[i]];
    }

    /* reverse sequence */
    for (i=0;i<t->Ned/2;i++) {
	swap(t->edBases[i],t->edBases[t->Ned-i-1],temp_char);
	swap(t->edPos[i],t->edPos[t->Ned-i-1],temp_int2);
	swap(t->edConf[i],t->edConf[t->Ned-i-1],temp_char);
    }

    /* screen position */
    t->disp_offset = t->read->NPoints - t->disp_offset - t->disp_width;

    t->comp ^= 1;

    /* Calculate trace positions */
    trace_init_pos(t);
}
