#include "IO.h"
#include "misc.h"
#include "io-reg.h"
#include "fort.h"
#include "dis_readings.h"

/*
 * This is now just an interface to disassemble readings to move all readings
 * from a specified reading number and onwards to a new contig.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int break_contig(GapIO *io, int rnum) {
    int *reads = NULL;
    int nreads, ret;

    /* Worst case memory usage. Temporary so we don't care greatly */
    if (!(reads = (int *)xmalloc((NumReadings(io)+1) * sizeof(int))))
	return -1;

    /* Produce reading list */
    for (nreads = 0; rnum; rnum = io_rnbr(io, rnum)) {
	reads[nreads++] = rnum;
    }

    /* Disassemble */
    ret = disassemble_readings(io, reads, nreads, 2 /* move */, 1 /* dup */);

    xfree(reads);
    return ret;
}
