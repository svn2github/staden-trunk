/*
 * As gap-thrash2, but we ensure that we re-read records periodically.
 * Note that a gel_read won't work as the data is cached.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "IO.h"

/*#define T_MAX_READ_LEN 4096 */
#define T_MAX_READ_LEN 66
#define NR 77

/*
 * Create a new reading of random length with a start of used data somewhere
 * within this reading, and an end somewhere between the start and the length.
 */
int t_add_reading(GapIO *io, int num) {
    char seq[T_MAX_READ_LEN+1];
    int1 conf[T_MAX_READ_LEN+1];
    int2 opos[T_MAX_READ_LEN+1];
    int i, len, cval;
    int2 length, start, end;
    
    len = random() % (T_MAX_READ_LEN-1)+1;
    cval = random()%100;
    for (i=0; i<len; i++) {
	seq[i] = "ACGT"[random() % 4];
	conf[i] = cval;
	opos[i] = i;
    }
    length = len;
    start = random() % len;
    end = (random() % (length - start)) + start;

    return io_write_seq(io, num, &length, &start, &end, seq, conf, opos);
}

/*
 * Delete a reading
 */
void t_del_reading(GapIO *io, int num) {
    char gel[4096+1];
    GReadings r;

    io_deallocate_reading(io, num);
    gel_read(io, NumReadings(io), r);
    gel_write(io, num, r);
    NumReadings(io)--;
    Nreadings(io)--; /* GBUG: shouldn't be needed */
    deallocate(io, arr(GCardinal, io->readings, Nreadings(io)));
    /* YUK! should write database record back */
    GT_Write(io,GR_Database,&io->db,sizeof(io->db),GT_Database);
    /* write array */
    ArrayWrite(io, io->db.readings, io->db.Nreadings, io->readings);
}

/*
 * Read a reading. We read everything that we've written, including
 * sequence, conf, orig pos, and GReadings structure. We need to make sure
 * that we use GT_Read for the GReadings struct as gel_read simply uses
 * a memory cache.
 */
void t_read_reading(GapIO *io, int num) {
    char seq[T_MAX_READ_LEN+1];
    int1 conf[T_MAX_READ_LEN+1];
    int2 opos[T_MAX_READ_LEN+1];
    int2 length, start, end;
    GReadings r;
    int i, cval;;

    /*
     * Read the reading and check that it's identical to the memory cache
     */
    GT_Read(io, arr(GCardinal, io->readings, num-1),
	    &r, sizeof(r), GT_Readings);
    if (memcmp(&r, arrp(GReadings, io->reading, num-1), sizeof(r)) != 0) {
	printf("Reading %d is not identical to memory cache\n", num);
	abort();
    }
    
    io_read_seq(io, num, &length, &start, &end, seq, conf, opos);

    /*
     * Check obvious size constraints.
     */
    if (length > T_MAX_READ_LEN + 1 || length < 0) {
	printf("Invalid length\n");
	abort();
    }
    if (start < 0 || start > length) {
	printf("Invalid start\n");
	abort();
    }
    if (end < 0 || end > length || end < start) {
	printf("Invalid end\n");
	abort();
    }

    /*
     * The sequence should consist purely of A, C, G and Ts
     */
    for (i=0; i<length; i++) {
	switch(seq[i]) {
	case 'A':
	case 'C':
	case 'G':
	case 'T':
	    break;

	default:
	    printf("Unknown sequence characters\n");
	    abort();
	}
    }

    /*
     * The confidence is a single value throughout
     */
    cval = conf[0];
    for (i=1; i<length; i++) {
	if (cval != conf[0]) {
	    printf("Invalid confidence value\n");
	    abort();
	}
    }

    /*
     * The original positions are in ascending order
     */
    for (i=0; i<length; i++) {
	if (opos[i] != i) {
	    printf("Invalid original position\n");
	    abort();
	}
    }
}

void thrash(GapIO *io) {
    int j;
    int changed = 0;
    int deleted = 0;
    int flushed = 0;
    int reads = 0;
    int z;

    for (z=0;z<1000;z++) {
	int x;

	if ((x = NumReadings(io)+(NR/10)) > NR)
	    x = NR;
	j = random()%x+1;
	if (random()%10 > 3) {
	    printf(">>> mod/add_reading %d\n",j);
	    t_add_reading(io, j);
	    changed++;
	} else {
	    if (j > NumReadings(io))
		continue;

	    printf(">>> deallocate reading %d\n", j);
	    /*io_deallocate_reading(io, j);*/
	    t_del_reading(io, j);
	    deleted++;
	}
	if (NumReadings(io)) {
#if 1
	    j = random()%NumReadings(io) + 1;
	    printf(">>> read reading %d\n", j);
	    t_read_reading(io, j);
	    reads++;
#else
	    for (j=1; j<=NumReadings(io); j++)
		t_read_reading(io,j);
	    reads += NumReadings(io);
#endif
	}
	if (random()%5 == 0) {
	    printf(">>> flush\n");
	    flush2t(io);
	    flushed++;

	    printf("read=%d, changed=%d, deleted=%d, flushed=%d\n",
		   reads, changed, deleted, flushed);
	}
	/*dmalloc_verify(0L);*/
    }
}

int main(int argc, char **argv) {
    GapIO *io;
    int status;
    time_t seed = time(NULL);

    if (argc == 2) {
	seed = atoi(argv[1]);
    }
    printf("Seed=%ld\n", seed);
    srandom(seed);

    system("/bin/rm thrash.0*");
    io = open_db("thrash", "0", &status, 1, 0);
    thrash(io);

    close_db(io);
}
