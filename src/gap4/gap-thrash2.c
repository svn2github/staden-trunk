#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "IO.h"

#define T_MAX_READ_LEN 4096

/*
 * Create a new reading of random length with a start of used data somewhere
 * within this reading, and an end somewhere between the start and the length.
 */
int t_add_reading(GapIO *io, int num) {
    char seq[T_MAX_READ_LEN+1];
    int1 conf[T_MAX_READ_LEN+1];
    int2 opos[T_MAX_READ_LEN+1];
    int i, len;
    int2 length, start, end;
    
    len = random() % (T_MAX_READ_LEN-1)+1;
    for (i=0; i<len; i++) {
	seq[i] = "ACGT"[random() % 4];
	conf[i] = random()%100;
	opos[i] = i;
    }
    length = len;
    start = random() % len;
    end = (random() % (length - start)) + start;

    return io_write_seq(io, num, &length, &start, &end, seq, conf, opos);
}

t_del_reading(GapIO *io, int num) {
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

#define NR 1000

void thrash(GapIO *io) {
    int j;
    int changed = 0;
    int deleted = 0;
    int flushed = 0;
    
    for (;;) {
	int x;

	if ((x = NumReadings(io)+50) > NR)
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
	if (random()%100 == 0) {
	    printf(">>> flush\n");
	    flush2t(io);
	    flushed++;

	    printf("changed=%d, deleted=%d, flushed=%d\n",
		   changed, deleted, flushed);
	}
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

    system("/bin/rm thrash2.0*");
    io = open_db("thrash2", "0", &status, 1, 0);
    thrash(io);

    close_db(io);
}
