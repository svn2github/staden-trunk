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
}

/*
 * The following thrash code, when called with a seed of 666, causes the G
 * library to corrupt it's data structures. Somehow the assert of
 * cache->refs being 1 is broken (it's 0). The code below is
 * incorrect 'gap' level code, but ought to be fine at the g level.
 *
 * The cause is due to deallocating or using a view that's already
 * been deallocated.
 */
void thrash(GapIO *io) {
    int i, j;
    int buf[20];
    
    for (i=0; i<20; i++)
	buf[i] = 0;

    for (i=0; i<1000000; i++) {
	j = random()%10+1;
	if (random()%10 > 3) {
	    if (buf[j])
		continue;
	    printf(">>> add_reading %d\n",j);
	    t_add_reading(io, j);
	    buf[j] = 1;
	} else {
	    if (buf[j]) {
		printf(">>> deallocate reading %d\n", j);
		/*io_deallocate_reading(io, j);*/
		t_del_reading(io, j);
		buf[j] = 0;
	    }
	}
	if (random()%5 == 0 || 1) {
	    printf(">>> flush\n");
	    flush2t(io);
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
