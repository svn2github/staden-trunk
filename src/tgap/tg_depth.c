/*
 * Author:         James Bonfield, Feb 2007
 *                 Wellcome Trust Sanger Institute
 *
 * g_view:
 * Loads and displays the contents of a g-library format database as created
 * by the g_index application.
 */

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <curses.h>
#include <signal.h>
#include <ctype.h>
#include <string.h>

#include "array.h"
#include "misc.h"
#include "tg_gio.h"

static void test_mode2(GapIO *io, contig_t **c, int x1, int x2, int bpv) {
    int i;
    track_t *t;

    t = contig_get_track(io, c, x1, x2, TRACK_READ_DEPTH, bpv);
    for (i = 0; i < (x2-x1+1)/bpv; i++) {
	printf("%d\t%d\n", x1+i*bpv, arr(int, t->data, i));
    }

    cache_flush(io);
    gio_close(io);
    exit(0);
}

/* ------------------------------------------------------------------------ */
int main(int argc, char **argv) {
    GapIO *io;
    contig_t *c;
    int cnum = 0;
    int read_only = 0;
    int x1 = 0, x2 = 1024, bpv = 1;

    if (argc > 2) x1  = atoi(argv[2]);
    if (argc > 3) x2  = atoi(argv[3]);
    if (argc > 4) bpv = atoi(argv[4]);

    if (NULL == (io = gio_open(argv[1], read_only, 0))) {
	fprintf(stderr, "Unable to open db: %s\n", argv[1]);
	return 1;
    }

    io->contig_num = cnum;
    gio_read_contig(io, cnum, &c); 
    cache_incr(io, c);

    test_mode2(io, &c, x1, x2, bpv);


    return 0;
}
