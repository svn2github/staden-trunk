#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO.h"
#include "template.h"
#include "xalloc.h"

/*
 * Stub functions to avoid link errors from IO.c.
 * This is the 'note execution' code, which we wish to avoid anyway.
 */
void execute_database_notes(GapIO *io, char *type) {}
void process_rawdata_note(GapIO *io) {}
void fix_notes(GapIO *io) {}

/* Maximum length of a legitimate sequence */
#define MAXREADLEN 2000

/*
 * Produce a quality value profile based on all sequences of a specific
 * chemistry. Chemistry==0 implies all sequences regardless of chemistry.
 * 
 * Returns: a count of how many samples have been added into total_conf;
 */
int read_conf_dist(GapIO *io, int chemistry, int *total_conf) {
    int i, j;
    GReadings r;
    int1 conf[MAXREADLEN];
    int count = 0;


    for (i = 1; i <= NumReadings(io); i++) {
	gel_read(io, i, r);
	if (!r.confidence)
	    continue;

	if (r.length >= MAXREADLEN)
	    continue;

	if (chemistry && r.chemistry != chemistry)
	    continue;

	if (DataRead(io, r.confidence, conf, r.length * sizeof(int1),
		     sizeof(int1)))
	    continue;

	/* Always work on top strand data */
	if (r.sense == 1) {
	    int left, right;
	    for (left = 0, right = r.length-1; left < right; left++, right--) {
		int tmp = conf[left];
		conf[left] = conf[right];
		conf[right] = tmp;
	    }
	}


	/* All < conf 10 - probably not a real */
	for (j = 0; j < r.length; j++) {
	    if (conf[j] >= 10)
		break;
	}
	if (j == r.length)
	    continue;

	/* All > conf 40 - probably not a real sequence */
	for (j = 0; j < r.length; j++) {
	    if (conf[j] < 40)
		break;
	}
	if (j == r.length)
	    continue;

	/* All same value - probably not a real sequence */
	for (j = 1; j < r.length; j++) {
	    if (conf[j] != conf[0])
		break;
	}
	if (j == r.length)
	    continue;

	/* Accumulate totals and counts */
	for (j = 0; j < r.length; j++) {
	    total_conf[j] += conf[j];
	}

	count++;
    }

    return count;
}

static void usage(void) {
    fprintf(stderr, "read_conf_dist chemistry GAP4_D.B ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    GapIO *io;
    char *project, *version;
    int status;
    int argn;
    int chemistry;
    int total_conf[MAXREADLEN];
    double total_dconf[MAXREADLEN];
    int count = 0;
    int i;

    if (argc < 3) {
	usage();
    }

    chemistry = atoi(argv[1]);

    memset(total_conf, 0, MAXREADLEN * sizeof(int));

    for (argn = 2; argn < argc; argn++) {
	project = argv[argn];
	version = strchr(project, '.');
	if (!version) {
	    fprintf(stderr, "Malformed project name\n");
	    return 1;
	}
	*version++ = 0;

	fprintf(stderr, "Opening project '%s' version '%s'\n",
		project, version);
	io = open_db(project, version, &status, 0, 1);
	if (!io) {
	    fprintf(stderr, "Couldn't open database\n");
	    return 1;
	}

	count += read_conf_dist(io, chemistry, total_conf);
	
	close_db(io);
    }

    fprintf(stderr, "Count = %d\n", count);
    for (i = 0; i < MAXREADLEN; i++) {
	total_dconf[i] = total_conf[i] / (double)count;
    }

    for (i = 0; i < MAXREADLEN; i++) {
	printf("%4d %f\n", i, total_dconf[i]);
    }

    return 0;
}
