#include <stdio.h>
#include <stdlib.h>

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

void read_size_dist(GapIO *io) {
    int i;
    GReadings r;

    for (i = 1; i <= NumReadings(io); i++) {
	gel_read(io, i, r);
	printf("%d %d\n", r.sequence_length, r.chemistry);
    }
}

static void usage(void) {
    fprintf(stderr, "insert_size_dist GAP4_D.B ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    GapIO *io;
    char *project, *version;
    int status;
    int argn;

    if (argc < 2) {
	usage();
    }

    for (argn = 1; argn < argc; argn++) {
	project = argv[argn];
	version = strchr(project, '.');
	if (!version) {
	    fprintf(stderr, "Malformed project name\n");
	    return 1;
	}
	*version++ = 0;

	printf("Opening project '%s' version '%s'\n", project, version);
	io = open_db(project, version, &status, 0, 1);
	if (!io) {
	    fprintf(stderr, "Couldn't open database\n");
	    return 1;
	}

	read_size_dist(io);
	
	close_db(io);
    }

    return 0;
}
