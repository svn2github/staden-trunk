#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "copy_db.h"
#include "misc.h"
#include "licence.h"

static void usage(void) {
    fprintf(stderr, "copy_db [-v] [-f] source ... destination\n");
    exit(1);
}

int main(int argc, char **argv) {
    char *pf, *pt, *to;
    GapIO *iof, *iot;
    int status, verbose = 0;

    extern int gap_fatal_errors;
    extern int maxdb;

    check_licence();

    maxdb = 20000;

    for (;argc > 1 && *argv[1] == '-'; argv++, argc--) {
	if (strcmp(argv[1], "-v") == 0)
	    verbose = 1;
	else if (strcmp(argv[1], "-f") == 0)
	    gap_fatal_errors = 0;
	else
	    usage();
    }

    if (argc < 3) {
	usage();
    }

    to = argv[argc-1];
    pt = strrchr(to, '.'); pt[2] = 0;

    if (file_exists(to)) {
	char ans, buf[100];

	printf("Database %s already exists.\n"
	       "Do you wish to overwrite it? [n] ", to);
	scanf("%c", &ans);
	if (ans != 'y' && ans != 'Y')
	    return 0;
	strcpy(buf, to);
	remove(buf);
	strcpy(buf+strlen(to), ".aux");
	remove(buf);
	strcpy(buf+strlen(to), ".BUSY");
	remove(buf);
    }

    if (NULL == pt) {
	fprintf(stderr, "Malformed database name. Should be \"PROJ.V\"\n");
	return 2;
    }
    *pt = 0;
    if (NULL == (iot = open_db(to, pt+1, &status, 1, 0))) {
	fprintf(stderr, "Couldn't create database\n");
	return 5;
    }

    for (; argc >= 3; argc--, argv++) {
	pf = strrchr(argv[1], '.'); pf[2] = 0;
	if (strcmp(argv[1], to) == 0) {
	    fprintf(stderr, "copy_db: %s and %s are identical (not copied).\n",
		    to, to);
	    return 7;
	}

	if (NULL == pf) {
	    fprintf(stderr, "Malformed database name. Should be \"PROJ.V\"\n");
	    return 2;
	}

	*pf = 0;
	
	if (verbose)
	    printf("Copying database %s.%s to database %s.%s\n",
		   argv[1], pf+1, to, pt+1);

	if (NULL == (iof = open_db(argv[1], pf+1, &status, 0, 1))) {
	    fprintf(stderr, "Couldn't read database\n");
	    return 4;
	}

	if (-1 == copy_database(iof, iot, verbose, gap_fatal_errors)) {
	    fprintf(stderr, "Couldn't copy database\n");
	    return 6;
	}
	
	close_db(iof);
    }

    close_db(iot);

    return 0;
}

/*
 * Stub functions to avoid link errors from IO.c.
 * This is the 'note execution' code, which we wish to avoid anyway.
 */
void execute_database_notes(GapIO *io, char *type) {}
void process_rawdata_note(GapIO *io) {}
void fix_notes(GapIO *io) {}
