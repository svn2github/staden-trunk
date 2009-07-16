#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "copy_db.h"
#include "misc.h"
#include "gap-create.h"

#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif
 
extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;

static void usage(void) {
    fprintf(stderr, "copy_db [-v] [-f] [-b 32/64] [-T] source ... destination\n");
    fprintf(stderr, "    -v     verbose\n");
    fprintf(stderr, "    -f     fix problems\n");
    fprintf(stderr, "    -b 32  32-bit format\n");
    fprintf(stderr, "    -b 64  64-bit format\n");
    fprintf(stderr, "    -T     remove tags (annotations)\n");
    exit(1);
}

int main(int argc, char **argv) {
    char *pf, *pt, *to;
    GapIO *iof, *iot;
    int status, verbose = 0;
    int bitsize = G_32BIT;
    int c;
    int notags = 0;

    extern int gap_fatal_errors;
    extern int maxdb;

    maxdb = 20000;

    while ((c = getopt(argc, argv, "vfb:T")) != -1) {
	switch (c) {
	case 'v':
	    verbose = 1;
	    break;

	case 'f':
	    gap_fatal_errors = 0;
	    break;

	case 'b':
	    bitsize = atoi(optarg);
	    if (bitsize != 32 && bitsize != 64)
		usage();
	    else
		bitsize = (bitsize == 32) ? G_32BIT : G_64BIT;
	    break;

	case 'T':
	    notags=1;
	    break;

	default:
	    usage();
	}
    }

    if (argc - optind < 2) {
	usage();
    }

    to = argv[argc-1];
    pt = strrchr(to, '.');
    if (NULL == pt) {
	fprintf(stderr, "Malformed database name. Should be \"PROJ.V\"\n");
	return 2;
    }
    pt[2] = 0;

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

    *pt = 0;
    set_db_bitsize(bitsize);
    if (NULL == (iot = open_db(to, pt+1, &status, 1, 0))) {
	fprintf(stderr, "Couldn't create database\n");
	return 5;
    }

    for (; optind < argc-1; optind++) {
	char *fn = argv[optind];
	pf = strrchr(fn, '.'); pf[2] = 0;
	if (strcmp(fn, to) == 0) {
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
		   fn, pf+1, to, pt+1);

	if (NULL == (iof = open_db(fn, pf+1, &status, 0, 1))) {
	    fprintf(stderr, "Couldn't read database\n");
	    return 4;
	}

	if (-1 == copy_database(iof, iot, verbose, gap_fatal_errors, notags)) {
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
