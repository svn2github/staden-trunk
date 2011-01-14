/*
 * Fetches all or several specific comment(s) from the SCF comment section.
 * Comments are referenced by a Field-ID.
 *
 * Options:
 *        -s 	Silent - no error messages produced (just error codes).
 *        -q	Query - returns 0 or 3 depending on whether the Field-ID
 *              exists. No messages sent to stdout.
 *        -c	Supresses display of Field-ID.
 * 
 * Returns 0 for success.
 *         1 for usage failure
 *         2 for file IO failure
 *         3 for unknown Field-ID
 */

#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#include <os.h>
#include <misc.h>
#include <io_lib/mFILE.h>
#include <io_lib/scf.h>
#include <io_lib/xalloc.h>
#include <io_lib/compress.h>

/* 6/1/99 johnt - need to implicitly import globals from DLLs with Visual C++ */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif


extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;

#define MAX_COMMENT 1024

static int silent = 0, query = 0, suppress = 0;
static int do_it(char *file, int argc, char **argv);

__NORETURN__ static void usage(void) {
    fprintf(stderr, "Usage: get_scf_field [-sqc] scf_file\n");
    fprintf(stderr, "or     get_scf_field [-sqc] scf_file Field-ID ...\n");

    exit(1);
}


int main(int argc, char **argv) {
    char *file;
    int c;

    optarg = NULL;
    while ((c = getopt(argc, argv, "sqc")) != -1) {
	switch (c) {
	case 's':
	    silent = 1;
	    break;
	case 'q':
	    query = 1;
	    break;
	case 'c':
	    suppress = 1;
	    break;
	case '?':
	    usage();
	}
    }

    if (optind == argc)
	usage();
    else
	file = argv[optind++];

    return do_it(file, argc-optind, &argv[optind]);
}


static int do_it(char *file, int argc, char **argv) {
    mFILE *fp;
    Header h;
    char *c, *p, *p2;
    int all = 0;
    int found = 0;

    if (NULL == (fp = fopen_compressed(file, NULL))) {
	if (!silent)
	    perror(file);
	return 2;
    }

    /* Read the SCF header */

    if (-1 == read_scf_header(fp, &h)) {
	if (!silent)
	    fprintf(stderr, "Failed to read header\n");
	return 2;
    }

    /* Read our comments */
    if (NULL == (c = (char *)xmalloc(h.comments_size + 1)))
	return 2;
    
    if (-1 == mfseek(fp, h.comments_offset, SEEK_SET)) {
	if (!silent)
	    perror("seek()");
	return 2;
    }

    if (-1 == read_scf_comment(fp, c, h.comments_size)) {
	if (!silent)
	    perror("read() of comments");
	return 2;
    }

    if (argc == 0)
	all = 1;

    p = c;

    for (; p2 = strstr(p, "\n"); p = p2+1) {
	char *q;

	*p2 = '\0';
	if (NULL == (q = strchr(p, '=')))
	    continue;

	*q = '\0';
	
	if (!all) {
	    int i;

	    for (i=0; i<argc; i++)
		if (strcmp(argv[i], p) == 0)
		    break;

	    if (i == argc)
		continue;
	}

	if (suppress)
	    printf("%s\n", q+1);
	else
	    printf("%s=%s\n", p, q+1);

	found = 1;
    }

    return found ? 0 : 3;
}
