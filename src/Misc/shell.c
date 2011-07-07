#include "misc.h"
#include <stdio.h>

void shell_call(char *command, char *output, int len)
{
    FILE *pipe;
    char *a;

    output[0]='\0';
    /* some systems don't prototype this */
    pipe = (FILE *)popen(command,"r");
    if (!fgets(output,len,pipe))
	*output = 0;
    pclose(pipe);

    /* clobber last new line */
    for (a=output;*a && *a != '\n'; a++);
    *a = '\0';
}

/* Getopt for windows */
#ifdef _MSC_VER

/*
 * An implementation of the getopt() function to be used on platforms which
 * do not have a system version (eg. MS Windows)
 */
#include <stdio.h>
#include <string.h>

char *optarg = NULL;
int optind=1;
int opterr=1;

int getopt(int argc, char * const argv[], const char *optstring) {
    static char *arg = NULL;
    char *optcode;
    int opt;

    /* Initialisation */
    if (!arg) {
	arg = "";
	optind = 1;
    }

    /* End of this argument (or starting up), move on to next. */
    if (!*arg) {
	if (optind >= argc) {
	    arg = NULL;
	    return -1;
	}

	arg = argv[optind++];

	/* Arguments end when we see a non '-' sign */
	if (*arg != '-') {
	    optind--;
	    arg = NULL;
	    return -1;
	}
	arg++;
    }

    /* Get the option. '-' by itself as an option also ends the parsing */
    if (!(opt = *arg++)) {
	optind--;
	arg = NULL;
	return -1;
    }

    /* Find the index in the optstring */
    if (NULL == (optcode = strchr(optstring, opt))) {
	if (opterr)
	    fprintf(stderr, "%s: illegal option -- %c\n", argv[0], opt);
	return '?';
    }

    /* If next char is ':', set optarg to point to the argument */
    if (*++optcode == ':') {
	if (*arg) {
	    optarg = arg;
	    arg = "";
	} else {
	    if (optind < argc) {
		optarg = argv[optind++];
	    } else {
		if (opterr)
		    fprintf(stderr, "%s: option requires an argument -- %c\n",
			    argv[0], opt);
		optarg = "?";
	    }
	}
    }

    return opt;
}

#endif

