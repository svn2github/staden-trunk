#include <staden_config.h>

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "renz_utils.h"
#include "sequence_formats.h"
#include "misc.h"
#include "getfile.h"
#include "dna_utils.h"
#include "genetic_code.h"

extern char *optarg;
extern int optind;

/*
 * Converts a restriction enzyme name into an index into the enzymefile.
 *
 * This is needed because the enzyme code in seq_utils doesn't work on enzyme
 * names, but on indexes to the enzyme files.
 *
 * Returns line number for success (starting from zero)
 *         -1 for error.
 */
#define MAXLINE 1024
static int renz2index(char *name) {
    char path[FILENAME_MAX];
    char line[MAXLINE+1];
    char *cp;
    int lineno;
    FILE *fp;
    
    expandpath("$STADENROOT/tables/RENZYM.ALL", path);
    if (NULL == (fp = fopen(path, "r"))){
	return -1;
    }

    lineno = 0;
    while (fgets(line, MAXLINE, fp) != NULL) {
	if (cp = strchr(line, '/'))
	    *cp = 0;
	if (strcasecmp(line, name) == 0) {
	    fclose(fp);
	    return lineno;
	}
	lineno++;
    }

    fclose(fp);

    return -1;
}

/*
 * Prints out the position of an enzyme cut site.
 * If vp is set to 1, we also print up the vector-primer file format
 * information.
 *
 * Returns 0 for success
 *        -1 for failure
 */
static int find_renz(char *in_file, char *enz,
		     int vp, int size_left, int size_right) {
    R_Enz *renzyme = NULL;
    int nrenzyme = 0;
    R_Match *rmatch;
    int nrmatch = 0;
    char strline[1024];
    int i, j;
    char *vector = NULL;
    int vector_len;

    if (0 != get_seq(&vector, 100000, &vector_len, in_file, NULL)) {
	fprintf(stderr, "Failed to read file '%s'\n", in_file);
	return -1;
    }

    if (NULL == (rmatch = (R_Match*)xcalloc(MAXMATCHES, sizeof(R_Match))))
	return -1;

    sprintf(strline, "%d", i = renz2index(enz));
    if (-1 == i) {
	goto tidyup;
    }

    open_renz_file("$STADTABL/RENZYM.ALL", strline, 1, &renzyme, &nrenzyme);

    if (1 != FindMatches(renzyme, nrenzyme, vector, vector_len,
			 1 /* circular */, &rmatch, &nrmatch)) {
	fprintf(stderr, "Failed in FindMatches()\n");
	nrmatch = 0;
	goto tidyup;
    }

    if (nrmatch == 1) {
	if (vp) {
	    char name[1024], *cp;
	    char *seq, *seqp;

	    /* Needs to be circular, so we tripple the sequence and
	     * set a pointer to the middle copy.
	     */
	    seq = (char *)xmalloc(3*vector_len+1);
	    sprintf(seq, "%.*s%.*s%.*s",
		    vector_len, vector,
		    vector_len, vector,
		    vector_len, vector);
	    seqp = seq + vector_len;

	    /* Extract a sequence name from the filename */
	    if (cp = strrchr(in_file, '/'))
		cp++;
	    else
		cp = in_file;
	    strcpy(name, cp);

	    if (cp = strchr(name, '.'))
		*cp = 0;
	    
	    strcat(name, "/");
	    strcat(name, enz);
	    
	    if (cp = strrchr(in_file, '/'))
		cp++;
	    else
		cp = in_file;

	    printf("%s\t\t%.*s\t%.*s\t%s\n", name,
		   size_left, seqp + rmatch[0].cut_pos - size_left - 1,
		   size_right, seqp + rmatch[0].cut_pos - 1, cp);

	    xfree(seq);
	    
	} else {
	    printf("%d\n", rmatch[0].cut_pos);
	}
    } else {
	if (nrmatch > 1)
	    fprintf(stderr, "Found more than one match\n");
	else
	    fprintf(stderr, "Enzyme not found in sequence\n");
    }

 tidyup:
    if (vector)
	xfree(vector);
    xfree(rmatch);
    if (renzyme) {
	for (i = 0; i < nrenzyme; i++) {
	    xfree(renzyme[i].name);
	    for (j = 0; j < renzyme[i].num_seq; j++) {
		xfree(renzyme[i].seq[j]);
	    }
	    xfree(renzyme[i].seq);
	    xfree(renzyme[i].cut_site);
	}
	xfree(renzyme);
    }

    return nrmatch == 1 ? 0 : -1;
}

static void usage(void) {
    fprintf(stderr, "Usage: find_renz [-vp] enzyme filename ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    int ret = 0;
    int vp = 0;
    char *enz;
    int size_left = 50;
    int size_right = 50;

    set_char_set(1);    /* 1 == DNA */
    set_dna_lookup();   /* general lookup and complementing */
    set_iubc_lookup();  /* iubc codes for restriction enzymes */
    init_genetic_code();

    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (strcasecmp(*argv, "-vp") == 0) {
            vp = 1;
	} else if (**argv != '-') {
	    break;
        } else {
            usage();
        }
    }

    if (argc == 0)
	usage();

    enz = *argv++;
    argc--;

    for (;argc > 0; argc--, argv++)
	ret |= find_renz(*argv, enz, vp, size_left, size_right);
    
    return ret ? 1 : 0;
}

