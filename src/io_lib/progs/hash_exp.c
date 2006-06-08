#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "hash_table.h"
#include "os.h"
#include "mFILE.h"

HashFile *build_index(FILE *fp) {
    char line[8192];
    char rname[8192];
    size_t pos = 0, last = 0;
    char *c = "magic";
    HashFile *hf;

    /* Create the hash table */
    hf = HashFileCreate(0, HASH_DYNAMIC_SIZE);
    hf->headers = (HashFileSection *)malloc(sizeof(*hf->headers));
    hf->nheaders = 0;

    for (*rname = 0; c;) {
	c = fgets(line, 8192, fp);

	if (c == NULL || strncmp("ID   ", line, 5) == 0) {
	    /* Add this entry; it extends from 'last' to 'pos' */
	    if (*rname) {
		HashData hd;
		HashFileItem *hfi = (HashFileItem *)calloc(1, sizeof(*hfi));
		
		hfi->header = 0;
		hfi->footer = 0;
		hfi->pos = last;
		hfi->size = pos - last;
		hd.p = hfi;
		
		HashTableAdd(hf->h, rname, strlen(rname), hd, NULL);
	    }

	    /* Remember this ID line for when we meet the next */
	    if (c) {
		char *nl = strchr(c, '\n');
		if (nl)
		    *nl = 0;

		strcpy(rname, c+5);
	    }

	    last = pos;
	}
	pos = ftell(fp);
    }

    HashTableStats(hf->h, stderr);

    return hf;
}

int main(int argc, char **argv) {
    HashFile *hf;
    FILE *fp;

    if (argc != 2) {
	fprintf(stderr, "Usage: hash_exp exp_file_ball > exp.hash\n");
	return 1;
    }
    if (NULL == (fp = fopen(argv[1], "r"))) {
	perror(argv[1]);
	return 1;
    }

    hf = build_index(fp);
    hf->archive = strdup(argv[1]);
    HashFileSave(hf, stdout, 0);
    return 0;
}
