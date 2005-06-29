#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"

/* Copies a single named file to stdout */
void extract(HashFile *hf, char *file) {
    size_t len;
    char *data;

    if (data = HashFileExtract(hf, file, &len)) {
	fwrite(data, len, 1, stdout);
	free(data);
    }
    return;
}

int main(int argc, char **argv) {
    char *fofn = NULL;
    char *hash;
    HashFile *hf;

    /* process command line arguments of the form -arg */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-I") == 0) {
	    argv++;
	    fofn = *argv;
	    argc--;
	}
    }

    if (argc < 2) {
	fprintf(stderr, "Usage: hash_extract [-I fofn] hashfile [name ...]\n");
	return 1;
    }
    hash = argv[0];
    argc--;
    argv++;

    if (NULL == (hf = HashFileOpen(hash))) {
	perror(hash);
	return 1;
    }

    if (fofn) {
	FILE *fofnfp;
	char file[256];

	if (strcmp(fofn, "-") == 0) {
	    fofnfp = stdin;
	} else {
	    if (NULL == (fofnfp = fopen(fofn, "r"))) {
		perror(fofn);
		return 1;
	    }
	}

	while (fgets(file, 255, fofnfp)) {
	    char *c;
	    if (c = strchr(file, '\n'))
		*c = 0;

	    extract(hf, file);
	}

	fclose(fofnfp);
    }

    for (; argc; argc--, argv++) {
	extract(hf, *argv);
    }

    HashFileDestroy(hf);

    return 0;
}
