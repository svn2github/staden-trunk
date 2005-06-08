#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"

/* Copies a single named file to stdout */
void extract(HashFile *hf, char *file) {
    uint64_t pos;
    uint32_t size;
    char *data;

    // printf("==== %s ====\n", file);

    /* Obtain the pos+size */
    if (-1 == HashFileQuery(hf, (uint8_t *)file, strlen(file), &pos, &size))
	return;

    /* Copy the data to stdout */
    data = (char *)malloc(size);
    fseek(hf->afp, pos, SEEK_SET);
    fread(data, size, 1, hf->afp);
    fwrite(data, size, 1, stdout);
    free(data);
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

    HashFileClose(hf);

    return 0;
}
