#include <stdio.h>
#include <errno.h>
#include <stdlib.h>

#include "seqIOABI.h"

int main(int argc, char **argv) {
    FILE *fp;
    uint_4 indexO;
    char str[256];
    uint_4 label;
    int count = 1;

    if (argc < 2 || argc > 4) {
	fprintf(stderr, "Usage: %s [ident [count]] file\n", argv[0]);
	return 1;
    }
    
    if (NULL == (fp = fopen(argv[argc-1], "rb"))) {
	perror(argv[argc-1]);
	return 2;
    }

    if (-1 == getABIIndexOffset(fp, &indexO)) {
	fprintf(stderr, "Couldn't find ABI file index\n");
	return 3;
    }

    if (argc == 2) {
	return dump_labels(fp, (off_t)indexO);
    }

    label = ((((argv[1][0]<<8) + argv[1][1])<<8) + argv[1][2]<<8) + argv[1][3];

    if (argc == 4)
	count = atoi(argv[2]);

    if (-1 != getABIString(fp, (off_t)indexO, label, count, str))
	puts(str);
    else {
	fprintf(stderr, "Couldn't find that identifier\n");
	return 4;
    }

    return 0;
}
