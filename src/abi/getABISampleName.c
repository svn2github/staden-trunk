#include <stdio.h>
#include <errno.h>
#include <sys/types.h>

#include "seqIOABI.h"

int main(int argc, char **argv) {
    FILE *fp;
    uint_4 indexO;
    char name[256];

    if (argc != 2) {
	fprintf(stderr, "Usage: %s file\n", argv[0]);
	return 1;
    }
    
    if (NULL == (fp = fopen(argv[1], "rb"))) {
	perror(argv[1]);
	return 2;
    }

    if (-1 == getABIIndexOffset(fp, &indexO)) {
	fprintf(stderr, "Couldn't find ABI file index\n");
	return 3;
    }

    if (-1 != getABIString(fp, (off_t)indexO, SMPLLabel, 1, name))
	puts(name);

    return 0;
}
