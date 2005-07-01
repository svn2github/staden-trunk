#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "abiIO.h"

int main(int argc, char **argv) {
    abi_t *abi;
    abi_index_t *ind1, *ind2;
    uint32_t rund, runt;
    int verbose = 0;
    char *prog = argv[0];

    if (argc >= 2 && 0 == strcmp(argv[1], "-v")) {
	verbose=1;
	argc--;
	argv++;
    }
    if (argc != 2) {
	fprintf(stderr, "Usage: %s file\n", prog);
	return 1;
    }
    
    if (NULL == (abi = read_abi(argv[1]))) {
	perror(argv[1]);
	return 2;
    }

    if (NULL == (ind1 = find_abi_index(abi, ABI_LABEL("RUND"), 1)) ||
	NULL == (ind2 = find_abi_index(abi, ABI_LABEL("RUNT"), 1))) {
	fprintf(stderr, "Couldn't find RUND/RUNT labels\n");
	return 3;
    }

    rund = *((uint32_t *)ind1->data);
    runt = *((uint32_t *)ind2->data);

    printf((verbose ?
	    "year %04d, month %02d, day %02d, time %02d:%02d:%02d\n" :
	    "%04d%02d%02d.%02d%02d%02d\n"),
	   ((rund & 0xff) << 8) | ((rund >> 8) & 0xff),
	   (rund >> 16) & 0xff,
	   (rund >> 24) & 0xff,
	   runt & 0xff,
	   (runt >> 8) & 0xff,
	   (runt >> 16) & 0xff);

    del_abi_t(abi);

    return 0;
}

