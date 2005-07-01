#include <stdio.h>
#include <errno.h>

#include "abiIO.h"

int main(int argc, char **argv) {
    abi_t *abi;
    abi_index_t *ind;

    if (argc != 2) {
	fprintf(stderr, "Usage: %s file\n", argv[0]);
	return 1;
    }
    
    if (NULL == (abi = read_abi(argv[1]))) {
	perror(argv[1]);
	return 2;
    }

    if (NULL == (ind = find_abi_index(abi, ABI_LABEL("CMNT"), 1))) {
	fprintf(stderr, "Couldn't find CMNT label\n");
	return 3;
    }

    printf("%.*s\n", *ind->data, ind->data + 1);

    del_abi_t(abi);

    return 0;
}
