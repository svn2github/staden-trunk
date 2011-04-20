/*
 * This file implements various atoi/strtol implementations for integers
 * of known sizes, rather than variable sized "long" vs "long long"
 * implementations.
 */

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

#include "misc.h"


int32_t strtol32(const char *nptr, char **endptr, int base) {
    return strtol(nptr, endptr, base);
}

int64_t strtol64(const char *nptr, char **endptr, int base) {
    if (sizeof(int) >= 8 || sizeof(long) >= 8) {
	return strtol(nptr, endptr, base);
    } else if (sizeof(long long) >= 8) {
	return strtoll(nptr, endptr, base);
    } else {
	fprintf(stderr, "Unknown strto* function for 64-bit types\n");
	return strtoll(nptr, endptr, base);
    }
}

int32_t atoi32(const char *nptr) {
    return strtol32(nptr, (char **)NULL, 10);
}

int32_t atoi64(const char *nptr) {
    return strtol64(nptr, (char **)NULL, 10);
}
