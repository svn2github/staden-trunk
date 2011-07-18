#include <string.h>
#include "misc.h"
#include "tg_gio.h"

/*
 * Converts an unsigned integer value into a 7-bit encoded format.
 * We store the bottom 7 bits of value with either 0 or 1 for the top-bit
 * depending on whether any bits are left. We keep repeating this until
 * all significant bits of value have been used.
 *
 * Ie 15551 = hex 3cbf = 0011 1100 1011 1111 becomes:
 *
 *     111 1001   011 1111 (0x3cbf input)
 *    0111 1001  1011 1111 (0x79bf output)
 *
 * Takes an unsigned 32-bit integer and stores in out.
 * Returns the number of bytes written to 'out'
 */
int int2u7(uint32_t in, unsigned char *out) {
    unsigned char *cp = out;

    while (in >= 128) {
	*cp++ = (in & 0x7f) | 0x80;
	in >>= 7;
    }
    *cp++ = in;

    return cp-out;
}

/*
 * Takes a 7-bit encoded value in 'u7' and stores in a
 * 32-bit unsigned int pointed to by 'out'.
 *
 * Returns the number of bytes read from u7.
 */
int u72int(unsigned char *u7, uint32_t *out) {
    uint32_t ret = 0;
    int b = 0, c = 0;

    ret = *u7 & 0x7f;
    while (*u7++ & 0x80) {
	c++;
	ret |= ((uint64_t)(*u7 & 0x7f)) << (b += 7);
    }

    *out = ret;
    return c+1;
}

/*
 * Converts an signed value into a 7-bit encoded format.
 * We store the bottom 7 bits of value with either 0 or 1 for the top-bit
 * depending on whether any bits are left. We keep repeating this until
 * all significant bits of value have been used.
 *
 * Ie 15551 = hex 3cbf = 0011 1100 1011 1111 becomes:
 *
 *     111 1001   011 1111 (0x3cbf input)
 *    0111 1001  1011 1111 (0x79bf output)
 *
 * Takes an unsigned 32-bit integer and stores in out.
 * Returns the number of bytes written to 'out'
 */
int int2s7(int32_t in, unsigned char *out) {
    unsigned char *cp = out;
    uint32_t u = (ABS(in) << 1) | (in < 0);

    while (u >= 128) {
	*cp++ = (u & 0x7f) | 0x80;
	u >>= 7;
    }
    *cp++ = u;

    return cp-out;
}

/*
 * Takes a 7-bit encoded value in 'u7' and stores in a
 * 32-bit signed int pointed to by 'out'.
 *
 * Returns the number of bytes read from u7.
 */
int s72int(unsigned char *u7, int32_t *out) {
    uint32_t ret;
    int b = 0, c = 0;

    ret = *u7 & 0x7f;
    while (*u7++ & 0x80) {
	c++;
	ret |= (*u7 & 0x7f) << (b += 7);
    }

    *out = (ret & 1) ? -(ret >> 1) : (ret >> 1);
    return c+1;
}


/*
 * Converts an unsigned 64-bit integer value into a 7-bit encoded format.
 * We store the bottom 7 bits of value with either 0 or 1 for the top-bit
 * depending on whether any bits are left. We keep repeating this until
 * all significant bits of value have been used.
 *
 * Ie 15551 = hex 3cbf = 0011 1100 1011 1111 becomes:
 *
 *     111 1001   011 1111 (0x3cbf input)
 *    0111 1001  1011 1111 (0x79bf output)
 *
 * Takes an unsigned 64-bit integer and stores in out.
 * Returns the number of bytes written to 'out'
 */
int intw2u7(uint64_t in, unsigned char *out) {
    unsigned char *cp = out;

    while (in >= 128) {
	*cp++ = (in & 0x7f) | 0x80;
	in >>= 7;
    }
    *cp++ = in;

    return cp-out;
}

/*
 * Takes a 7-bit encoded value in 'u7' and stores in a
 * 64-bit unsigned int pointed to by 'out'.
 *
 * Returns the number of bytes read from u7.
 */
int u72intw(unsigned char *u7, uint64_t *out) {
    uint64_t ret = 0;
    int b = 0, c = 0;

    ret = *u7 & 0x7f;
    while (*u7++ & 0x80) {
	c++;
	ret |= (*u7 & 0x7f) << (b += 7);
    }

    *out = ret;
    return c+1;
}

/*
 * Converts an signed value into a 7-bit encoded format.
 * We store the bottom 7 bits of value with either 0 or 1 for the top-bit
 * depending on whether any bits are left. We keep repeating this until
 * all significant bits of value have been used.
 *
 * Ie 15551 = hex 3cbf = 0011 1100 1011 1111 becomes:
 *
 *     111 1001   011 1111 (0x3cbf input)
 *    0111 1001  1011 1111 (0x79bf output)
 *
 * Takes an unsigned 32-bit integer and stores in out.
 * Returns the number of bytes written to 'out'
 */
int intw2s7(int64_t in, unsigned char *out) {
    unsigned char *cp = out;
    uint64_t u = (ABS(in) << 1) | (in < 0);

    while (u >= 128) {
	*cp++ = (u & 0x7f) | 0x80;
	u >>= 7;
    }
    *cp++ = u;

    return cp-out;
}

/*
 * Takes a 7-bit encoded value in 'u7' and stores in a
 * 32-bit signed int pointed to by 'out'.
 *
 * Returns the number of bytes read from u7.
 */
int s72intw(unsigned char *u7, int64_t *out) {
    uint64_t ret;
    int b = 0, c = 0;

    ret = *u7 & 0x7f;
    while (*u7++ & 0x80) {
	c++;
	ret |= ((uint64_t)(*u7 & 0x7f)) << (b += 7);
    }

    *out = (ret & 1) ? -(ret >> 1) : (ret >> 1);
    return c+1;
}


/* Like atoi() but for tg_rec, which is typically atol */
tg_rec atorec(const char *str) {
    if (sizeof(tg_rec) == sizeof(int)) {
	return atoi(str);
    } else if (sizeof(tg_rec) == sizeof(long)) {
	return atol(str);
    } else if (sizeof(tg_rec) == sizeof(long long)) {
	return atoll(str);
    } else {
	fprintf(stderr, "no suitable ato*() func for tg_rec\n");
	return -1;
    }
}

