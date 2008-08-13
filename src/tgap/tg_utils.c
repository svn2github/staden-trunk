#include <string.h>
#include "misc.h"
#include "tg_gio.h"

/*
 * Converts an integer value into a 7-bit encoded format.
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
    int n = 0;

    out[0] = in & 0x7f;
    while (in >= 128) {
	out[n++] |= 128;
	in >>= 7;
	out[n] = in & 0x7f;
    }

    return n+1;
}

/*
 * Takes a 7-bit encoded value in 'u7' and stores in a
 * 32-bit unsigned int pointed to by 'out'.
 *
 * Returns the number of bytes read from u7.
 */
int u72int(unsigned char *u7, uint32_t *out) {
    uint32_t ret = 0;
    int b = 0;

    ret |= *u7 & 0x7f;
    while (*u7 & 0x80) {
	u7++;
	b += 7;
	ret |= (*u7 & 0x7f) << b;
    }

    *out = ret;
    return b/7+1;
}
