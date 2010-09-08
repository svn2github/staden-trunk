#ifndef _TG_UTILS_H_
#define _TG_UTILS_H_

#include "tg_struct.h"

/*
 * Packing and unpacking of unsigned and signed integers of variable size.
 * These fit values from 0-127 (or -63 to 63) in 1 byte, 0-16383 in 2 bytes,
 * etc.
 *
 * They return the number of bytes written to or read from 'out'.
 */
int int2u7(uint32_t in, unsigned char *out);
int u72int(unsigned char *u7, uint32_t *out);
int int2s7(int32_t in, unsigned char *out);
int s72int(unsigned char *u7, int32_t *out);

/*
 * As above, but 64-bit versions. Sorry about the names but this avoids
 * oddities like int642s7.
 */
int intw2u7(uint64_t in, unsigned char *out);
int u72intw(unsigned char *u7, uint64_t *out);
int intw2s7(int64_t in, unsigned char *out);
int s72intw(unsigned char *u7, int64_t *out);

/*
 * Macro versions of above - about 35% faster.
 *
 * They modify the variables directly. Ie instead of:
 *   cp += s72int(cp, &i);
 * use
 *   S72INT(cp, i);   //increments cp too
 *
 * WARNING: not all behave well with complex expressions as arguments.
 */

#define INT2U7(in,cp) do {			\
    uint32_t _i = (in);				\
    while (_i >= 128) {				\
	*(cp)++ = (_i & 0x7f) | 0x80;		\
	_i >>= 7;				\
    }						\
    *(cp)++ = _i;				\
} while (0)

#define U72INT(u7,o) do {			\
    int _b = 0;					\
    o = *(u7) & 0x7f;				\
    while (*(u7)++ & 0x80) {			\
	o |= (*(u7) & 0x7f) << (_b += 7);	\
    }						\
} while(0)

#define INT2S7(in,cp) do {			\
    uint32_t _i = (ABS(in) << 1) | ((in) < 0);	\
    while (_i >= 128) {				\
	*(cp)++ = (_i & 0x7f) | 0x80;		\
	_i >>= 7;				\
    }						\
    *(cp)++ = _i;				\
} while (0)

#define S72INT(u7,o) do {			\
    int _b = 0;					\
    o = *(u7) & 0x7f;				\
    while (*(u7)++ & 0x80) {			\
	o |= (*(u7) & 0x7f) << (_b += 7);	\
    }						\
    o = (o & 1) ? -(o >> 1) : (o >> 1);		\
} while(0)

#define INTW2U7(in,cp) do {			\
    uint64_t _i = (in);				\
    while (_i >= 128) {				\
	*(cp)++ = (_i & 0x7f) | 0x80;		\
	_i >>= 7;				\
    }						\
    *(cp)++ = _i;				\
} while (0)

#define U72INTW(u7,o) U72INT((u7),(o))

#define INTW2S7(in,cp) do {			\
    uint64_t _i = (ABS(in) << 1) | ((in) < 0);	\
    while (_i >= 128) {				\
	*(cp)++ = (_i & 0x7f) | 0x80;		\
	_i >>= 7;				\
    }						\
    *(cp)++ = _i;				\
} while (0)

#define S72INTW(u7,o) S72INT((u7),(o))

/* Like atoi() but for tg_rec, which is typically atol */
tg_rec atorec(const char *str);

#endif /* _TG_UTILS_H_ */
