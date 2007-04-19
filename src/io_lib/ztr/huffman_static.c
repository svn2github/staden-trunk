#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "huffman_static.h"
#include "ztr.h"

typedef struct {
    signed int symbol; /* 0-255 character, 256 = exception code, 257 = EOF */
    int nbits;
    unsigned int code;
} huffman_code_t;

#define SYM_ANY 256
#define SYM_EOF 257

/* For DNA */
static huffman_code_t codes_dna[] = {
    {SYM_ANY, 5,}, {SYM_EOF, 5,},
    {'A',  2,}, {'C',  3,}, {'G',  2,}, {'T',  2,}, {'-',  4,}
};

/* For Traces */
static huffman_code_t codes_traces[] = {
    {SYM_ANY, 13,}, {SYM_EOF, 13,},
    {  0,  1,}, {  1,  4,}, {  2,  6,}, {  3,  6,}, {  4,  7,}, {  5,  7,},
    {  6,  7,}, {  7,  7,}, {  8,  7,}, {  9,  7,}, { 10,  8,}, { 11,  8,},
    { 12,  8,}, { 13,  8,}, { 14,  8,}, { 15,  8,}, { 16,  8,}, { 17,  8,},
    { 18,  8,}, { 19,  8,}, { 20,  8,}, { 21,  8,}, { 22,  8,}, { 23,  8,},
    { 24,  8,}, { 25,  8,}, { 26,  8,}, { 27,  8,}, { 28,  8,}, { 29,  8,},
    { 30,  8,}, { 31,  9,}, { 32,  9,}, { 33,  9,}, { 34,  9,}, { 35,  9,},
    { 36,  9,}, { 37,  9,}, { 38,  9,}, { 39,  9,}, { 40,  9,}, { 41,  9,},
    { 42,  9,}, { 43,  9,}, { 44,  9,}, { 45,  9,}, { 46,  9,}, { 47,  9,},
    { 48,  9,}, { 49,  9,}, { 50,  9,}, { 51,  9,}, { 52,  9,}, { 53,  9,},
    { 54,  9,}, { 55,  9,}, { 56,  9,}, { 57,  9,}, { 58,  9,}, { 59,  9,},
    { 60,  9,}, { 61,  9,}, { 62,  9,}, { 63,  9,}, { 64,  9,}, { 65,  9,},
    { 66,  9,}, { 67,  9,}, { 68,  9,}, { 69,  9,}, { 70,  9,}, { 71,  9,},
    { 72,  9,}, { 73,  9,}, { 74,  9,}, { 75,  9,}, { 76,  9,}, { 77,  9,},
    { 78,  9,}, { 79,  9,}, { 80,  9,}, { 81,  9,}, { 82,  9,}, { 83,  9,},
    { 84,  9,}, { 85,  9,}, { 86,  9,}, { 87,  9,}, { 88,  9,}, { 89,  9,},
    { 90,  9,}, { 91,  9,}, { 92, 10,}, { 93, 10,}, { 94, 10,}, { 95, 10,},
    { 96, 10,}, { 97, 10,}, { 98, 10,}, { 99, 10,}, {100, 10,}, {101, 10,},
    {102, 10,}, {103, 10,}, {104, 10,}, {105, 10,}, {106, 10,}, {107, 10,},
    {108, 10,}, {109, 10,}, {110, 10,}, {111, 10,}, {112, 10,}, {113, 10,},
    {114, 10,}, {115, 10,}, {116, 10,}, {117, 10,}, {118, 10,}, {119, 10,},
    {120, 10,}, {121, 10,}, {122, 10,}, {123, 10,}, {124, 10,}, {125, 10,},
    {126, 10,}, {127, 10,}, {128, 10,}, {129, 10,}, {130, 10,}, {131, 10,},
    {132, 10,}, {133, 10,}, {134, 10,}, {135, 10,}, {136, 10,}, {137, 10,},
    {138, 10,}, {139, 10,}, {140, 10,}, {141, 10,}, {142, 10,}, {143, 10,},
    {144, 10,}, {145, 10,}, {146, 10,}, {147, 10,}, {148, 10,}, {149, 10,},
    {150, 10,}, {151, 10,}, {152, 10,}, {153, 10,}, {154, 10,}, {155, 10,},
    {156, 10,}, {157, 10,}, {158, 10,}, {159, 10,}, {160, 10,}, {161, 10,},
    {162, 10,}, {163, 10,}, {164, 10,}, {165, 10,}, {166, 10,}, {167, 10,},
    {168, 10,}, {169, 10,}, {170, 10,}, {171, 10,}, {172, 10,}, {173, 10,},
    {174, 10,}, {175, 10,}, {176, 10,}, {177, 10,}, {178, 10,}, {179, 10,},
    {180, 10,}, {181, 10,}, {182, 10,}, {183, 10,}, {184, 10,}, {185, 10,},
    {186, 10,}, {187, 10,}, {188, 10,}, {189, 10,}, {190, 10,}, {191, 10,},
    {192, 10,}, {193, 10,}, {194, 10,}, {195, 10,}, {196, 10,}, {197, 10,},
    {198, 10,}, {199, 10,}, {200, 10,}, {201, 10,}, {202, 10,}, {203, 10,},
    {204, 10,}, {205, 10,}, {206, 10,}, {207, 10,}, {208, 10,}, {209, 10,},
    {210, 10,}, {211, 10,}, {212, 10,}, {213, 10,}, {214, 10,}, {215, 10,},
    {216, 10,}, {217, 10,}, {218, 10,}, {219, 10,}, {220, 10,}, {221, 10,},
    {222, 10,}, {223, 10,}, {224, 10,}, {225, 10,}, {226, 10,}, {227, 10,},
    {228, 10,}, {229, 10,}, {230, 10,}, {231, 10,}, {232, 10,}, {233, 10,},
    {234, 10,}, {235, 10,}, {236, 10,}, {237, 10,}, {238, 10,}, {239, 10,},
    {240, 10,}, {241, 10,}, {242, 10,}, {243, 10,}, {244, 10,}, {245, 10,},
    {246, 10,}, {247, 10,}, {248, 10,}, {249, 10,}, {250, 11,}, {251, 10,},
    {252, 11,}, {253, 11,}, {254, 12,}, {255, 10,}
};

/* For 4x uncompressed confidence values */
static huffman_code_t codes_conf[] = {
    {SYM_ANY, 13,}, {SYM_EOF, 13,},
    {  0,  6,}, {  1,  6,}, {  2,  7,}, {  3,  7,}, {  4,  7,}, {  5,  7,},
    {  6,  7,}, {  7,  7,}, {  8,  8,}, {  9,  8,}, { 10,  8,}, { 11,  9,},
    { 12,  9,}, { 13,  9,}, { 14,  9,}, { 15,  9,}, { 16, 10,}, { 17,  9,},
    { 18, 10,}, { 19, 12,}, { 20,  9,}, { 21, 11,}, { 22, 10,}, { 23, 10,},
    { 24, 11,}, { 25, 11,}, { 27,  9,}, { 30,  3,}, {226,  1,}, {229,  7,},
    {231,  7,}, {232,  8,}, {233,  8,}, {234,  7,}, {235,  9,}, {236,  7,},
    {237,  8,}, {238,  7,}, {239,  7,}, {240,  7,}, {241,  7,}, {242,  7,},
    {243,  7,}, {244,  7,}, {245,  7,}, {246,  6,}, {247,  6,}, {248,  6,},
    {249,  6,}, {250,  6,}, {251,  6,}, {252,  6,}, {253,  6,}, {254,  6,},
    {255,  6,}
};

/* For RLE compressed solexa confidence values */
static huffman_code_t codes_conf_rle[] = {
    {SYM_ANY, 20,}, {SYM_EOF, 20,},
    {  0,  3,}, {  1,  5,}, {  2,  6,}, {  3,  6,}, {  4,  6,}, {  5,  6,},
    {  6,  6,}, {  7,  7,}, {  8,  7,}, {  9,  7,}, { 10,  7,}, { 11,  8,},
    { 12,  8,}, { 13,  8,}, { 14,  7,}, { 15,  8,}, { 16,  8,}, { 17,  8,},
    { 18,  9,}, { 19, 10,}, { 20,  8,}, { 21,  9,}, { 22,  9,}, { 23,  8,},
    { 24, 10,}, { 25,  9,}, { 26, 10,}, { 27,  8,}, { 28, 13,}, { 29, 12,},
    { 30,  4,}, { 31, 13,}, { 32, 11,}, { 33, 12,}, { 34, 13,}, { 35, 11,},
    { 36, 11,}, { 37, 13,}, { 38, 12,}, { 39, 10,}, { 40, 12,}, { 41, 13,},
    { 42, 10,}, { 43, 11,}, { 44, 12,}, { 45, 10,}, { 46, 13,}, { 47, 14,},
    { 48, 11,}, { 49, 13,}, { 50, 15,}, { 51, 12,}, { 52, 17,}, { 53, 15,},
    { 54, 12,}, { 55, 16,}, { 56, 16,}, { 57, 12,}, { 58, 15,}, { 59, 15,},
    { 60, 13,}, { 61, 15,}, { 62, 15,}, { 63, 12,}, { 64, 14,}, { 65, 16,},
    { 66, 12,}, { 67, 15,}, { 68, 16,}, { 69, 12,}, { 70, 13,}, { 71, 15,},
    { 72, 12,}, { 73, 14,}, { 74, 15,}, { 75, 12,}, { 76, 13,}, { 77,  4,},
    { 78, 12,}, { 79, 13,}, { 80, 15,}, { 81,  9,}, { 83, 18,}, { 95, 18,},
    { 97, 18,}, { 99, 17,}, {100, 18,}, {101, 18,}, {102, 18,}, {109,  6,},
    {114, 18,}, {115, 17,}, {116, 18,}, {117, 18,}, {122, 19,}, {226,  3,},
    {229,  6,}, {231,  7,}, {232,  7,}, {233,  7,}, {234,  7,}, {235,  8,},
    {236,  6,}, {237,  7,}, {238,  7,}, {239,  7,}, {240,  6,}, {241,  6,},
    {242,  6,}, {243,  6,}, {244,  6,}, {245,  6,}, {246,  6,}, {247,  6,},
    {248,  6,}, {249,  6,}, {250,  6,}, {251,  5,}, {252,  5,}, {253,  5,},
    {254,  6,}, {255,  5,}
};

static huffman_code_t *codes = NULL;
static int ncodes = 0;
static int code_set = 0;

static huffman_code_t *lookup[258];

typedef struct block {
    unsigned char *data;
    size_t alloc;
    size_t byte;
    int bit;
} block_t;

typedef struct {
    unsigned int code_start;
    int ncodes;
    unsigned int *symbols;
} decode_t;
static decode_t decode_tab[32];

void print_bits(unsigned int val, int nbits) {
    unsigned int mask = 1 << (nbits-1);
    do {
	printf("%d", (val & mask) ? 1 : 0);
	mask >>= 1;
    } while(--nbits);
}

/*
 * Sort huffman_code_t by their code bit-lengths
 */
int sort_func(const void *p1, const void *p2) {
    const huffman_code_t *c1 = (const huffman_code_t *)p1;
    const huffman_code_t *c2 = (const huffman_code_t *)p2;
    return c1->nbits - c2->nbits;
}

/*
 * Generates canonical huffman codes based on their bit-lengths.
 * This has useful encoding/decoding properties, but also means we only
 * need to store the bit-lengths to regenerate the same code set.
 */
void canonical_codes(huffman_code_t *codes, int ncodes) {
    int i;
    unsigned int code, last_len;
    huffman_code_t *exception = NULL;

    /* Sort by bit-length */
    qsort(codes, ncodes, sizeof(*codes), sort_func);

    /* Initialise code length lookup tables */
    for (i = 0; i <= 31; i++) {
	decode_tab[i].code_start = 0;
	decode_tab[i].ncodes = 0;
	decode_tab[i].symbols = NULL;
    }

    /* Generate codes */
    for (i = 0; i < ncodes; i++) {
	int nbits = codes[i].nbits;

	if (i == 0) {
	    code = 0;
	    last_len = nbits;
	    decode_tab[nbits].symbols = (unsigned int *)malloc(ncodes * 
							       sizeof(int));
	    decode_tab[nbits].code_start = code;
	} else {
	    code++;
	}
	if (nbits > last_len) {
	    code <<= (nbits - last_len);
	    last_len = nbits;
	    decode_tab[nbits].symbols = (unsigned int *)malloc(ncodes * 
							       sizeof(int));
	    decode_tab[nbits].code_start = code;
	}
	codes[i].code = code;
	/*
	 * printf("CODE %c %d ", codes[i].symbol, codes[i].nbits, code);
	 * print_bits(codes[i].code, codes[i].nbits);
	 * printf("\n");
	 */

	decode_tab[nbits].symbols[decode_tab[nbits].ncodes++]
	    = codes[i].symbol;
    }

    /* Produce fast lookup table */
    for (i = 0; i < 256; i++) {
	lookup[i] = NULL;
    }
    for (i = 0; i < ncodes; i++) {
        lookup[codes[i].symbol] = &codes[i];
	if (codes[i].symbol == SYM_ANY)
	    exception = &codes[i];
    }
    for (i = 0; i < 256; i++) {
	if (!lookup[i])
	    lookup[i] = exception;
    }
}

void store_bits(block_t *block, unsigned int val, int nbits) {
    /* Slow, but simple */
    unsigned int mask = 1 << (nbits-1);
    do {
	int bit = 1 << block->bit;
	if (block->byte >= block->alloc) {
	    block->alloc += 8192;
	    block->data = (unsigned char *)realloc(block->data, block->alloc);
	}

	if (val & mask)
	    block->data[block->byte] |= bit;
	else
	    block->data[block->byte] &= ~bit;

	if (++block->bit == 8) {
	    block->bit = 0;
	    block->byte++;
	}
	mask >>= 1;
    } while(--nbits);
}

/* Returns bit (0 or 1) on success
 *         -1 on failure or EOF.
 */
int get_bit(block_t *block) {
    int b;

    if (block->byte >= block->alloc)
	return -1;

    b = block->data[block->byte] & (1 << block->bit);
    if (8 == ++block->bit) {
	block->bit = 0;
	block->byte++;
    }

    return b ? 1 : 0;
}

static int huffman_init(int code) {
    if (code_set != code) {
	switch(code) {
	case CODE_DNA:
	    /* fprintf(stderr, "using DNA\n"); */
	    codes = codes_dna;
	    ncodes = sizeof(codes_dna)/sizeof(*codes);
	    break;

	case CODE_TRACES:
	    /* fprintf(stderr, "using TRACES\n"); */
	    codes = codes_traces;
	    ncodes = sizeof(codes_traces)/sizeof(*codes);
	    break;

	case CODE_CONF:
	    /* fprintf(stderr, "using CONF\n"); */
	    codes = codes_conf;
	    ncodes = sizeof(codes_conf)/sizeof(*codes);
	    break;

	case CODE_CONF_RLE:
	    /* fprintf(stderr, "using CONF_RLE\n"); */
	    codes = codes_conf_rle;
	    ncodes = sizeof(codes_conf_rle)/sizeof(*codes);
	    break;

	default:
	    fprintf(stderr, "Unknown huffman code set '%d'\n",
		    code_set);
	    return -1;
	}

	canonical_codes(codes, ncodes);
	code_set = code;
    }

    return 0;
}

unsigned char *encode_memory(unsigned char *data, int len, int *nbytes,
			     int code_set) {
    int i;
    block_t out;

    huffman_init(code_set);

    out.data = malloc(len+2);
    out.alloc = len+2;
    out.byte = 0;
    out.data[out.byte++] = ZTR_FORM_STHUFF;
    out.data[out.byte++] = code_set;
    out.bit = 0;

    for (i = 0; i < len; i++) {
	/* print_bits(lookup[data[i]]->code, lookup[data[i]]->nbits); */
	if (lookup[data[i]]) {
	    store_bits(&out, lookup[data[i]]->code, lookup[data[i]]->nbits);
	    if (lookup[data[i]]->symbol == SYM_ANY) {
		store_bits(&out, data[i], 8);
	    }
	}
	/* printf("%c %d %d\n", data[i], lookup[data[i]]->nbits, lookup[data[i]]->code); */
    }
    store_bits(&out, lookup[SYM_EOF]->code, lookup[SYM_EOF]->nbits);

    *nbytes = out.byte + (out.bit != 0);
    out.data = realloc(out.data, out.byte);
    return out.data;
}

unsigned char *decode_memory(unsigned char *data, int len, int *nbytes) {
    block_t in, out;
    int b, nbits = 0;
    unsigned int val = 0;

    in.data = data;
    in.alloc = len;
    in.byte = 2;
    in.bit = 0;

    out.data = NULL;
    out.alloc = 0;
    out.byte = 0;
    out.bit = 0;

    huffman_init(in.data[1]);

    while((b = get_bit(&in)) != -1) {
	int sym;
	val = (val << 1) + b;
	nbits++;
	if (val >= decode_tab[nbits].code_start &&
	    val <  decode_tab[nbits].code_start + decode_tab[nbits].ncodes) {
	    sym = decode_tab[nbits].symbols[val-decode_tab[nbits].code_start];
	    if (sym == SYM_ANY) {
		sym = (get_bit(&in) << 7)
		    + (get_bit(&in) << 6)
		    + (get_bit(&in) << 5)
		    + (get_bit(&in) << 4)
		    + (get_bit(&in) << 3)
		    + (get_bit(&in) << 2)
		    + (get_bit(&in) << 1)
		    + (get_bit(&in) << 0);
	    } else if (sym == SYM_EOF) {
		break;
	    }

	    if (out.byte >= out.alloc) {
		out.alloc += 8192;
		out.data = (unsigned char *)realloc(out.data, out.alloc);
	    }
	    out.data[out.byte++] = sym;
	    nbits = val = 0;
	}
    }

    *nbytes = out.byte;
    return out.data;
}

unsigned char *load(int *lenp) {
    unsigned char *data = NULL;
    int dsize = 0;
    int dcurr = 0, len;

    do {
	if (dsize - dcurr < 8192) {
	    dsize = dsize ? dsize * 2 : 8192;
	    data = realloc(data, dsize);
	}

	len = read(0, data + dcurr, 8192);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }

    *lenp = dcurr;
    return data;
}

void dump_decode_table(void) {
    int i, j;
    for (i = 0; i < 31; i++) {
	if (!decode_tab[i].ncodes)
	    continue;
	
	printf("=== len %d ==\n", i);
	for (j = 0; j < decode_tab[i].ncodes; j++) {
	    printf("%c ",
		   isprint(decode_tab[i].symbols[j])
		   ? decode_tab[i].symbols[j]
		   : '?');
	    print_bits(decode_tab[i].code_start + j, i);
	    printf("\n");
	}
    }
}

#if 0
int main(int argc, char **argv) {
    unsigned char *data, *out;
    int len, out_len;
    int decode = 0;

    if (argc > 1 && strcmp(argv[1], "-d") == 0) {
	decode = 1;
    }
    data = load(&len);

    if (!decode) {
	out = encode_memory(data, len, &out_len, CODE_DNA);
    } else {
	out = decode_memory(data, len, &out_len);
    }
    
    write(1, out, out_len);

    return 0;
}
#endif
