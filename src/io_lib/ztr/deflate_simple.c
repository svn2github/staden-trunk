/*
 * This code is a simplistic implementation of the Deflate algorithm.
 * See http://www.ietf.org/rfc/rfc1951.txt for details on this.
 *
 * The simplification is that it ONLY supports the huffman encoding
 * step and does not attempt to do any LZ-style string matching to generate
 * distance codes. (These generally do not improve data compression for our
 * desired use.)
 *
 * It has been written here, instead of using zlib, so that we can separate
 * out the encoding of the huffman tree from the compression of the data
 * stream into separate memory sections with the intent to optimise
 * compression of very small blocks of data by sharing one set of frequency
 * tables (ie huffman tree) with multiple sets of compressed data blocks.
 *
 * Old decompression routes (huffman_decode) can be found here too, but
 * the restore_codes() function will need a rewrite before this can be used.
 * Instead for now you are encouraged to use zlib for decompression.
 *
 * James Bonfield, 2007
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <unistd.h>

#include "deflate_simple.h"


/* #define TIME_IT */
/* #define TEST_MAIN */

#ifdef TIME_IT
#  include <sys/time.h>
#endif

/*
 * ---------------------------------------------------------------------------
 * Local structs & defines
 */

/* Used in tree construction only */
typedef struct node {
    int count;
    int sym; /* char or SYM_EOF */
    struct node *parent;
    struct node *next;
} node_t;

#define SYM_EOF 256

static void output_code_set(FILE *fp, huffman_code_t *codes, int ncodes);
static void output_code_set2(FILE *fp, huffman_code_t *codes, int ncodes);

/*
 * ---------------------------------------------------------------------------
 * Our standard precomputed tables, for DNA, English text, etc.
 */

/* DNA */
/*
 * A   00
 * C   01
 * G   110
 * T   10
 * N   1110
 * EOF 11110
 * ?   11111*
 */
static huffman_code_t codes_dna[] = {
    {'A',  2}, {'C',  2}, {'T',  2}, {'G',  3}, {'N',  4}, {  0,  5},
    {SYM_EOF,  6},
    {  1, 13}, {  2, 13}, {  3, 13}, {  4, 13}, {  5, 13}, {  6, 13},
    {  7, 14}, {  8, 14}, {  9, 14}, { 10, 14}, { 11, 14}, { 12, 14},
    { 13, 14}, { 14, 14}, { 15, 14}, { 16, 14}, { 17, 14}, { 18, 14},
    { 19, 14}, { 20, 14}, { 21, 14}, { 22, 14}, { 23, 14}, { 24, 14},
    { 25, 14}, { 26, 14}, { 27, 14}, { 28, 14}, { 29, 14}, { 30, 14},
    { 31, 14}, { 32, 14}, { 33, 14}, { 34, 14}, { 35, 14}, { 36, 14},
    { 37, 14}, { 38, 14}, { 39, 14}, { 40, 14}, { 41, 14}, { 42, 14},
    { 43, 14}, { 44, 14}, { 45, 14}, { 46, 14}, { 47, 14}, {'0', 14},
    {'1', 14}, {'2', 14}, {'3', 14}, {'4', 14}, {'5', 14}, {'6', 14},
    {'7', 14}, {'8', 14}, {'9', 14}, { 58, 14}, { 59, 14}, { 60, 14},
    { 61, 14}, { 62, 14}, { 63, 14}, { 64, 14}, {'B', 14}, {'D', 14},
    {'E', 14}, {'F', 14}, {'H', 14}, {'I', 14}, {'J', 14}, {'K', 14},
    {'L', 14}, {'M', 14}, {'O', 14}, {'P', 14}, {'Q', 14}, {'R', 14},
    {'S', 14}, {'U', 14}, {'V', 14}, {'W', 14}, {'X', 14}, {'Y', 14},
    {'Z', 14}, { 91, 14}, { 92, 14}, { 93, 14}, { 94, 14}, { 95, 14},
    { 96, 14}, {'a', 14}, {'b', 14}, {'c', 14}, {'d', 14}, {'e', 14},
    {'f', 14}, {'g', 14}, {'h', 14}, {'i', 14}, {'j', 14}, {'k', 14},
    {'l', 14}, {'m', 14}, {'n', 14}, {'o', 14}, {'p', 14}, {'q', 14},
    {'r', 14}, {'s', 14}, {'t', 14}, {'u', 14}, {'v', 14}, {'w', 14},
    {'x', 14}, {'y', 14}, {'z', 14}, {123, 14}, {124, 14}, {125, 14},
    {126, 14}, {127, 14}, {128, 14}, {129, 14}, {130, 14}, {131, 14},
    {132, 14}, {133, 14}, {134, 14}, {135, 14}, {136, 14}, {137, 14},
    {138, 14}, {139, 14}, {140, 14}, {141, 14}, {142, 14}, {143, 14},
    {144, 14}, {145, 14}, {146, 14}, {147, 14}, {148, 14}, {149, 14},
    {150, 14}, {151, 14}, {152, 14}, {153, 14}, {154, 14}, {155, 14},
    {156, 14}, {157, 14}, {158, 14}, {159, 14}, {160, 14}, {161, 14},
    {162, 14}, {163, 14}, {164, 14}, {165, 14}, {166, 14}, {167, 14},
    {168, 14}, {169, 14}, {170, 14}, {171, 14}, {172, 14}, {173, 14},
    {174, 14}, {175, 14}, {176, 14}, {177, 14}, {178, 14}, {179, 14},
    {180, 14}, {181, 14}, {182, 14}, {183, 14}, {184, 14}, {185, 14},
    {186, 14}, {187, 14}, {188, 14}, {189, 14}, {190, 14}, {191, 14},
    {192, 14}, {193, 14}, {194, 14}, {195, 14}, {196, 14}, {197, 14},
    {198, 14}, {199, 14}, {200, 14}, {201, 14}, {202, 14}, {203, 14},
    {204, 14}, {205, 14}, {206, 14}, {207, 14}, {208, 14}, {209, 14},
    {210, 14}, {211, 14}, {212, 14}, {213, 14}, {214, 14}, {215, 14},
    {216, 14}, {217, 14}, {218, 14}, {219, 14}, {220, 14}, {221, 14},
    {222, 14}, {223, 14}, {224, 14}, {225, 14}, {226, 14}, {227, 14},
    {228, 14}, {229, 14}, {230, 14}, {231, 14}, {232, 14}, {233, 14},
    {234, 14}, {235, 14}, {236, 14}, {237, 14}, {238, 14}, {239, 14},
    {240, 14}, {241, 14}, {242, 14}, {243, 14}, {244, 14}, {245, 14},
    {246, 14}, {247, 14}, {248, 14}, {249, 14}, {250, 14}, {251, 14},
    {252, 14}, {253, 14}, {254, 14}, {255, 14},
};

/* DNA with a few ambiguity codes */
static huffman_code_t codes_dna_ambig[] = {
    {'A',  2}, {'C',  2}, {'T',  2}, {'G',  3}, {'N',  4}, {  0,  7},
    { 45,  7}, {'B',  8}, {'D',  8}, {'H',  8}, {'K',  8}, {'M',  8},
    {'R',  8}, {'S',  8}, {'V',  8}, {'W',  8}, {'Y',  8}, {SYM_EOF, 11},
    {226, 14}, {  1, 15}, {  2, 15}, {  3, 15}, {  4, 15}, {  5, 15},
    {  6, 15}, {  7, 15}, {  8, 15}, {  9, 15}, { 10, 15}, { 11, 15},
    { 12, 15}, { 13, 15}, { 14, 15}, { 15, 15}, { 16, 15}, { 17, 15},
    { 18, 15}, { 19, 15}, { 20, 15}, { 21, 15}, { 22, 15}, { 23, 15},
    { 24, 15}, { 25, 15}, { 26, 15}, { 27, 15}, { 28, 15}, { 29, 15},
    { 30, 15}, { 31, 15}, { 32, 15}, { 33, 15}, { 34, 15}, { 35, 15},
    { 36, 15}, { 37, 15}, { 38, 15}, { 39, 15}, { 40, 15}, { 41, 15},
    { 42, 15}, { 43, 15}, { 44, 15}, { 46, 15}, { 47, 15}, {'0', 15},
    {'1', 15}, {'2', 15}, {'3', 15}, {'4', 15}, {'5', 15}, {'6', 15},
    {'7', 15}, {'8', 15}, {'9', 15}, { 58, 15}, { 59, 15}, { 60, 15},
    { 61, 15}, { 62, 15}, { 63, 15}, { 64, 15}, {'E', 15}, {'F', 15},
    {'I', 15}, {'J', 15}, {'L', 15}, {'O', 15}, {'P', 15}, {'Q', 15},
    {'U', 15}, {'X', 15}, {'Z', 15}, { 91, 15}, { 92, 15}, { 93, 15},
    { 94, 15}, { 95, 15}, { 96, 15}, {'a', 15}, {'b', 15}, {'c', 15},
    {'d', 15}, {'e', 15}, {'f', 15}, {'g', 15}, {'h', 15}, {'i', 15},
    {'j', 15}, {'k', 15}, {'l', 15}, {'m', 15}, {'n', 15}, {'o', 15},
    {'p', 15}, {'q', 15}, {'r', 15}, {'s', 15}, {'t', 15}, {'u', 15},
    {'v', 15}, {'w', 15}, {'x', 15}, {'y', 15}, {'z', 15}, {123, 15},
    {124, 15}, {125, 15}, {126, 15}, {127, 15}, {128, 15}, {129, 15},
    {130, 15}, {131, 15}, {132, 15}, {133, 15}, {134, 15}, {135, 15},
    {136, 15}, {137, 15}, {138, 15}, {139, 15}, {140, 15}, {141, 15},
    {142, 15}, {143, 15}, {144, 15}, {145, 15}, {146, 15}, {147, 15},
    {148, 15}, {149, 15}, {150, 15}, {151, 15}, {152, 15}, {153, 15},
    {154, 15}, {155, 15}, {156, 15}, {157, 15}, {158, 15}, {159, 15},
    {160, 15}, {161, 15}, {162, 15}, {163, 15}, {164, 15}, {165, 15},
    {166, 15}, {167, 15}, {168, 15}, {169, 15}, {170, 15}, {171, 15},
    {172, 15}, {173, 15}, {174, 15}, {175, 15}, {176, 15}, {177, 15},
    {178, 15}, {179, 15}, {180, 15}, {181, 15}, {182, 15}, {183, 15},
    {184, 15}, {185, 15}, {186, 15}, {187, 15}, {188, 15}, {189, 15},
    {190, 15}, {191, 15}, {192, 15}, {193, 15}, {194, 15}, {195, 15},
    {196, 15}, {197, 15}, {198, 15}, {199, 15}, {200, 15}, {201, 15},
    {202, 15}, {203, 15}, {204, 15}, {205, 15}, {206, 15}, {207, 15},
    {208, 15}, {209, 15}, {210, 15}, {211, 15}, {212, 15}, {213, 15},
    {214, 15}, {215, 15}, {216, 15}, {217, 15}, {218, 15}, {219, 15},
    {220, 15}, {221, 15}, {222, 15}, {223, 15}, {224, 15}, {225, 15},
    {227, 15}, {228, 15}, {229, 15}, {230, 15}, {231, 15}, {232, 15},
    {233, 15}, {234, 15}, {235, 15}, {236, 15}, {237, 15}, {238, 15},
    {239, 15}, {240, 15}, {241, 15}, {242, 15}, {243, 15}, {244, 15},
    {245, 15}, {246, 15}, {247, 15}, {248, 15}, {249, 15}, {250, 15},
    {251, 15}, {252, 15}, {253, 15}, {254, 15}, {255, 15},
};

/* English text */
static huffman_code_t codes_english[] = {
    { 32,  3}, {'e',  3}, {'a',  4}, {'i',  4}, {'n',  4}, {'o',  4},
    {'s',  4}, {'t',  4}, {'d',  5}, {'h',  5}, {'l',  5}, {'r',  5},
    {'u',  5}, { 10,  6}, { 13,  6}, { 44,  6}, {'c',  6}, {'f',  6},
    {'g',  6}, {'m',  6}, {'p',  6}, {'w',  6}, {'y',  6}, { 46,  7},
    {'b',  7}, {'v',  7}, { 34,  8}, {'I',  8}, {'k',  8}, { 45,  9},
    {'A',  9}, {'N',  9}, {'T',  9}, { 39, 10}, { 59, 10}, { 63, 10},
    {'B', 10}, {'C', 10}, {'E', 10}, {'H', 10}, {'M', 10}, {'S', 10},
    {'W', 10}, {'x', 10}, { 33, 11}, {'0', 11}, {'1', 11}, {'F', 11},
    {'G', 11}, {  0, 15}, {  1, 15}, {  2, 15}, {  3, 15}, {  4, 15},
    {  5, 15}, {  6, 15}, {  7, 15}, {  8, 15}, {  9, 15}, { 11, 15},
    { 12, 15}, { 14, 15}, { 15, 15}, { 16, 15}, { 17, 15}, { 18, 15},
    { 19, 15}, { 20, 15}, { 21, 15}, { 22, 15}, { 23, 15}, { 24, 15},
    { 25, 15}, { 26, 15}, { 27, 15}, { 28, 15}, { 29, 15}, { 30, 15},
    { 31, 15}, { 35, 15}, { 36, 15}, { 37, 15}, { 38, 15}, { 40, 15},
    { 41, 15}, { 42, 15}, { 43, 15}, { 47, 15}, {'2', 15}, {'3', 15},
    {'4', 15}, {'5', 15}, {'6', 15}, {'7', 15}, {'8', 15}, {'9', 15},
    { 58, 15}, { 60, 15}, { 61, 15}, { 62, 15}, { 64, 15}, {'D', 15},
    {'J', 15}, {'K', 15}, {'L', 15}, {'O', 15}, {'P', 15}, {'Q', 15},
    {'R', 15}, {'U', 15}, {'V', 15}, {'X', 15}, {'Y', 15}, {'Z', 15},
    { 91, 15}, { 92, 15}, { 93, 15}, { 94, 15}, { 95, 15}, { 96, 15},
    {'j', 15}, {'q', 15}, {'z', 15}, {123, 15}, {124, 15}, {125, 15},
    {126, 15}, {127, 15}, {128, 15}, {129, 15}, {130, 15}, {131, 15},
    {132, 15}, {133, 15}, {134, 15}, {135, 15}, {136, 15}, {137, 15},
    {138, 15}, {139, 15}, {140, 15}, {141, 15}, {142, 15}, {143, 15},
    {144, 15}, {145, 15}, {146, 15}, {147, 15}, {148, 15}, {149, 15},
    {150, 15}, {151, 15}, {152, 15}, {153, 15}, {154, 15}, {155, 15},
    {156, 15}, {157, 15}, {158, 15}, {159, 15}, {160, 15}, {161, 15},
    {162, 15}, {163, 15}, {164, 15}, {165, 15}, {166, 15}, {167, 15},
    {168, 15}, {169, 15}, {170, 15}, {171, 15}, {172, 15}, {173, 15},
    {174, 15}, {175, 15}, {176, 15}, {177, 15}, {178, 15}, {179, 15},
    {180, 15}, {181, 15}, {182, 15}, {183, 15}, {184, 15}, {185, 15},
    {186, 15}, {187, 15}, {188, 15}, {189, 15}, {190, 15}, {191, 15},
    {192, 15}, {193, 15}, {194, 15}, {195, 15}, {196, 15}, {197, 15},
    {198, 15}, {199, 15}, {200, 15}, {201, 15}, {202, 15}, {203, 15},
    {204, 15}, {205, 15}, {206, 15}, {207, 15}, {208, 15}, {209, 15},
    {210, 15}, {211, 15}, {212, 15}, {213, 15}, {214, 15}, {215, 15},
    {216, 15}, {217, 15}, {218, 15}, {219, 15}, {220, 15}, {221, 15},
    {222, 15}, {223, 15}, {224, 15}, {225, 15}, {226, 15}, {227, 15},
    {228, 15}, {229, 15}, {230, 15}, {231, 15}, {232, 15}, {233, 15},
    {234, 15}, {235, 15}, {236, 15}, {237, 15}, {238, 15}, {239, 15},
    {240, 15}, {241, 15}, {242, 15}, {243, 15}, {244, 15}, {245, 15},
    {246, 15}, {247, 15}, {248, 15}, {249, 15}, {250, 15}, {251, 15},
    {252, 15}, {253, 15}, {254, 15}, {255, 15}, {SYM_EOF, 15},
};

/*
 * ---------------------------------------------------------------------------
 * Block_t structure support
 */

/*
 * Allocates and returns a new block_t struct of a specified default size.
 * Size maybe zero to defer allocation.
 *
 * Returns newly created block_t* on success
 *         NULL on failure
 */
block_t *block_create(size_t size) {
    block_t *b = (block_t *)malloc(sizeof(*b));
    if (!b)
	return NULL;

    b->data = NULL;
    b->alloc = size;
    b->byte = 0;
    b->bit = 0;

    if (size && NULL == (b->data = calloc(size, 1))) {
	free(b);
	return NULL;
    }

    return b;
}

/*
 * Deallocates memory created by block_create().
 * keep_data is a boolean which if true requests that the data held within
 * the block should not be deallocated as it is in use elsewhere.
 */
void block_destroy(block_t *blk, int keep_data) {
    if (!blk)
	return;

    if (!keep_data && blk->data)
	free(blk->data);
    
    free(blk);
}

/*
 * Ensures a block_t holds at least 'size' bytes.
 * Newly allocated data is initialised to zero.
 *
 * Returns 0 on success
 *        -1 on failure, leaving block pointing to the existing data
 */
int block_resize(block_t *blk, size_t size) {
    unsigned char *newp = NULL;

    if (!blk)
	return -1;

    /* Grow size to next power of 2, if we're growing */
    if (size > blk->alloc) {
	size--;
	size |= size >> 1;
	size |= size >> 2;
	size |= size >> 4;
	size |= size >> 8;
	size |= size >> 16;
	size++;
    }

    if (NULL == (newp = realloc(blk->data, size)))
	return -1;
    else
	blk->data = newp;

    if (size > blk->alloc)
	memset(&blk->data[blk->alloc], 0, size - blk->alloc);
    blk->alloc = size;

    return 0;
}


/*
 * ---------------------------------------------------------------------------
 * Tree building and code generation functions
 */

/*
 * Reverses the order of bits in the bottom nbits of val.
 * Returns the bit-reverse value.
 */
unsigned int bit_reverse(unsigned int val, int nbits) {
    unsigned int new = 0, i;

    for (i = 0; i < nbits; i++) {
	new = (new << 1) | (val & 1);
	val >>= 1;
    }

    return new;
}


/*
 * Generates canonical huffman codes given a set of symbol bit lengths.
 * The results are stored within the supplied huffman_codes_t struct.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int canonical_codes(huffman_codes_t *c) {
    int i, j;
    unsigned int code, last_len;
    int clens[33];
    int offs[33];
    huffman_code_t ctmp[258];
    signed int symtab[258];

    /* Sort by bit-length, subfield symbol - much faster than qsort() */
    for (i = 0; i < 258; i++)
	symtab[i] = -1;
    for (i = 0; i < c->ncodes; i++)
	symtab[c->codes[i].symbol] = i;
    for (i = 0; i <= 32; i++)
	offs[i] = clens[i] = 0;
    for (i = 0; i < c->ncodes; i++)
	clens[c->codes[i].nbits]++;
    for (i = 1; i <= 32; i++)
	offs[i] = offs[i-1] + clens[i-1];
    for (i = 0; i < 258; i++) {
	if (symtab[i] != -1)
	    ctmp[offs[c->codes[symtab[i]].nbits]++] = c->codes[symtab[i]];
    }
    memcpy(c->codes, ctmp, c->ncodes * sizeof(huffman_code_t));

    /*
     * Force all codes to be <= max_code_len. This is needed due to the
     * 15-bit length limitation of Deflate literal codes and the 7-bit 
     * limit of the code bit-length table.
     */
    /* Find first point of failure */
    for (i = 0; i < c->ncodes; i++) {
	if (c->codes[i].nbits > c->max_code_len)
	    break;
    }
    /*
     * From here on we shrink the length of the current code by increasing
     * the length of an earlier symbol, at last_code.
     */
    if (i != c->ncodes) {
	int delta = 0;

	/*
	fprintf(stderr, "=== REORDERING %d ===\n", c->code_set);
	output_code_set(stderr, c->codes, c->ncodes);
	output_code_set2(stderr, c->codes, c->ncodes);
	*/

	for (; i < c->ncodes; i++) {
	    int k, cur_len;

	    c->codes[i].nbits -= delta;
	    if (c->codes[i].nbits <= c->max_code_len)
		continue;

	    for (j = i; j >= 0 && c->codes[j].nbits >= c->max_code_len; j--)
		;
	    if (j < 0) {
		fprintf(stderr,
			"Too many symbols to fit in bit-length requirements\n");
		fprintf(stderr, "=== FAILING AT %d ===\n", c->code_set);
		output_code_set(stderr, c->codes, c->ncodes);
		output_code_set2(stderr, c->codes, c->ncodes);
		abort();
	    }

	    /*
	    fprintf(stderr, "Changing code %d/%d to len %d\n",
		    c->codes[i].symbol, c->codes[j].symbol,
		    c->codes[j].nbits+1);
	    */
	    cur_len = c->codes[i].nbits;
	    c->codes[i].nbits = ++c->codes[j].nbits;

	    /*
	     * Shrink the next code by one, or if none at that bit-length
	     * the next 2, and so on
	     */
	    delta = 1;
	    for (k = i+1; delta && k < c->ncodes; k++) {
		while (c->codes[k].nbits > cur_len) {
		    delta *= 2;
		    cur_len++;
		}
		c->codes[k].nbits--;
		delta--;
	    }
	    assert(delta == 0);
	}

	/*
	fprintf(stderr, "=== REORDERED TO %d ===\n", c->code_set);
	output_code_set(stderr, c->codes, c->ncodes);
	output_code_set2(stderr, c->codes, c->ncodes);
	*/

	/* Ordering is shot - regenerate via brute force way */
	return canonical_codes(c);
    }


    /* Generate codes */
    code = last_len = 0; /* stop warning */
    for (i = 0; i < c->ncodes; i++) {
	int nbits = c->codes[i].nbits;

	if (i == 0) {
	    code = 0;
	    last_len = nbits;
	} else {
	    code++;
	}
	if (nbits > last_len) {
	    code <<= (nbits - last_len);
	    last_len = nbits;
	}
	c->codes[i].code = bit_reverse(code, nbits);
    }

    /* Reindex so the symbol is the primary index into codes */
    for (i = 0; i <= 257; i++) {
	c->lookup[i].nbits = 0;
    }
    for (i = 0; i < c->ncodes; i++) {
        c->lookup[c->codes[i].symbol] = c->codes[i];
    }

    return 0;
}

static int node_compar(const void *vp1, const void *vp2) {
    const node_t *n1 = (const node_t *)vp1;
    const node_t *n2 = (const node_t *)vp2;

    /*
     * The sort order is vital here. This needs to return the same collating
     * order on all systems so that differing qsort() functions will not
     * swap around symbols with the same bit lengths, hence we sort by both
     * fields to force a unique stable ordering.
     */
    if (n1->count != n2->count)
	return n1->count - n2->count;
    else
	return n2->sym - n1->sym;
}

/*
 * Computes the huffman bit-lengths for a data set. We don't care
 * about the actual tree, just how deep the symbols end up.
 *
 * Huffman trees are constructed by constructing a set of nodes
 * initially containing the symbol and it's frequency. We then merge
 * the two least used nodes to produce a new node with a combined
 * frequency. Repeat until one root node is left.
 *
 * data/len is the input data to analyse.
 *
 * 'eof' is a boolean to indicate whether the EOF symbol should be included
 * in the symbols produced.
 *
 * all_codes is a boolean to indicate whether we should include symbols not
 * found in the input data set. (This was used to create the static lookup
 * tables.)
 *
 * Returns 0 on success
 *         -1 on failure
 */
int calc_bit_lengths(huffman_codes_t *c, unsigned char *data, int len,
		     int eof, int all_codes) {
    int i, ncodes, node_start;
    node_t nodes[258+257], *head, *new = &nodes[258];
    int map[258];

    /*
     * Initialise nodes. We build a map of ASCII character code to node
     * number. (By default it's a simple 1:1 mapping unless legal_chars is
     * defined.)
     */
    ncodes = 0;
    if (eof) {
	nodes[ncodes].sym = SYM_EOF;
	nodes[ncodes].count = 1;
	nodes[ncodes].parent = NULL;
	map[SYM_EOF] = ncodes++;
    }

    /* All 256 chars existing at a minimal level */
    for (i = 0; i < 256; i++) {
	nodes[ncodes].sym = i;
	nodes[ncodes].count = 0;
	nodes[ncodes].parent = NULL;
	map[i] = ncodes++;
    }

    /* Calc freqs */
    for (i = 0; i < len; i++) {
	nodes[map[data[i]]].count++;
    }

    /* Sort by counts, smallest first and form a sorted linked list */
    qsort(nodes, ncodes, sizeof(*nodes), node_compar);
    for (i = 0; i < ncodes; i++) {
	nodes[i].next = i+1 < ncodes ? &nodes[i+1] : NULL;
    }

    /* Skip symbols that do not occur, unless all_codes is true */
    node_start = 0;
    if (!all_codes) {
	while (node_start < ncodes && nodes[node_start].count == 0)
	    node_start++;
	ncodes -= node_start;
    }

    /* Repeatedly merge two smallest values */
    head = &nodes[node_start];
    while (head && head->next) {
	node_t *after = head->next, *n;
	int sum = head->count + head->next->count;
	
	for (n = head->next->next; n; after = n, n = n->next) {
	    if (sum < n->count)
		break;
	}

	/* Produce a new summation node and link it in place */
	after->next = new;
	new->next = n;
	new->sym = '?';
	new->count = sum;
	new->parent = NULL;
	head->parent = new;
	head->next->parent = new;
	head = head->next->next;

	new++;
    }

    /* Walk up tree computing the bit-lengths for our symbols */
    c->ncodes = ncodes;
    c->codes = (huffman_code_t *)malloc(c->ncodes * sizeof(*c->codes));
    if (NULL == c->codes) {
	return -1;
    }

    for (i = 0; i < ncodes; i++) {
	int len = 0;
	node_t *n;
	for (n = nodes[i + node_start].parent; n; n = n->parent) {
	    len++;
	}

	c->codes[i].symbol = nodes[i + node_start].sym;
	c->codes[i].freq   = nodes[i + node_start].count;
	c->codes[i].nbits  = len;
    }

    return 0;
}

/*
 * Initialises and returns a huffman_codes_t struct from a specified code_set.
 * If code_set is not one of the standard predefined values then the
 * input data is analysed using calc_bit_lengths() above to produce the
 * optimal set of huffman codes, otherwise we return predefined values.
 *
 * 'eof' is a boolean to indicate whether the EOF symbol should be included
 * in the symbols produced.
 *
 * all_codes is a boolean to indicate whether we should include symbols not
 * found in the input data set. (This was used to create the static lookup
 * tables.)
 *
 * Returns huffman_codes_t* on success; free using huffman_codes_destroy().
 *         NULL on failure.
 */
huffman_codes_t *generate_code_set(int code_set, unsigned char *data, int len,
				   int eof, int max_code_len, int all_codes) {
    huffman_codes_t *c;

    if (NULL == (c = (huffman_codes_t *)malloc(sizeof(*c))))
	return NULL;

    c->code_set = code_set;
    c->max_code_len = max_code_len;
    
    if (code_set >= 128 || code_set == CODE_INLINE) { 
	if (-1 == calc_bit_lengths(c, data, len, eof, all_codes)) {
	    free(c);
	    return NULL;
	}
    } else {
	switch(code_set) {
	case CODE_DNA:
	    c->codes = codes_dna;
	    c->ncodes = sizeof(codes_dna)/sizeof(*c->codes);
	    break;

	case CODE_DNA_AMBIG:
	    c->codes = codes_dna_ambig;
	    c->ncodes = sizeof(codes_dna_ambig)/sizeof(*c->codes);
	    break;

	case CODE_ENGLISH:
	    c->codes = codes_english;
	    c->ncodes = sizeof(codes_english)/sizeof(*c->codes);
	    break;

	default:
	    fprintf(stderr, "Unknown huffman code set '%d'\n", code_set);
	    return NULL;
	}
    }

    canonical_codes(c);

    return c;
}

huffman_codes_t *get_code_set(int code_set) {
    return code_set < CODE_USER
	? generate_code_set(code_set, NULL, 0, 1, MAX_CODE_LEN, 0)
	: NULL;
}

void huffman_codes_destroy(huffman_codes_t *c) {
    if (!c)
	return;

    switch(c->code_set) {
    case CODE_DNA:
    case CODE_DNA_AMBIG:
    case CODE_ENGLISH:
	/* All these use static tables */
	break;

    default:
	if (c->codes)
	    free(c->codes);
    }

    free(c);
}

/*
 * ---------------------------------------------------------------------------
 * Encoding and decoding related functions
 */

/*
 * Can store up to 24-bits worth of data encoded in an integer value
 * Possibly we'd want to have a less optimal store_bits function when dealing
 * with nbits > 24, but for now we assume the codes generated are never
 * that big. (Given this is only possible with 121392 or more
 * characters with exactly the correct frequency distribution we check
 * for it elsewhere.)
 */
static void store_bits(block_t *block, unsigned int val, int nbits) {
    /* fprintf(stderr, " store_bits: %02x %d\n", val, nbits); */

#if 1
    {
	unsigned int curr = block->data[block->byte];
	curr |= (val & ((1 << nbits)-1)) << block->bit;
	block->bit += nbits;
	while (block->bit >= 8) {
	    block->data[block->byte++] = curr & 0xff;
	    curr >>= 8;
	    block->bit -= 8;
	}
	block->data[block->byte] = curr & 0xff;
    }
    return;

#else

    {
      /* Slow, but simple */
      unsigned int mask = 1;
      int bit = 1 << block->bit;
      do {
	if (val & mask)
	    block->data[block->byte] |= bit;
	/*
	 * Data should be zeroed anyway, so this is not needed.
	 *
	 * else
	 *    block->data[block->byte] &= ~bit;
	*/

	if (++block->bit == 8) {
	    block->bit = 0;
	    block->byte++;
	    block->data[block->byte] = 0;
	    bit = 1;
	} else {
	    bit <<= 1;
	}
	mask <<= 1;
      } while(--nbits);
    }

#endif
}

/* stores nbytes bytes, padding to align on the next byte boundary */
void store_bytes(block_t *block, unsigned char *val, int nbytes) {
    /* Align */
    if (block->bit != 0) {
	block->byte++;
	block->bit = 0;
    }

    /* Resize */
    block_resize(block, block->byte + nbytes);

    /* Store */
    memcpy(&block->data[block->byte], val, nbytes);
    block->byte += nbytes;
}


/*
 * Encodes the huffman symbol bit-lengths as a serialised block of data
 * suitable for storing in a ZTR "ZLBH" chunk. This uses the Deflate
 * storage format defined in RFC1951.
 *
 * Returns: 0 on success
 *          -1 on failure
 */
#ifndef MAX
#    define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
int store_codes(block_t *out, huffman_codes_t *c, int last_block) {
    int i;
    unsigned char bl_code[257]; /* bit-length codes and for codes 16-18 */
    unsigned char bl_opt[257];  /*     the operand to the blcode */
    unsigned char sorted_codes[258];
    int bl_freq[19]; /* frequency of bit-length codes produced */
    int bl_count;
    huffman_codes_t *bl_cds = NULL;
    int hclen_order[] = {
	16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    };
    int hlit, hdist, hclen, hcmap[19];

    if (out->alloc < out->byte + 1000) {
	out->alloc = out->byte + 1000;
	if (NULL == (out->data = realloc(out->data, out->alloc)))
	    return -1;
    }

    /* Header details */
    store_bits(out, last_block != 0,  1); /* last block */
    store_bits(out, 2,  2); /* dynamic huffman */

    /*
     *-----------------------------------------------------------------
     * Reformat the dynamic code bit-lengths into an alphabet of 19 
     * "code length" symbols as defined in RFC1951.
     */
    memset(sorted_codes, 0, 258);
    for (i = 0; i < c->ncodes; i++) {
	sorted_codes[c->codes[i].symbol] = c->codes[i].nbits;
    }
    for (i = 0; i < 19; i++)
	bl_freq[i] = 0;

    bl_count = 0;
    for (i = 0; i < 257; ) {
	int j = i+1, n;
	int v = sorted_codes[i];
	while (j < 257 && sorted_codes[j] == v)
	    j++;

	n = j-i; /* n = run-length */
	/* fprintf(stderr, "value=%d, run_len=%d\n", v, n); */
	if (v == 0) {
	    /* bit-len zero => no code and uses code 17/18 for run len */
	    while (n > 0) {
		while (n >= 11) {
		    bl_freq[18]++;
		    bl_code[bl_count] = 18;
		    bl_opt[bl_count] = MIN(n, 138)-11;
		    n -= bl_opt[bl_count]+11;
		    bl_count++;
		}

		while (n >= 3) {
		    bl_freq[17]++;
		    bl_code[bl_count] = 17;
		    bl_opt[bl_count] = MIN(n, 10)-3;
		    n -= bl_opt[bl_count]+3;
		    bl_count++;
		}

		switch (n) {
		case 2: bl_code[bl_count++] = 0; bl_freq[0]++; n--;
		case 1: bl_code[bl_count++] = 0; bl_freq[0]++; n--;
		}
	    }
	} else if (v <= 15) {
	    /* non-zero code, uses code 16 for run-len */
	    if (n >= 4) {
		bl_freq[v]++;
		bl_code[bl_count++] = v;
		n--;
		while (n >= 3) {
		    bl_freq[16]++;
		    bl_code[bl_count] = 16;
		    bl_opt[bl_count] = MIN(n, 6)-3;
		    n -= bl_opt[bl_count]+3;
		    bl_count++;
		}
	    }

	    switch(n) {
	    case 3: bl_code[bl_count++] = v; bl_freq[v]++; n--;
	    case 2: bl_code[bl_count++] = v; bl_freq[v]++; n--;
	    case 1: bl_code[bl_count++] = v; bl_freq[v]++; n--;
	    }
	} else {
	    fprintf(stderr, "Unsupported code length: %d\n", v);
	}

	i = j;
    }
    hlit = 257;

    /* Add a single distance code of zero bits. This means that there
     * are no distance codes used at all.
     */
    bl_code[bl_count++] = 0;
    bl_freq[0]++;
    hdist = 1;

    /* Produce new huffman codes for our code-length symbols. */
    bl_cds = generate_code_set(CODE_INLINE, bl_code, bl_count, 0, 7, 0);

    /*
     *-----------------------------------------------------------------
     * Output the "code length" bit-lengths, 3 bits at a time.
     *
     * Compute how many HCLEN code length values we need, using the
     * predefined order in the RFC.
     */
    for (hclen = 19; hclen > 0; hclen--) {
	if (bl_freq[hclen_order[hclen-1]])
	    break;
    }

    store_bits(out, hlit-257,  5);
    store_bits(out, hdist-1, 5);
    store_bits(out, hclen-4, 4);

    for (i = 0; i < 19; i++)
	hcmap[i] = -1;
    for (i = 0; i < bl_cds->ncodes; i++)
	hcmap[bl_cds->codes[i].symbol] = i;
    for (i = 0; i < hclen; i++) {
	if (hcmap[hclen_order[i]] >= 0) {
	    store_bits(out, bl_cds->codes[hcmap[hclen_order[i]]].nbits, 3);
	} else {
	    store_bits(out, 0, 3);
	}
    }


    /*
     *----------------------------------------------------------------
     * Finally output the original bit-lengths using the code-len codes.
     */
    for (i = 0; i < bl_count; i++) {
	huffman_code_t *c = &bl_cds->codes[hcmap[bl_code[i]]];
	store_bits(out, c->code, c->nbits); 
	/*fprintf(stderr, "bl_code %d (opt %d), code %d/%d\n",
	  bl_code[i], bl_opt[i], c->code, c->nbits); */
	switch(bl_code[i]) {
	case 18:
	    store_bits(out, bl_opt[i], 7);
	    break;
	case 17:
	    store_bits(out, bl_opt[i], 3);
	    break;
	case 16:
	    store_bits(out, bl_opt[i], 2);
	    break;
	}
    }

    huffman_codes_destroy(bl_cds);

    return 0;
}

/*
 * This is the opposite of the store_codes() function. It loads generates
 * huffman_codes_t structs from the a serialised data stream as presented
 * in the above format.
 *
 * The input data is the data-string. On return the number of bytes
 * consumed will be returned in *len_used (if non NULL).
 * This is to allow stripping off of the huffman codes from a longer
 * array of data (ie probably followed by the STHUFF encoded chunk
 * itself).
 *
 * Returns: malloced huffman_codes_t structure on success.
 *          NULL on failure.
 */

/* FIXME: no longer valid as we need to handle deflate code ordering */
huffman_codes_t *restore_codes(unsigned char *data, int *len_used) {
    fprintf(stderr, "restore_codes() not yet implemented for Deflate\n");
    return NULL;
}

/*
 * Given a set of huffman codes and a block of data this compresses
 * and returns the data block.
 *
 * Returns: 0 on success
 *          -1 on failure
 */
int huffman_encode(block_t *blk, huffman_codes_t *c, int code_set,
		   unsigned char *data, int len) {
    int i;
    huffman_code_t *lookup;

#ifdef TIME_IT
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
#endif    

    /*
     * The maximum size to encode len bytes is < 9 bits per symbol
     * (not quite 8 due to an EOF symbol) plus the overhead of the bit-length
     * tree. That in turn, with alternating 8/9 bit-lengths would max out
     * as 258*8 + 5+5+4 + 19*3 + 258*5 bits (429 bytes), but in practice
     * I'm not even sure if it's possible to construct such a set of code
     * lengths that would compress that poor.
     *
     * This of course assumes we're using appropriate compression codes.
     * Given a user may give a completely inappropriate code we have to
     * assume every symbol is actually 15 bits instead of < 9 on average.
     *
     * We ensure blk here is large enough for the worst case scenario so we
     * don't incur overheads in store_bits().
     */
    if (blk->alloc <= 429 + 2*len + 2 + blk->byte) {
	blk->alloc  = 429 + 2*len + 2 + blk->byte;
	blk->data = realloc(blk->data, blk->alloc);
	if (!blk->data)
	    return -1;
    }

    if (!c) {
	/* No code set, so derive our own */
	c = generate_code_set(code_set, data, len, 1, MAX_CODE_LEN, 0);
    }

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "code production took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
    gettimeofday(&t1, NULL);
#endif    

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "alloc + storing codes took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
    gettimeofday(&t1, NULL);
#endif    

    lookup = c->lookup;
    for (i = 0; i < len; i++) {
	store_bits(blk, lookup[data[i]].code,
		   lookup[data[i]].nbits);
    }
    store_bits(blk, c->lookup[SYM_EOF].code, c->lookup[SYM_EOF].nbits);

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "encoding took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
    gettimeofday(&t1, NULL);
#endif    

    assert(blk->alloc > blk->byte);
    /* We probably massively overallocated, so return some of it back */
    blk->data = realloc(blk->data, blk->byte+1);
    blk->alloc = blk->byte+1;

    return 0;
}

typedef struct {
    /* Graph construction */
    unsigned short c[2]; /* child node */
      signed short l[2]; /* symbol to emit on transition. -1 => none */
} htree_t;

typedef struct {
    /* Byte-wise jumping table */
    unsigned short jump;
    unsigned char symbol[4];
    unsigned char nsymbols;
    unsigned char top_bit;   /* bit 9 of symbol[] */
} h_jump4_t;

/*
 * The opposite of huffman_encode().
 *
 * Returns: malloced uncompressed data on success, of length *nbytes.
 *          NULL on failure.
 *
 * Method 1
 * --------
 *
 * At any node in our tree we can precompute a lookup table so that upon
 * reading the next 'k' bits we know the new node we'd end up in and what
 * symbols to export.
 * Then decoding simply works in fixed sets of k bits at a time.
 *
 * We use k=4 for efficient table space (they fit neatly in cache) and ease
 * of decoding 4-bits at a time. k=8 is about 20% faster as reading the input
 * byte by byte is easy, but the setup time is substantially longer
 * (16x at a guess) and the lookup tables no longer fit in the L1 cache.
 */
unsigned char *huffman_decode(huffman_codes_t *c,
			      unsigned char *data, int len, int *nbytes) {
    block_t *in, *out;
    htree_t t[513]; /* number of internal nodes */
    int i, j, n;
    int private_codes = 0;
    int new_node;
    h_jump4_t J4[513][16];

#ifdef TIME_IT
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
#endif    

    if (NULL == (in = block_create(len)))
	return NULL;
    in->byte = 1;

    if (NULL == (out = block_create(8*len))) {
	block_destroy(in, 0);
	return NULL;
    }
    if (in->data[in->byte] == 0) {
	int used;
	unsigned char *data = &in->data[in->byte];

	/* inline code-set, so read */
	c = restore_codes(data, &used);
	private_codes = 1;
	in->byte += used;
    } else if (in->data[in->byte] < CODE_USER) {
	/* static predefined codes */
	c = get_code_set(in->data[in->byte]);
	private_codes = 1;
	in->byte++;
    } else {
	/* Use passed in supplied code_set */
	in->byte++;
    }

    /*
    output_code_set (stderr, c->codes, c->ncodes);
    output_code_set2(stderr, c->codes, c->ncodes);
    */

    /* Construct the tree from the codes */
    new_node = 1;
    t[0].l[0] = t[0].l[1] = -1;
    t[0].c[0] = t[0].c[1] = 0;
    for (i = 0; i < c->ncodes; i++) {
	int n = 0;
	unsigned int v = c->codes[i].code;

	for (j = 0; j < c->codes[i].nbits-1; j++) {
	    int b = v & 1;
	    if (t[n].c[b]) {
		n = t[n].c[b];
	    } else {
		n = (t[n].c[b] = new_node++);
		t[n].c[0] = t[n].c[1] = 0;
		t[n].l[0] = t[n].l[1] = -1;
	    }
	    v >>= 1;
	}
	/* last bit */
	t[n].l[v & 1] = c->codes[i].symbol;
    }

    /* Build the 16 wide lookup table per node */
    for (n = 0; n < new_node; n++) {
	for (j = 0; j < 16; j++) {
	    unsigned int v = j;
	    int n2 = n;
	    h_jump4_t *hj = &J4[n][j];
	    hj->nsymbols = 0;
	    hj->top_bit = 0;
	    for (i = 0; i < 4; i++) {
		int b = v & 1;
		if (t[n2].l[b] >= 0) {
		    hj->symbol[hj->nsymbols++] = t[n2].l[b];
		    if (t[n2].l[b] == SYM_EOF)
			if (!hj->top_bit)
			    hj->top_bit |= 1 << (hj->nsymbols-1);
		}
		n2 = t[n2].c[b];
		v >>= 1;
	    }
	    hj->jump = n2;
	}
    }

#if 0
    /* Debug output */
    for (n = 0; n < new_node; n++) {
	printf("Node %d, c[]={%d,%d}, l[]={%d,%d}\n",
	       n, t[n].c[0], t[n].c[1], t[n].l[0], t[n].l[1]);
	for (i = 0; i < 256; i++) {
	    printf("\t%02x %s =>%02d, ", i, print_8rev(i), J4[n][i].jump);
	    for (k = 0; k < J4[n][i].nsymbols; k++) {
		if (isprint(J4[n][i].symbol[k]))
		    printf(" '%c'", J4[n][i].symbol[k]);
		else
		    printf(" %03d", J4[n][i].symbol[k]);
	    }
	    printf("\n");
	}
    }
#endif


#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "setup2 took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
    gettimeofday(&t1, NULL);
#endif    

    /*
     * Decoding.
     * We handle data nibble by nibble using the nibble to get an
     * h_jump4_t lookup from the J4[] table.
     * If top_bit is clear then we know we have no funny business (SYM_EOF)
     * so we use a fast decoding technique, otherwise we have to do a slower
     * loop with a check.
     */
     {
	 int node_num = 0, last_node = 0;
	 unsigned char *cp = out->data;
	 unsigned char *last_cp = cp;
	 h_jump4_t *x = &J4[node_num][data[i] & 0x0f];
	 int l = x->nsymbols;

	 /*
	  * This is the tight loop, so we over-optimise here by ignoring EOF
	  * and relying on knowing 'len' - the input data stream length.
	  * This allows us to ignore the 9-bit data and only operate on
	  * the basic 0-255 symbols, glossing over the minor issue that EOF
	  * will look like an ordinary symbol.
	  */
	 for (i = in->byte; i < len; i++) {
	     last_cp = cp;
	     last_node = node_num;

	     x = &J4[node_num][data[i] & 0x0f];
	     l = x->nsymbols;

	     for (j = 0; j < l; j++) {	
		 *cp++ = x->symbol[j];
	     }
	     node_num = x->jump;

	     x = &J4[node_num][(data[i] >> 4) & 0x0f];
	     l = x->nsymbols;
	     
	     for (j = 0; j < l; j++) {	
		 *cp++ = x->symbol[j];
	     }
	     node_num = x->jump;
	 }

	 /*
	  * The above optimisation has unfortunately added EOF to our data
	  * along with whatever else was packed in the last byte after the
	  * EOF symbol. So we rewind one byte and finish off decoding
	  * taking appropriate care and attention to get it right this time.
	  */
	 cp = last_cp;
	 node_num = last_node;
	 i--;
	 x = &J4[node_num][data[i] & 0x0f];
	 l = x->nsymbols;

	 for (j = 0; j < l; j++) {	
	     unsigned short sym = x->symbol[j] |
		 (((x->top_bit & (1<<j)) != 0) << 8);
	     if (sym == SYM_EOF)
		 goto term;
	     else
		 *cp++ = sym;
	 }
	 
	 node_num = x->jump;
	 x = &J4[node_num][(data[i] >> 4) & 0x0f];
	 l = x->nsymbols;
	     
	 for (j = 0; j < l; j++) {	
	     unsigned short sym = x->symbol[j] |
		 (((x->top_bit & (1<<j)) != 0) << 8);
	     if (sym == SYM_EOF)
		 goto term;
	     else
		 *cp++ = sym;
	 }

     term:
	 out->byte = cp - out->data;
	 out->data = (unsigned char *)realloc(out->data, out->byte+1);
     }

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "decoding2 took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
#endif    

    if (private_codes)
	huffman_codes_destroy(c);

    *nbytes = out->byte;
    return out->data; /* memory leak - 'out' */
}


/*
 * ---------------------------------------------------------------------------
 * Debug code. This turns the library into a stand-alone program for
 * easy debugging.x
 */
static void output_code_set(FILE *fp, huffman_code_t *codes, int ncodes) {
    int i, j;
    int nbits = 0;

    fprintf(fp, "static huffman_code_t codes_FIXME[] = {\n");
    for (i = j = 0; i < ncodes; i++) {
	nbits += codes[i].nbits * codes[i].freq;
	if (j == 0)
	    fprintf(fp, "    ");
	if (codes[i].symbol == SYM_EOF) {
	    fprintf(fp, "{SYM_EOF,%3d}, ", codes[i].nbits);
	    j = 10;
	} else {
	    if (isalnum(codes[i].symbol)) {
		fprintf(fp, "{'%c',%3d}, ", codes[i].symbol, codes[i].nbits);
	    } else {
		fprintf(fp, "{%3d,%3d}, ", codes[i].symbol, codes[i].nbits);
	    }
	}
	j++;
	
	if (j >= 6) {
	    fputc('\n', fp);
	    j = 0;
	}
    }
    if (j)
	fputc('\n', fp);
    fprintf(fp, "};\n");
    fprintf(fp, "/* Tested on %d bits of compressed data */\n", nbits);
}

static void output_code_set2(FILE *fp, huffman_code_t *codes, int ncodes) {
    int i;

    fprintf(fp, "huffman_code_t = {\n");
    for (i = 0; i < ncodes; i++) {
	fprintf(fp, "\t%d:\t%3d %c %2d %04x %d\n",
		i,codes[i].symbol, codes[i].symbol,
		codes[i].nbits, codes[i].code,
		codes[i].freq);
    }
    fprintf(fp, "};\n");
}


/*
 * --------------------------------------------------------------------------
 * A test main() to create an application capable of compressing and
 * uncompressing stdin.
 */

#ifdef TEST_MAIN
#include <zlib.h>

/*
 * Slurps the entirety of stdin into a malloced buffer and returns a pointer
 * to it.
 *
 * Returns: malloced buffer on success, *lenp equal to length
 *          NULL on failure
 */
static unsigned char *load(int *lenp) {
    unsigned char *data = NULL;
    int dsize = 0;
    int dcurr = 0, len;

    do {
	if (dsize - dcurr < 8192) {
	    dsize = dsize ? dsize * 2 : 8192;
	    if (NULL == (data = realloc(data, dsize)))
		return NULL;
	}

	len = read(0, data + dcurr, 8192);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
	return NULL;
    }

    *lenp = dcurr;
    return data;
}

/*
 * Returns 0 for success
 *        -1 for failure.
 */
int decode_main(unsigned char *data, int len, int code_set) {
    huffman_codes_t *cds;
    z_stream zstr;
    block_t *blk_in;
    unsigned char *out;
    int err, out_len;

    /*
     * FIXME: decoding not working in this code now - uses zlib instead
     *
     * out = huffman_decode(NULL, data, len, &out_len);
     * write(1, out, out_len);
     * free(out);
     */
    blk_in = block_create(1000 + len);

    if (code_set != 0) {
	cds = generate_code_set(code_set, NULL, 0, 1, MAX_CODE_LEN, 0);
	store_codes(blk_in, cds, 1);
    }

    if (blk_in->bit != 0) {
	blk_in->data[blk_in->byte] |= data[0];
	memcpy(&blk_in->data[blk_in->byte+1], data+1, len-1);
    } else {
	memcpy(&blk_in->data[blk_in->byte], data, len);
    }

    /* simple hack as a proof of concept only */
    zstr.zalloc = (alloc_func)0;
    zstr.zfree = (free_func)0;
    zstr.opaque = (voidpf)0;
    zstr.next_in = blk_in->data;
    zstr.avail_in = blk_in->byte + len;
    out = malloc(out_len = 65536);

    if ((err = inflateInit2(&zstr, -15)) != Z_OK) {
	fprintf(stderr, "zlib errror in inflateInit2(): %d/%s\n",
		err, zstr.msg);
	return -1;
    }

    do {
	zstr.next_out = out;
	zstr.avail_out = out_len;

	err = inflate(&zstr, Z_FINISH);
	if (err == Z_STREAM_END || err == Z_OK || err == Z_BUF_ERROR) {
	    write(1, out, zstr.total_out);
	    zstr.total_out = 0;
	} else {
	    fprintf(stderr, "zlib errror in inflate(): %d/%s\n",
		    err, zstr.msg);
	    return 1;
	}
    } while (err != Z_FINISH && err != Z_STREAM_END);

    write(1, out, zstr.total_out);
    inflateEnd(&zstr);
    block_destroy(blk_in, 0);
    free(out);

    if (code_set != 0)
	huffman_codes_destroy(cds);

    return 0;
}

/*
 * Returns 0 for success
 *        -1 for failure.
 */
int encode_main(unsigned char *data, int len, int code_set,
		int blk_size, int dump_tree, int exit_after_tree) {
    /* Encoding */
    unsigned char *d2 = data;
    block_t *blk;
    huffman_codes_t *cds;

    blk = block_create(8192);

#if 0
    /* Zlib header - see http://www.ietf.org/rfc/rfc1950.txt */
    blk->data[blk->byte++] = 0x78;
    blk->data[blk->byte++] = 0x01;
    /* Would need adler32(data,len) footer too */
#endif

    fprintf(stderr, "Input %d bytes\n", len);

    do {
	int l2 = len > blk_size ? blk_size : len;
	if (code_set != 0)
	    l2 = len; /* predefined code-sets have final-block bit set */
	cds = generate_code_set(code_set, d2, l2, 1, MAX_CODE_LEN, 0);
	if (!cds)
	    return -1;

	if (dump_tree) {
	    output_code_set(stdout, cds->codes, cds->ncodes);
	    output_code_set2(stdout, cds->codes, cds->ncodes);
	    if (exit_after_tree)
		return 0;
	}

	store_codes(blk, cds, l2 == len);
	if (code_set != 0)
	    blk->data[blk->byte = 0] = 0;  /* starting bit no. preseved */
	if (exit_after_tree) {
	    write(1, blk->data, blk->byte + (blk->bit != 0));
	    return 0;
	}
	huffman_encode(blk, cds, code_set, d2, l2);
	huffman_codes_destroy(cds);
	len -= l2;
	d2  += l2;

    } while (len > 0);
    write(1, blk->data, blk->byte + (blk->bit != 0));

    fprintf(stderr, "output %ld bytes\n", blk->byte + (blk->bit != 0));
    block_destroy(blk, 0);

    return 0;
}

int main(int argc, char **argv) {
    unsigned char *data;
    int len;
    int decode = 0;
    int dump_tree = 0;
    int exit_after_tree = 0;
    int code_set = CODE_INLINE;
    int c;
    int blk_size = 0x7fff;
    int r;

    while ((c = getopt(argc, argv, "c:detxl:b:h")) != -1) {
	switch (c) {
	case 'b':
	    blk_size = atoi(optarg);
	    break;

	case 'c':
	    code_set = atoi(optarg);
	    break;

	case 'd':
	    decode = 1;
	    break;

	case 'e':
	    decode = 0;
	    break;

	case 't':
	    dump_tree = 1;
	    break;

	case 'x':
	    exit_after_tree = 1;
	    break;

	default:
	    fprintf(stderr, "Usage: huffman_static [options] < stdin > stdout\n");
	    fprintf(stderr, "    Decoding options\n");
	    fprintf(stderr, "        -d\tdecode flag\n");
	    fprintf(stderr, "    Encoding options\n");
	    fprintf(stderr, "        -e\tencode flag\n");
	    fprintf(stderr, "        -b size\tspecify the block-size\n");
	    fprintf(stderr, "        -c code\tspecify code-set. 0 => inline\n");
	    fprintf(stderr, "        -t\tpretty-print the code-set used\n");
	    fprintf(stderr, "        -x\texit after outputting code-set\n");
	    exit(c != 'h');
	}
    }

    data = load(&len);

    r = decode
	? decode_main(data, len, code_set)
	: encode_main(data, len, code_set,
		      blk_size, dump_tree, exit_after_tree);

    free(data);

    return r == 0 ? 0 : 1;
}
#endif
