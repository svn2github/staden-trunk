#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <unistd.h>

#include "huffman_static.h"

/*
 * ---------------------------------------------------------------------------
 * Local structs & defines
 */

/* Use for store_bits() and get_bits() */
typedef struct block {
    unsigned char *data;
    size_t alloc;
    size_t byte;
    int bit;
} block_t;

/* Used in tree construction only */
typedef struct node {
    int count;
    int sym; /* char, SYM_EOF or SYM_ANY */
    struct node *parent;
    struct node *next;
} node_t;

#define SYM_EOF 256
#define SYM_ANY 257

/*
 * ---------------------------------------------------------------------------
 * Our standard precomputed tables, for DNA, Solexa traces, confidence values,
 * etc.
 */

/* For DNA; 4 choices depending on rarest base type */
static huffman_code_t codes_dna_a[] = {
    {SYM_ANY, 5,}, {SYM_EOF, 5,},
    {'A',  3,}, {'C',  2,}, {'G',  2,}, {'T',  2,}, {'-',  4,}
};

static huffman_code_t codes_dna_c[] = {
    {SYM_ANY, 5,}, {SYM_EOF, 5,},
    {'A',  2,}, {'C',  3,}, {'G',  2,}, {'T',  2,}, {'-',  4,}
};

static huffman_code_t codes_dna_g[] = {
    {SYM_ANY, 5,}, {SYM_EOF, 5,},
    {'A',  2,}, {'C',  2,}, {'G',  3,}, {'T',  2,}, {'-',  4,}
};

static huffman_code_t codes_dna_t[] = {
    {SYM_ANY, 5,}, {SYM_EOF, 5,},
    {'A',  2,}, {'C',  2,}, {'G',  2,}, {'T',  3,}, {'-',  4,}
};

/* For Solexa Traces */
static huffman_code_t codes_solexa_trace[] = {
    {  0,  2}, {  1,  3}, {  2,  4}, {  3,  5}, {  6,  6}, {  5,  6},
    {  4,  6}, { 11,  7}, { 10,  7}, {  9,  7}, {  8,  7}, {  7,  7},
    { 20,  8}, { 21,  8}, { 19,  8}, { 18,  8}, { 17,  8}, { 16,  8},
    { 15,  8}, { 14,  8}, { 13,  8}, { 12,  8}, {224,  9}, {211,  9},
    {188,  9}, {182,  9}, {184,  9}, {200,  9}, {217,  9}, {199,  9},
    {160,  9}, {173,  9}, {204,  9}, {198,  9}, {206,  9}, {176,  9},
    {186,  9}, {205,  9}, {192,  9}, {185,  9}, {169,  9}, {196,  9},
    {197,  9}, {202,  9}, {177,  9}, {194,  9}, {179,  9}, {183,  9},
    {187,  9}, {168,  9}, {181,  9}, {174,  9}, {150,  9}, {139,  9},
    {172,  9}, {178,  9}, {170,  9}, {154,  9}, {144,  9}, {167,  9},
    {189,  9}, {149,  9}, {155,  9}, {171,  9}, {195,  9}, {175,  9},
    {145,  9}, {159,  9}, {157,  9}, {165,  9}, {152,  9}, {164,  9},
    {162,  9}, {156,  9}, {163,  9}, {134,  9}, {143,  9}, {148,  9},
    {138,  9}, {161,  9}, {'s',  9}, {137,  9}, {125,  9}, {146,  9},
    {166,  9}, {151,  9}, {158,  9}, {'x',  9}, {124,  9}, {'i',  9},
    {123,  9}, {126,  9}, {'y',  9}, {136,  9}, {129,  9}, {127,  9},
    {132,  9}, {141,  9}, {153,  9}, {'z',  9}, {133,  9}, {140,  9},
    {130,  9}, {'v',  9}, {135,  9}, {142,  9}, {'k',  9}, {'t',  9},
    {131,  9}, { 93,  9}, {'n',  9}, {'w',  9}, {'e',  9}, {'m',  9},
    {147,  9}, {'u',  9}, {128,  9}, {'f',  9}, {'V',  9}, {'r',  9},
    {'g',  9}, {'o',  9}, {'c',  9}, {'j',  9}, {'U',  9}, {'P',  9},
    {'q',  9}, { 95,  9}, {'l',  9}, {'b',  9}, {'d',  9}, {'R',  9},
    {'N',  9}, { 96,  9}, {'O',  9}, {'C',  9}, {'J',  9}, {'p',  9},
    {'h',  9}, {'a',  9}, {'Y',  9}, {'W',  9}, { 91,  9}, {'Z',  9},
    { 94,  9}, {'X',  9}, {'S',  9}, { 92,  9}, {'K',  9}, {'M',  9},
    {'F',  9}, {'L',  9}, {'T',  9}, { 62,  9}, {'G',  9}, {'D',  9},
    {'H',  9}, {'I',  9}, {'B',  9}, {'E',  9}, { 60,  9}, { 47,  9},
    { 61,  9}, { 64,  9}, {'Q',  9}, { 63,  9}, {'6',  9}, {'A',  9},
    { 59,  9}, {'2',  9}, { 45,  9}, { 40,  9}, {'0',  9}, {'5',  9},
    {'1',  9}, {'9',  9}, { 58,  9}, {'7',  9}, { 43,  9}, {'4',  9},
    { 33,  9}, { 37,  9}, { 38,  9}, { 42,  9}, { 44,  9}, {'8',  9},
    {'3',  9}, { 35,  9}, { 46,  9}, { 41,  9}, { 36,  9}, { 32,  9},
    { 39,  9}, { 34,  9}, { 31,  9}, { 29,  9}, { 30,  9}, { 28,  9},
    { 26,  9}, { 27,  9}, { 25,  9}, { 24,  9}, { 23,  9}, { 22,  9},
    {240, 10}, {245, 10}, {213, 10}, {248, 10}, {255, 10}, {247, 10},
    {252, 10}, {254, 10}, {230, 10}, {222, 10}, {232, 10}, {244, 10},
    {235, 10}, {215, 10}, {251, 10}, {246, 10}, {218, 10}, {242, 10},
    {253, 10}, {214, 10}, {226, 10}, {238, 10}, {243, 10}, {237, 10},
    {249, 10}, {241, 10}, {225, 10}, {234, 10}, {228, 10}, {208, 10},
    {239, 10}, {233, 10}, {201, 10}, {220, 10}, {193, 10}, {231, 10},
    {236, 10}, {219, 10}, {212, 10}, {229, 10}, {221, 10}, {210, 10},
    {190, 10}, {223, 10}, {180, 10}, {216, 10}, {203, 10}, {209, 10},
    {227, 10}, {191, 10}, {207, 10}, {SYM_EOF, 11},
    {250, 11},
};

/* For 4x uncompressed confidence values */
static huffman_code_t codes_conf[] = {
    {216,  1}, { 40,  3}, {251,  4}, {  0,  6}, {239,  7}, {241,  7},
    {240,  7}, {242,  7}, {243,  7}, {246,  7}, {244,  7}, {245,  7},
    {247,  7}, {248,  7}, {250,  7}, {249,  7}, {253,  7}, {252,  7},
    {254,  7}, {255,  7}, {218,  8}, {220,  8}, {221,  8}, { 12,  8},
    {219,  8}, {222,  8}, {  9,  8}, {223,  8}, {  8,  8}, { 11,  8},
    { 10,  8}, {224,  8}, {229,  8}, {225,  8}, {227,  8}, {  7,  8},
    {228,  8}, {226,  8}, {  6,  8}, {231,  8}, {  5,  8}, {230,  8},
    {232,  8}, {  3,  8}, {  4,  8}, {234,  8}, {  2,  8}, {235,  8},
    {  1,  8}, {233,  8}, {237,  8}, {236,  8}, {238,  8}, { 23,  9},
    { 28,  9}, { 26,  9}, { 22,  9}, { 24,  9}, { 25,  9}, { 18,  9},
    { 20,  9}, { 21,  9}, { 17,  9}, { 19,  9}, { 15,  9}, { 14,  9},
    { 16,  9}, {217,  9}, { 13,  9}, { 38, 10}, { 36, 10}, { 35, 10},
    { 34, 10}, { 33, 10}, { 37, 10}, { 27, 10}, { 29, 10}, { 31, 10},
    { 32, 10}, { 30, 10}, {SYM_EOF, 11},
    { 39, 11}
};

/* For RLE compressed solexa confidence values */
static huffman_code_t codes_conf_rle[] = {
    {216,  2}, {  0,  4}, {'M',  4}, { 40,  4}, {251,  5}, {245,  6},
    {247,  6}, {  5,  6}, {248,  6}, {250,  6}, {249,  6}, {253,  6},
    {  4,  6}, {252,  6}, {254,  6}, {255,  6}, {  1,  6}, {226,  7},
    { 12,  7}, {231,  7}, {230,  7}, {232,  7}, { 10,  7}, {  3,  7},
    {  9,  7}, {234,  7}, { 11,  7}, {  2,  7}, {235,  7}, {233,  7},
    {237,  7}, {  8,  7}, {236,  7}, {238,  7}, {  7,  7}, {239,  7},
    {241,  7}, {240,  7}, {  6,  7}, {242,  7}, {149,  7}, {243,  7},
    {246,  7}, {244,  7}, { 25,  8}, { 22,  8}, { 23,  8}, { 24,  8},
    { 20,  8}, { 18,  8}, { 21,  8}, {217,  8}, { 19,  8}, {218,  8},
    { 17,  8}, {220,  8}, {221,  8}, {219,  8}, {222,  8}, { 15,  8},
    {223,  8}, { 16,  8}, { 13,  8}, {224,  8}, {229,  8}, { 14,  8},
    {225,  8}, {227,  8}, {228,  8}, { 38,  9}, { 39,  9}, { 34,  9},
    { 35,  9}, { 37,  9}, { 36,  9}, { 33,  9}, { 32,  9}, { 31,  9},
    { 27,  9}, { 29,  9}, { 28,  9}, { 30,  9}, { 26,  9}, { 42, 11},
    {'H', 11}, {'B', 11}, { 45, 11}, { 63, 11}, {'9', 11}, {'0', 11},
    {'5', 12}, {'K', 12}, { 60, 12}, {'N', 12}, {'E', 12}, {'3', 12},
    {'6', 12}, {'I', 13}, {'Z', 13}, {'1', 13}, {'O', 13}, {'T', 13},
    { 64, 13}, {'i', 13}, { 41, 13}, {'W', 13}, {'Q', 13}, { 43, 14},
    { 62, 14}, {'L', 14}, {'R', 14}, {'S', 14}, {'U', 14}, {'f', 14},
    {'8', 14}, {'G', 14}, { 93, 14}, {'c', 14}, { 46, 14}, { 61, 14},
    { 44, 14}, { 58, 14}, { 59, 14}, {'A', 15}, {'F', 15}, {'X', 15},
    { 95, 15}, { 47, 15}, {'2', 15}, {'7', 15}, {'D', 15}, { 94, 15},
    { 96, 15}, {'4', 16}, {'C', 16}, {'J', 16}, {'V', 16}, {'Y', 16},
    {'a', 16}, {'o', 16}, {'P', 16}, { 92, 16}, {SYM_EOF, 17},
    { 91, 17}, {'b', 17}, {'g', 17}, {'n', 17}, {'h', 18}, {'l', 18}
};


/*
 * ---------------------------------------------------------------------------
 * Tree building and code generation functions
 */
static void print_bits(unsigned int val, int nbits) {
    unsigned int mask = 1 << (nbits-1);
    do {
	printf("%d", (val & mask) ? 1 : 0);
	mask >>= 1;
    } while(--nbits);
}

/*
 * Sort huffman_code_t by their code bit-lengths
 */
static int sort_func(const void *p1, const void *p2) {
    const huffman_code_t *c1 = (const huffman_code_t *)p1;
    const huffman_code_t *c2 = (const huffman_code_t *)p2;
    return c1->nbits - c2->nbits;
}

/*
 * Generates canonical huffman codes based on their bit-lengths.
 * This has useful encoding/decoding properties, but also means we only
 * need to store the bit-lengths to regenerate the same code set.
 */
static void canonical_codes(huffman_codes_t *c) {
    int i;
    unsigned int code, last_len;
    huffman_code_t *exception = NULL;

    /* Sort by bit-length */
    qsort(c->codes, c->ncodes, sizeof(*c->codes), sort_func);

    /* Initialise code length lookup tables */
    for (i = 0; i < 258; i++) {
	c->d_tab[i].code_start = 0;
	c->d_tab[i].ncodes = 0;
	c->d_tab[i].symbols = NULL;
    }

    /* Generate codes */
    for (i = 0; i < c->ncodes; i++) {
	int nbits = c->codes[i].nbits;

	if (i == 0) {
	    code = 0;
	    last_len = nbits;
	    c->d_tab[nbits].symbols = (unsigned int *)malloc(c->ncodes * 
							     sizeof(int));
	    c->d_tab[nbits].code_start = code;
	} else {
	    code++;
	}
	if (nbits > last_len) {
	    code <<= (nbits - last_len);
	    last_len = nbits;
	    c->d_tab[nbits].symbols = (unsigned int *)malloc(c->ncodes * 
							     sizeof(int));
	    c->d_tab[nbits].code_start = code;
	}
	c->codes[i].code = code;

	c->d_tab[nbits].symbols[c->d_tab[nbits].ncodes++]
	    = c->codes[i].symbol;
    }

    /* Produce fast lookup table */
    for (i = 0; i < 256; i++) {
	c->lookup[i] = NULL;
    }
    for (i = 0; i < c->ncodes; i++) {
        c->lookup[c->codes[i].symbol] = &c->codes[i];
	if (c->codes[i].symbol == SYM_ANY)
	    exception = &c->codes[i];
    }
    for (i = 0; i < 256; i++) {
	if (!c->lookup[i])
	    c->lookup[i] = exception;
    }
}

static int node_compar(const void *vp1, const void *vp2) {
    return ((const node_t *)vp1)->count - ((const node_t *)vp2)->count;
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
 * data/len is the input data to analyse. legal_chars is a way to restrict the
 * input data down to a specific set of legal characters. We allow for storing
 * out of bounds characters "as-is" by having a SYM_ANY symbol as an escape
 * mechanism. legal_chars may be passed in as NULL in which case codes for all
 * 256 character values are generated.
 *
 * Returns 0 on success
 *         -1 on failure
 */
int calc_bit_lengths(huffman_codes_t *c, unsigned char *data, int len,
		     unsigned char *legal) {
    int i, ncodes, node_start;
    node_t nodes[258+257], *head, *new = &nodes[258];
    int map[258];

    /*
     * Initialise nodes. We build a map of ASCII character code to node
     * number. (By default it's a simple 1:1 mapping unless legal_chars is
     * defined.)
     */
    nodes[0].sym = SYM_EOF;
    nodes[0].count = 1;
    nodes[0].parent = NULL;
    map[SYM_EOF] = 0;

    nodes[1].sym = SYM_ANY;
    nodes[1].count = 1;
    nodes[1].parent = NULL;
    map[SYM_ANY] = 1;

    if (legal) {
	unsigned char *cp;
	int i;
	for (i = 2; i < 258; i++)
	    map[i] = map[SYM_ANY];

	for (i = 2, cp = legal; *cp; cp++, i++) {
	    nodes[i].sym = *cp;
	    nodes[i].count = 0;
	    nodes[i].parent = NULL;
	    map[*cp] = i;
	}
	ncodes = i;
    } else {
	for (i = 0; i < 256; i++) {
	    nodes[i+1].sym = i;
	    nodes[i+1].count = 0;
	    nodes[i+1].parent = NULL;
	    map[i] = i+1;
	}
	ncodes = 257;
    }

    /* Calc freqs */
    for (i = 0; i < len; i++) {
	nodes[map[data[i]]].count++;
    }

    /* Protect against having 258 nodes (see store_codes() for why) */
    if (ncodes == 258) {
	for (i = 0; i < ncodes; i++) {
	    /* Remove SYM_ANY as it's not needed when all 256 chars found */
	    if (nodes[i].sym == SYM_ANY) {
		nodes[i].count = 0;
		break;
	    }
	}
    }

    /* Sort by counts, smallest first and form a sorted linked list */
    qsort(nodes, ncodes, sizeof(*nodes), node_compar);
    for (i = 0; i < ncodes; i++) {
	nodes[i].next = i+1 < ncodes ? &nodes[i+1] : NULL;
	/*
	printf("%p: %d (%c) = %d\n",
	       &nodes[i], nodes[i].sym, nodes[i].sym, nodes[i].count);
	*/
    }

    /* Skip codes with no frequency. Should we do this? */
    for (node_start = 0; node_start < ncodes && nodes[node_start].count == 0;)
	node_start++;
    ncodes -= node_start;

    /* Repeatedly merge two smallest values */
    head = &nodes[node_start];
    while (head && head->next) {
	node_t *after = head->next, *n;
	int sum = head->count + head->next->count;
	
	/*
	printf("Merge node %p (%c/%d) with %p (%c/%d) => %p (%d)\n",
	       head, head->sym, head->count,
	       head->next, head->next->sym, head->next->count,
	       new, sum);
	*/

	for (n = head->next->next; n; after = n, after = n, n = n->next) {
	    if (sum <= n->count)
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
	c->codes[i].nbits = len;
	/*
	print(stderr, "%d(%c) len %d\n",
	      nodes[i+node_start].sym, nodes[i+node_start].sym, len);
	*/
    }

    return 0;
}

/*
 * Initialises and returns a huffman_codes_t struct from a specified code_set.
 * If code_set is not one of the standard predefined values then the
 * input data is analysed using calc_bit_lengths() above to produce the
 * optimal set of huffman codes, otherwise we return predefined values.
 *
 * Returns huffman_codes_t* on success; free using huffman_codes_destroy().
 *         NULL on failure.
 */
huffman_codes_t *generate_code_set(int code_set, unsigned char *data, int len,
				   unsigned char *legal) {
    huffman_codes_t *c;

    if (NULL == (c = (huffman_codes_t *)malloc(sizeof(*c))))
	return NULL;

    c->code_set = code_set;
    
    if (code_set >= 128 || code_set == CODE_INLINE) { 
	if (-1 == calc_bit_lengths(c, data, len, legal)) {
	    free(c);
	    return NULL;
	}
    } else {
	switch(code_set) {
	case CODE_DNA: {
	    int F[256], i;
	    /* Pick optimal DNA code set */
	    memset(F, 0, 256*sizeof(*F));
	    for (i = 0; i < len; i++)
		F[data[i]]++;
	    if (F['A'] < F['C']) {
		if (F['A'] < F['G']) {
		    if (F['A'] < F['T'])
			c->codes = codes_dna_a, code_set = CODE_DNA_A;
		    else
			c->codes = codes_dna_t, code_set = CODE_DNA_T;
		} else {
		    if (F['G'] < F['T'])
			c->codes = codes_dna_g, code_set = CODE_DNA_G;
		    else
			c->codes = codes_dna_t, code_set = CODE_DNA_T;
		}
	    } else {
		if (F['C'] < F['G']) {
		    if (F['C'] < F['T'])
			c->codes = codes_dna_c, code_set = CODE_DNA_C;
		    else
			c->codes = codes_dna_t, code_set = CODE_DNA_T;
		} else {
		    if (F['G'] < F['T'])
			c->codes = codes_dna_g, code_set = CODE_DNA_G;
		    else
			c->codes = codes_dna_t, code_set = CODE_DNA_T;
		}
	    }
	    c->ncodes = sizeof(codes_dna_a)/sizeof(*c->codes);
	    break;
	}

	case CODE_DNA_A:
	    c->codes = codes_dna_a;
	    c->ncodes = sizeof(codes_dna_a)/sizeof(*c->codes);
	    break;

	case CODE_DNA_C:
	    c->codes = codes_dna_c;
	    c->ncodes = sizeof(codes_dna_c)/sizeof(*c->codes);
	    break;

	case CODE_DNA_G:
	    c->codes = codes_dna_g;
	    c->ncodes = sizeof(codes_dna_g)/sizeof(*c->codes);
	    break;

	case CODE_DNA_T:
	    c->codes = codes_dna_t;
	    c->ncodes = sizeof(codes_dna_t)/sizeof(*c->codes);
	    break;

	case CODE_SOLEXA_TRACE:
	    c->codes = codes_solexa_trace;
	    c->ncodes = sizeof(codes_solexa_trace)/sizeof(*c->codes);
	    break;

	case CODE_CONF:
	    c->codes = codes_conf;
	    c->ncodes = sizeof(codes_conf)/sizeof(*c->codes);
	    break;

	case CODE_CONF_RLE:
	    c->codes = codes_conf_rle;
	    c->ncodes = sizeof(codes_conf_rle)/sizeof(*c->codes);
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
	? generate_code_set(code_set, NULL, 0, NULL)
	: NULL;
}

void huffman_codes_destroy(huffman_codes_t *c) {
    int i;

    if (!c)
	return;

    switch(c->code_set) {
    case CODE_DNA:
    case CODE_SOLEXA_TRACE:
    case CODE_CONF:
    case CODE_CONF_RLE:
	/* All these use static tables */
	break;

    default:
	if (c->codes)
	    free(c->codes);
    }

    for (i = 0; i < 258; i++) {
	if (c->d_tab[i].symbols)
	    free(c->d_tab[i].symbols);
    }

    free(c);
}

/*
 * ---------------------------------------------------------------------------
 * Encoding and decoding related functions
 */

/* Can store up to 32-bits worth of data encoded in an integer value */
static void store_bits(block_t *block, unsigned int val, int nbits) {
    /* Slow, but simple */
    unsigned int mask = 1 << (nbits-1);

    /* Ensure enough room */
    while (block->byte + 3 + nbits/8 >= block->alloc) {
	block->alloc += 8192;
	block->data = (unsigned char *)realloc(block->data, block->alloc);
    }

    do {
	int bit = 1 << block->bit;

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
	}
	mask >>= 1;
    } while(--nbits);
}

/* stores nbytes bytes, padding to align on the next byte boundary */
static void store_bytes(block_t *block, unsigned char *val, int nbytes) {
    /* Align */
    if (block->bit != 0) {
	block->byte++;
	block->bit = 0;
    }

    /* Resize */
    while (block->byte + nbytes >= block->alloc) {
	block->alloc += 8192;
	block->data = (unsigned char *)realloc(block->data, block->alloc);
    }

    /* Store */
    memcpy(&block->data[block->byte], val, nbytes);
    block->byte += nbytes;
}

/* Returns bit (0 or 1) on success
 *         -1 on failure or EOF.
 */
static int get_bit(block_t *block) {
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

/*
 * Encodes the huffman symbol bit-lengths as a serialised block of data
 * suitable for storing in a ZTR "HUFF" chunk.
 *
 * Format:
 * 1 byte   = 0 (always zero => RAW encoding format)
 * 1 byte   = SYM_ANY length (0 if SYM_ANY not used).
 * 1 byte   = SYM_EOF length
 * 1 bytes  = total number of remaining symbols (minus SYM_ANY/SYM_EOF)
 * 1 byte   = no. symbols with code len 1 (N1)
 * N1 bytes = symbols
 * 1 byte   = no. symbols with code len 2 (N2)
 * N2 bytes = symbols
 * ...
 * data ends when all symbols are accounted for.
 *
 * Returns: malloced char array of size *comp_len on success
 *          NULL on failure.
 */
unsigned char *store_codes(huffman_codes_t *c, unsigned int *comp_len) {
    int i, j, clen, max_clen;
    unsigned char *str, *cp;

    /*
    fprintf(stderr, "=== CODE SET %d ===\n", c->code_set);
    output_code_set(stderr, c->codes, c->ncodes);
    */

    *comp_len = c->ncodes*2 + 5;
    if (NULL == (str = (unsigned char *)malloc(*comp_len)))
	return NULL;

    assert(c->ncodes-2 >= 0 && c->ncodes-2 <= 255);

    /*
     * Maximum 258 symbols (256+SYM_ANY and SYM_EOF).
     * However if we have all 256 input values measured then the SYM_ANY
     * frequency must be zero. Hence we our maximum number of codes is really
     * 257.
     *
     * Therefore the maximum code length (in the case of frequencies
     * escalating in a fibbonacci sequence) is 256. (ncodes-1).
     * => use bit length 0 to mean bit length 256.
     *
     * Addendum: I just realised this is all horribly hypothetical. In
     * order to achieve a bit-length of 256 we'd need a file of
     * length 6*10^53. Realistically in maxxed out 64-bit file size
     * can have at most bit-lengths of 92 and 4Gb max file gives
     * bit-lengths 45.
     */

    cp = str;
    *cp++ = 0; /* RAW encoding format */
    *cp++ = c->code_set;
    *cp = 0;   /* Default SYM_ANY bit length */

    /* Slow but sure method of encoding. */
    for (i = j = 0; i < c->ncodes; i++) {
	if (c->codes[i].symbol == SYM_ANY) {
	    cp[0] = c->codes[i].nbits;
	    j++;
	} else if (c->codes[i].symbol == SYM_EOF) {
	    cp[1] = c->codes[i].nbits;
	    j++;
	}
    }
    cp += 2;
    *cp++ = c->ncodes - j; /* Number of remaining codes */

    /* Find Maximum symbol length */
    for (max_clen = i = 0; i < c->ncodes; i++) {
	if (max_clen < c->codes[i].nbits) {
	    if (c->codes[i].symbol != SYM_ANY && c->codes[i].symbol != SYM_EOF)
		max_clen = c->codes[i].nbits;
	}
    }

    /*
     * FIXME: Inefficient; faster is we presort or do one-pass to count
     * bit-length frequencies and a 2nd pass to fill out the data stream.
     */
    for (clen = 1; clen <= max_clen; clen++) {
	int count = 0;
	for (i = 0; i < c->ncodes; i++) {
	    if (c->codes[i].symbol == SYM_ANY || c->codes[i].symbol == SYM_EOF)
		continue;
	    if (c->codes[i].nbits == clen)
		count++;
	}
	*cp++ = count;
	for (i = 0; i < c->ncodes; i++) {
	    if (c->codes[i].symbol == SYM_ANY || c->codes[i].symbol == SYM_EOF)
		continue;
	    if (c->codes[i].nbits == clen) {
		*cp++ = c->codes[i].symbol;
	    }
	}
    }

    *comp_len = cp-str;
    return str;
}

/*
 * As store_codes() above, but instead of storing just the symbols used
 * we encode all 258 possible symbols. The benefit is that when most symbols
 * are present anyway then reordering and writing out every single one can
 * produce a smaller output string. (As we know all are present we
 * only then need to store the bit-lengths rather than the symbols
 * themselves.)
 *
 * Output format:
 * 1 byte      zero (ZTR "RAW" format byte)
 * 1? byte     SYM_ANY length
 * 1? byte     SYM_EOF length
 * 1 byte      N symbols (N+1 infact)
 * N? bytes    bit-lengths for symbols 0 to N inclusive.
 *
 * The "1? byte" is because most bit-lengths are one byte except for the
 * rarely needed bit-length 255 and 256 which need 2 byte
 * encodings. (See source.)
 *
 * NB: unused at present, but ideally when ncodes is high (> 100?) we should
 * use this representation instead and then compress the result using
 * huffman compression and the old store_codes() storage mechanism.
 * Experimentally on a solexa SMP4 segments this reduces the code string
 * size from 259 bytes to 79 bytes.
 */
static unsigned char *store_codes_all(huffman_codes_t *c, unsigned int *length)
{
    int i;
    unsigned char *str, *cp;
    int blen[258];

    /*
     * 257 possible codes => 256 possible bit lengths.
     * Therefore we cannot distinguish between code not used and code
     * length 256. Instead we encode some values in 2 bytes.
     * unused: 0
     * 255:    255 255
     * 256:    256 256
     */
    *length = 300; /* enough! */
    if (NULL == (str = (unsigned char *)malloc(*length)))
	return NULL;

    /* Lookup table indexed on symbol, returning bit-length used */
    memset(blen, 0, 258 * sizeof(*blen));
    for (i = 0; i < c->ncodes; i++) {
	blen[(unsigned int)c->codes[i].symbol] = c->codes[i].nbits;
    }

    cp = str;
    *cp++ = 0; /* RAW encoding format */

    /* Output all 258 bit-lengths */
    for (i = 0; i < 258; i++) {
	switch (blen[i]) {
	case 255:
	    *cp++ = 255;
	    *cp++ = 255;
	    break;
	case 256:
	    *cp++ = 255;
	    *cp++ = 0;
	    break;
	default:
	    *cp++ = blen[i];
	}
    }

    *length = cp-str;
    return str;
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
huffman_codes_t *restore_codes(unsigned char *data, int *len_used) {
    unsigned char *cp = data;
    int nfound, clen;
    int i, j = 0;
    huffman_codes_t *c;

    /* Sanity checking and memory allocs */
    if (*cp++ != 0) /* RAW encoding byte */
    	return NULL;

    if (NULL == (c = (huffman_codes_t *)malloc(sizeof(*c))))
	return NULL;

    c->code_set = *cp++;
    if (c->code_set < CODE_USER && c->code_set != CODE_INLINE)
	/* 1-127 are reserved */
	return NULL;

    c->codes = (huffman_code_t *)malloc(258*sizeof(*c->codes));
    if (NULL == c->codes) {
	free(c->codes);
	return NULL;
    }

    /* User specified codes - load them up and initialise huffman system */
    c->codes[j].symbol = SYM_ANY;
    c->codes[j].nbits = *cp++;
    if (c->codes[j].nbits) /* SYM_ANY is optional */
	j++;
    c->codes[j].symbol = SYM_EOF;
    c->codes[j++].nbits = *cp++;

    c->ncodes = *cp++;
    if (c->ncodes == 0)
	c->ncodes = 256;
    for (clen = 1, nfound = 0; nfound < c->ncodes; clen++) {
	int run_len = *cp++;
	for (i = 0; i < run_len; i++) {
	    c->codes[j].symbol = *cp++;
	    c->codes[j++].nbits  = clen;
	    nfound++;
	}
    }

    if (len_used)
	*len_used = cp-data;
    c->ncodes = j; /* EOF and ANY */

    canonical_codes(c);

    /*
    fprintf(stderr, "=== CODE SET %d ===\n", c->code_set);
    output_code_set(stderr, c->codes, c->ncodes);
    */

    return c;
}

/*
 * Given a set of huffman codes and a block of data this compresses
 * and returns the data block.
 *
 * Returns: malloced compressed data on success, of length *nbytes.
 *          NULL on failure.
 */
unsigned char *huffman_encode(huffman_codes_t *c, int code_set,
			      unsigned char *data, int len, int *nbytes) {
    int i;
    block_t out;
    unsigned char *code_str = NULL;

    if (!c) {
	/* No code set, so derive our own */
	c = generate_code_set(code_set, data, len, NULL);
    }

    out.data = calloc(len+2, 1); /* avoids most realloc needs */
    out.alloc = len+2;
    out.byte = 0;
    out.data[out.byte++] = ZTR_FORM_STHUFF;
    out.data[out.byte++] = c->code_set;
    out.bit = 0;

    if (c->code_set == 0) {
	/* inlined code-set => output them here */
	unsigned int clen;
	code_str = store_codes(c, &clen);
	store_bytes(&out, code_str, clen);
	free(code_str);
    }

    for (i = 0; i < len; i++) {
	if (c->lookup[data[i]]) {
	    store_bits(&out, c->lookup[data[i]]->code,
		       c->lookup[data[i]]->nbits);
	    if (c->lookup[data[i]]->symbol == SYM_ANY) {
		store_bits(&out, data[i], 8);
	    }
	}
    }
    store_bits(&out, c->lookup[SYM_EOF]->code, c->lookup[SYM_EOF]->nbits);

    *nbytes = out.byte + (out.bit != 0);
    out.data = realloc(out.data, *nbytes); /* shrink to real size needed */
    return out.data;
}

/*
 * The opposite of huffman_encode().
 *
 * Returns: malloced uncompressed data on success, of length *nbytes.
 *          NULL on failure.
 */
unsigned char *huffman_decode(huffman_codes_t *c,
			      unsigned char *data, int len, int *nbytes) {
    block_t in, out;
    int b, nbits = 0;
    unsigned int val = 0;
    int private_codes = 0;

    in.data = data;
    in.alloc = len;
    in.byte = 1;
    in.bit = 0;

    out.data = NULL;
    out.alloc = 0;
    out.byte = 0;
    out.bit = 0;

    if (in.data[in.byte] == 0) {
	int used;
	unsigned char *data = &in.data[in.byte];

	/* inline code-set, so read */
	c = restore_codes(data, &used);
	private_codes = 1;
	in.byte += used;
    } else if (in.data[in.byte] < CODE_USER) {
	/* static predefined codes */
	c = get_code_set(in.data[in.byte]);
	private_codes = 1;
	in.byte++;
    } else {
	/* Use passed in supplied code_set */
	in.byte++;
    }

    while((b = get_bit(&in)) != -1) {
	int sym;
	val = (val << 1) + b;
	nbits++;
	if (val >= c->d_tab[nbits].code_start &&
	    val <  c->d_tab[nbits].code_start + c->d_tab[nbits].ncodes) {
	    sym = c->d_tab[nbits].symbols[val-c->d_tab[nbits].code_start];
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

    if (private_codes)
	huffman_codes_destroy(c);

    *nbytes = out.byte;
    return out.data;
}


/*
 * ---------------------------------------------------------------------------
 * Debug code. This turns the library into a stand-alone program for
 * easy debugging.x
 */

static void dump_decode_table(huffman_codes_t *c) {
    int i, j;
    for (i = 0; i < 258; i++) {
	if (!c->d_tab[i].ncodes)
	    continue;
	
	printf("=== len %d ==\n", i);
	for (j = 0; j < c->d_tab[i].ncodes; j++) {
	    printf("%c ",
		   isprint(c->d_tab[i].symbols[j])
		   ? c->d_tab[i].symbols[j]
		   : '?');
	    print_bits(c->d_tab[i].code_start + j, i);
	    printf("\n");
	}
    }
}

static void output_code_set(FILE *fp, huffman_code_t *codes, int ncodes) {
    int i, j;

    fprintf(fp, "static huffman_code_t codes_FIXME[] = {\n");
    for (i = j = 0; i < ncodes; i++) {
	if (j == 0)
	    fprintf(fp, "    ");
	if (codes[i].symbol == SYM_ANY) {
	    fprintf(fp, "{SYM_ANY,%3d}, ", codes[i].nbits);
	    j = 10;
	} else if (codes[i].symbol == SYM_EOF) {
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
}

#if 0
static unsigned char *load(int *lenp) {
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

int main(int argc, char **argv) {
    unsigned char *data, *out;
    int len, out_len;
    int decode = 0;
    int dump_tree = 0, exit_after_tree = 0;
    int code_set = CODE_INLINE;
    char *legal_chars = NULL;
    int c;
    huffman_codes_t *cds;

    while ((c = getopt(argc, argv, "c:detxl:")) != -1) {
	switch (c) {
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

	case 'l':
	    legal_chars = optarg;
	    break;

	default:
	    fprintf(stderr, "Usage: huffman_static [options] < stdin > stdout\n");
	    fprintf(stderr, "    Decoding options\n");
	    fprintf(stderr, "        -d\tdecode flag\n");
	    fprintf(stderr, "    Encoding options\n");
	    fprintf(stderr, "        -e\tencode flag\n");
	    fprintf(stderr, "        -c code\tspecify code-set. 0 => inline\n");
	    fprintf(stderr, "        -l chars\tspecify legal characters to use in auto code-set\n");
	    fprintf(stderr, "        -t\tpretty-print the code-set used\n");
	    fprintf(stderr, "        -x\texit after outputting code-set\n");
	    exit(1);
	}
    }

    data = load(&len);

    if (decode) {
	out = huffman_decode(NULL, data, len, &out_len);
    } else {
	/* Encoding */
	cds = generate_code_set(code_set, data, len, legal_chars);

	if (dump_tree) {
	    output_code_set(stdout, cds->codes, cds->ncodes);
	    if (exit_after_tree)
		return 0;
	}

	out = huffman_encode(cds, data, len, &out_len);
	huffman_codes_destroy(cds);
    }

    write(1, out, out_len);

    free(data);
    free(out);

    return 0;
}
#endif
