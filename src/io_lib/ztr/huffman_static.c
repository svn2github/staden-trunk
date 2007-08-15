#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <unistd.h>

#include "huffman_static.h"


/* #define TIME_IT */
/* #define TEST_MAIN */

#ifdef TIME_IT
#  include <sys/time.h>
#endif

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

#ifdef SLOW_DECODE
static void dump_decode_table(huffman_codes_t *c);
#endif
static void output_code_set(FILE *fp, huffman_code_t *codes, int ncodes);
static void output_code_set2(FILE *fp, huffman_code_t *codes, int ncodes);

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

/* English text */
static huffman_code_t codes_english[] = {
    {'e',  3}, { 32,  3}, {'s',  4}, {'i',  4}, {'n',  4}, {'o',  4},
    {'a',  4}, {'t',  4}, {'u',  5}, {'l',  5}, {'d',  5}, {'r',  5},
    {'h',  5}, {'y',  6}, { 44,  6}, {'p',  6}, {'g',  6}, {'w',  6},
    {'f',  6}, {'m',  6}, { 13,  6}, {'c',  6}, { 10,  6}, {'v',  7},
    { 46,  7}, {'b',  7}, { 34,  8}, {'I',  8}, {'k',  8}, {'A',  9},
    {'N',  9}, {'T',  9}, { 45,  9}, {'B', 10}, {'E', 10}, {'M', 10},
    {'H', 10}, { 39, 10}, {'S', 10}, {'W', 10}, { 59, 10}, {'x', 10},
    {'C', 10}, {'z', 11}, {'G', 11}, {'R', 11}, {'F', 11}, {'0', 11},
    {'O', 11}, {'1', 11}, {'j', 11}, {'P', 11}, {'q', 11}, { 33, 11},
    {'L', 11}, { 63, 11}, {'5', 12}, {'U', 12}, { 58, 12}, {'3', 12},
    {'2', 12}, {'Y', 12}, {'D', 12}, {  8, 13}, { 41, 13}, { 40, 13},
    {'9', 13}, {'V', 13}, {'J', 13}, {'8', 13}, {'4', 13}, {'7', 13},
    {'6', 13}, {'X', 14}, { 42, 14}, { 64, 14}, { 47, 14}, {'K', 14},
    { 95, 14}, {'Z', 16}, { 92, 16}, { 96, 16}, {'Q', 16}, { 35, 17},
    { 36, 17}, {124, 17}, { 91, 17}, { 93, 17}, { 37, 18}, { 43, 19},
    { 60, 19}, { 61, 19}, { 62, 19}, { 94, 19}, {123, 19}, {125, 19},
    {126, 19}, { 38, 19}, {SYM_EOF, 20}, {SYM_ANY, 20}
};

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
 */
static void canonical_codes(huffman_codes_t *c) {
    int i;
    unsigned int code, last_len;
    huffman_code_t *exception = NULL;
    int clens[33];
    int offs[33];
    huffman_code_t ctmp[258];

    /* Sort by bit-length - much faster than qsort() */
    for (i = 0; i <= 32; i++)
	offs[i] = clens[i] = 0;
    for (i = 0; i < c->ncodes; i++)
	clens[c->codes[i].nbits]++;
    for (i = 1; i <= 32; i++)
	offs[i] = offs[i-1] + clens[i-1];
    for (i = 0; i < c->ncodes; i++)
	ctmp[offs[c->codes[i].nbits]++] = c->codes[i];
    memcpy(c->codes, ctmp, c->ncodes * sizeof(huffman_code_t));

    /* Ensure SYM_ANY and SYM_EOF are the first of their bit-length */
    for (i = 0; i < c->ncodes; i++) {
	if (c->codes[i].symbol == SYM_ANY || c->codes[i].symbol == SYM_EOF) {
	    int j;
	    for (j = i-1; j >= 0; j--) {
		if (c->codes[j].nbits == c->codes[j+1].nbits) {
		    huffman_code_t tmp = c->codes[j];
		    c->codes[j] = c->codes[j+1];
		    c->codes[j+1] = tmp;
		}
	    }
	}
    }

    /* Initialise code length lookup tables */
#ifdef SLOW_DECODE
    for (i = 0; i < 258; i++) {
    	c->d_tab[i].code_start = 0;
    	c->d_tab[i].ncodes = 0;
    	c->d_tab[i].symbols = NULL;
    }
#endif

    /* Generate codes */
    code = last_len = 0; /* stop warning */
    for (i = 0; i < c->ncodes; i++) {
	int nbits = c->codes[i].nbits;

	if (i == 0) {
	    code = 0;
	    last_len = nbits;
#ifdef SLOW_DECODE
	    c->d_tab[nbits].symbols = (unsigned int *)malloc(c->ncodes * 
							     sizeof(int));
	    c->d_tab[nbits].code_start = code;
#endif
	} else {
	    code++;
	}
	if (nbits > last_len) {
	    code <<= (nbits - last_len);
	    last_len = nbits;
#ifdef SLOW_DECODE
	    c->d_tab[nbits].symbols = (unsigned int *)malloc(c->ncodes * 
							     sizeof(int));
	    c->d_tab[nbits].code_start = code;
#endif
	}
	c->codes[i].code = bit_reverse(code, nbits);

#ifdef SLOW_DECODE
	c->d_tab[nbits].symbols[c->d_tab[nbits].ncodes++]
	    = c->codes[i].symbol;
#endif
    }

    /* Reindex so the symbol is the primary index into codes */
    for (i = 0; i < 256; i++) {
	c->lookup[i].nbits = 0;
    }
    for (i = 0; i < c->ncodes; i++) {
        c->lookup[c->codes[i].symbol] = c->codes[i];
	if (c->codes[i].symbol == SYM_ANY)
	    exception = &c->codes[i];
    }
    if (exception) {
	for (i = 0; i < 256; i++) {
	    if (!c->lookup[i].nbits)
		c->lookup[i] = *exception;
	}
    }
}

static int node_compar(const void *vp1, const void *vp2) {
    const node_t *n1 = (const node_t *)vp1;
    const node_t *n2 = (const node_t *)vp2;

    /*
     * The sort order is vital here. This needs to return the same collating
     * order on all systems so that differing qsort() functions will not
     * swap around symbols with the same bit lengths, hence we sort by both
     * fields to force a unique stable ordering.
     *
     * Secondly it is important that we store high symbols first as this
     * is exploited by the store and restore codes functions. The SYM_ANY
     * and SYM_EOF codes are stored out of order, so we need to make sure
     * after resorting they end up in the correct location after
     * canonical codes runs. Hence we put high values first.
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
	if (len > 24) {
	    fprintf(stderr, "Bit lengths too large for current optimisations\n");
	    return -1;
	}
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
			c->codes = codes_dna_a, c->code_set = CODE_DNA_A;
		    else
			c->codes = codes_dna_t, c->code_set = CODE_DNA_T;
		} else {
		    if (F['G'] < F['T'])
			c->codes = codes_dna_g, c->code_set = CODE_DNA_G;
		    else
			c->codes = codes_dna_t, c->code_set = CODE_DNA_T;
		}
	    } else {
		if (F['C'] < F['G']) {
		    if (F['C'] < F['T'])
			c->codes = codes_dna_c, c->code_set = CODE_DNA_C;
		    else
			c->codes = codes_dna_t, c->code_set = CODE_DNA_T;
		} else {
		    if (F['G'] < F['T'])
			c->codes = codes_dna_g, c->code_set = CODE_DNA_G;
		    else
			c->codes = codes_dna_t, c->code_set = CODE_DNA_T;
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
	? generate_code_set(code_set, NULL, 0, NULL)
	: NULL;
}

void huffman_codes_destroy(huffman_codes_t *c) {
    if (!c)
	return;

    switch(c->code_set) {
    case CODE_DNA:
    case CODE_DNA_A:
    case CODE_DNA_C:
    case CODE_DNA_G:
    case CODE_DNA_T:
    case CODE_ENGLISH:
	/* All these use static tables */
	break;

    default:
	if (c->codes)
	    free(c->codes);
    }

#ifdef SLOW_DECODE
    for (i = 0; i < 258; i++) {
    	if (c->d_tab[i].symbols)
	    free(c->d_tab[i].symbols);
    }
#endif

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
#if 0
    /* NB: huffman_encode ensures we don't need this check now */
    /* Ensure enough room */
    while (block->byte + 3 + nbits/8 >= block->alloc) {
	block->alloc += 8192;
	block->data = (unsigned char *)realloc(block->data, block->alloc);
    }
#endif

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


/*
 * Encodes the huffman symbol bit-lengths as a serialised block of data
 * suitable for storing in a ZTR "HUFF" chunk.
 *
 * Format:
 * 1 byte   = 0 (always zero => RAW encoding format)
 * 1 byte   = code_set number
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
#if 1
unsigned char *store_codes(huffman_codes_t *c, unsigned int *comp_len) {
    int i, j;
    unsigned char *str;
    int offsets[257];
    int lenfreq[257], maxlen;

    /*
    fprintf(stderr, "=== CODE SET %d ===\n", c->code_set);
    output_code_set(stderr, c->codes, c->ncodes);
    */

    *comp_len = c->ncodes*2 + 5;
    if (NULL == (str = (unsigned char *)malloc(*comp_len)))
	return NULL;

    assert(c->ncodes-2 >= 0 && c->ncodes-2 <= 255);

    /*
     * Algorithm:
     * Count the number of occurrences of each bit-length.
     * From this compute the starting point for each bit-length.
     * Then loop through symbols storing in their correct location.
     */

    /* Count the number of occurrences of each bit-length. */
    for (i = 0; i < 257; i++)
	lenfreq[i] = offsets[i] = 0;
    for (maxlen = i = 0; i < c->ncodes; i++) {
	if (c->codes[i].symbol == SYM_ANY)
	    continue;
	if (c->codes[i].symbol == SYM_EOF)
	    continue;
	lenfreq[c->codes[i].nbits]++;
	if (maxlen < c->codes[i].nbits)
	    maxlen = c->codes[i].nbits;
    }

    /* Compute the starting point for each bit-length. */
    offsets[1] = 5; /* code_set, SYM_ANY len, SYM_EOF len, nsymbols */
    str[offsets[1]++] = lenfreq[1];
    for (i = 2; i <= maxlen; i++) {
	offsets[i] = offsets[i-1] + lenfreq[i-1];
	str[offsets[i]++] = lenfreq[i];
    }
    *comp_len = offsets[maxlen] + lenfreq[maxlen];

    /* Store symbols in their correct location */
    for (i = 0; i < c->ncodes; i++) {
	if (c->codes[i].symbol == SYM_ANY)
	    continue;
	if (c->codes[i].symbol == SYM_EOF)
	    continue;

#if 0
	fprintf(stderr, "str[%d]=%02x (nbits=%d offsets=%d)\n",
		offsets[c->codes[i].nbits],
		c->codes[i].symbol,
		c->codes[i].nbits,
		offsets[c->codes[i].nbits]);
#endif
	str[offsets[c->codes[i].nbits]++] = c->codes[i].symbol;
    }

    /* Add initial formatting bits and SYM_ANY/EOF values*/
    str[0] = 0; /* RAW format */
    str[1] = c->code_set;
    str[2] = 0; /* default SYM_ANY */
    for (i = j = 0; i < c->ncodes; i++) {
	if (c->codes[i].symbol == SYM_ANY) {
	    str[2] = c->codes[i].nbits;
	    j++;
	} else if (c->codes[i].symbol == SYM_EOF) {
	    str[3] = c->codes[i].nbits;
	    j++;
	}
    }
    str[4] = c->ncodes - j; /* Number of remaining codes */

    return str;
}

#else

/* Slow old version */

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
     *
     * If we limit the code to fitting a code in a 32-bit integer then
     * this means the input data size can be at max 5.7Mb long.
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
#endif

#if 0
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
#endif

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

    //fprintf(stderr, "=== CODE SET %d ===\n", c->code_set);
    //output_code_set(stderr, c->codes, c->ncodes);
    //dump_decode_table(c);

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
    huffman_code_t *lookup;

#ifdef TIME_IT
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
#endif    

    if (!c) {
	/* No code set, so derive our own */
	c = generate_code_set(code_set, data, len, NULL);
    }

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "code production took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
    gettimeofday(&t1, NULL);
#endif    

    out.alloc = len*4+10;
    out.data = calloc(out.alloc, 1); /* worst case, avoids realloc */
    out.byte = 0;
    out.data[out.byte++] = ZTR_FORM_STHUFF;
    out.bit = 0;

    if (c->code_set == 0) {
	/* inlined code-set => output them here */
	unsigned int clen;
	code_str = store_codes(c, &clen);
	store_bytes(&out, code_str, clen);
	free(code_str);
    } else {
	out.data[out.byte++] = c->code_set;
    }

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "alloc + storing codes took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
    gettimeofday(&t1, NULL);
#endif    

    lookup = c->lookup;
    if (lookup[SYM_ANY].nbits) {
	for (i = 0; i < len; i++) {
	    store_bits(&out, lookup[data[i]].code,
		       lookup[data[i]].nbits);
	    if (lookup[data[i]].symbol == SYM_ANY) {
		store_bits(&out, data[i], 8);
	    }
	}
    } else {
	for (i = 0; i < len; i++) {
	    store_bits(&out, lookup[data[i]].code,
		       lookup[data[i]].nbits);
	}
    }
    store_bits(&out, c->lookup[SYM_EOF].code, c->lookup[SYM_EOF].nbits);

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "encoding took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
    gettimeofday(&t1, NULL);
#endif    

    *nbytes = out.byte + (out.bit != 0);
    out.data = realloc(out.data, *nbytes); /* shrink to real size needed */
    return out.data;
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
 *
 * Included in this is some code to simplify the SYM_ANY handling too
 * (not including in the k=8 case below this function), by adding an
 * additional 256 entries to the code set and pretending SYM_ANY 
 * is an internal node.
 */
unsigned char *huffman_decode(huffman_codes_t *c,
			      unsigned char *data, int len, int *nbytes) {
    block_t in, out;
    htree_t t[513]; /* number of internal nodes */
    int i, j, n;
    int private_codes = 0;
    int new_node;
    h_jump4_t J4[513][16];
    int need_any = 0;

#ifdef TIME_IT
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
#endif    

    in.data = data;
    in.alloc = len;
    in.byte = 1;
    in.bit = 0;

    out.alloc = 8 * len; /* max length */
    out.data = (unsigned char *)malloc(out.alloc);
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

	if (c->codes[i].symbol == SYM_ANY) {
	    need_any = i+1;
	    continue;
	}

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

    /*
     * To handle SYM_ANY we basically just add on another 256 symbols of
     * length len(SYM_ANY)+8.
     */
    if (need_any) {
	int any_code  = c->codes[need_any-1].code;
	int any_nbits = c->codes[need_any-1].nbits;
	for (i = 0; i < 256; i++) {
	    int n = 0;
	    int v = any_code + (i << any_nbits);

	    for (j = 0; j < any_nbits+7; j++) {
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
	    t[n].l[v & 1] = i;
	}
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
	 unsigned char *cp = out.data;
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
	 for (i = in.byte; i < len; i++) {
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
	 out.byte = cp - out.data;
	 out.data = (unsigned char *)realloc(out.data, out.byte+1);
     }

#ifdef TIME_IT
    gettimeofday(&t2, NULL);
    fprintf(stderr, "decoding2 took %ld microsec\n",
	    (t2.tv_sec -  t1.tv_sec) * 1000000 +
	    (t2.tv_usec - t1.tv_usec));
#endif    

    if (private_codes)
	huffman_codes_destroy(c);

    *nbytes = out.byte;
    return out.data;
}


#if 0
/*
 * Method 1b
 * ---------
 * As method 1 above, but here we use k=8 so we can read an entire
 * byte at a time.
 */

/*
 * Store 9 bit data as 8-bits + top-bit. Top-bit only needed for the
 * special case of SYM_EOF or SYM_ANY, where upon 'special' is true.
 */
typedef struct {
    /* Byte-wise jumping table */
    unsigned short jump;
    unsigned char symbol[8];
    unsigned char nsymbols:4;
    unsigned char special:4; /* set to >0 if SYM_EOF or SYM_ANY in symbol[] */
    unsigned char top_bit;   /* bit 9 of symbol[], valid if special != 0 */
} h_jump8_t;

unsigned char *huffman_decode(huffman_codes_t *c,
			      unsigned char *data, int len, int *nbytes) {
    block_t in, out;
    htree_t t[257]; /* number of internal nodes */
    int i, j, k, n, b, nbits = 0;
    unsigned int val = 0;
    int private_codes = 0;
    int new_node;
    h_jump8_t J[257][256]; /* jump table, indexed by node and byte */
    h_jump8_t B[7][256]; /* for restarting from 1-7 bits; jumps back into J[] */
    int need_any = 0;

    in.data = data;
    in.alloc = len;
    in.byte = 1;
    in.bit = 0;

    out.alloc = 8 * len; /* max length */
    out.data = (unsigned char *)malloc(out.alloc);
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
    
    /* Construct the tree from the codes */
    new_node = 1;
    t[0].l[0] = t[0].l[1] = -1;
    t[0].c[0] = t[0].c[1] = 0;
    for (i = 0; i < c->ncodes; i++) {
	int n = 0;
	unsigned int v = c->codes[i].code;

	if (c->codes[i].symbol == SYM_ANY)
	    need_any = 1;

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

#if 0
    for (n = 0; n < new_node; n++) {
	printf("Node %d, c[]={%d,%d}, l[]={%d,%d}\n",
	       n, t[n].c[0], t[n].c[1], t[n].l[0], t[n].l[1]);
    }
    printf("\n");
#endif

    /* Build the 256 wide lookup table per node */
    for (n = 0; n < new_node; n++) {
	for (j = 0; j < 256; j++) {
	    unsigned int v = j;
	    int n2 = n;
	    h_jump8_t *hj = &J[n][j];
	    hj->nsymbols = 0;
	    hj->special = 0;
	    hj->top_bit = 0;
	    for (i = 0; i < 8; i++) {
		int b = v & 1;
		if (t[n2].l[b] >= 0) {
		    hj->symbol[hj->nsymbols++] = t[n2].l[b] & 0xff;
		    if (t[n2].l[b] == SYM_EOF || t[n2].l[b] == SYM_ANY)
			if (!hj->special) {
			    hj->special = i+1; /* next bit number */
			    hj->top_bit |= 1 << (hj->nsymbols-1);
			}
		}
		n2 = t[n2].c[b];
		v >>= 1;
	    }
	    hj->jump = n2;
	}
    }

    /*
     * Build 7x256 lookup tables to act as jump points coming from SYM_ANY
     * decoding. Here we just want the last N bits (7 to 1 inclusive, as all
     * 8 is just a matter of using the standard lookup tables).
     * Each time the starting state will be node 0.
     */
    if (need_any) {
	for (n = 0; n < 7; n++) {
	    for (j = 0; j < 256; j++) {
		unsigned int v = j;
		int n2 = 0;
		h_jump8_t *hj = &B[n][j];
		hj->nsymbols = 0;
		hj->special = 0;
		hj->top_bit = 0;
		v >>= (7-n); /* last n bits only */
		for (i = 0; i <= n; i++) {
		    int b = v & 1;
		    if (t[n2].l[b] >= 0) {
			hj->symbol[hj->nsymbols++] = t[n2].l[b];
			if (t[n2].l[b] == SYM_EOF || t[n2].l[b] == SYM_ANY)
			    if (!hj->special) {
				hj->special = (7-n)+i+1; /* next bit number */
				hj->top_bit |= 1 << (hj->nsymbols-1);
			    }
		    }
		    n2 = t[n2].c[b];
		    v >>= 1;
		}
		hj->jump = n2;
	    }
	}
    }

    /* Fast byte-wise lookup */
    {
	int node_num = 0;
	unsigned char *cp = out.data;

	for (i = in.byte; i < len; i++) {
	    h_jump8_t *x = &J[node_num][data[i]];
	    int l = x->nsymbols;

	restart: /* Sorry! */
	    if (!x->special) {
		for (j = 0; j < l; j++) {	
		    *cp++ = x->symbol[j];
		}
	    } else {
		for (j = 0; j < l; j++) {
		    unsigned short sym = x->symbol[j] |
			(((x->top_bit & (1<<j)) != 0) << 8);
		    if (sym == SYM_EOF) {
			goto term; /* I know, I know... */
		    } else if (sym == SYM_ANY) {
			int b, v, byte = data[i];
			/*
			 * This is nasty. SYM_ANY means the next 8-bits
			 * encode a complete uncompressed symbol.
			 * So we need to break out of the byte-by-byte
			 * decoding method and manually extract the 8-bits
			 * to obtain the escaped symbol. Ugly, but the
			 * fastest way to deal with exceptions without adding
			 * all possible 256 SYM_ANY values into the graph.
			 *
			 * We need to know the bit-number that the SYM_ANY
			 * occured. We know we can only get a maximum of 1
			 * per byte, so we store the bit_number+1 (ie
			 * next bit) in the 'special' flag.
			 */
			for (b = v = 0; b < 8; b++, x->special++) {
			    if (x->special >= 8) {
				byte = data[++i];
				x->special = 0;
			    }
			    v = (v << 1) | ((byte & (1<<x->special)) != 0);
			}
			*cp++ = bit_reverse(v, 8);

			/*
			 * Now we need to manually restart. We use the 7
			 * B[] lookup tables to restart from, or J[] is we're
			 * on an exact byte boundary.
			 */
			x = (x->special == 8)
			    ? &J[0][data[++i]]
			    : &B[7 - x->special][byte];
			l = x->nsymbols;
			goto restart;
		    } else {
			*cp++ = x->symbol[j];
		    }
		}
	    }

	    node_num = x->jump;
	}
    term:
	out.byte = cp - out.data;

	out.data = (unsigned char *)realloc(out.data, out.byte+1);
    }

    if (private_codes)
	huffman_codes_destroy(c);

    *nbytes = out.byte;
    return out.data;
}
#endif

#ifdef SLOW_DECODE
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

unsigned char *huffman_decode(huffman_codes_t *c,
			      unsigned char *data, int len, int *nbytes) {
    block_t in, out;
    int i, j, k, n, b, nbits = 0;
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

    /* Slow but easy to understand method */
    while((b = get_bit(&in)) != -1) {
	int sym;
	val = (val << 1) + b;
	nbits++;
	if (val >= c->d_tab[nbits].code_start &&
	    val <  c->d_tab[nbits].code_start + c->d_tab[nbits].ncodes) {
	    sym = c->d_tab[nbits].symbols[val-c->d_tab[nbits].code_start];
	    if (sym == SYM_ANY) {
		/* Cannot use one statement due to no sequence-points */
		sym  = (get_bit(&in) << 7);
		sym += (get_bit(&in) << 6);
		sym += (get_bit(&in) << 5);
		sym += (get_bit(&in) << 4);
		sym += (get_bit(&in) << 3);
		sym += (get_bit(&in) << 2);
		sym += (get_bit(&in) << 1);
		sym += (get_bit(&in) << 0);
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
#endif


/*
 * ---------------------------------------------------------------------------
 * Debug code. This turns the library into a stand-alone program for
 * easy debugging.x
 */
#ifdef SLOW_DECODE
static void print_bits(unsigned int val, int nbits) {
    unsigned int mask = 1 << (nbits-1);
    do {
	printf("%d", (val & mask) ? 1 : 0);
	mask >>= 1;
    } while(--nbits);
}

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
#endif

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

static void output_code_set2(FILE *fp, huffman_code_t *codes, int ncodes) {
    int i;

    fprintf(fp, "huffman_code_t = {\n");
    for (i = 0; i < ncodes; i++) {
	fprintf(fp, "\t%d:\t%d %2d %04x\n",
		i,codes[i].symbol, codes[i].nbits, codes[i].code);
    }
    fprintf(fp, "};\n");
}


/*
 * --------------------------------------------------------------------------
 * A test main() to create an application capable of compressing and
 * uncompressing stdin.
 */

#ifdef TEST_MAIN
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
	cds = generate_code_set(code_set, data, len,
				(unsigned char *)legal_chars);
	if (!cds)
	    return 1;

	if (dump_tree) {
	    output_code_set(stdout, cds->codes, cds->ncodes);
	    output_code_set2(stdout, cds->codes, cds->ncodes);
	    if (exit_after_tree)
		return 0;
	}

	out = huffman_encode(cds, code_set, data, len, &out_len);
	huffman_codes_destroy(cds);
    }
    write(1, out, out_len);

    free(data);
    free(out);

    return 0;
}
#endif
