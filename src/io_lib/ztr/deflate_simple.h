#ifndef _DEFLATE_SIMPLE_H_
#define _DEFLATE_SIMPLE_H_

#ifdef __cplusplus
extern "C" {
#endif

/* inlined codes */
#define CODE_INLINE	  0   

/* predefined codes */
#define CODE_DNA	  1   /* DNA, uppercase only */
#define CODE_DNA_AMBIG	  2   /* DNA, uc with ambiguity codes */
#define CODE_ENGLISH      3   /* English text */

/* predefined elsewhere in HUFF chunks, 128 onwards */
#define CODE_USER	  128 

#define MAX_CODE_LEN	  15 /* maximum allowed by RFC 1951 */

#ifndef ZTR_FORM_STHUFF
#  define ZTR_FORM_STHUFF  77
#endif

/* A single symbol and it's encoding */
typedef struct {
    signed int symbol; /* 0-255 character, 256 = exception code, 257 = EOF */
    int nbits;
    unsigned int code;
    int freq;
} huffman_code_t;

/* Used for canonical codes. Maps symbols for the same bit-length to codes */
typedef struct {
    unsigned int code_start;
    int ncodes;
    unsigned int *symbols;
} decode_t;

/* A collection of huffman_code_t along with decoding optimisations */
typedef struct {
    huffman_code_t *codes;
    int ncodes;
    int code_set; /* (128-255) The user specified number for this encoding */
    huffman_code_t lookup[258]; /* Mapping of symbol character to code */
    int max_code_len;
} huffman_codes_t;

/* Use for store_bits() and get_bits() */
typedef struct block {
    unsigned char *data;
    size_t alloc;
    size_t byte;
    int bit;
} block_t;

block_t *block_create(size_t size);
void block_destroy(block_t *blk, int keep_data);
int block_resize(block_t *blk, size_t size);


int huffman_encode(block_t *blk, huffman_codes_t *c, int code_set,
		   unsigned char *data, int len);

unsigned char *huffman_decode(huffman_codes_t *c,
			      unsigned char *data, int len, int *nbytes);

huffman_codes_t *generate_code_set(int code_set, unsigned char *data, int len,
				   int eof, int max_code_len, int all_codes);
int store_codes(block_t *out,
		huffman_codes_t *c,
		int last_block);
huffman_codes_t *restore_codes(unsigned char *data, int *len_used);
void huffman_codes_destroy(huffman_codes_t *c);

#ifdef __cplusplus
}
#endif

#endif /* _DEFLATE_SIMPLE_H_ */
