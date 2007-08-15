#ifndef _HUFFMAN_STATIC_H_
#define _HUFFMAN_STATIC_H_

#ifdef __cplusplus
extern "C" {
#endif

/* inlined codes */
#define CODE_INLINE	  0   

/* predefined codes */
#define CODE_DNA	  1   /* auto-pick from following 4 */
#define CODE_DNA_A	  2   
#define CODE_DNA_C	  3   
#define CODE_DNA_G	  4   
#define CODE_DNA_T	  5   
#define CODE_ENGLISH      6

/* predefined elsewhere in HUFF chunks, 128 onwards */
#define CODE_USER	  128 

#ifndef ZTR_FORM_STHUFF
#  define ZTR_FORM_STHUFF  77
#endif

/* A single symbol and it's encoding */
typedef struct {
    signed int symbol; /* 0-255 character, 256 = exception code, 257 = EOF */
    int nbits;
    unsigned int code;
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
    decode_t d_tab[258]; /* indexed by symbol encoding bit-lengths */
    huffman_code_t lookup[258]; /* Mapping of symbol character to code */
} huffman_codes_t;

unsigned char *huffman_encode(huffman_codes_t *c, int code_set,
			      unsigned char *data, int len, int *nbytes);

unsigned char *huffman_decode(huffman_codes_t *c,
			      unsigned char *data, int len, int *nbytes);

huffman_codes_t *generate_code_set(int code_set, unsigned char *data, int len,
				   unsigned char *legal);
unsigned char *store_codes(huffman_codes_t *c, unsigned int *comp_len);
huffman_codes_t *restore_codes(unsigned char *data, int *len_used);
void huffman_codes_destroy(huffman_codes_t *c);

#ifdef __cplusplus
}
#endif

#endif /* _HUFFMAN_STATIC_H_ */
