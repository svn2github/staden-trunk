#ifndef _HUFFMAN_STATIC_H_
#define _HUFFMAN_STATIC_H_

#ifdef __cplusplus
extern "C" {
#endif

#define CODE_UNSET	0
#define CODE_DNA	1
#define CODE_TRACES	2
#define CODE_CONF	3
#define CODE_CONF_RLE	4


unsigned char *encode_memory(unsigned char *data, int len, int *nbytes,
			     int code_set);

unsigned char *decode_memory(unsigned char *data, int len, int *nbytes);

#ifdef __cplusplus
}
#endif

#endif /* _HUFFMAN_STATIC_H_ */
