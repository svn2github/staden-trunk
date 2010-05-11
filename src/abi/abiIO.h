#ifndef _ABI_H_
#define _ABI_H_

#include <os.h>
#include <io_lib/os.h>
#include <io_lib/fpoint.h>

/*
 * The ABI magic number - "ABIF"
 */
#define ABI_MAGIC	((int_4) ((((('A'<<8)+'B')<<8)+'I')<<8)+'F')

#define ABI_LABEL(a) ((int_4) ((((((a)[0]<<8)+(a)[1])<<8)+(a)[2])<<8)+(a)[3])

/*
 * ---------------------------------------------------------------------------
 * Fetches values from character array 'd' which is stored in big-endian 
 * format and not necessarily word aligned.
 *
 * The returned values are in machine native endianess.
 */

/* signed 2-byte shorts */
#define get_be_int2(d) \
    ((((unsigned char *)d)[0]<<8) + \
     (((unsigned char *)d)[1]))

/* signed 4-byte ints */
#define get_be_int4(d) \
    ((((unsigned char *)d)[0]<<24) + \
     (((unsigned char *)d)[1]<<16) + \
     (((unsigned char *)d)[2]<< 8) + \
     (((unsigned char *)d)[3]))

/* 4-byte float */
#define get_be_float(d) (int_to_float(get_be_int4(d)))

typedef struct {
    uint_4 magic;        /* "ABIF" - magic number */
    char stuff1[12];      /* don't care what this is */
    int  ilen;            /* Length of an index entry */
    int  nindex;          /* Number of index entries */
    int  isize;		  /* Total size of the index directory */
    int  ioff;		  /* Location of the directory */
    char spare[98];       /* Unused? */
} abi_header_t;

typedef struct {
    uint_4 id;       /* DATA */
    uint_4 idv;      /* eg 8,9,10,11 */
    uint_2 format;   /* 4=short, etc */
    uint_2 isize;    /* size of format, eg 2 */
    uint_4 icount;   /* count of items */
    uint_4 size;     /* icount * isize */
    uint_4 offset;   /* file offset (or data if size <= 4) */
    uint_4 junk_ptr; /* junk_ptr in ABI file */
    char  *data;     /* in memory copy */
} abi_index_t;

typedef struct {
    abi_header_t header;  /* File header */
    abi_index_t *index;   /* Pointer to the index */
    char        *data;    /* in memory copy of the entire ABI file */
    size_t       data_sz; /* size of 'data' */
} abi_t;

double get_be_double(char *d);
abi_t *new_abi_t(char *data, size_t dsize);
void del_abi_t(abi_t *abi);
int decode_abi_header(abi_t *abi);
int decode_abi_index(abi_t *abi);
abi_t *decode_abi_file(char *data, size_t data_sz);
abi_t *read_abi(char *fn);
abi_index_t *find_abi_index(abi_t *abi, uint_4 id, uint_4 idv);
int dump_index(abi_t *abi);

#endif /* _ABI_H_ */
