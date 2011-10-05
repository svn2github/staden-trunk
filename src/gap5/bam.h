#ifndef _BAM_H_
#define _BAM_H_

#include <inttypes.h>
#include <zlib.h>
#include <io_lib/hash_table.h>

/* BAM header structs */
typedef struct tag_list {
    char *value; /* NULL => end of tags */
    int key;
    int length;  
} tag_list_t;

typedef struct {
    char *name;
    uint32_t len;
} bam_ref_t;


/* The main bam sequence struct */
typedef struct {
    uint32_t alloc; /* total size of this struct + 'data' onwards */
    uint32_t blk_size;

    /* Unpacked copy of cigar_len */
    uint32_t cigar_len;

    /* The raw bam block follows, in same order as on the disk */
    int32_t  ref;
    int32_t  pos;
    uint32_t bin_mq_nl;
    //uint32_t bin:16, map_qual:8, name_len:8;
    uint32_t flag_nc;
    //uint32_t flag:16, cigar_len:16;
    int32_t  len;
    int32_t  mate_ref;
    int32_t  mate_pos;
    int32_t  ins_size;

    /* Followed by arbitrary bytes of packed data */
    unsigned char data; /* unknown size */
} bam_seq_t;


/* Auxillary field handling */
typedef union {
    char  *s;
    int    i;
    float  f;
    double d;
} bam_aux_t;


/*
 * Our bam stream consists of a zlib gzFile stream and a buffer for it to
 * output to. This allows us to call small bam_read requests while
 * translating these to fewer, larger, gzread calls. The overhead is
 * therefore minimal.
 */
#define Z_BUFF_SIZE 65536 /* Max size of a BGZF block */
typedef struct {
    int fd;
    z_stream s;

    unsigned char in[Z_BUFF_SIZE];
    unsigned char *in_p;
    size_t in_sz;

    unsigned char out[Z_BUFF_SIZE];
    unsigned char *out_p;
    size_t out_sz;

    /* BAM specifics */
    int32_t next_len;

    uint32_t header_len;
    char *header;

    HashTable *ref_hash;  /* keyed on SN, value = numeric id */
    HashTable *rg_hash;   /* keyed on ID, value = taglist ptr */
    uint32_t nref;
    bam_ref_t *ref;

    /* Cached bam_seq_t, to avoid excessive mallocs */
    bam_seq_t *bs;
    int bs_size;

    /* Boolean to indicate if we've finished the most recent z stream */
    int z_finish;

    /* Indicates whether gzipped, and if so with bgzf extra fields */
    int gzip;
    int bgzf; 

    /* Whether BAM or SAM format */
    int bam;

    /* If true, skip auxillary field parsing while reading SAM */
    int no_aux;

    /* line number (when in SAM mode) */
    int line;
} bam_file_t;

/* Decoding the above struct */
#define bam_map_qual(b)  (((b)->bin_mq_nl >> 8) & 0xff)
#define bam_bin(b)       ((b)->bin_mq_nl >> 16)
#define bam_name_len(b)  (((b)->bin_mq_nl) & 0xff)
#define bam_set_map_qual(b,v) \
    ((b)->bin_mq_nl = ((b)->bin_mq_nl & 0xffff00ff) | (((v) & 0xff)<<8))
#define bam_set_bin(b,v) \
    ((b)->bin_mq_nl = ((b)->bin_mq_nl & 0x0000ffff) | (((v) & 0xffff)<<16))
#define bam_set_name_len(b,v) \
    ((b)->bin_mq_nl = ((b)->bin_mq_nl & 0xffffff00) | ((v) & 0xff))
//#define bam_map_qual(b)  ((b)->map_qual)
//#define bam_bin(b)       ((b)->bin)
//#define bam_name_len(b)  ((b)->name_len)

#define bam_flag(b)      ((b)->flag_nc >> 16)
#define bam_set_flag(b,v) \
    ((b)->flag_nc = ((b)->flag_nc & 0x0000ffff) | (((v) & 0xffff)<<16))

#if 0
#  define bam_cigar_len(b) ((b)->flag_nc & 0xffff)
#  define bam_set_cigar_len(b,v) \
    ((b)->flag_nc = ((b)->flag_nc & 0xffff0000) | ((v) & 0xffff))
#else
#  define bam_cigar_len(b)        ((b)->cigar_len)
#  define bam_set_cigar_len(b, v) ((b)->cigar_len = (v))
#endif

#define bam_strand(b)    ((bam_flag((b)) & BAM_FREVERSE) != 0)
#define bam_name(b)      ((char *)(&(b)->data))
#define bam_cigar(b)     ((uint32_t *)(bam_name((b)) + bam_name_len((b))))
#define bam_seq(b)       (((char *)bam_cigar((b))) + 4*bam_cigar_len(b))
#define bam_qual(b)      (bam_seq(b) + (int)(((b)->len+1)/2))
#define bam_aux(b)       (bam_qual(b) + (b)->len)

/* Rounds up to the next multiple of 4 */
#define round4(v) (((v-1)&~3)+4)

/* BAM flags */
#define BAM_FPAIRED        1
#define BAM_FPROPER_PAIR   2
#define BAM_FUNMAP         4
#define BAM_FMUNMAP        8
#define BAM_FREVERSE      16
#define BAM_FMREVERSE     32
#define BAM_FREAD1        64
#define BAM_FREAD2       128
#define BAM_FSECONDARY   256
#define BAM_FQCFAIL      512
#define BAM_FDUP        1024

/* CIGAR operations, taken from samtools bam.h */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

enum cigar_op {
    BAM_CMATCH=0,
    BAM_CINS=1,
    BAM_CDEL=2,
    BAM_CREF_SKIP=3,
    BAM_CSOFT_CLIP=4,
    BAM_CHARD_CLIP=5,
    BAM_CPAD=6,
    BAM_CBASE_MATCH=7,
    BAM_CBASE_MISMATCH=8
};

/* ----------------------------------------------------------------------
 * Function prototypes
 * We only support reading, so basically we have open, read, close along
 * with some utility functions for querying aux records.
 */
bam_file_t *bam_open(char *fn, char *mode);
void bam_close(bam_file_t *b);
int bam_next_seq(bam_file_t *b, bam_seq_t **bsp);
char *bam_aux_find(bam_seq_t *b, char *key);
tag_list_t *bam_find_rg(bam_file_t *b, char *id);
int bam_aux_iter(bam_seq_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val);

/* Taken from samtools/bam.h */
#define bam_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam_nt16_rev_table "=ACMGRSVTWYHKDBN"

#endif /* _BAM_H_ */
