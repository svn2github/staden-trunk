#ifndef _SRF_H_
#define _SRF_H_

#include <inttypes.h>
#include "mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SRF_MAGIC		"SSRF"
#define SRF_VERSION             "1.1"

#define SRFB_CONTAINER 		'S'
#define SRFB_TRACE_HEADER	'H'
#define SRFB_TRACE_BODY		'R'
#define SRFB_INDEX		'I'

/*--- Public structures */

/* Container header - several per file */
typedef struct {
    uint32_t next_block_offset;
    char block_type;
    char version[256];
    char container_type;
    char base_caller[256];
    char base_caller_version[256];
} srf_cont_hdr_t;

/* Trace header - several per container */
typedef struct {
    char block_type;
    char id_prefix[256];
    uint32_t trace_hdr_size;
    unsigned char *trace_hdr;
} srf_trace_hdr_t;

/* Trace body - several per trace header */
typedef struct {
    uint32_t next_block_offset;
    char block_type;
    char read_id[256];
    unsigned char flags;
    uint32_t trace_size;
    unsigned char *trace;
} srf_trace_body_t;

#define SRF_READ_FLAG_BAD (1<<0)

/* Index - one per file */
typedef struct {
    uint32_t next_block_offset;
    char block_type;
} srf_index_t;

/* Master SRF object */
typedef struct {
    FILE *fp;
    srf_cont_hdr_t    ch; /* Cached copies of the most recent container, */
    srf_trace_hdr_t   th; /*   trace header */
    srf_trace_body_t  tb; /*   and trace body blocks */
} srf_t;

/* Indexing */
typedef struct {
    char     magic[4];
    char     version[4];
    uint64_t size;
    uint32_t n_container;
    uint32_t n_data_block_hdr;
    uint32_t n_buckets;
    uint8_t  hash_func;
} srf_index_hdr_t;

#define SRF_INDEX_HDR_SIZE 29 /* size of the above struct when on disk */
#define SRF_INDEX_MAGIC    "Ihsh"
#define SRF_INDEX_VERSION  "1.00"


/*--- Initialisation */
srf_t *srf_create(FILE *fp);
srf_t *srf_open(char *fn, char *mode);
void srf_destroy(srf_t *srf, int auto_close);

/*--- Base type I/O methods */

int srf_write_pstring(srf_t *srf, char *str);
int srf_read_pstring(srf_t *srf, char *str);

int srf_read_uint32(srf_t *srf, uint32_t *val);
int srf_write_uint32(srf_t *srf, uint32_t val);

int srf_read_uint64(srf_t *srf, uint64_t *val);
int srf_write_uint64(srf_t *srf, uint64_t val);

/*--- Mid level I/O - srf block */
srf_cont_hdr_t *srf_construct_cont_hdr(srf_cont_hdr_t *ch,
				       char *bc,
				       char *bc_version);
void srf_destroy_cont_hdr(srf_cont_hdr_t *ch);
int srf_read_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch);
int srf_write_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch);

srf_trace_hdr_t *srf_construct_trace_hdr(srf_trace_hdr_t *th,
					 char *prefix,
					 unsigned char *header,
					 uint32_t header_sz);
void srf_destroy_trace_hdr(srf_trace_hdr_t *th);
int srf_read_trace_hdr(srf_t *srf, srf_trace_hdr_t *th);
int srf_write_trace_hdr(srf_t *srf, srf_trace_hdr_t *th);

srf_trace_body_t *srf_construct_trace_body(srf_trace_body_t *th,
					   char *suffix,
					   unsigned char *body,
					   uint32_t body_size);
void srf_destroy_trace_body(srf_trace_body_t *th);
int srf_write_trace_body(srf_t *srf, srf_trace_body_t *th);
int srf_read_trace_body(srf_t *srf, srf_trace_body_t *th, int no_trace);

int srf_read_index_hdr(srf_t *srf, srf_index_hdr_t *hdr);
int srf_write_index_hdr(srf_t *srf, srf_index_hdr_t *hdr);

/*--- Higher level I/O functions */
mFILE *srf_next_trace(srf_t *srf, char *name);

int srf_next_block_type(srf_t *srf); /* peek ahead */
int srf_next_block_details(srf_t *srf, uint64_t *pos, char *name);

int srf_find_trace(srf_t *srf, char *trace,
		   uint64_t *cpos, uint64_t *hpos, uint64_t *dpos);

#ifdef __cplusplus
}
#endif

#endif /* _SRF_H_ */
