#ifndef MAQMAP_H_
#define MAQMAP_H_

//#define MAX_READLEN 128
//#define MAX_READLEN 64


#define MAX_NAMELEN 36

#define MAQMAP_FORMAT_OLD 0
#define MAQMAP_FORMAT_NEW -1

#define PAIRFLAG_FF      0x01
#define PAIRFLAG_FR      0x02
#define PAIRFLAG_RF      0x04
#define PAIRFLAG_RR      0x08
#define PAIRFLAG_PAIRED  0x10
#define PAIRFLAG_DIFFCHR 0x20
#define PAIRFLAG_NOMATCH 0x40

#include <string.h>
#include <zlib.h>

typedef unsigned char bit8_t;
typedef unsigned bit32_t;
typedef unsigned long long bit64_t;

/*
  name: read name
  size: the length of the read
  seq: read sequence (see also below)
  seq[MAX_READLEN-1]: single end mapping quality (equals to map_qual if not paired)
  map_qual: the final mapping quality
  alt_qual: the lower quality of the two ends (equals to map_qual if not paired)
  flag: status of the pair
  dist: offset of the mate (zero if not paired)
  i1: mismatch and rough sum of errors of the best hit
  i2: mismatch and rough sum of errors of the second best hit
  c[2]: count of all 0- and 1-mismatch hits on the reference
 */
typedef struct
{
	bit8_t seq[64]; /* the last base is the single-end mapping quality. */
	bit8_t size, map_qual, i1, i2, c[2], flag, alt_qual;
	bit32_t seqid, pos;
	int dist;
	char name[MAX_NAMELEN];
} maqmap64_t;

typedef struct
{
	bit8_t seq[128]; /* the last base is the single-end mapping quality. */
	bit8_t size, map_qual, i1, i2, c[2], flag, alt_qual;
	bit32_t seqid, pos;
	int dist;
	char name[MAX_NAMELEN];
} maqmap128_t;

typedef struct
{
	int format, n_ref;
	char **ref_name;
	bit64_t n_mapped_reads;
} maqmap_t;

#ifdef __cplusplus
extern "C" {
#endif
	maqmap_t *maq_new_maqmap(void);
	void maq_delete_maqmap(maqmap_t *mm);
	void maqmap_write_header(gzFile fp, const maqmap_t *mm);
	maqmap_t *maqmap_read_header(gzFile fp);
        int maq_detect_size(gzFile fp);

#ifdef __cplusplus
}
#endif

#endif
