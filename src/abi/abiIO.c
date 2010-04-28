#include <staden_config.h>

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "abiIO.h"

/* 8-byte doubles */
#ifdef SP_BIG_ENDIAN
double get_be_double(char *d) {
    int i;
    union {
	double d;
	char c[8];
    } cvt;

    for (i = 0; i < 8; i++) {
	cvt.c[7-i] = d[i];
    }

    return cvt.d;
}
#else
double get_be_double(char *d) {
    int i;
    union {
	double d;
	char c[8];
    } cvt;

    for (i = 0; i < 8; i++) {
	cvt.c[i] = d[i];
    }

    return cvt.d;
}
#endif

/*
 * Data will have been malloced by the caller. Passing it here gives ownership
 * of data to the abi_t structure such that the caller should then not free
 * or otherwise modify data.
 */
abi_t *new_abi_t(char *data, size_t dsize) {
    abi_t *abi;
    if (NULL == (abi = (abi_t *)calloc(1, sizeof(abi_t))))
	return NULL;

    abi->data = data;
    abi->data_sz = dsize;

    return abi;
}

void del_abi_t(abi_t *abi) {
    if (!abi)
	return;

    if (abi->index)
	free(abi->index);

    if (abi->data)
	free(abi->data);

    free(abi);
}

int decode_abi_header(abi_t *abi) {
    uint_2 i2;
    uint_4 i4;

    if (!abi->data || abi->data_sz < 128)
	return -1;

    memcpy(&i4, &abi->data[0x00], 4);
    abi->header.magic = be_int4(i4);
    if (abi->header.magic != ABI_MAGIC)
	return -1;

    memcpy(abi->header.stuff1, &abi->data[0x04], 12);
    memcpy(abi->header.spare,  &abi->data[0x1e], 98);

    memcpy(&i2, &abi->data[0x10], 2);
    abi->header.ilen = be_int2(i2);

    memcpy(&i2, &abi->data[0x14], 2);
    abi->header.nindex = be_int2(i2);

    memcpy(&i4, &abi->data[0x16], 4);
    abi->header.isize = be_int4(i4);

    memcpy(&i4, &abi->data[0x1a], 4);
    abi->header.ioff = be_int4(i4);

    /*
    printf("ABI: ilen=%d, nindex=%d, isize=%d, ioff=%d\n",
	   abi->header.ilen,
	   abi->header.nindex,
	   abi->header.isize,
	   abi->header.ioff);
    */
    return 0;
}

int decode_abi_index(abi_t *abi) {
    int i;
    char *p;

    if (abi->header.ioff + abi->header.ilen * abi->header.nindex >
	abi->data_sz)
	return -1;

    abi->index = (abi_index_t *)calloc(abi->header.nindex,
				       sizeof(abi_index_t));
    for (i = 0; i < abi->header.nindex; i++) {
	p = &abi->data[abi->header.ioff + abi->header.ilen * i];
	abi->index[i].id       = get_be_int4(&p[0]);
	abi->index[i].idv      = get_be_int4(&p[4]);
	abi->index[i].format   = get_be_int2(&p[8]);
	abi->index[i].isize    = get_be_int2(&p[10]);
	abi->index[i].icount   = get_be_int4(&p[12]);
	abi->index[i].size     = get_be_int4(&p[16]);
	abi->index[i].offset   = get_be_int4(&p[20]);
	abi->index[i].junk_ptr = get_be_int4(&p[24]);

	if (abi->index[i].size <= 4)
	    abi->index[i].data = &p[20];
	else
	    abi->index[i].data = &abi->data[abi->index[i].offset];
    }

    /*
    for (i = 0; i < abi->header.nindex; i++) {
	char c1, c2, c3, c4;
	c1 = (abi->index[i].id >> 24) & 0xff;
	c2 = (abi->index[i].id >> 16) & 0xff;
	c3 = (abi->index[i].id >>  8) & 0xff;
	c4 = (abi->index[i].id >>  0) & 0xff;
	
	printf("%d: %c%c%c%c %d %d %d %d %d %d %d\n",
	       i,
	       c1,c2,c3,c4,
	       abi->index[i].idv,
	       abi->index[i].format,
	       abi->index[i].isize,
	       abi->index[i].icount,
	       abi->index[i].size,
	       abi->index[i].offset,
	       abi->index[i].junk_ptr);
    }
    */

    return 0;
}

abi_t *decode_abi_file(char *data, size_t data_sz) {
    abi_t *abi = new_abi_t(data, data_sz);
    decode_abi_header(abi);
    decode_abi_index(abi);

    return abi;
}

abi_t *read_abi(char *fn) {
    FILE *fp;
    size_t sz = 0, len;
    char *data = NULL;
    char block[8192];
    
    if (strcmp(fn, "-") == 0) {
	fp = stdin;
    } else {
	fp = fopen(fn, "rb");
	if (!fp)
	    return NULL;
    }

    while (len = fread(block, 1, 8192, fp)) {
	data = (char *)realloc(data, sz+len);
	memcpy(&data[sz], block, len);
	sz += len;
    }

    if (fp != stdin)
	fclose(fp);
    
    return decode_abi_file(data, sz);
}

abi_index_t *find_abi_index(abi_t *abi, uint_4 id, uint_4 idv) {
    int i;
    for (i = 0; i < abi->header.nindex; i++) {
	if (abi->index[i].id == id && abi->index[i].idv == idv)
	    return &abi->index[i];
    }

    return NULL;
}


int dump_index(abi_t *abi) {
    int i;
    unsigned char c1, c2, c3, c4;

    for (i = 0; i < abi->header.nindex; i++) {
	abi_index_t *ind = &abi->index[i];

	c1 = (ind->id >> 24) & 0xff;
	c2 = (ind->id >> 16) & 0xff;
	c3 = (ind->id >>  8) & 0xff;
	c4 = (ind->id >>  0) & 0xff;

	printf("%c%c%c%c %3d %04x %04x %08x %08x %08x %08x\n",
	       c1, c2, c3, c4,
	       ind->idv,
	       ind->format,
	       ind->isize,
	       ind->icount,
	       ind->size,
	       ind->offset,
	       ind->junk_ptr);
    }

    return 0;
}
