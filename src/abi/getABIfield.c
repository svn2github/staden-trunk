#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "seqIOABI.h"
#include "mach-io.h"
#include "fpoint.h"

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


abi_t *new_abi_t(char *data, size_t dsize) {
    abi_t *abi;
    if (NULL == (abi = (abi_t *)calloc(1, sizeof(abi_t))))
	return NULL;

    abi->data = data;
    abi->data_sz = dsize;

    return abi;
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

void dump_raw(char *data, size_t sz) {
    fwrite(data, 1, sz, stdout);
}

void dump_int1(char *data, size_t sz, int hex, int separator) {
    size_t i;
    for (i = 0; i < sz;) {
	if (hex)
	    printf("%02x", data[i] & 0xff);
	else
	    printf("%d", (signed char)(data[i] & 0xff));
	if (++i < sz)
	    putchar(separator);
    }
}

void dump_int2(char *data, size_t sz, int hex, int separator) {
    size_t i;

    sz /= 2;
    for (i = 0; i < sz;) {
	int_2 i2 = get_be_int2(&data[i*2]);
	if (hex)
	    printf("%04x", i2 &0xffff);
	else
	    printf("%d", i2);
	if (++i < sz)
	    putchar(separator);
    }
}

void dump_int4(char *data, size_t sz, int hex, int separator) {
    size_t i;

    sz /= 4;
    for (i = 0; i < sz;) {
	int_4 i4 = get_be_int4(&data[i*4]);
	printf(hex ? "%08x" : "%d", i4);
	if (++i < sz)
	    putchar(separator);
    }
}

void dump_pstr(char *data, size_t sz, int hex, int separator) {
    if (hex)
	dump_int1(data, sz, 1, separator);
    else
	printf("%.*s", *data, data+1);
}

void dump_cstr(char *data, size_t sz, int hex, int separator) {
    if (hex)
	dump_int1(data, sz, 1, separator);
    else
	printf("%.*s", (int)sz, data);
}

void dump_float(char *data, size_t sz, int hex, int separator) {
    size_t i;

    if (hex) {
	dump_int4(data, sz, hex, separator);
	return;
    }

    for (i = 0; i < sz; i+=4) {
	printf("%f", get_be_float(&data[i]));
	if (i+4 < sz)
	    putchar(separator);
    }
}

void dump_double(char *data, size_t sz, int hex, int separator) {
    size_t i;

    if (hex) {
	dump_int4(data, sz, hex, separator);
	return;
    }

    for (i = 0; i < sz; i+=8) {
	printf("%f", get_be_double(&data[i]));
	if (i+8 < sz)
	    putchar(separator);
    }
}

void dump_date(char *data, size_t sz, int hex, int separator, int date_sep) {
    size_t i;
    unsigned char *p = (unsigned char *)data;

    if (hex) {
	dump_int4(data, sz, hex, separator);
	return;
    }

    for (i = 0; i < sz; i+=4) {
	printf("%04d%c%02d%c%02d",
	       p[i]*256 + p[i+1], date_sep ? date_sep : '/',
	       p[i+2], date_sep ? date_sep : '/',
	       p[i+3]);
	if (i+4 < sz)
	    putchar(separator);
    }
}

void dump_time(char *data, size_t sz, int hex, int separator, int date_sep) {
    size_t i;
    unsigned char *p = (unsigned char *)data;

    if (hex) {
	dump_int4(data, sz, hex, separator);
	return;
    }

    for (i = 0; i < sz; i+=4) {
	printf("%02d%c%02d%c%02d%c%02d",
	       p[i], date_sep ? date_sep : ':',
	       p[i+1], date_sep ? date_sep : ':',
	       p[i+2], date_sep ? date_sep : '.',
	       p[i+3]);
	if (i+4 < sz)
	    putchar(separator);
    }
}

/*
 * Style: 0 = auto
 *       -1 = hex
 *       -2 = raw
 *       +? = that format number
 */
int dumpABIBlock(abi_index_t *ind, int style, int separator, int date_sep) {
    if (style == -2) {
	dump_raw(ind->data, ind->size);
    } else {
	int hex = (style == -1);

	/*
	 * Guess some style. PBAS, APrX, MODL etc are type 2 (1-byte char),
	 * but are really strings. Conversly PCON is also type 2 and is a
	 * series of numbers.
	 */
	if (style == 0 && ind->format == 2) {
	    int i;
	    for (i = 0; i < ind->size; i++)
		if (!(isprint(ind->data[i]) || isspace(ind->data[i])))
		    break;
	    if (i == ind->size)
		style = 19; /* C-string */
	}


	if (style <= 0)
	    style = ind->format;

	switch(style) {
	case 2:
	    dump_int1(ind->data, ind->size, hex, separator);
	    break;
	    
	case 4:
	    dump_int2(ind->data, ind->size, hex, separator);
	    break;
	    
	case 5:
	    dump_int4(ind->data, ind->size, hex, separator);
	    break;
	    
	case 7:
	    dump_float(ind->data, ind->size, hex, separator);
	    break;
	    
	case 8:
	    dump_double(ind->data, ind->size, hex, separator);
	    break;
	    
	case 10:
	    dump_date(ind->data, ind->size, hex, separator, date_sep);
	    break;
	    
	case 11:
	    dump_time(ind->data, ind->size, hex, separator, date_sep);
	    break;
	    
	case 18:
	    dump_pstr(ind->data, ind->size, hex, separator);
	    break;
	    
	case 19:
	    dump_cstr(ind->data, ind->size, hex, separator);
	    break;
	    
	default:
	    /* Default to hex dump for unknown ones */
	    dump_int1(ind->data, ind->size, 1, separator);
	    break;
	}
    }

    if (separator != '\n')
	puts("");

    return 0;
}

int dump_abi(char *fn, abi_t *abi, int style, int separator, int date_sep,
	     int tagged) {
    int i;
    for (i = 0; i < abi->header.nindex; i++) {
	abi_index_t *ind = &abi->index[i];
	if (tagged) {
	    printf("%s %c%c%c%c %d ",
		   fn,
		   (ind->id >> 24) & 0xff,
		   (ind->id >> 16) & 0xff,
		   (ind->id >>  8) & 0xff,
		   (ind->id >>  0) & 0xff,
		   ind->idv);
	}

	dumpABIBlock(ind, style, separator, date_sep);
    }

    return 0;
}

int doit(char *fn, int argc, char **argv,
	 int style, int separator, int date_sep, int tagged, int all) {
    abi_t *abi;
    int retcode = 0;
    int arg = 0;

    if (NULL == (abi = read_abi(fn))) {
	perror(fn);
	return -1;
    }

    if (all) {
	retcode = dump_abi(fn, abi, style, separator, date_sep, tagged);
    } else if (arg == argc) {
	retcode = dump_index(abi);
    } else {
	for (; arg < argc; arg++) {
	    uint_4 label;
	    int count;
	    char *lstr = argv[arg];
	    abi_index_t *ind;
	    
	    label = ((((lstr[0]<<8) + lstr[1])<<8) +
		     lstr[2]<<8) + lstr[3];
	    
	    count = (++arg < argc) ? atoi(argv[arg]) : 1;
	    
	    if (tagged)
		printf("%s %.4s %d ", fn, lstr, count);

	    if (NULL == (ind = find_abi_index(abi, label, count))) {
		retcode = -1;
		continue;
	    }

	    if (-1 == dumpABIBlock(ind, style, separator, date_sep))
		retcode = -1;
	}
    }

    return retcode;
}

void usage(void) {
    fprintf(stderr, "Usage: getABIfield [optons] [--] file [NAME [COUNT]] ... \n");
    fprintf(stderr, "    -n x           display name 'x'. List index is not set.\n");
    fprintf(stderr, "    -l             set the output field separator to newline\n");
    fprintf(stderr, "    -h             hex mode\n");
    fprintf(stderr, "    -r             raw mode\n");
    fprintf(stderr, "    -t             tagged format output\n");
    fprintf(stderr, "    -a             dump all blocks\n");
    fprintf(stderr, "    -f x           reformat (1=int1, 4=int2, 5=int4, 7=float, 8=double,\n");
    fprintf(stderr, "                             10=date, 11=time, 18=P-string, 19=C-string)\n");
    fprintf(stderr, "    -I fofn        file of filenames, use \"-I -\" to read fofn from stdin\n");
    fprintf(stderr, "    -F sep         separator between elements in a block\n");
    fprintf(stderr, "    -D sep         separator within date and time format items\n");
    fprintf(stderr, "    --             force end of command line options\n");
    fprintf(stderr, "    file           ABI filename to read. Use \"-\" for stdin\n");
    exit(1);
}

int main(int argc, char **argv) {
    int style = 0, separator = ' ', date_sep = '\0';
    int arg, retcode = 0, tagged = 0, all = 0;
    char *fofn = NULL, file[1024];
    FILE *fofn_fp;

    /* Check args */
    for (arg = 1; arg < argc; arg++) {
	if (strcmp(argv[arg], "-f") == 0) {
	    if (++arg == argc)
		usage();
	    style = atoi(argv[arg]);

	} else if (strcmp(argv[arg], "-h") == 0) {
	    style=-1;

	} else if (strcmp(argv[arg], "-r") == 0) {
	    style=-2;

	} else if (strcmp(argv[arg], "-a") == 0) {
	    all = 1;

	} else if (strcmp(argv[arg], "-F") == 0) {
	    if (++arg == argc)
		usage();
	    separator=*argv[arg];

	} else if (strcmp(argv[arg], "-D") == 0) {
	    if (++arg == argc)
		usage();
	    date_sep=*argv[arg];

	} else if (strcmp(argv[arg], "-l") == 0) {
	    separator='\n';

	} else if (strcmp(argv[arg], "-t") == 0) {
	    tagged = 1;

	} else if (strcmp(argv[arg], "-I") == 0) {
	    if (++arg == argc)
		usage();
	    fofn=argv[arg];

	} else if (strcmp(argv[arg], "--") == 0) {
	    arg++;
	    break;

	} else if (argv[arg][0] == '-' && argv[arg][1] != '\0') {
	    usage();

	} else {
	    break;
	}
    }

    /* Single file */
    if (!fofn) {
	if (argc-arg < 1)
	    usage();
	return doit(argv[arg], argc-arg-1, &argv[arg+1],
		    style, separator, date_sep, tagged, all);
    }

    /* File of filenames */
    if (strcmp(fofn, "-") == 0) {
	fofn_fp = stdin;
    } else {
	if (NULL == (fofn_fp = fopen(fofn, "r"))) {
	    perror(fofn);
	    return 1;
	}
    }

    while (fgets(file, 1024, fofn_fp)) {
	char *cp;
	if (cp = strrchr(file, '\n'))
	    *cp = 0;

	retcode |= doit(file, argc-arg, &argv[arg],
			style, separator, date_sep, tagged, all);
    }

    if (fofn_fp != stdin) {
	fclose(fofn_fp);
    }

    return retcode;
}
