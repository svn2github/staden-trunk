#include <staden_config.h>

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "abiIO.h"
#include <io_lib/mach-io.h> /* get_be_float, get_be_double */

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
int dumpABIBlock(abi_index_t *ind, int style, int separator, int date_sep,
		 int line_sep) {
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
	putchar(line_sep);

    return 0;
}

int dump_abi(char *fn, abi_t *abi, int style, int separator, int date_sep,
	     int line_sep, int tagged) {
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

	dumpABIBlock(ind, style, separator, date_sep, line_sep);
    }

    if (line_sep != '\n')
	puts("");

    return 0;
}

int doit(char *fn, int argc, char **argv,
	 int style, int separator, int date_sep, int line_sep,
	 int tagged, int all) {
    abi_t *abi;
    int retcode = 0;
    int arg = 0;

    if (NULL == (abi = read_abi(fn))) {
	perror(fn);
	return -1;
    }

    if (all) {
	retcode = dump_abi(fn, abi, style, separator, date_sep, line_sep,
			   tagged);
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
	    
	    if ((arg == 1 || line_sep == '\n') && tagged)
		printf("%s %.4s %d ", fn, lstr, count);
	    else if (tagged)
		printf("%.4s %d ", lstr, count);

	    if (NULL == (ind = find_abi_index(abi, label, count))) {
		retcode = -1;
		continue;
	    }

	    if (-1 == dumpABIBlock(ind, style, separator, date_sep, line_sep))
		retcode = -1;
	}
    }

    if (line_sep != '\n')
	puts("");

    return retcode;
}

void usage(void) {
    fprintf(stderr, "Usage: getABIfield [optons] [--] file [NAME [COUNT]] ... \n");
    fprintf(stderr, "    -l             set the output field separator to newline\n");
    fprintf(stderr, "    -h             hex mode\n");
    fprintf(stderr, "    -r             raw mode\n");
    fprintf(stderr, "    -t             tagged format output\n");
    fprintf(stderr, "    -a             dump all blocks\n");
    fprintf(stderr, "    -f x           reformat (1=int1, 4=int2, 5=int4, 7=float, 8=double,\n");
    fprintf(stderr, "                             10=date, 11=time, 18=P-string, 19=C-string)\n");
    fprintf(stderr, "    -I fofn        file of filenames, use \"-I -\" to read fofn from stdin\n");
    fprintf(stderr, "    -F sep         separator between elements in a block\n");
    fprintf(stderr, "    -L sep         separator between blocks listed in a file\n");
    fprintf(stderr, "    -D sep         separator within date and time format items\n");
    fprintf(stderr, "    --             force end of command line options\n");
    fprintf(stderr, "    file           ABI filename to read. Use \"-\" for stdin\n");
    exit(1);
}

int main(int argc, char **argv) {
    int style = 0, separator = ' ', date_sep = '\0', line_sep = '\n';
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

	} else if (strcmp(argv[arg], "-L") == 0) {
	    if (++arg == argc)
		usage();
	    line_sep=*argv[arg];

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
		    style, separator, date_sep, line_sep, tagged, all);
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
			style, separator, date_sep, line_sep, tagged, all);
    }

    if (fofn_fp != stdin) {
	fclose(fofn_fp);
    }

    return retcode;
}
