#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>

#include "seqIOABI.h"
#include "mach-io.h"

int dump_labels2(FILE *fp, off_t indexO) {
    off_t entryNum = -1;
    uint_4 entryLabel, entryLw1, entryLw2, entryLw3;
    uint_4 entryLw4, entryLw5, entryLw6;

    do {
	entryNum++;

	if (fseek(fp, indexO+(entryNum*IndexEntryLength), 0) != 0)
	    return 0;

	if (!be_read_int_4(fp, &entryLabel))
	    return 0;

	if (!be_read_int_4(fp, &entryLw1))
	    return 0;

	if (!be_read_int_4(fp, &entryLw2))
	    return 0;

	if (!be_read_int_4(fp, &entryLw3))
	    return 0;

	if (!be_read_int_4(fp, &entryLw4))
	    return 0;

	if (!be_read_int_4(fp, &entryLw5))
	    return 0;

	if (!be_read_int_4(fp, &entryLw6))
	    return 0;

	if (entryLabel) {
	    unsigned char c1, c2, c3, c4;

	    c1 = (entryLabel >> 24) & 0xff;
	    c2 = (entryLabel >> 16) & 0xff;
	    c3 = (entryLabel >>  8) & 0xff;
	    c4 = (entryLabel >>  0) & 0xff;

	    if (!isprint(c1))
		break;

	    printf("%c%c%c%c %3d %08x %08x %08x %08x %08x\n",
		   c1, c2, c3, c4, entryLw1,
		   entryLw2, entryLw3, entryLw4, entryLw5, entryLw6);
	}
    } while (entryLabel);

    return 0;
}

int dumpABIraw(FILE *fp, off_t indexO, uint_4 label, uint_4 count) {
    uint_4 off;
    uint_4 len;
    uint_4 i;
    unsigned char data;

    if (off = getABIIndexEntryLW(fp, indexO, label, count, 4, &len)) {
	if (!len)
	    return 0;

	/* Determine offset */
	if (len <= 4)
	    off += 20;
	else
	    getABIIndexEntryLW(fp, indexO, label, count, 5, &off);

	getABIIndexEntryLW(fp, indexO, label, count, 4, &len);

	/* Read length byte */
	fseek(fp, off, 0);

	/* Read data */
	for (i = 0; i < len; i++) {
	    fread(&data, 1, 1, fp);
	    printf("%c", data);
	}

	return 0;
    } else
	return -1;
}

int main(int argc, char **argv) {
    FILE *fp;
    uint_4 indexO;
    uint_4 label;
    int count = 1;

    if (argc < 2 || argc > 4) {
	fprintf(stderr, "Usage: %s [ident [count]] file\n", argv[0]);
	return 1;
    }
    
    if (NULL == (fp = fopen(argv[argc-1], "rb"))) {
	perror(argv[argc-1]);
	return 2;
    }

    if (-1 == getABIIndexOffset(fp, &indexO)) {
	fprintf(stderr, "Couldn't find ABI file index\n");
	return 3;
    }

    if (argc == 2) {
	return dump_labels2(fp, (off_t)indexO);
    }

    label = ((((argv[1][0]<<8) + argv[1][1])<<8) + argv[1][2]<<8) + argv[1][3];

    if (argc == 4)
	count = atoi(argv[2]);

    if (-1 == dumpABIraw(fp, (off_t)indexO, label, count)) {
	fprintf(stderr, "Couldn't find that identifier\n");
	return 4;
    }

    return 0;
}
