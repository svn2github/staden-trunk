#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "hash_table.h"
#include "sff.h"
#include "os.h"
#include "mFILE.h"

int main(int argc, char **argv) {
    mFILE *mf;
    HashFile *hf;
    sff_common_header *ch;
    sff_read_header *rh;
    int i;
    char *sff;
    char hdr[31];
    uint64_t hsize;
    FILE *fp;
    
    if (argc != 2) {
	fprintf(stderr, "Usage: hash_sff sff_file\n");
	return 1;
    }

    /* open (and read) the entire sff file */
    sff = argv[1];
    if (NULL == (mf = mfopen(sff, "rb"))) {
	perror(argv[0]);
	return 1;
    }

    /* Read the common header */
    ch = read_sff_common_header(mf);

    /* Create the hash table */
    hf = HashFileCreate(0, HASH_DYNAMIC_SIZE);

    /* Add the SFF common header as a hash file-header */
    hf->headers = (HashFileSection *)malloc(sizeof(*hf->headers));
    hf->nheaders = 1;
    hf->headers[0].pos = 0;
    hf->headers[0].size = ch->header_len;
    hf->headers[0].cached_data = NULL;

    /* Read the index items, adding to the hash */
    for (i = 0; i < ch->nreads; i++) {
	int dlen;
	uint32_t offset;
	HashData hd;
	HashFileItem *hfi = (HashFileItem *)calloc(1, sizeof(*hfi));

	offset = mftell(mf);
	rh = read_sff_read_header(mf);
	dlen = (2*ch->flow_len + 3*rh->nbases + 7) & ~7;
	mfseek(mf, dlen, SEEK_CUR);
	
	hfi->header = hf->nheaders;
	hfi->footer = 0;
	hfi->pos = offset;
	hfi->size = mftell(mf) - offset;
	hd.p = hfi;

	HashTableAdd(hf->h, rh->name, rh->name_len, hd, NULL);
    }

    mfclose(mf);

    HashTableStats(hf->h, stderr);

    /* Open the sff file and add the header */
    fp = fopen(sff, "rb+");
    fread(hdr, 1, 31, fp);
    fseek(fp, 0, SEEK_END);
    //hsize = HashFileSave(hf, stdout, 0);
    hsize = HashFileSave(hf, fp, 0);
    HashFileDestroy(hf);
    *(uint64_t *)(hdr+8)  = be_int8(mftell(mf));
    *(uint32_t *)(hdr+16) = be_int4(hsize);
    fprintf(stderr, "hsize=%ld\n", (long)hsize);
    fseek(fp, 0, SEEK_SET);
    fwrite(hdr, 1, 31, fp);
    fclose(fp);

    return 0;
}
