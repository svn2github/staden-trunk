/*
 * File: abi.c
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created: 20/9/92 jkb - hacked from abi.c
 * Updated:
 *
 */


#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mach-io.h"
    
typedef struct {
    char id[12];
    int_2 s1;
    int_2 s2;
    int_2 s3;
    int_2 s4;
    int_2 recs;			/* number of record entries in index */
    int_4 l1;
    int_4 offset;		/* offset of index */
    int_4 l2;
} Header;

typedef struct {
    char id[4];
    int_4 index;
    int_2 type;
    int_2 size;
    int_4 N;
    int_4 length;
    int_4 ptr1;
    int_4 ptr2;
} Rec;

/*
 * Meaningful ids: (Explanations c/o jes)
 *
 * GMBF       Gel type
 * DATA.1-4   Raw data block
 * DATA.5-8   Central (rather featureless) data block
 * DATA.9-12  Processed data block
 * FWO_       Base order (see seqIOABI.c)
 * PBAS.1-2   Base sequence
 * PLOC.1-2   Base positions
 * S/N%       Signal strengths (array of 4 int2s)
 * SMPL       Sample name
 * SPAC       Base spacing (float)
 *
 * Explanations c/o lfw
 * LIMC      Array of four int2s
 *           { Base Call Start, Base Call Base 1, 0 (always?), Base Call End }
 *
 * It appears that if the value for an id occupies less than or equal four
 * bytes, it is shoe-horned into the ptr1 field. Otherwise, ptr1 holds
 * the byte offset in the file where the data can be found.
 * 
 *
 *
 *
 */


int read_header(FILE *fp, Header *h)
{
    fseek(fp,0,0);
    if (fread(&h->id[0],sizeof(h->id),1,fp)==0) return 1;
    if (be_read_int_2(fp,(uint_2 *)&h->s1)==0) return 1;
    if (be_read_int_2(fp,(uint_2 *)&h->s2)==0) return 1;
    if (be_read_int_2(fp,(uint_2 *)&h->s3)==0) return 1;
    if (be_read_int_2(fp,(uint_2 *)&h->s4)==0) return 1;
    if (be_read_int_2(fp,(uint_2 *)&h->recs)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&h->l1)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&h->offset)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&h->l2)==0) return 1;
    return 0;
}

int read_rec(FILE *fp, Rec *r)
{
    if (fread(&r->id[0],sizeof(r->id),1,fp)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&r->index)==0) return 1;
    if (be_read_int_2(fp,(uint_2 *)&r->type)==0) return 1;
    if (be_read_int_2(fp,(uint_2 *)&r->size)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&r->N)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&r->length)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&r->ptr1)==0) return 1;
    if (be_read_int_4(fp,(uint_4 *)&r->ptr2)==0) return 1;
    return 0;
}


int main(int argc, char *argv[]) {
  FILE *fp;
  Rec *recs;
  Header header;
  int i, verbose=0;
  int rund=0, runt=0;

  if (argc >= 2 && strcmp(argv[1], "-v") == 0) {
    argc--;
    argv++;
    verbose++;
  }

  if (argc != 2) {
    fprintf(stderr, "Usage: getABIdate filename\n");
    return 1;
  }

  if (NULL ==  (fp=fopen(argv[1],"rb"))) {
    perror(argv[1]);
    return 1;
  }

  if (read_header(fp, &header)) {
    fclose(fp);
    fprintf(stderr, "Failed to read header for %s\n", argv[1]);
    return 2;
  }

  if (NULL == (recs = (Rec *) malloc(header.recs * sizeof(Header)))) {
    perror("malloc");
    return 3;
  }

  fseek(fp,header.offset,0);
  for (i=0; i<header.recs; i++) {
    if (read_rec(fp,&recs[i])) {
      fclose(fp);
      free(recs);
      fprintf(stderr, "Failed to read header data for %s\n", argv[1]);
      return 4;
    }

    if (recs[i].index == 1 && strncmp(recs[i].id, "RUND",4) == 0) {
      /* RUND, being a 4 byte value, fits in the ptr1 field and so is
       * stored there.
       */
      rund = recs[i].ptr1;
    }

    if (recs[i].index == 1 && strncmp(recs[i].id, "RUNT",4) == 0) {
      runt = recs[i].ptr1;
    }
  }

  printf((verbose ?
	 "year %04d, month %02d, day %02d, time %02d:%02d:%02d\n" :
	 "%04d%02d%02d.%02d%02d%02d\n"),
	 rund >> 16, (rund >> 8) & 0xff, rund & 0xff,
	 runt >> 24, (runt >> 16) & 0xff, (runt >> 8) & 0xff);

  return 0;
}

