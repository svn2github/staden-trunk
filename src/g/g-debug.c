/*
 * File: g-debug.c
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: Debug routines for the g library
 * 
 * Created: 11 June 1996
 */

#include <stdio.h>
#include <unistd.h> /* IMPORT: lseek */
#include <fcntl.h> /* IMPORT: O_RDWR */
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "misc.h"
#include "array.h"
#include "freetree.h"

#include "g-files.h"
#include "g-os.h"
#include "g-error.h"
#include "g-io.h" /* IMPORT: low_level_vector */
#include "g-db.h" /* IMPORT: panic_shutdown() */
#include "g-defs.h" /* IMPORT: G_AUX_SUFFIX */
#include "xalloc.h"

#define MAX_BUF 100

static int g_read_aux_header(int fdaux, AuxHeader *header)
/*
 * Read the header from the aux file
 */
{
    /*
     * NOTE -
     * for generality we avoid using machine independant IO here
     */

    /* LOW LEVEL IO HERE */
    /* return (read(fdaux,header,sizeof(AuxHeader))!=sizeof(AuxHeader)); */
    return (low_level_vector[GOP_READ_AUX_HEADER])(fdaux,header,1);

}

static int g_read_aux_index(int fdaux, AuxIndex *idx)
/*
 * Read a record from the index of the aux file
 */
{
    /*
     * NOTE -
     * for generality we avoid using machine independant IO here
     */

    /* LOW LEVEL IO HERE */
    /* return (read(fdaux,idx,sizeof(AuxIndex))!=sizeof(AuxIndex)); */
    return (low_level_vector[GOP_READ_AUX_INDEX])(fdaux,idx,1);

}




static GToggle g_toggle_state(GTimeStamp time, AuxIndex *idx)
/*
 * return the most recent image for a record from index
 * which was created on or before time
 */
{
    GToggle r,i;
    GTimeStamp t;
    r = G_NO_TOGGLE;
    t = 0; /* assumes signed type */
    for (i=0;i<2;i++) {
	if (idx->time[i] <= time &&
	    idx->time[i] >= t) {
	    t = idx->time[i];
	    r = i;
	}
    }

    return r;
}


/*
 * Dumps a G database in an ASCII readable form. Modelled around the start
 * of the g_open_file() routine.
 */
void g_dump_file(char *fn) {
    GFile *gfile = NULL;
    char fnaux[1024];
    AuxIndex aux_ind;
    int i;
    
#define ABORT(E)\
    {\
	 g_free_gfile(gfile); \
	 gfile = NULL; \
	 (void)gerr_set(E); \
	 perror("ABORT"); \
	 return; \
    }
    
    /* check file name isn't too long */
    if (strlen(fn) + strlen(G_AUX_SUFFIX) >= sizeof(fnaux))
	ABORT(GERR_NAME_TOO_LONG);
    strcpy(fnaux, fn);
    strcat(fnaux, G_AUX_SUFFIX);

    /* allocate new data structure - GFile */
    gfile = g_new_gfile();
    if (gfile == NULL)
	ABORT(GERR_OUT_OF_MEMORY);

    /* set file name */
    if ((gfile->fname = (char *)xmalloc(strlen(fn)+1)) != NULL)
	strcpy(gfile->fname, fn);
    
    /* open file and its aux */
    /* LOW LEVEL IO HERE */
    if ((gfile->fd = open(fn, O_RDONLY)) == -1)
	ABORT(GERR_OPENING_FILE);
    /* LOW LEVEL IO HERE */
    if ((gfile->fdaux = open(fnaux, O_RDONLY)) == -1)
	ABORT(GERR_OPENING_FILE);

    /* LOW LEVEL IO HERE */
    if (-1 == lseek(gfile->fdaux, 0, 0))
	ABORT(GERR_SEEK_ERROR);
    if (g_read_aux_header(gfile->fdaux, &gfile->header))
	ABORT(GERR_READ_ERROR);

    printf("** \n");
    printf("** Opening file %s\n",fn);
    printf("**    file_size = %"PRIGImage"\n", gfile->header.file_size);
    printf("**   block_size = %"PRIGCardinal"\n", gfile->header.block_size);
    printf("**  num_records = %"PRIGCardinal"\n", gfile->header.num_records);
    printf("**  max_records = %"PRIGCardinal"\n", gfile->header.max_records);
    printf("**    last_time = %"PRIGCardinal"\n", gfile->header.last_time);
    printf("**        flags = %"PRIGHFlags"\n", gfile->header.flags);
    printf("** \n");

    /* allocate index */
    gfile->Nidx = gfile->header.num_records;
    if ((gfile->idx = ArrayCreate(sizeof(Index), gfile->Nidx)) == NULL )
	ABORT(GERR_OUT_OF_MEMORY);

    (void) ArrayRef(gfile->idx, gfile->Nidx-1);
    for(i = 0; i < gfile->Nidx; i++)
	arr(Index, gfile->idx, i).flags = G_INDEX_NEW;

    /* read aux index and initialise */
    /* LOW LEVEL IO HERE */
    if (-1 == lseek(gfile->fdaux, sizeof(AuxHeader), 0))
	ABORT(GERR_SEEK_ERROR);

    /* force Array.max field to be updated */
    (void)ArrayRef(gfile->idx, gfile->header.num_records-1);

    printf("global_time %08x\n", gfile->header.last_time);

    for (i = 0; i < gfile->header.num_records; i++) {
	char buf[MAX_BUF];
	int toggle, len, len_r;

	/* Load index for this record */
	if (g_read_aux_index(gfile->fdaux, &aux_ind))
	    ABORT(GERR_READ_ERROR);

	/* Compute toggle */
	toggle = g_toggle_state(gfile->header.last_time, &aux_ind);

	/* LOW LEVEL IO HERE */
	if (-1 == lseek(gfile->fd, aux_ind.image[toggle], 0))
	    ABORT(GERR_SEEK_ERROR);

	len = MIN(aux_ind.used[toggle], MAX_BUF);
	/* LOW LEVEL IO HERE */
	if (-1 == (len_r = read(gfile->fd, buf, len)))
	    ABORT(GERR_READ_ERROR);
	if (len_r != len) {
	    fprintf(stderr, "WARNING: Read too short. Requested %d, got %d\n",
		    len, len_r);
	}

	printf("record %05d pos %020"PRIGImage" len %08d : %08x",
	       i, aux_ind.image[toggle], aux_ind.used[toggle],
	       (((((buf[0] << 8) + buf[1]) << 8) + buf[2]) << 8) + buf[3]);
	if (len > 4)
	    printf(" %c%c%c%c%c%c%c%c\n",
		   isprint(buf[4])?buf[4]:'.',
		   isprint(buf[5])?buf[5]:'.',
		   isprint(buf[6])?buf[6]:'.',
		   isprint(buf[7])?buf[7]:'.',
		   isprint(buf[8])?buf[8]:'.',
		   isprint(buf[9])?buf[9]:'.',
		   isprint(buf[10])?buf[10]:'.',
		   isprint(buf[11])?buf[11]:'.');
	else
	    putchar('\n');
    }

#undef ABORT

    g_free_gfile(gfile);
    return;
}
