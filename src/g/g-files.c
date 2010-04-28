/*
 * File: g-files.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: routines for gap server file i/o
 *
 * Created: prior to 18 September 1992
 * Updated:
 *
 */

#include <staden_config.h>

#include <stdio.h>
#include <unistd.h> /* IMPORT: lseek */
#include <fcntl.h> /* IMPORT: O_RDWR */
#include <string.h>
#include <errno.h>
/*#include <malloc.h>*/

/*
 * mmap() support is just an experiment to allow benchmarking. It has been
 * tested on DEC Alphas for read-access (write-access via mmap doesn't
 * exist yet).
 * The results were small (10%) when the data is cached, but in the more
 * normal case (no cached data) the speed difference is negligible.
 *
 * Tested on linux on a LARGE (EIMER1: approx 900Mb) it was taking
 * 8mins real time to open and 25secs CPU time without MMAP
 * and 4.5mins/12sec with MMAP.
 *
 * This test database was largely fragmented though. After a copy_db the times
 * changed to 1m50s/18s w/o MMAP, ~2min/13s with MMAP.
 * Hence it seems the MMAP benefit is better solved by using copy_db to
 * "defrag" the database.
 */
#ifdef USE_MMAP
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "array.h"
#include "freetree.h"

#include "g-files.h"
#include "g-os.h"
#include "g-error.h"
#include "g-io.h" /* IMPORT: low_level_vector */
#include "g-db.h" /* IMPORT: panic_shutdown() */
#include "g-defs.h" /* IMPORT: G_AUX_SUFFIX */
#include "xalloc.h"

/* johnt 6/1/99 - O_BINARY required to open files in binary mode for MS Windows */
#ifndef O_BINARY
#define O_BINARY 0
#endif
 

static int g_read_aux_header(GFile *gfile, AuxHeader *header)
/*
 * Read the header from the aux file
 */
{
    int err;
    errno = 0;
    err = gfile->low_level_vector[GOP_READ_AUX_HEADER](gfile->fdaux,header,1);
    
    /* Auto-sense low-level vector set use based on the bit-size */
    if (!err) {
	if (header->format == G_32BIT) {
	    gfile->low_level_vector = gfile->swapped
		? low_level_vectors_swapped32
		: low_level_vectors32;
	} else {
	    gfile->low_level_vector = gfile->swapped
		? low_level_vectors_swapped64
		: low_level_vectors64;
	}
    }

    return err;
}

static int g_read_aux_index(GFile *gfile, AuxIndex *idx, int num)
/*
 * Read records from the index of the aux file
 */
{
    /*
     * NOTE -
     * for generality we avoid using machine independant IO here
     */
    return (gfile->low_level_vector[GOP_READ_AUX_INDEX])(gfile->fdaux,idx,num);
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
 * EXTERNAL ROUTINES
 */    


GFile *g_open_file(char *fn, int read_only)
/*
 * Open a file and its associated index
 * On Error:
 *	Returns NULL
 *      Sets xerrno
 */
{
    GFile *gfile;
    char fnaux[1024];
    GCardinal i;
    int tree_init;
    AuxIndex *idx_arr = NULL;
    off_t recsize;

    /* g_dump_file(fn); */

    gfile = NULL;

#define ABORT(E)\
    {\
         if (idx_arr) \
             xfree(idx_arr); \
	 g_free_gfile(gfile); \
	 gfile = NULL; \
	 (void)gerr_set(E); \
	 return NULL; \
    }

    /* check file name isn't too long */
    if ( strlen(fn) + strlen(G_AUX_SUFFIX) >= sizeof(fnaux) )
	ABORT(GERR_NAME_TOO_LONG);
    strcpy(fnaux,fn);
    strcat(fnaux,G_AUX_SUFFIX);

    /* allocate new data structure - GFile */
    gfile = g_new_gfile(G_32BIT);
    if (gfile == NULL)
	ABORT(GERR_OUT_OF_MEMORY);

    /* check access privilages */
    /* YUK! - to do */

    /* set file name */
    if ( (gfile->fname = (char *)xmalloc(strlen(fn)+1)) != NULL )
	strcpy(gfile->fname,fn);
    
    /* open file and its aux */
    /* LOW LEVEL IO HERE */
    errno = 0;
    if (read_only || (gfile->fd = open(fn,O_RDWR|O_BINARY)) == -1 )
	if ( !read_only || (gfile->fd = open(fn,O_RDONLY|O_BINARY)) == -1 )
	    ABORT(GERR_OPENING_FILE);
    /* LOW LEVEL IO HERE */
    errno = 0;
    if (read_only || (gfile->fdaux = open(fnaux,O_RDWR|O_BINARY)) == -1 )
	if ( !read_only || (gfile->fdaux = open(fnaux,O_RDONLY|O_BINARY)) == -1 )
	    ABORT(GERR_OPENING_FILE);

#if 0
    /*
     * WARNING: See warning in gap4/actf.c for why this is a problem and
     * has been commented out.
     */
    if (lockf(gfile->fd, F_TEST, 0)) {
	fprintf(stderr, "*** File %s is locked (inuse) by another process\n",
		fn);
	if (!read_only) {
	    close(gfile->fdaux);
	    ABORT(GERR_OPENING_FILE);
	}
    }
    if (!read_only) {
	lockf(gfile->fd, F_LOCK, 0);
    }
#endif


#ifdef USE_MMAP
    {
	struct stat sb;

	stat(fn, &sb);
	gfile->fdmap = (char *)mmap(NULL, sb.st_size, PROT_READ,
				    MAP_FILE | /* MAP_FIXED | */ MAP_SHARED,
				    gfile->fd, 0);
    }
#endif

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1==lseek(gfile->fdaux,0,0))
	ABORT(GERR_SEEK_ERROR);

    if (g_read_aux_header(gfile, &gfile->header))
	ABORT(GERR_READ_ERROR);

/*
    fprintf(stderr,"** \n");
    fprintf(stderr,"** Opening file %s\n",fn);
    fprintf(stderr,"**    file_size = %d\n",  gfile->header.file_size);
    fprintf(stderr,"**   block_size = %d\n", gfile->header.block_size);
    fprintf(stderr,"**  num_records = %d\n",gfile->header.num_records);
    fprintf(stderr,"**  max_records = %d\n",gfile->header.max_records);
    fprintf(stderr,"**    last_time = %d\n",  gfile->header.last_time);
    fprintf(stderr,"**        flags = %d\n",      gfile->header.flags);
    fprintf(stderr,"**    free_time = %d\n",  gfile->header.free_time);
    fprintf(stderr,"** \n");
*/

    /* allocate index */
    gfile->Nidx = gfile->header.num_records;
    if ( (gfile->idx = ArrayCreate(sizeof(Index),gfile->Nidx)) == NULL )
	ABORT(GERR_OUT_OF_MEMORY);

    (void) ArrayRef(gfile->idx,gfile->Nidx-1);
    for(i=0;i<gfile->Nidx;i++) arr(Index,gfile->idx,i).flags = G_INDEX_NEW;

    /* force Array.max field to be updated */
    (void)ArrayRef(gfile->idx,gfile->header.num_records-1);

    /* read cached freetree */
    /* LOW LEVEL IO HERE */
    errno = 0;
    recsize = (gfile->header.format == G_32BIT)
	? sizeof(AuxIndex32)
	: sizeof(AuxIndex);
    lseek(gfile->fdaux, sizeof(AuxHeader) +
	  gfile->header.num_records * recsize, SEEK_SET);
#ifdef CACHE_FREETREE
    gfile->freetree = (gfile->header.format == G_32BIT)
	? freetree_load_int4(gfile->fdaux, gfile->header.last_time)
	: freetree_load_int8(gfile->fdaux, gfile->header.last_time);
    tree_init = gfile->freetree ? 1 : 0;
#else
    gfile->freetree = NULL;
    tree_init = 0;
#endif

    /* allocate freetree, if not loaded */
    if (!tree_init) {
	gfile->freetree = freetree_create(0,MAX_GImage);
	if (gfile->freetree==NULL)
	    ABORT(GERR_OUT_OF_MEMORY);
    }

    /* read aux index and initialise */
    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1==lseek(gfile->fdaux,sizeof(AuxHeader),0))
	ABORT(GERR_SEEK_ERROR);

    /* Temporarily load the entire index in one step */
    idx_arr = (AuxIndex *)xmalloc(gfile->header.num_records *
				  sizeof(AuxIndex) + 1);
    if (g_read_aux_index(gfile, idx_arr, gfile->header.num_records))
	ABORT(GERR_READ_ERROR);
    
    for (i=0;i<gfile->header.num_records;i++) {
	AuxIndex *aidx;
	GToggle toggle;

	/* aidx points to the current processed index record */
	aidx = &idx_arr[i];

	/* set toggle to be current index of fruity chews */
	toggle = g_toggle_state(gfile->header.last_time, aidx);
	    
	if (toggle != G_NO_TOGGLE) {
	    if (gfile->header.flags & G_BLOCK_SIZE_CHANGED) {
		arr(Index,gfile->idx,i).aux_allocated = aidx->used[toggle];
	    } else {
		arr(Index,gfile->idx,i).aux_allocated =
		    FREETREE_BLOCK(aidx->used[toggle],
				   gfile->header.block_size);
	    }
	    arr(Index,gfile->idx,i).aux_image = aidx->image[toggle];
	    arr(Index,gfile->idx,i).aux_time = aidx->time[toggle];
	    arr(Index,gfile->idx,i).aux_used = aidx->used[toggle];
	} else {
	    printf("No toggle for record %d\n", i);
	}

	/* update freetree */
	if (toggle != G_NO_TOGGLE) {
	    GImage image;
	    image = aidx->image[toggle];
	    if (image!=G_NO_IMAGE) {
		int err;
		/*
		 * NOTE: allocated versus used space in file
		 * We optimise on file space usage since we registered only
		 * the used (arr(Index,gfile->idx,i).aux.used[t]) space, rather
		 * than the whole of the allocated space. To do this we
		 * would have to set the allocated space to equal the used
		 * space. It is not necessary to modify index records in file.
		 */
		if (!tree_init && !read_only) {
		    err = freetree_register(gfile->freetree, image,
					    arr(Index,gfile->idx,i)
					       .aux_allocated);
		    if (err == -1) {
			if (err == TREE_OVERLAP) {
			    fprintf(stderr,
				    "** SERIOUS CORRUPTION DETECTED "
				    "- file %s\n",
				    g_filename(gfile));
			    fprintf(stderr,"** record %d: "
				    "image overlaps that of another record\n",
				    i);
			}
			ABORT(err);
		    }
		}

		/* flag as in use (G_INDEX_USED isn't needed or ever used) */
		arr(Index,gfile->idx,i).flags = G_INDEX_NONE;
	    } else {
		/* 08/10/96
		 * This is not an error, but some debugging info to see
		 * whether or not we had a problem before when we set all
		 * records to G_INDEX_NONE, rather than just the ones with
		 * images (which stops them being initialised later).
		 */
		/* printf("Leaving record %d as G_INDEX_NEW\n", i); */
	    }
	}

    }

    xfree(idx_arr);
    
#undef ABORT

    /*
      tree_print(gfile->freetree);
      tree_print_lists(gfile->freetree);
    */

    return gfile;
}


void g_close_file(GFile *g)
/*
 * Close a file and its associated index
 */
{
    /* tidy up files if necessary */

    /* free data structures */
    g_free_gfile(g);
}

/* 6/1/99 johnt - Don't have fcntl with Visual C++ */
#if !defined(_MSC_VER) && !defined(__APPLE__) && !defined(__MINGW32__)
int g_sync_aux_on(GFile *gfile)
/*
 * Force a flush of all prior data and then all subsequent data is to be
 * written in SYNC mode.
 */
{
    int fdaux = gfile->fdaux;
#if 0
    char c;
#endif

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1 == fcntl(fdaux, F_SETFL, O_RDWR | O_SYNC))
	return gerr_set(GERR_SYNC);

#if 0
    /* Force a SYNC write to flush all prior data. */
    if (-1 == lseek(fdaux, 0, 0)) return gerr_set(GERR_SEEK_ERROR);
    if (-1 == read(fdaux, &c, 1)) return gerr_set(GERR_READ_ERROR);
    lseek(fdaux, 0, 0);
    if (-1 == write(fdaux, &c, 1)) return gerr_set(GERR_WRITE_ERROR);
#endif
    if (-1 == fsync(fdaux)) return gerr_set(GERR_SYNC);

    return 0;
}

int g_sync_aux_off(GFile *gfile)
/*
 * All subsequent data should not be written in SYNC mode.
 */
{
    int fdaux = gfile->fdaux;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1 == fcntl(fdaux, F_SETFL, O_RDWR))
	return gerr_set(GERR_SYNC);

    return 0;
}


int g_sync_on(GFile *gfile)
/*
 * Force a flush of all prior data and then all subsequent data is to be
 * written in SYNC mode.
 */
{
    int fd = gfile->fd;
    char c;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1 == fcntl(fd, F_SETFL, O_RDWR | O_SYNC))
	return gerr_set(GERR_SYNC);

#if 0
    /* Sync all data to the disk */
    fsync(gfile->fd);
    fsync(gfile->fdaux);
#endif

    /*
     * Force a SYNC write to flush all prior data.
     *
     * Faster than the fsync() method (0.6x time), but still much slower
     * than no sync at at (~2.0x time).
     */
    if (-1 == lseek(fd, 0, 0)) return gerr_set(GERR_SEEK_ERROR);
    if (-1 == read(fd, &c, 1)) return gerr_set(GERR_READ_ERROR);
    lseek(fd, 0, 0);
    if (-1 == write(fd, &c, 1)) return gerr_set(GERR_WRITE_ERROR);

    return 0;
}

int g_sync_off(GFile *gfile)
/*
 * All subsequent data should not be written in SYNC mode.
 */
{
    int fd = gfile->fd;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1 == fcntl(fd, F_SETFL, O_RDWR))
	return gerr_set(GERR_SYNC);

    return 0;
}

#endif


int g_write_aux_header(GFile *gfile)
/*
 * Write the header of the aux file
 */
{
    int fdaux = gfile->fdaux;
    AuxHeader *header = &gfile->header;
    int check;

    /* check arguments */
    if (gfile==NULL) return gerr_set(GERR_INVALID_ARGUMENTS);

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1==lseek(fdaux,0,0)) return gerr_set(GERR_SEEK_ERROR);
    /* LOW LEVEL IO HERE */
    errno = 0;
    /* check = write(fdaux,header,sizeof(AuxHeader)); */
    check = (gfile->low_level_vector[GOP_WRITE_AUX_HEADER])(fdaux,header,1);
    if (check)
	return gerr_set(GERR_WRITE_ERROR);
    else
	return 0;
}



int g_write_aux_index(GFile *gfile, GCardinal rec)
/*
 * Read a record from the index of the aux file
 */
{
    AuxIndex idx;
    int fdaux = gfile->fdaux;
    int check;
    off_t recsize;
    
    idx.image[0] = arr(Index,gfile->idx,rec).aux_image;
    idx.time [0] = arr(Index,gfile->idx,rec).aux_time;
    idx.used [0] = arr(Index,gfile->idx,rec).aux_used;
    idx.image[1] = 0;
    idx.time [1] = 0;
    idx.used [1] = 0;

    /* check arguments */
    if (gfile==NULL) return gerr_set(GERR_INVALID_ARGUMENTS);

    /* LOW LEVEL IO HERE */
    errno = 0;
    recsize = (gfile->header.format == G_32BIT)
	? sizeof(AuxIndex32)
	: sizeof(AuxIndex);
    if (-1==lseek(fdaux,sizeof(AuxHeader)+rec*recsize,0))
	return gerr_set(GERR_SEEK_ERROR);

    /* LOW LEVEL IO HERE */
    /* check = write(fdaux,idx,sizeof(AuxIndex)); */
    errno = 0;
    check = (gfile->low_level_vector[GOP_WRITE_AUX_INDEX])(fdaux,&idx,1);
    if (check)
	return gerr_set(GERR_WRITE_ERROR);
    else
	return 0;
}





int g_remove_client(GFile *gfile, GClient client)
/*
 * Remove all locks in this file for this client
 */
{

    /* check arguments */
    if (gfile==NULL) return gerr_set(GERR_INVALID_ARGUMENTS);

    /*
     * if this client has locked the file, remove the lock
     */
    if (gfile->flock_client == client &&
	gfile->flock_status == G_FLOCK_LOCKED) {

	/* reset flock variables */
	gfile->flock_status = G_FLOCK_NONE;
	gfile->flock_client = 0;
	gfile->flock_view = -1;

    }

    return 0;
}





char *g_filename(GFile *gfile)
/*
 * return file name of open file
 */
{
    if (gfile->fname == NULL)
	return "(unknown)";
    else
	return gfile->fname;

}

int g_check_header(GFile *gfile)
/*
 * Checks whether the on-disk copy matches the in-memory copy of the
 * header.
 *
 * The rationale behind this is that users still manually remove the BUSY
 * file because "they know best". They then act all suprised when their
 * database becomes corrupt. Our alternative here is to check the master
 * time-stamp in the Aux header and to simply abort Gap4 (losing all unsaved
 * edits) if it has been changed external to this process.
 *
 * For efficiencies sake, we do not want to call this function too often.
 * Hence it is planned to be called only on the first disk update after a
 * flush. This still leaves room for race conditions (multiple people editing
 * simultaenously), but it covers the more common case of a gap4 session
 * left open for ages and then going back to it later.
 */
{
    AuxHeader diskheader;

    if (gfile == NULL)
	return gerr_set(GERR_INVALID_ARGUMENTS);

    /* Re-read from disk */
    if (-1==lseek(gfile->fdaux,0,0))
	return gerr_set(GERR_SEEK_ERROR);

    g_read_aux_header(gfile, &diskheader);

    if (diskheader.last_time != gfile->header.last_time) {
	fprintf(stderr, "** SERIOUS PROBLEM - file %s\n",
		g_filename(gfile));
	fprintf(stderr, "** Time stamp modified by another process\n");
	fprintf(stderr, "** Hint: DO NOT REMOVE LOCK FILES!\n**\n");
	fprintf(stderr, "** The '%s.log' file contains information on\n",
		g_filename(gfile));
	fprintf(stderr, "** who else has the database open.\n");
	panic_shutdown();
    }

    return 0;
}


