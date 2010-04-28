/*
 * File: g-requests.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: gap server requests
 *
 * Created: 18-Sep-1992
 * Updated:
 *
 */

#include <staden_config.h>

#include <stdio.h>		/* IMPORT: NULL */
/*#include <malloc.h>*/
#include <unistd.h>		/* IMPORT: lseek */
#include <string.h>		/* IMPORT: memset() */
#include <errno.h>

#include "array.h"
#include "freetree.h"

#include "g-request.h"
#include "g-files.h"
#include "g-defs.h"		/* IMPORT: G_LOCKs */
#include "g-db.h"		/* IMPORT: panic_shutdown() */
#include "g-error.h"
#include "g-misc.h"		/* IMPORT: G_Assert */

#include "xalloc.h"


/*************************************************************
 * Utility routines
 ************************************************************/

/*
 * return the current image for this record
 */
#define get_current_image(gfile,rec) \
    (arr(Index,(gfile)->idx,(rec)).aux_image)





#define GERR_ERROR 1

static int initialise_record(GFile *gfile, GCardinal rec)
/*
 * initialise a new record
 */
{

    /* don't allow a old record to be initialised */
    if (! (arr(Index,gfile->idx,rec).flags & G_INDEX_NEW) ) return GERR_ERROR;


    /* initialise */
    arr(Index,gfile->idx,rec).aux_image = G_NO_IMAGE;
    arr(Index,gfile->idx,rec).aux_time = G_YEAR_DOT;
    arr(Index,gfile->idx,rec).aux_allocated = (GCardinal)0;
    arr(Index,gfile->idx,rec).aux_used = (GCardinal)0;
    arr(Index,gfile->idx,rec).flags = (GFlags)0;

    return 0;
}



/* ARGSUSED */
void init_cache(GDB *gdb, GFile *gfile, GCardinal rec, GLock  lock, GView view)
/*
 * Find a cache entry suitable for this rec-lock combination
 * If an existing one already exists, use it
 * otherwise create a new cache entry
 */
{
    GImage image;
    Cache *cache;

    /*
     * get default image from rec
     * if this is an unitialised record we will need to
     *   initialise and allocate a record
     */
    if (arr(Index,gfile->idx,rec).flags & G_INDEX_NEW)
	/* an error only means the record has already been initialised */
	(void)initialise_record(gfile,rec);

    image = get_current_image(gfile,rec);

    /* find a cache for this image */
    cache = &arr(View,gdb->view,view).lcache;

    /* initialise */
    cache->rec   = rec;
    cache->image = image;
    cache->allocated = arr(Index,gfile->idx,rec).aux_allocated;
    cache->used = arr(Index,gfile->idx,rec).aux_used;
}






/*************************************************************
 * checking routines
 *************************************************************/

#define IDX_INCREMENT 10

static int check_GIOVec(GIOVec *v, GCardinal vcnt, int *len)
/*
 *
 */
{
    int i;
    if (v==NULL || vcnt<0) return 1;

    *len = 0;
    for(i=0;i<vcnt;i++) {
	if(v[i].len <= 0 ||
	   v[i].buf == NULL) return 1;
	*len+=v[i].len;
    }

    return 0;

}



static int check_record(GFile *gfile, GCardinal rec)
/*
 * Check that rec is a valid record.
 * Don't allow:
 *    1. rec < 0
 * Allow anything else, but:
 *    1. if rec >= gfile->Nidx, need to expand array gfile->idx;
 *    2. if rec >= gfile->header.num_records, need to initialise entry in gfile->idx;
 */
{
    if ( rec >= gfile->Nidx ) {
	/*
	 * need to make room in gfile->idx array - realloc
	 */
	GCardinal i;
	GCardinal new_Nidx;

	new_Nidx = (rec+1)+IDX_INCREMENT;

	/* for extension of array */
	(void) ArrayRef(gfile->idx,new_Nidx-1);

	/*
	 * Set entries [Nidx..newNidx] to NEW
	 * We will delay the actual initialisation till when we
	 * need to access each individual entry
	 */
	for (i=gfile->Nidx;i<new_Nidx;i++)
	    arr(Index,gfile->idx,i).flags = G_INDEX_NEW;
	gfile->Nidx = new_Nidx;

    }

    /*
     * IMPORTANT NOTE:
     * We postpone the initialisation of new records in the auxilliary index,
     * adjusting the number of images, and writing to the aux index until we
     * do an g_unlock(). This is so that we do not have to do unnecessary
     * work until it is really needed.
     */

    return 0;
}

#define check_mode(mode) ((mode) > G_LOCK_EX)

#define check_client_mode(gdb,c,mode) \
    (arr(Client,(gdb)->client,(c)).max_lock < (mode))

#define check_file_lock(gfile, client) \
    ((gfile)->flock_status == G_FLOCK_LOCKED && \
     (gfile)->flock_client != (client))

#define check_view(gdb, v) \
    ((v)<0 || (v)>=(gdb)->Nview || \
     (arr(View,gdb->view,(v)).flags & G_VIEW_FREE))

#define check_client(gdb, c) \
    ((c)<0 || (c)>=(gdb)->Nclient)



/*************************************************************
 * low level io follows
 *************************************************************/



static int read_image_(int fd, GImage image, GCardinal used, void *buf, GCardinal len)
/*
 * Read `len' characters from image (where there are `used' bytes available).
 * Zero pad buf with 0 bytes
 */
{
    int in;

    if (image == G_NO_IMAGE) {
	/*
	 * don't want to do anything here really. Shouldn't rely on reading from
	 * uninitialised data anyway
	 */
	in = 0;
    } else {
	int check;
	in = (len>used)?used:len;
	/* LOW LEVEL IO HERE */
	errno = 0;
	if (-1==lseek(fd, (off_t)image, 0))
	    return gerr_set(GERR_SEEK_ERROR);

	/* LOW LEVEL IO HERE */
	errno = 0;
	check = read(fd, buf, in);
	if (check !=  in)
	    return gerr_set(GERR_READ_ERROR);
    }
    
    /* pad buf with 0 bytes */
    memset((char*)buf+in,'\0', len-in);
    
    return 0;
}

#ifdef USE_MMAP
static int mmap_read_image_(char *addr, GImage image, GCardinal used, void *buf, GCardinal len)
/*
 * Read `len' characters from image (where there are `used' bytes available).
 * Zero pad buf with 0 bytes
 */
{
    int in;

    if (image == G_NO_IMAGE) {
	/*
	 * don't want to do anything here really. Shouldn't rely on reading from
	 * uninitialised data anyway
	 */
	in = 0;
    } else {
	in = (len>used)?used:len;
	memcpy(buf, addr+image, in);
    }
    
    /* pad buf with 0 bytes */
    memset((char*)buf+in,'\0', len-in);
    
    return 0;
}
#endif







static int readv_image_(int fd, GImage image, GCardinal used, GIOVec *v, GCardinal vcnt)
/*
 * This mocks readv(). I originally wrote it using readv, but it was more
 * concise using read().
 */
{
    int parti, partj, count;
    int check;

    parti = 0;
    partj = 0;
    if (image != G_NO_IMAGE) {

	/* LOW LEVEL IO HERE */
	errno = 0;
	if (-1==lseek(fd, (off_t)image, 0))
	    return gerr_set(GERR_SEEK_ERROR);

	if (used > 0 && vcnt > 0) {
	    /*
	     * Read what we can
	     */
	    count = 0;
	    for(parti=0;parti<vcnt && count<used;parti++) {
		
		partj = v[parti].len < used-count ? v[parti].len : used-count;
		
		/* LOW LEVEL IO HERE */
		errno = 0;
		check = read(fd, v[parti].buf, partj);
		if (check !=  partj)
		    return gerr_set(GERR_READ_ERROR);
		count += partj;
	    }
	    parti--;
	}
    }

    
    /* pad buffers with 0 bytes */
    if(parti<vcnt) {
	/* fill a partial block */
	if(v[parti].len > partj)
	    memset((char *)v[parti].buf+partj,'\0',v[parti].len-partj);
	/* fill remaining full blocks */
	for (parti++;parti<vcnt;parti++)
	    memset(v[parti].buf,'\0',v[parti].len);
    }
    
    return 0;
}

#ifdef USE_MMAP
static int mmap_readv_image_(char *addr, GImage image, GCardinal used, GIOVec *v, GCardinal vcnt)
/*
 * This mocks readv(). I originally wrote it using readv, but it was more
 * concise using read().
 */
{
    int parti, partj, count;

    parti = 0;
    partj = 0;
    if (image != G_NO_IMAGE) {

	if (used > 0 && vcnt > 0) {
	    /*
	     * Read what we can
	     */
	    count = 0;
	    for(parti=0;parti<vcnt && count<used;parti++) {
		
		partj = v[parti].len < used-count ? v[parti].len : used-count;
		
		memcpy(v[parti].buf, addr+image+count, partj);
		count += partj;
	    }
	    parti--;
	}
    }

    
    /* pad buffers with 0 bytes */
    if(parti<vcnt) {
	/* fill a partial block */
	if(v[parti].len > partj)
	    memset((char *)v[parti].buf+partj,'\0',v[parti].len-partj);
	/* fill remaining full blocks */
	for (parti++;parti<vcnt;parti++)
	    memset(v[parti].buf,'\0',v[parti].len);
    }
    
    return 0;
}
#endif

static int write_zeros(int fd, int i)
/*
 * Write out 'i' zero bytes
 */
{
    static long zbuf[4] = {0L,0L,0L,0L};
    int check;

    for(; (size_t)i>=sizeof(zbuf); i-=sizeof(zbuf)) {
	/* LOW LEVEL IO HERE */
	errno = 0;
	check = write(fd, zbuf, sizeof(zbuf));
	if (check !=  sizeof(zbuf))
	    return gerr_set(GERR_WRITE_ERROR);
    }
    /* LOW LEVEL IO HERE */
    if (i) {
	errno = 0;
	check = write(fd, zbuf, i);
	if (check !=  i)
	    return gerr_set(GERR_WRITE_ERROR);
    }

    return 0;
}



static int write_image_(int fd, GImage image, GCardinal allocated, void *buf, GCardinal len)
/*
 * Write `len' characters to image (where there are `allocated' bytes available).
 * Zero pad buf with 0 bytes
 */
{
    int check;
    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1==lseek(fd, (off_t)image, 0))
	return gerr_set(GERR_SEEK_ERROR);

    /* LOW LEVEL IO HERE */
    errno = 0;
    check = write(fd, buf, (int)len);
    if (check !=  len)
	return gerr_set(GERR_WRITE_ERROR);

    /*
     * Fill to allocated size with 0 bytes ????
     * Do it anyway!
     */
    return (allocated-len) > 0 ? write_zeros(fd,allocated-len) : 0;
}

static int writev_image_(int fd, GImage image, GCardinal allocated, GIOVec *v, GCardinal vcnt)
/*
 * This mocks writev().
 */
{
    int parti, partj, count;
    int check;

    count = 0;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (-1==lseek(fd, (off_t)image, 0))
	return gerr_set(GERR_SEEK_ERROR);

    if (allocated > 0 && vcnt > 0) {
	/*
	 * Write what we can
	 */
	for(parti=0;parti<vcnt && count<allocated;parti++) {

	    partj = v[parti].len < allocated-count ? v[parti].len : allocated-count;

	    /* LOW LEVEL IO HERE */
	    errno = 0;
	    check = write(fd, v[parti].buf, partj);
	    if (check !=  partj)
		return gerr_set(GERR_WRITE_ERROR);
	    count += partj;
	}
    }

    return (allocated-count) > 0 ? write_zeros(fd,allocated-count) : 0;
}






/*************************************************************
 * external routines start here
 *************************************************************/



/* ARGSUSED */
GView g_lock_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GLock lock)
/*
 *
 */
{
    GFile *gfile;
    GView view;

    /* check arguments */
    if(gdb==NULL || check_client(gdb,c)) {
	(void) gerr_set(GERR_INVALID_ARGUMENTS);
	return -1;
    }
    gfile = gdb->gfile;

    if (check_record(gfile,rec))
	return -1;


    /*
     * Get a new view.
     * We do this before getting a cache entry because
     * it is a lot easier to destroy a view if we have
     * problems with getting the cache.
     */
    if ( (view = g_new_view(gdb)) == -1) {
	(void)gerr_set(GERR_OUT_OF_MEMORY);
	return -1;
    }

    /* find an existing cache entry for this record and lock */
    init_cache(gdb,gfile,rec,lock,view);

    /*
     * do our duty with this view
     */
    arr(View,gdb->view,view).client = c;
    arr(View,gdb->view,view).flags = G_VIEW_USED;

    return view;
}




/* ARGSUSED */
int g_upgrade_(GDB *gdb, GClient c, GView v, GLock lock)
/*
 * Upgrade (or downgrade) a lock
 */
{
    return 0;
}



static GTimeStamp next_edtime(GFile *gfile)
/*
 * 
 */
{
    GTimeStamp edtime;

    edtime = gfile->header.last_time;
    edtime++;
    if (edtime==0) {
	/*
	 * Panic! Time stamps are about to wrap around
	 */
	fprintf(stderr,"** SERIOUS PROBLEM - file %s\n",g_filename(gfile));
	fprintf(stderr,"** time stamp wrap around\n");
	panic_shutdown();
    }

    return edtime;
}


static void initialise_records(GFile *gfile, GCardinal rec)
/*
 * The file `gfile' now has `rec' records in it.
 * Initialiase new records if need be.
 */
{
    if (rec > gfile->header.num_records) {
	GCardinal i; /* loop variable */
	int err;
	for(i=gfile->header.num_records; i<rec; i++) {
	    err = initialise_record(gfile,i);
	    /*
	     * An error here means we're trying to initialise an already
	     * initialised record. Ignore!
	     */
	    err = g_write_aux_index(gfile,i);
	    /*
	     * If we get an error here, something is seriously wrong!
	     */
	    if (err) {
		fprintf(stderr,"** SERIOUS PROBLEM\n");
		fprintf(stderr,"** record %d: failed to write to index.\n",i);
		panic_shutdown();
	    }
	}
	gfile->header.num_records = rec;
    }

}

static void update_record(GFile *gfile, GCardinal rec, GImage image,
			  GCardinal allocated, GCardinal used,
			  GTimeStamp edtime)
{
    GImage old_image;
    GCardinal old_allocated; /* space allocated for previous image */
    int err;
    
    /*
     * we may need to initialise other uninitialised records
     */
    initialise_records(gfile,rec+1);
    
    /*
     * we may need to adjust the size of the file
     */
    if (gfile->header.file_size < image + (GCardinal)allocated) {
	gfile->header.file_size = image + (GCardinal)allocated;
    }
    
    /*
     *
     */
    old_image = arr(Index,gfile->idx,rec).aux_image;
    old_allocated = arr(Index,gfile->idx,rec).aux_allocated;
    
    /* need to update idx image + time */
    arr(Index,gfile->idx,rec).aux_image = image;
    arr(Index,gfile->idx,rec).aux_allocated = allocated;
    
    arr(Index,gfile->idx,rec).aux_used = used;
    arr(Index,gfile->idx,rec).aux_time = edtime;

    if (image == G_NO_IMAGE) {
	arr(Index,gfile->idx,rec).flags = G_INDEX_NEW;	
    }
    
    err = g_write_aux_index(gfile,rec);
    /*
     * If we get an error here, something is seriously wrong!
     */
    if (err) {
	fprintf(stderr,"** SERIOUS PROBLEM - file %s\n",g_filename(gfile));
	fprintf(stderr,"** record %d: failed to write to index.\n",rec);
	panic_shutdown();
    }
    
    /* check through the cache list */
    if (old_image != G_NO_IMAGE) {
	err = freetree_unregister(gfile->freetree,old_image,old_allocated);
	if (err) {
	    gerr_set(err);
	    fprintf(stderr,"** SERIOUS PROBLEM - file %s\n",
		    g_filename(gfile));
	    fprintf(stderr,"** In update_record(): "
		    "freetree_unregister returned error code %d.\n",err);
	    panic_shutdown();
	}
    }
}







static void update_header(GFile *gfile, GTimeStamp edtime)
/*
 *
 */
{
    
    int err;

/*
    fsync(gfile->fd);
    fsync(gfile->fdaux);
*/
    /*
     * updates have been made - write the header
     */
    gfile->header.last_time = edtime;
    err = g_write_aux_header(gfile);
    /*
     * If we get an error here, something is seriously wrong!
     */
    if (err) {
	fprintf(stderr,"** SERIOUS PROBLEM - file %s\n",g_filename(gfile));
	fprintf(stderr,"** failed to write to file header.\n");
	panic_shutdown();
    }
/*
    fsync(gfile->fd);
    fsync(gfile->fdaux);
*/
}









static int g_unlock_views(GDB *gdb, GView v)
/*
 *
 */
{
    GFile *gfile;
    GTimeStamp edtime;
    int updates;
    GView nextv;

    /* check arguments */
    if (gdb==NULL) return gerr_set(GERR_INVALID_ARGUMENTS);

    /*
     * Check v
     */
    if (-1==v) return 0;
    
    /*
     * assume all views are for the same gfile
     */
    gfile = gdb->gfile;

    /*
     * Determine edit time for all updates
     */
    edtime = next_edtime(gfile);

    /* initialise number of updates */
    updates = 0;

    /*
     * unlock each view in turn
     */
    for(; v!=-1; v=nextv) {

	Cache *cache;
	GImage image; /* default image for record after unlock */

	nextv = arr(View,gdb->view,v).next;
	cache = &arr(View,gdb->view,v).lcache;

	/*
	 * Has there been an update?
	 */
	if (arr(View,gdb->view,v).flags & G_VIEW_UPDATED &&
	    !(arr(View,gdb->view,v).flags & G_VIEW_ABANDONED) ) {
	    /*
	     * Yes - we will have to write to the aux file
	     */

	    /*
	     * Check that the allocated size is still the same. If we can,
	     * shrink it.
	     */
	    {
		int nalloc = FREETREE_BLOCK(cache->used,
					    gfile->header.block_size);
		if (cache->allocated > nalloc) {
		    freetree_unregister(gfile->freetree,
					cache->image + nalloc,
					cache->allocated - nalloc);
		    cache->allocated = nalloc;
		}
	    }

	    update_record(gfile, cache->rec, cache->image, cache->allocated,
			  cache->used, edtime);
	    /* count updates to file */
	    updates++;
    
	}


	/*
	 * Are we just flushing?
	 */
	if ( (arr(View,gdb->view,v).flags & G_VIEW_FLUSHED)  &&
	    !(arr(View,gdb->view,v).flags & G_VIEW_UNLOCKED)  &&
	    !(arr(View,gdb->view,v).flags & G_VIEW_ABANDONED) ) {
	    /*
	     * Yes - reset flags on view
	     */
	    arr(View,gdb->view,v).flags = G_VIEW_USED;
	    arr(View,gdb->view,v).next = -1;
	} else {
	    
	    
	    /*
	     * Restoring the lock on the record.
	     * There may be a better way of doing this.
	     */
	    image = get_current_image(gfile,cache->rec);	    
	    
	    /* free view */
	    g_free_view(gdb,v);
	    
	    /*
	     * remove reference from cache
	     * if refs==0 and its image is different from default for record
	     *    freetree_unregister();
	     */
	    {
		/* reclaim image if no references and it is no longer the default */
		if ( cache->image != image && cache->image != G_NO_IMAGE ) {
		    int err;
		    err = freetree_unregister(gfile->freetree,
					      cache->image,cache->allocated);
		    if (err) {
			gerr_set(err);
			fprintf(stderr,"** SERIOUS PROBLEM - file %s\n",
				g_filename(gfile));
			fprintf(stderr,"** In g_unlock_views(): "
				"freetree_unregister returned error code %d.\n"
				,err);
			panic_shutdown();
		    }
		}
	    }
	}
    }

    /*
     * Sync disabled due to speed issues - we haven't knowingly seen any
     * cases of this anyway. It's only protection against when networks
     * or machines crash, in which case users blame their sysadmins rather
     * than us!
     */
    /* g_sync_on(gfile); */
    if (updates) update_header(gfile,edtime);
    /* g_sync_off(gfile); */

    gfile->check_header = 1;

    return 0;
}





/* ARGSUSED */
static int g_unlock_view(GDB *gdb, GClient c, GView v, GFlags mode)
/*
 *
 */
{
    GFile *gfile;

    /* check arguments */
    if (gdb==NULL || check_client(gdb,c) || check_view(gdb,v))
	return gerr_set(GERR_INVALID_ARGUMENTS);
    gfile = gdb->gfile;
    if (gfile==NULL)
	return gerr_set(GERR_INVALID_ARGUMENTS);

    /*
     * Special action if the file is locked
     */

    /*
     * 
     * Check the file isn't locked, and if it is then we are the person who
     * owns the lock.  If we can't perform this operation because we're
     * locked out we may want to place the client on a queue.  This isn't
     * supported yet.
     */
    if (check_file_lock(gfile,arr(View,gdb->view,v).client))
	return gerr_set(GERR_WOULD_BLOCK);

    if (gfile->flock_status == G_FLOCK_LOCKED) {
	/*
	 * check to see if this view is already pending unlocked/flushed/abandoned
	 */
	if (!(arr(View,gdb->view,v).flags &
	      (G_VIEW_UNLOCKED|G_VIEW_ABANDONED|G_VIEW_FLUSHED|G_VIEW_FREE))) {
	    /* no it isn't - add it to pending list */
	    arr(View,gdb->view,v).next = gfile->flock_view;
	    gfile->flock_view = v;
	}
	arr(View,gdb->view,v).flags |= mode;

    } else {
	arr(View,gdb->view,v).next = -1;
	arr(View,gdb->view,v).flags |= mode;
	return g_unlock_views(gdb,v);
    }

    return 0;
}


/* ARGSUSED */
int g_unlock_(GDB *gdb, GClient c, GView v)
{
    return g_unlock_view(gdb,c,v,G_VIEW_UNLOCKED);
}



/* ARGSUSED */
int g_abandon_(GDB *gdb, GClient c, GView v)
{
    return g_unlock_view(gdb,c,v,G_VIEW_ABANDONED);
}



/* ARGSUSED */
int g_flush_(GDB *gdb, GClient c, GView v)
{
    return g_unlock_view(gdb,c,v,G_VIEW_FLUSHED);
}




/* ARGSUSED */
int g_read_(GDB *gdb, GClient c, GView v, void *buf, GCardinal len)
/*
 * read an image
 */
{
    Cache *cache;

    /* check arguments */
    if (gdb==NULL || buf==NULL || len<=0 || check_client(gdb,c) || check_view(gdb,v))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    cache = &arr(View,gdb->view,v).lcache;
#ifdef USE_MMAP
    return mmap_read_image_(gdb->gfile->fdmap,
			    cache->image,
			    cache->used,
			    buf,
			    len);
#else
    return read_image_(gdb->gfile->fd,
		       cache->image,
		       cache->used,
		       buf,
		       len);
#endif
}





/* ARGSUSED */
int g_readv_(GDB *gdb, GClient c, GView v, GIOVec *vec, GCardinal vcnt)
/*
 * read an image
 */
{
    Cache *cache;
    int len;

    /* check arguments */
    if (gdb==NULL || check_GIOVec(vec,vcnt,&len) || check_client(gdb,c) || check_view(gdb,v))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    cache = &arr(View,gdb->view,v).lcache;

#ifdef USE_MMAP
    return mmap_readv_image_(gdb->gfile->fdmap,
			     cache->image,
			     cache->used,
			     vec,
			     vcnt);
#else
    return readv_image_(gdb->gfile->fd,
			cache->image,
			cache->used,
			vec,
			vcnt);
#endif
}






/* ARGSUSED */
static int update_cache_for_write(GDB *gdb, GClient c, GView v, GCardinal len,
				  int remove, Cache **cache_ret)
/*
 *
 */
{
    Cache *cache;
    View *view;
    int err;

    /* check arguments */
    /* assume already done */

    /* set cache variable */
    view = arrp(View,gdb->view,v);
    cache = &view->lcache;

    /* does the view correspond to this user? */
    /* YUK! need to do something here */

    if (! (view->flags & G_VIEW_UPDATED)) {
	/*
	 * special action on first write
	 */
	GImage image;
	GCardinal allocate;

	if (remove) {
	    image = G_NO_IMAGE;
	    allocate = 0;
	} else {
	    /* need to allocate a new image */
	    allocate = FREETREE_BLOCK(len,gdb->gfile->header.block_size);
	    if (-1 == (image = freetree_allocate(gdb->gfile->freetree,
						 allocate)))
		return get_xerrnum();
	}

	cache->image     = image;
	cache->allocated = allocate;
	cache->used      = len;

	view->flags |= G_VIEW_UPDATED;
    } else if (view->lcache.allocated < len || remove) {
	/*
	 * The allocated space is too small for the data to be written!
	 * We will need to reallocate. Note that we could use
	 * freetree_realloc here but it appears to have virtually no impact
	 * on performance, and this code is tried and tested (unlike the
	 * realloc).
	 *
	 * OR
	 *
	 * We're removing a record that we've already updated. In this case
	 * it'll have been allocated its own space already when first updated,
	 * so we need to deallocate this before resetting the cache to
	 * G_NO_IMAGE.
	 */
	GImage image;
	GCardinal allocate;

	if (cache->image != G_NO_IMAGE) {
	    /* free old image */
	    err = freetree_unregister(gdb->gfile->freetree,
				      cache->image,cache->allocated);
	    if (err) {
		gerr_set(err);
		fprintf(stderr,"** SERIOUS PROBLEM - file %s\n",
			g_filename(gdb->gfile));
		fprintf(stderr,"** In g_write_(): "
			"freetree_unregister returned error code %d.\n",err);
		panic_shutdown();
	    }
	} else {
	    /* 08/10/96 - jkb
	     * This isn't an error, but rather a check to see whether a
	     * fix has actually solved some real problems. My guess is that
	     * we never triggered the bug.
	     */
	    printf("Reusing(%d) a deleted record (corrected) - please mail jkb@mrc-lmb.cam.ac.uk\n", remove);
	}

	/* need to allocate a new image */
	if (!remove) {
	    allocate = FREETREE_BLOCK(len,gdb->gfile->header.block_size);
	    if (-1 == (image = freetree_allocate(gdb->gfile->freetree,
						 allocate)))
		return get_xerrnum();
	} else {
	    allocate = 0;
	    image = G_NO_IMAGE;
	}

	cache->image     = image;
	cache->allocated = allocate;
	cache->used      = len;
    } else {
	/*
	 * NB: This may mean that the original amount allocated here is no
	 * longer necessary. This isn't a performance impact as it'll be
	 * reused if the record grows again (whilst continuing this run), or
	 * flagged as only using as much as needed if we quit and reopen.
	 */
	cache->used = len;
    }

    /*
     * NOTE: An idea to optimise space allocation.
     * If len is significantly less than the allocated space we could
     * deallocate then reallocate, thus freeing a substantial part of the
     * previously allocated but unused space.
     *
     * Better still we could unregister the end portion only.
     */


    *cache_ret = cache;
    return 0;
}




/* ARGSUSED */
int g_write_(GDB *gdb, GClient c, GView v, void *buf, GCardinal len)
/*
 *
 */
{
    Cache *cache;
    int err;

    /* check arguments */
    if (gdb==NULL || buf==NULL || len<=0 || check_client(gdb,c) || check_view(gdb,v))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    if (gdb->gfile->check_header) {
	g_check_header(gdb->gfile);
	gdb->gfile->check_header = 0;
    }

    if (err = update_cache_for_write(gdb, c, v, len, 0, &cache))
	return err;

    /* do the write */
    return write_image_(gdb->gfile->fd,
			cache->image,
			cache->allocated,
			buf,
			len);
}



/* ARGSUSED */
int g_writev_(GDB *gdb, GClient c, GView v, GIOVec *vec, GCardinal vcnt)
/*
 *
 */
{
    Cache *cache;
    int err;
    int len;

    /* check arguments */
    if (gdb==NULL || check_GIOVec(vec,vcnt,&len) || check_client(gdb,c) || check_view(gdb,v))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    if (gdb->gfile->check_header) {
	g_check_header(gdb->gfile);
	gdb->gfile->check_header = 0;
    }

    if (err = update_cache_for_write(gdb, c, v, len, 0, &cache))
	return err;

    /* do the write */
    return writev_image_(gdb->gfile->fd,
			 cache->image,
			 cache->allocated,
			 vec,
			 vcnt);
}


/*
 * Remove the contents of a view. This marks the view as having G_NO_IMAGE
 * so that upon the next unlock it'll be unregistered and marked as a
 * G_NEW_IMAGE.
 */
/* ARGSUSED */
int g_remove_(GDB *gdb, GClient c, GView v) {
    int err;
    Cache *cache;

    /* check arguments */
    if (gdb==NULL || check_client(gdb,c) || check_view(gdb,v))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    if (gdb->gfile->check_header) {
	g_check_header(gdb->gfile);
	gdb->gfile->check_header = 0;
    }

    /*
     * Find a cache for this 'write' request. This is either our own cache
     * or a new (and now the default) one if our own was shared.
     */
    if (err = update_cache_for_write(gdb, c, v, 0, 1, &cache))
	return err;

    return 0;
}


/* ARGSUSED */
int g_fast_read_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, void *buf, GCardinal len)
/*
 *
 */
{
    int err;
    GFile *gfile;

    /* check arguments */
    if (gdb==NULL || buf==NULL || len<=0 || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    gfile = gdb->gfile;

    if (err=check_record(gfile,rec)) return err;


    /*
     * if this is an unitinitalised record we will need to
     *   initialise and allocate a record
     */
    if (arr(Index,gfile->idx,rec).flags & G_INDEX_NEW)
	/* an error only means the record has already been initialised */
	(void)initialise_record(gfile,rec);

#ifdef USE_MMAP
    return mmap_read_image_(gfile->fdmap,
			    arr(Index,gfile->idx,rec).aux_image,
			    arr(Index,gfile->idx,rec).aux_used,
			    buf,
			    len);
#else
    return read_image_(gfile->fd,
		       arr(Index,gfile->idx,rec).aux_image,
		       arr(Index,gfile->idx,rec).aux_used,
		       buf,
		       len);
#endif
}




/* ARGSUSED */
int g_fast_readv_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt)
/*
 * There is high code redundancy between this routine and g_fast_read_N_().
 * The only lines that differ are flagged with FOLLOWING LINE DIFFERS.
 */
{
    int err;
    GFile *gfile;
    /*FOLLOWING LINE DIFFERS*/
    int len;

    /* check arguments */
    /*FOLLOWING LINE DIFFERS*/
    if (gdb==NULL || check_GIOVec(vec,vcnt, &len) || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    gfile = gdb->gfile;

    if (err=check_record(gfile,rec)) return err;


    /*
     * if this is an unitinitalised record we will need to
     *   initialise and allocate a record
     */
    if (arr(Index,gfile->idx,rec).flags & G_INDEX_NEW)
	/* an error only means the record has already been initialised */
	(void)initialise_record(gfile,rec);

#ifdef USE_MMAP
    return mmap_readv_image_(gfile->fdmap,
			     arr(Index,gfile->idx,rec).aux_image,
			     arr(Index,gfile->idx,rec).aux_used,
    /*FOLLOWING LINE DIFFERS*/
			     vec,
    /*FOLLOWING LINE DIFFERS*/
			     vcnt);
#else
    return readv_image_(gfile->fd,
			arr(Index,gfile->idx,rec).aux_image,
			arr(Index,gfile->idx,rec).aux_used,
    /*FOLLOWING LINE DIFFERS*/
			vec,
    /*FOLLOWING LINE DIFFERS*/
			vcnt);
#endif
   
}










/* ARGSUSED */
int g_fast_write_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, void *buf, GCardinal len)
/*
 *
 */
{
    int err;
    GTimeStamp edtime;
    GImage image;
    GCardinal allocate;
    GFile *gfile;

    if (gdb==NULL || buf==NULL || len<= 0 || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    gfile = gdb->gfile;
    if (err=check_record(gfile,rec)) return err;


    /*
     * if this is an unitinitalised record we will need to
     *   initialise and allocate a record
     */
    if (arr(Index,gfile->idx,rec).flags & G_INDEX_NEW)
	(void)initialise_record(gfile,rec);


    /* get next edtime */
    edtime = next_edtime(gfile);


    /*
     * NOTE -
     * We don't have to recover the old image here (when no one else is using it)
     * because it is done in update_record().
     */


    /* need to allocate a new image */
    allocate = FREETREE_BLOCK(len,gfile->header.block_size);
    image = freetree_allocate(gfile->freetree,allocate);
    if (image==-1) return gerr_set(GERR_FILE_FULL);

    /* write to image */
    if (err = write_image_(gfile->fd, image, allocate, buf, len)) return err;


    /* update record */
    update_record(gfile, rec, image, allocate, len, edtime);
    
    /* update header */
    update_header(gfile, edtime);
    
    return 0;

}







/* ARGSUSED */
int g_fast_writev_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt)
/*
 * There is high code redundancy between this routine and g_fast_write_N_().
 * The only lines that differ are flagged with FOLLOWING LINE DIFFERS.
 */
{
    int err;
    GTimeStamp edtime;
    GImage image;
    GCardinal allocate;
    GFile *gfile;
    /*FOLLOWING LINE DIFFERS*/
    int len;

    /*FOLLOWING LINE DIFFERS*/
    if (gdb==NULL || check_GIOVec(vec,vcnt,&len) || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);

    gfile = gdb->gfile;
    if (err=check_record(gfile,rec)) return err;


    /*
     * if this is an unitinitalised record we will need to
     *   initialise and allocate a record
     */
    if (arr(Index,gfile->idx,rec).flags & G_INDEX_NEW)
	(void)initialise_record(gfile,rec);


    /* get next edtime */
    edtime = next_edtime(gfile);


    /*
     * NOTE -
     * We don't have to recover the old image here (when no one else is using it)
     * because it is done in update_record().
     */


    /* need to allocate a new image */
    allocate = FREETREE_BLOCK(len,gfile->header.block_size);
    image = freetree_allocate(gfile->freetree,allocate);
    if (image==-1) return gerr_set(GERR_FILE_FULL);

    /* write to image */
    /*FOLLOWING LINE DIFFERS*/
    if (err = writev_image_(gfile->fd, image, allocate, vec, vcnt)) return err;


    /* update record */
    update_record(gfile, rec, image, allocate, len, edtime);
    
    /* update header */
    update_header(gfile, edtime);
    
    return 0;

}





int g_lock_file_N_(GDB *gdb, GClient c, GFileN file_N)
/*
 * IMPORTANT NOTE:
 * The reason for locking the file is to guarentee atomicity
 * during locking and unlocking operations. Under no circumstances
 * should it be assumed that since we have a file lock we have exclusive
 * access to the file.
 */
{
    GFile *gfile;
    Client *client;

    /* check arguments */
    if (gdb==NULL || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);
    gfile = gdb->gfile;
    client = arrp(Client,gdb->client,c);

    /* check there is no outstanding lock on this file */
    if (gfile->flock_status == G_FLOCK_LOCKED) return gerr_set(GERR_WOULD_BLOCK);

    /* check we can lock this file */
    /* if (client->max_lock < G_LOCK_RW) return gerr_set(GERR_PERMISSION); */

    /* so we can lock! - initialise things */
    gfile->flock_status = G_FLOCK_LOCKED;
    gfile->flock_client = c;
    gfile->flock_view = -1;

    return 0;
}





int g_unlock_file_N_(GDB *gdb, GClient c, GFileN file_N)
/*
 *
 */
{
    int err;
    GFile *gfile;

    /* check arguments */
    if (gdb==NULL || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);
    gfile = gdb->gfile;

    /* check that we have a lock on this file */
    if (gfile->flock_client != c ||
	gfile->flock_status != G_FLOCK_LOCKED) return gerr_set(GERR_INVALID_ARGUMENTS);

    /* save error */
    err =  g_unlock_views(gdb,gfile->flock_view);
    /*
     * If the error is serious enough to do something about,
     * it would have been handled in the g_unlock_views() routine itself!!
     */

    /* reset flock variables */
    gfile->flock_status = G_FLOCK_NONE;
    gfile->flock_client = 0;
    gfile->flock_view = -1;


    /*
     * This should force propagation over NFS to be faster. We only do this
     * here as this function is called by Gap4's flush2t(io) call. We're
     * only interested in keeping the remote end consistent when Gap4
     * itself is committing a set of consistent changes to the DB.
     */

    /* LOW LEVEL IO HERE */
    fsync(gfile->fd);
    fsync(gfile->fdaux);

    return 0;
}



/* ARGSUSED */
int g_view_info_(GDB *gdb, GClient c, GView v, GViewInfo *info)
/*
 *
 */
{
    View *view;
    /* check arguments */
    if (gdb==NULL || info==NULL || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);
    view = arrp(View,gdb->view,v);

    info->image = view->lcache.image;
    info->allocated = view->lcache.allocated;
    info->used = view->lcache.used;
    info->lock = 0;

    return 0;
}


/* ARGSUSED */
int g_rec_info_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GRecInfo *info)
/*
 *
 */
{
    int err;
    GFile *gfile;
    /* check arguments */
    if (gdb==NULL || info==NULL || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);
    gfile = gdb->gfile;

    /* is this a valid record? */
    if (err=check_record(gfile,rec)) return gerr_set(err);

    if (arr(Index,gfile->idx,rec).flags & G_INDEX_NEW)
	(void)initialise_record(gfile,rec);

    info->image     = arr(Index,gfile->idx,rec).aux_image;
    info->allocated = arr(Index,gfile->idx,rec).aux_allocated;
    info->used      = arr(Index,gfile->idx,rec).aux_used;
    info->lock      = 0;

    return 0;
}



/* ARGSUSED */
int g_header_info_(GDB *gdb, GClient c, GFileN file_N, GHeaderInfo *info)
/*
 *
 */
{
    GFile *gfile;
    /* check arguments */
    if (gdb==NULL || info==NULL || check_client(gdb,c))
	return gerr_set(GERR_INVALID_ARGUMENTS);
    gfile = gdb->gfile;


    info->file_size = gfile->header.file_size;
    info->block_size = gfile->header.block_size;
    info->num_records = gfile->header.num_records;
    info->max_records = gfile->header.max_records;

    return 0;
}
