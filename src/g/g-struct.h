/*
 * File: g-struct.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for alloc and dealloc routines
 *
 * Created: 14 October 1992
 * Updated:
 *
 */
#ifndef _G_STRUCT_H_
#define _G_STRUCT_H_

#include "freetree.h"

#include "g-filedefs.h"
#include "g-os.h"
#include "array.h"


/***************************************************************/

/*
 * internal memory data structures
 */



/*
 * Structures pertaining to individual files in the database
 */


#define G_INDEX_NONE (GFlags) 0
#define G_INDEX_NEW  (GFlags) 1
#define G_INDEX_USED (GFlags) 2





typedef struct _cache_ {
    /*8*/ GImage image;		/* file offset */
    /*4*/ GCardinal allocated;	/* bytes allocated */
    /*4*/ GCardinal used;	/* bytes used */
    /*4*/ GCardinal rec;	/* index entry for this cache entry */
} Cache;



typedef struct {
#if 0
    /*4*/ GImage     loc_image;
    /*4*/ GTimeStamp loc_time;
    /*4*/ GCardinal  loc_used;
    /*4*/ GCardinal  loc_allocated;
#endif

    /*8*/ GImage     aux_image;
    /*4*/ GTimeStamp aux_time;
    /*4*/ GCardinal  aux_used;
    /*4*/ GCardinal  aux_allocated;	/* bytes allocated to record*/

    /*1*/ GFlags flags;
} Index;




typedef struct {
    GImage image;
    GCardinal allocated;
    GCardinal used;
    GLock lock;
} GViewInfo, GRecInfo;



typedef struct {
    GImage    file_size;	/* size of file in bytes */
    GCardinal block_size;       /* size of each block (for allocation only) */
    GCardinal num_records;	/* number of records in record file */
    GCardinal max_records;	/* max number of records in record file */
} GHeaderInfo;

typedef struct {
    void *buf;
    GCardinal len;
} GIOVec;			/* for g_readv, g_writev */





typedef struct _gfile_ {
    char *fname; /* name of file */
    int fd;
    int fdaux;
    AuxHeader header;
    /*
     * list of free space in file
     */
    free_tree *freetree;

    /*
     * Mapping of record numbers to images
     */
    GCardinal Nidx;
    Array idx;
    /*
     * For file locking
     */
    int flock_status; /* lock status; 0=no lock; 1=locked */
    GClient flock_client; /* client who set lock */
    GCardinal flock_view; /* views for existing file lock */

    /*
     * For checking external processes have no performed an update.
     */
    int check_header;

    /*
     * For mmaped IO
     */
    char *fdmap;
    char *fdauxmap;

    /* Mapping for which file reading functions to use */
    int (*(*low_level_vector))(int fd, void *x, int num);
    int swapped; /* true => byte-swapping is needed */
} GFile;


#define G_FLOCK_NONE 0
#define G_FLOCK_LOCKED 1


/*
 * Structures pertaining to whole database
 */

typedef struct _client_ {
    /* useful information pertaining to clients */
    int id;
    GLock max_lock; /* privilages */
} Client;


typedef struct _view_ {
    /*16*/ Cache lcache;
    /*4*/  GCardinal next;
    /*2*/  GClient client;
    /*1*/  GFlags flags;
} View;


typedef struct _gdb_ {
    /*
     * file data structure
     */
    GFile *gfile;
    /*
     * client data structure (YUK! should be Array)
     */
    Array client;
    GCardinal Nclient; /* maximum number of clients */
    /*
     * view data structure
     *
     * NOTE about Views
     *   Views are allocated as an array because we want the user
     *   to be able to reference them using an array index
     */
    Array view;
    GView Nview;
    /**/
    GView free_view;
    GCardinal ConnectedClients; /* number of connected clients */
} GDB;





/*
 * Macros for View.flags
 */
#define G_VIEW_NEW       (GFlags)(0)     /* initial state */
#define G_VIEW_USED      (GFlags)(1<<0)
#define G_VIEW_FREE      (GFlags)(1<<1)
#define G_VIEW_UPDATED   (GFlags)(1<<2)	/* view has been modified */
#define G_VIEW_ABANDONED (GFlags)(1<<3)	/* view has been requested to be abandoned */
#define G_VIEW_UNLOCKED  (GFlags)(1<<4) /* view has been requested to be unlocked  */
#define G_VIEW_FLUSHED   (GFlags)(1<<5) /* view has been requested to be flushed   */


    


/*************************************************************
 * function declarations
 *************************************************************/

extern GFile *g_new_gfile(int bitsize);
/*
 * create and initialise a new gfile structure
 */

extern void g_free_gfile(GFile *gfile);
/*
 * free gfile structure
 */



extern GDB *g_new_gdb(void);
/*
 * create and initialise a new gdb structure
 */


extern void g_free_gdb(GDB *gdb);
/*
 * free gdb structure
 */


extern GView g_new_view(GDB *gdb);
/*
 * allocate a new view
 */


extern void g_free_view(GDB *gdb, GView v);



#endif /*_G_STRUCT_H_*/
