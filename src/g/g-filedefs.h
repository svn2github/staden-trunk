/*
 * File: g-filedefs.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: file structures for server
 *
 * Created: prior to 18-Sep-1992
 * Updated:
 *
 */

#ifndef _G_FILEDEFS_H_
#define _G_FILEDEFS_H_

#include "g-os.h"

/* aux file header */



typedef struct {
    GCardinal file_size;	/* size of file in bytes */
    GCardinal block_size;       /* size of each block (for allocation only) */
    GCardinal num_records;	/* number of records in record file */
    GCardinal max_records;	/* max number of records in record file */
    GTimeStamp last_time;	/* last update time */
    GHFlags flags;		/* flags */
    GHFlags spare1;		/* a spare short */
    GTimeStamp free_time;	/* time of the last freetree save */
    int4 spare[9];		/* for later use */
} AuxHeader;	

#define G_BLOCK_SIZE_CHANGED (GFlags) 1

typedef struct {
    GImage     image[2];	/* offset in bytes */
    GTimeStamp time[2];		/* time last updated */
    GCardinal  used[2];		/* bytes used */
} AuxIndex;


#endif /*_G_FILEDEFS_H_*/




