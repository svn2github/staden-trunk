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

/*
 * Aux file header.
 * Note that both the 64-bit (AuxHeader) and the old 32-bit (AuxHeader32)
 * format have the last 4-byte word as "format". This is how we determine
 * the layout of the rest of structure.
 *
 * In memory, the g-library always usesthe 64-bit version, but on disk
 * we have to convert to and from the 32-bit version for backwards
 * compatibility.
 */
typedef struct {
    GImage    file_size;	/* size of file in bytes */
    GCardinal block_size;       /* size of each block (for allocation only) */
    GCardinal num_records;	/* number of records in record file */
    GCardinal max_records;	/* max number of records in record file */
    GTimeStamp last_time;	/* last update time */
    GHFlags flags;		/* flags */
    GHFlags spare1;		/* a spare short */
    GTimeStamp free_time;	/* time of the last freetree save */
    int4 spare[7];		/* for later use */
    int4 format;		/* 0=32-bit file_size, 1=64-bit file_size */
} AuxHeader;	

/* AuxHeader with 32-bit GImage size */
typedef struct {
    int4      file_size;	/* size of file in bytes */
    GCardinal block_size;       /* size of each block (for allocation only) */
    GCardinal num_records;	/* number of records in record file */
    GCardinal max_records;	/* max number of records in record file */
    GTimeStamp last_time;	/* last update time */
    GHFlags flags;		/* flags */
    GHFlags spare1;		/* a spare short */
    GTimeStamp free_time;	/* time of the last freetree save */
    int4 spare[8];		/* for later use */
    int4 format;		/* 0=32-bit file_size, 1=64-bit file_size */
} AuxHeader32;	

#define G_BLOCK_SIZE_CHANGED (GFlags) 1

/* AuxIndex with the default (64-bit) GImage size */
typedef struct {
    GImage     image[2];	/* offset in bytes */
    GTimeStamp time[2];		/* time last updated */
    GCardinal  used[2];		/* bytes used */
} AuxIndex;

/* AuxIndex with a 32-bit GImage size */
typedef struct {
    int4       image[2];	/* offset in bytes */
    GTimeStamp time[2];		/* time last updated */
    GCardinal  used[2];		/* bytes used */
} AuxIndex32;

#endif /*_G_FILEDEFS_H_*/




