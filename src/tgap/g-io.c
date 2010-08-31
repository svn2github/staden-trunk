/*
 * File: g-io.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: low level io
 *
 * Created: 10 December 1992
 * Updated:
 *
 */

#include <staden_config.h>

#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include "g-files.h"
#include "g-io.h"


#define WRITE(X) \
    do { \
    /* LOW LEVEL IO HERE */ \
    errno = 0; \
    if ( (write(fd,&X,sizeof(X))) != sizeof(X) ) return 1; \
    } while (0)

#define READ(X) \
    do { \
    /* LOW LEVEL IO HERE */ \
    errno = 0; \
    if ( (read(fd,&X,sizeof(X))) != sizeof(X) ) return 1; \
    } while (0)



/*
 ******************************************************************************
 * Conversion functions to swap between 32-bit and 64-bit memory structures
 * for the AuxHeader.
 * The 64-bit one is used in memory data structures, but on disk we may
 * use a 32-bit version for backwards compatibility.
 ******************************************************************************
 */
static void header_32to64(AuxHeader32 *rec32, AuxHeader *rec64) {
    rec64->file_size   = rec32->file_size;
    rec64->block_size  = rec32->block_size;
    rec64->num_records = rec32->num_records;
    rec64->max_records = rec32->max_records;
    rec64->last_time   = rec32->last_time;
    rec64->flags       = rec32->flags;
    rec64->spare1      = rec32->spare1;
    rec64->free_time   = rec32->free_time;
    rec64->free_record = rec32->free_record;
    rec64->spare[0]    = rec32->spare[0];
    rec64->spare[1]    = rec32->spare[1];
    rec64->spare[2]    = rec32->spare[2];
    rec64->spare[3]    = rec32->spare[3];
    rec64->spare[4]    = rec32->spare[4];
    rec64->spare[5]    = rec32->spare[5];
    rec64->format      = rec32->format;
}

static void header_64to32(AuxHeader *rec64, AuxHeader32 *rec32) {
    rec32->file_size   = rec64->file_size;
    rec32->block_size  = rec64->block_size;
    rec32->num_records = rec64->num_records;
    rec32->max_records = rec64->max_records;
    rec32->last_time   = rec64->last_time;
    rec32->flags       = rec64->flags;
    rec32->spare1      = rec64->spare1;
    rec32->free_time   = rec64->free_time;
    rec32->free_record = rec64->free_record;
    rec32->spare[0]    = rec64->spare[0];
    rec32->spare[1]    = rec64->spare[1];
    rec32->spare[2]    = rec64->spare[2];
    rec32->spare[3]    = rec64->spare[3];
    rec32->spare[4]    = rec64->spare[4];
    rec32->spare[5]    = rec64->spare[5];
    rec32->spare[6]    = 0;
    rec32->format      = rec64->format;
}


/*************************************************************
 * Read and write an aux header record
 *************************************************************/


int write_aux_header32_(int fd, void *recv, int num)
{
    AuxHeader *rec = recv;
    AuxHeader32 rec32;

    /* Convert 64-bit memory structure to 32-bit if needed */
    header_64to32(rec, &rec32);

    /* LOW LEVEL IO HERE */
    errno = 0;
    if ((write(fd, &rec32, sizeof(AuxHeader))) != sizeof(AuxHeader))
	return 1;

    return 0;
}


int write_aux_header64_(int fd, void *recv, int num)
{
    AuxHeader *rec = recv;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if ((write(fd, rec, sizeof(AuxHeader))) != sizeof(AuxHeader))
	return 1;

    return 0;
}


int write_aux_header_swapped32_(int fd, void *headerv, int num)
{
    AuxHeader *rec = headerv;
    AuxHeader32 rec32;

    if (rec->format != G_32BIT) {
	fprintf(stderr, "** Expected 32-bit file size data; not found\n");
	return 1;
    }

    /* Convert to 32-bit format */
    header_64to32(rec, &rec32);
    
    /* byte-swap */
    swap_int4(rec32.file_size,        rec32.file_size);
    swap_GCardinal(rec32.block_size,  rec32.block_size);
    swap_GCardinal(rec32.num_records, rec32.num_records);
    swap_GCardinal(rec32.max_records, rec32.max_records);
    swap_GTimeStamp(rec32.last_time,  rec32.last_time);
    swap_GHFlags(rec32.flags,         rec32.flags);
    swap_GHFlags(rec32.spare1,        rec32.spare1);
    swap_GTimeStamp(rec32.free_time,  rec32.free_time);
    swap_GCardinal(rec32.free_record, rec32.free_record);
    swap_int4(rec32.spare[0],         rec32.spare[0]);
    swap_int4(rec32.spare[1],         rec32.spare[1]);
    swap_int4(rec32.spare[2],         rec32.spare[2]);
    swap_int4(rec32.spare[3],         rec32.spare[3]);
    swap_int4(rec32.spare[4],         rec32.spare[4]);
    swap_int4(rec32.spare[5],         rec32.spare[5]);
    swap_int4(rec32.spare[6],         rec32.spare[6]);
    swap_int4(rec32.format,           rec32.format);

    /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, &rec32, sizeof(rec32)) != sizeof(rec32);
}


int write_aux_header_swapped64_(int fd, void *headerv, int num)
{
    AuxHeader *header = headerv;
    AuxHeader swapped;

    if (header->format != G_64BIT) {
	fprintf(stderr, "** Expected 64-bit file size data; not found\n");
	return 1;
    }

    swap_GImage(header->file_size,swapped.file_size);
    swap_GCardinal(header->block_size,swapped.block_size);
    swap_GCardinal(header->num_records,swapped.num_records);
    swap_GCardinal(header->max_records,swapped.max_records);
    swap_GTimeStamp(header->last_time,swapped.last_time);
    swap_GHFlags(header->flags,swapped.flags);
    swap_GHFlags(header->spare1,swapped.spare1);
    swap_GTimeStamp(header->free_time,swapped.free_time);
    swap_GCardinal(header->free_record, swapped.free_record);
    swap_int4(header->spare[0],swapped.spare[0]);
    swap_int4(header->spare[1],swapped.spare[1]);
    swap_int4(header->spare[2],swapped.spare[2]);
    swap_int4(header->spare[3],swapped.spare[3]);
    swap_int4(header->spare[4],swapped.spare[4]);
    swap_int4(header->spare[5],swapped.spare[5]);
    swap_int4(header->format,swapped.format);

    errno = 0;
    return write(fd, &swapped, sizeof(swapped)) != sizeof(swapped);
}

/*
 * Read the header in either 32-bit or 64-bit
 * Auto-senses 32-bit and 64-bit data to allow us to work out which format
 * g-database we are reading.
 */
int read_aux_header_(int fd, void *recv, int num)
{
    AuxHeader rec;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (sizeof(rec) != read(fd, &rec, sizeof(rec)))
	return 1;

    if (rec.format == G_32BIT) {
	/* Convert 32-bit layout to an in-memory 64-bit one */
	header_32to64((AuxHeader32 *)&rec, (AuxHeader *)recv);
    } else {
	/* Already a 64-bit data layout, so just copy */
	memcpy(recv, &rec, sizeof(rec));
    }

    return 0;
}


int read_aux_header_swapped_(int fd, void *headerv, int num)
{
    AuxHeader rec;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (sizeof(rec) != read(fd, &rec, sizeof(rec)))
	return 1;

    swap_int4(rec.format, rec.format);
    if (rec.format == G_32BIT) {
	/* byte-swap and remap to 64-bit at the same time */
	AuxHeader32 *rec32 = (AuxHeader32 *)&rec;
	AuxHeader32 swapped;
	int4 file_size;
	swap_int4(rec32->file_size, file_size);
	swapped.file_size = file_size;
	swap_GCardinal(rec32->block_size,  swapped.block_size);
	swap_GCardinal(rec32->num_records, swapped.num_records);
	swap_GCardinal(rec32->max_records, swapped.max_records);
	swap_GTimeStamp(rec32->last_time,  swapped.last_time);
	swap_GHFlags(rec32->flags,         swapped.flags);
	swap_GHFlags(rec32->spare1,        swapped.spare1);
	swap_GTimeStamp(rec32->free_time,  swapped.free_time);
	swap_GCardinal(rec32->free_record, swapped.free_record);
	swap_int4(rec32->spare[0],         swapped.spare[0]);
	swap_int4(rec32->spare[1],         swapped.spare[1]);
	swap_int4(rec32->spare[2],         swapped.spare[2]);
	swap_int4(rec32->spare[3],         swapped.spare[3]);
	swap_int4(rec32->spare[4],         swapped.spare[4]);
	swap_int4(rec32->spare[5],         swapped.spare[5]);
	swap_int4(rec32->spare[6],         swapped.spare[6]);
	swapped.format = G_32BIT;

	header_32to64(&swapped, (AuxHeader *)headerv);
    } else {
	swap_GImage(rec.file_size,      rec.file_size);
	swap_GCardinal(rec.block_size,  rec.block_size);
	swap_GCardinal(rec.num_records, rec.num_records);
	swap_GCardinal(rec.max_records, rec.max_records);
	swap_GTimeStamp(rec.last_time,  rec.last_time);
	swap_GHFlags(rec.flags,         rec.flags);
	swap_GHFlags(rec.spare1,        rec.spare1);
	swap_GTimeStamp(rec.free_time,  rec.free_time);
	swap_GCardinal(rec.free_record, rec.free_record);
	swap_int4(rec.spare[0],         rec.spare[0]);
	swap_int4(rec.spare[1],         rec.spare[1]);
	swap_int4(rec.spare[2],         rec.spare[2]);
	swap_int4(rec.spare[3],         rec.spare[3]);
	swap_int4(rec.spare[4],         rec.spare[4]);
	swap_int4(rec.spare[5],         rec.spare[5]);

	memcpy(headerv, &rec, sizeof(AuxHeader));
    }

    return 0;
}


/*************************************************************
 * Read and write an aux index entry
 *************************************************************/


int write_aux_index32_(int fd, void *recv, int num)
{
    AuxIndex *rec = recv;
    AuxIndex32 rec32;

    /* Convert from 64-bit to 32-bit based data structure */
    rec32.image[0] = rec->image[0];
    rec32.image[1] = rec->image[1];
    rec32.time[0] = rec->time[0];
    rec32.time[1] = rec->time[1];
    rec32.used[0] = rec->used[0];
    rec32.used[1] = rec->used[1];

   /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, &rec32, sizeof(rec32)) != sizeof(rec32);
}


int write_aux_index64_(int fd, void *recv, int num)
{
    AuxIndex *rec = recv;
   /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, rec, sizeof(*rec)) != sizeof(*rec);
}




int write_aux_index_swapped32_(int fd, void *idxv, int num)
{
    AuxIndex *idx = idxv;
    AuxIndex32 rec32;

    /* Convert from 64-bit to 32-bit based data structure */
    rec32.image[0] = idx->image[0];
    rec32.image[1] = idx->image[1];
    rec32.time[0]  = idx->time[0];
    rec32.time[1]  = idx->time[1];
    rec32.used[0]  = idx->used[0];
    rec32.used[1]  = idx->used[1];

    /* Byte-swap */
    swap_int4(rec32.image[0],      rec32.image[0]);
    swap_int4(rec32.image[1],      rec32.image[1]);
    swap_GTimeStamp(rec32.time[0], rec32.time[0]);
    swap_GTimeStamp(rec32.time[1], rec32.time[1]);
    swap_GCardinal(rec32.used[0],  rec32.used[0]);
    swap_GCardinal(rec32.used[1],  rec32.used[1]);

   /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, &rec32, sizeof(rec32)) != sizeof(rec32);
}

int write_aux_index_swapped64_(int fd, void *idxv, int num)
{
    AuxIndex *idx = idxv;
    AuxIndex swapped;

    /* Byte-swap */
    swap_GImage(idx->image[0],swapped.image[0]);
    swap_GImage(idx->image[1],swapped.image[1]);
    swap_GTimeStamp(idx->time[0],swapped.time[0]);
    swap_GTimeStamp(idx->time[1],swapped.time[1]);
    swap_GCardinal(idx->used[0],swapped.used[0]);
    swap_GCardinal(idx->used[1],swapped.used[1]);

    /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, &swapped, sizeof(swapped)) != sizeof(swapped);
}




int read_aux_index32_(int fd, void *recv, int num)
{
    AuxIndex *idx = recv;
    AuxIndex32 rec32;
    int i;

    for (i = 0; i < num; i++) {
	/* LOW LEVEL IO HERE */
	errno = 0;
	if (read(fd, &rec32, sizeof(rec32)) != (int)(sizeof(rec32)))
	    return i;

	/* Convert 32-bit struct to the 64-bit one used in memory */
	idx[i].image[0] = rec32.image[0];
	idx[i].image[1] = rec32.image[1];
	idx[i].time[0] = rec32.time[0];
	idx[i].time[1] = rec32.time[1];
	idx[i].used[0] = rec32.used[0];
	idx[i].used[1] = rec32.used[1];
    }

    return num;
}


int read_aux_index64_(int fd, void *recv, int num)
{
    AuxIndex *idx = recv;
    int ret;
    
    /* LOW LEVEL IO HERE */
    errno = 0;
    ret = read(fd, idx, sizeof(*idx)*num);
    return (int)(ret / sizeof(*idx));
}




int read_aux_index_swapped32_(int fd, void *idxv, int num)
{
    AuxIndex *idx = idxv;
    AuxIndex32 rec32;
    int i;

    for (i = 0; i < num; i++) {
	int4 image;

	/* LOW LEVEL IO HERE */
	errno = 0;
	if (read(fd, &rec32, sizeof(rec32)) != (int)(sizeof(rec32)))
	    return i;

	swap_int4(rec32.image[0], image);
	idx->image[0] = image;
	swap_int4(rec32.image[1], image);
	idx->image[1] = image;
	swap_GTimeStamp(rec32.time[0],idx->time[0]);
	swap_GTimeStamp(rec32.time[1],idx->time[1]);
	swap_GCardinal(rec32.used[0],idx->used[0]);
	swap_GCardinal(rec32.used[1],idx->used[1]);
	
	idx++;
    }

    return num;
}


int read_aux_index_swapped64_(int fd, void *idxv, int num)
{
    AuxIndex *idx = idxv;
    int i, ret;
    AuxIndex *swapped;

    /* LOW LEVEL IO HERE */
    errno = 0;
    ret = read(fd, idx, sizeof(*idx)*num);
    num = (int)(ret / sizeof(*idx));

    for (i = 0; i < num; i++) {
	swapped = &idx[i];

	swap_GImage(swapped->image[0],swapped->image[0]);
	swap_GImage(swapped->image[1],swapped->image[1]);
	swap_GTimeStamp(swapped->time[0],swapped->time[0]);
	swap_GTimeStamp(swapped->time[1],swapped->time[1]);
	swap_GCardinal(swapped->used[0],swapped->used[0]);
	swap_GCardinal(swapped->used[1],swapped->used[1]);
    }

    return num;
}


int seek_aux_index32_(int fd, void *dummy, int rec_num) {
    return lseek(fd, sizeof(AuxHeader) + rec_num * sizeof(AuxIndex32),
		 SEEK_SET);
}

int seek_aux_index64_(int fd, void *dummy, int rec_num) {
    return lseek(fd, sizeof(AuxHeader) + rec_num * sizeof(AuxIndex),
		 SEEK_SET);
}


/*************************************************************
 * Set swapped vector
 *************************************************************/

int (*low_level_vectors32[5])(int fd, void *x, int num) = {
    write_aux_header32_,
    write_aux_index32_,
    read_aux_header_,
    read_aux_index32_,
    seek_aux_index32_
    };

int (*low_level_vectors64[5])(int fd, void *x, int num) = {
    write_aux_header64_,
    write_aux_index64_,
    read_aux_header_,
    read_aux_index64_,
    seek_aux_index64_
    };

int (*low_level_vectors_swapped32[5])(int fd, void *x, int num) = {
    write_aux_header_swapped32_,
    write_aux_index_swapped32_,
    read_aux_header_swapped_,
    read_aux_index_swapped32_,
    seek_aux_index32_
    };

int (*low_level_vectors_swapped64[5])(int fd, void *x, int num) = {
    write_aux_header_swapped64_,
    write_aux_index_swapped64_,
    read_aux_header_swapped_,
    read_aux_index_swapped64_,
    seek_aux_index64_
    };
