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

#include <unistd.h>
#include <errno.h>
#include "g-files.h"
#include "g-io.h"


#define WRITE(X) \
    do { \
    /* LOW LEVEL IO HERE */ \
    errno = 0; \
    if ( (write(fd,&rec->X,sizeof(rec->X))) != sizeof(rec->X) ) return 1; \
    } while (0)

#define READ(X) \
    do { \
    /* LOW LEVEL IO HERE */ \
    errno = 0; \
    if ( (read(fd,&rec->X,sizeof(rec->X))) != sizeof(rec->X) ) return 1; \
    } while (0)



/*************************************************************
 * Read and write an aux header record
 *************************************************************/


int write_aux_header_(int fd, void *recv, int num)
/*
 *
 */
{
    AuxHeader *rec = recv;

    /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, rec, sizeof(*rec)) != sizeof(*rec);

#if 0
    WRITE(file_size);
    WRITE(block_size);
    WRITE(num_records);
    WRITE(max_records);
    WRITE(last_time);
    WRITE(flags);
    WRITE(spare1);
    WRITE(free_time);
    WRITE(spare[0]);
    WRITE(spare[1]);
    WRITE(spare[2]);
    WRITE(spare[3]);
    WRITE(spare[4]);
    WRITE(spare[5]);
    WRITE(spare[6]);
    WRITE(spare[7]);
    WRITE(spare[8]);

    return 0;
#endif
}




int write_aux_header_swapped_(int fd, void *headerv, int num)
/*
 *
 */
{
    AuxHeader *header = headerv;
    AuxHeader swapped;

    swap_GCardinal(header->file_size,swapped.file_size);
    swap_GCardinal(header->block_size,swapped.block_size);
    swap_GCardinal(header->num_records,swapped.num_records);
    swap_GCardinal(header->max_records,swapped.max_records);
    swap_GTimeStamp(header->last_time,swapped.last_time);
    swap_GHFlags(header->flags,swapped.flags);
    swap_GHFlags(header->spare1,swapped.spare1);
    swap_GTimeStamp(header->free_time,swapped.free_time);
    swap_int4(header->spare[0],swapped.spare[0]);
    swap_int4(header->spare[1],swapped.spare[1]);
    swap_int4(header->spare[2],swapped.spare[2]);
    swap_int4(header->spare[3],swapped.spare[3]);
    swap_int4(header->spare[4],swapped.spare[4]);
    swap_int4(header->spare[5],swapped.spare[5]);
    swap_int4(header->spare[6],swapped.spare[6]);
    swap_int4(header->spare[7],swapped.spare[7]);
    swap_int4(header->spare[8],swapped.spare[8]);

    return write_aux_header_(fd, &swapped, num);

}






int read_aux_header_(int fd, void *recv, int num)
/*
 *
 */
{
    AuxHeader *rec = recv;
    /* LOW LEVEL IO HERE */
    errno = 0;
    return read(fd, rec, sizeof(*rec)) != sizeof(*rec);

#if 0
    READ(file_size);
    READ(block_size);
    READ(num_records);
    READ(max_records);
    READ(last_time);
    READ(flags);
    READ(spare1);
    READ(free_time);
    READ(spare[0]);
    READ(spare[1]);
    READ(spare[2]);
    READ(spare[3]);
    READ(spare[4]);
    READ(spare[5]);
    READ(spare[6]);
    READ(spare[7]);
    READ(spare[8]);

    return 0;
#endif
}




int read_aux_header_swapped_(int fd, void *headerv, int num)
/*
 *
 */
{
    AuxHeader *header = headerv;
    int err;
    AuxHeader swapped;

    if (err=read_aux_header_(fd, &swapped, 1)) return err;

    swap_GCardinal(swapped.file_size,header->file_size);
    swap_GCardinal(swapped.block_size,header->block_size);
    swap_GCardinal(swapped.num_records,header->num_records);
    swap_GCardinal(swapped.max_records,header->max_records);
    swap_GTimeStamp(swapped.last_time,header->last_time);
    swap_GHFlags(swapped.flags,header->flags);
    swap_GHFlags(swapped.spare1,header->spare1);
    swap_GTimeStamp(swapped.free_time,header->free_time);
    swap_int4(swapped.spare[0],header->spare[0]);
    swap_int4(swapped.spare[1],header->spare[1]);
    swap_int4(swapped.spare[2],header->spare[2]);
    swap_int4(swapped.spare[3],header->spare[3]);
    swap_int4(swapped.spare[4],header->spare[4]);
    swap_int4(swapped.spare[5],header->spare[5]);
    swap_int4(swapped.spare[6],header->spare[6]);
    swap_int4(swapped.spare[7],header->spare[7]);
    swap_int4(swapped.spare[8],header->spare[8]);

    return 0;	

}



/*************************************************************
 * Read and write an aux index entry
 *************************************************************/


int write_aux_index_(int fd, void *recv, int num)
/*
 *
 */
{
    AuxIndex *rec = recv;
   /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, rec, sizeof(*rec)) != sizeof(*rec);

#if 0
    WRITE(image[0]);
    WRITE(image[1]);
    WRITE(time[0]);
    WRITE(time[1]);
    WRITE(used[0]);
    WRITE(used[1]);

    return 0;
#endif
}




int write_aux_index_swapped_(int fd, void *idxv, int num)
/*
 *
 */
{
    AuxIndex *idx = idxv;
    AuxIndex swapped;

    swap_GImage(idx->image[0],swapped.image[0]);
    swap_GImage(idx->image[1],swapped.image[1]);
    swap_GTimeStamp(idx->time[0],swapped.time[0]);
    swap_GTimeStamp(idx->time[1],swapped.time[1]);
    swap_GCardinal(idx->used[0],swapped.used[0]);
    swap_GCardinal(idx->used[1],swapped.used[1]);

/*    return write_aux_index_(fd, &swapped); */
    /* LOW LEVEL IO HERE */
    errno = 0;
    return write(fd, &swapped, sizeof(swapped)) != sizeof(swapped);
}




int read_aux_index_(int fd, void *recv, int num)
/*
 *
 */
{
    AuxIndex *idx = recv;
    /* LOW LEVEL IO HERE */
    errno = 0;
    return read(fd, idx, sizeof(*idx)*num) != (int)(sizeof(*idx)*num);
}




int read_aux_index_swapped_(int fd, void *idxv, int num)
/*
 *
 */
{
    AuxIndex *idx = idxv;
    int err;

    if (num == 1) {
	AuxIndex swapped;

	/* LOW LEVEL IO HERE */
	errno = 0;
	if (err = (read(fd, &swapped, sizeof(swapped)) != sizeof(swapped)))
	    return err;

	/* if (err=read_aux_index_(fd, &swapped, num)) return err; */
	
	swap_GImage(swapped.image[0],idx->image[0]);
	swap_GImage(swapped.image[1],idx->image[1]);
	swap_GTimeStamp(swapped.time[0],idx->time[0]);
	swap_GTimeStamp(swapped.time[1],idx->time[1]);
	swap_GCardinal(swapped.used[0],idx->used[0]);
	swap_GCardinal(swapped.used[1],idx->used[1]);
    } else {
	int i;
	AuxIndex *swapped;

	/* LOW LEVEL IO HERE */
	errno = 0;
	if (err = (read(fd, idx, sizeof(*idx)*num) != (int)(sizeof(*idx)*num)))
	    return err;

	for (i = 0; i < num; i++) {
	    swapped = &idx[i];

	    swap_GImage(swapped->image[0],swapped->image[0]);
	    swap_GImage(swapped->image[1],swapped->image[1]);
	    swap_GTimeStamp(swapped->time[0],swapped->time[0]);
	    swap_GTimeStamp(swapped->time[1],swapped->time[1]);
	    swap_GCardinal(swapped->used[0],swapped->used[0]);
	    swap_GCardinal(swapped->used[1],swapped->used[1]);
	}
    }

    return 0;

}




/*************************************************************
 * Set swapped vector
 *************************************************************/

int (*low_level_vectors[4])(int fd, void *x, int num) = {
    write_aux_header_,
    write_aux_index_,
    read_aux_header_,
    read_aux_index_
    };

int (*low_level_vectors_swapped[4])(int fd, void *x, int num) = {
    write_aux_header_swapped_,
    write_aux_index_swapped_,
    read_aux_header_swapped_,
    read_aux_index_swapped_
    };

/* 6/1/99 johnt - need to assign globals to force export under Visual C++ */
int (*(*low_level_vector))(int fd, void *x, int num) = NULL;



int (*(*set_low_level_vector(void)))(int fd, void *x, int num)
{
    int indian = 1;
    if ( *(char *)&indian )
	low_level_vector = low_level_vectors_swapped;
    else
	low_level_vector = low_level_vectors;
    
    return low_level_vector;
}

