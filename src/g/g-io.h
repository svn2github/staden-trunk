/*
 * File: g-io.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for low level io
 *
 * Created:
 * Updated:
 *
 */

#ifndef _G_IO_H_
#define _G_IO_H_

#include "g-filedefs.h"

#ifdef _MSC_VER
#  ifdef BUILDING_G_DLL
#    define G_EXPORT __declspec(dllexport)
#  else
#    define G_EXPORT __declspec(dllimport)
#  endif
#else
#  define G_EXPORT
#endif



extern int write_aux_header_(int fd, void *rec, int num);
extern int write_aux_header_swapped_(int fd, void *rec, int num);
extern int read_aux_header_(int fd, void *header, int num);
extern int read_aux_header_swapped_(int fd, void *header, int num);
extern int write_aux_index_(int fd, void *rec, int num);
extern int write_aux_index_swapped_(int fd, void *idx, int num);
extern int read_aux_index_(int fd, void *rec, int num);
extern int read_aux_index_swapped_(int fd, void *idx, int num);


extern G_EXPORT int (*low_level_vectors[4])(int fd, void *x, int num);
extern G_EXPORT int (*low_level_vectors_swapped[4])(int fd, void *x, int num);
extern G_EXPORT int (*(*low_level_vector))(int fd, void *x, int num);
extern G_EXPORT int (*(*set_low_level_vector(void)))(int fd, void *x, int num);

#define GOP_WRITE_AUX_HEADER 0
#define GOP_WRITE_AUX_INDEX  1
#define GOP_READ_AUX_HEADER  2
#define GOP_READ_AUX_INDEX   3

#endif /*_G_IO_H_*/
