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



extern int write_aux_header32_(int fd, void *rec, int num);
extern int write_aux_header_swapped32_(int fd, void *rec, int num);
extern int read_aux_header32_(int fd, void *header, int num);
extern int read_aux_header_swapped32_(int fd, void *header, int num);
extern int write_aux_index32_(int fd, void *rec, int num);
extern int write_aux_index_swapped32_(int fd, void *idx, int num);
extern int read_aux_index32_(int fd, void *rec, int num);
extern int read_aux_index_swapped32_(int fd, void *idx, int num);

extern int write_aux_header64_(int fd, void *rec, int num);
extern int write_aux_header_swapped64_(int fd, void *rec, int num);
extern int read_aux_header64_(int fd, void *header, int num);
extern int read_aux_header_swapped64_(int fd, void *header, int num);
extern int write_aux_index64_(int fd, void *rec, int num);
extern int write_aux_index_swapped64_(int fd, void *idx, int num);
extern int read_aux_index64_(int fd, void *rec, int num);
extern int read_aux_index_swapped64_(int fd, void *idx, int num);


extern G_EXPORT int (*low_level_vectors_swapped64[4])(int fd, void *x, int num);
extern G_EXPORT int (*low_level_vectors_swapped32[4])(int fd, void *x, int num);
extern G_EXPORT int (*low_level_vectors64[4])(int fd, void *x, int num);
extern G_EXPORT int (*low_level_vectors32[4])(int fd, void *x, int num);

#define GOP_WRITE_AUX_HEADER 0
#define GOP_WRITE_AUX_INDEX  1
#define GOP_READ_AUX_HEADER  2
#define GOP_READ_AUX_INDEX   3

#endif /*_G_IO_H_*/
