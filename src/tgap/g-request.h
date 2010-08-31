/*
 * File: g_request.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for g_request.c
 *
 * Created: 18-Sep-1992
 * Updated:	   
 *
 */

#ifndef _G_REQUEST_H_
#define _G_REQUEST_H_

#include "g-struct.h"
#include "g-os.h"

/* file locking */
extern int g_lock_file_N_(GDB *gdb, GClient c, GFileN file_N);
extern int g_unlock_file_N_(GDB *gdb, GClient c, GFileN file_N);
/* record locking */
extern GView g_lock_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GLock lock);
extern int g_upgrade_(GDB *gdb, GClient c, GView v,GLock lock);
extern int g_unlock_(GDB *gdb, GClient c, GView v);
extern int g_abandon_(GDB *gdb, GClient c, GView v);
extern int g_flush_(GDB *gdb, GClient c, GView v);
/* reading */
extern int g_read_(GDB *gdb, GClient c, GView v, void *buf,GCardinal len);
extern int g_readv_(GDB *gdb, GClient c, GView v, GIOVec *vec, GCardinal vcnt);
extern int g_fast_read_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, void *buf, GCardinal len);
extern int g_fast_readv_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt);
/* writing */
extern int g_write_(GDB *gdb, GClient c, GView v, void *buf,GCardinal len);
extern int g_writev_(GDB *gdb, GClient c, GView v, GIOVec *vec, GCardinal vcnt);
extern int g_fast_write_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, void *buf, GCardinal len);
extern int g_fast_writev_N_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt);
/* information */
extern int g_remove_(GDB *gdb, GClient c, GView v);
extern int g_view_info_(GDB *gdb, GClient c, GView v, GViewInfo *info);
extern int g_rec_info_(GDB *gdb, GClient c, GFileN file_N, GCardinal rec, GRecInfo *info);
extern int g_header_info_(GDB *gdb, GClient c, GFileN file_N, GHeaderInfo *info);
extern int g_free_rec_(GDB *gdb, GClient c, GFileN file_N);


#endif /*_G_REQUEST_H_*/

