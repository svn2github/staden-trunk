/*
 * File: gap-remote.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: gap server interface routine prototypes
 *
 * Created: 9 February 1993
 * Updated:
 *
 */

#ifndef _GAP_REMOTE_H_
#define _GAP_REMOTE_H_

#include "gap-if.h"

extern int remote_g_lock_file_N(GapClient *s, GCardinal file_N);
extern int remote_g_unlock_file_N(GapClient *s, GCardinal file_N);
extern GView remote_g_lock_N(GapClient *s, GCardinal file_N, GCardinal rec, GLock lock);
extern int remote_g_upgrade(GapClient *s, GView v, GLock lock);
extern int remote_g_unlock(GapClient *s, GView v);
extern int remote_g_abandon(GapClient *s, GView v);
extern int remote_g_flush(GapClient *s, GView v);
extern int remote_g_read(GapClient *s, GView v, void *buf, GCardinal len);
extern int remote_g_readv(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt);
extern int remote_g_write(GapClient *s, GView v, void *buf, GCardinal len);
extern int remote_g_writev(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt);
extern int remote_g_remove(GapClient *s, GView v);
extern int remote_g_view_info(GapClient *s, GView v, GViewInfo *info);
extern int remote_g_rec_info(GapClient *s, GCardinal file_N, GCardinal rec, GRecInfo *info);
extern int remote_g_header_info(GapClient *s, GCardinal file_N, GHeaderInfo *info);
extern int remote_g_fast_read_N(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len);
extern int remote_g_fast_readv_N(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt);
extern int remote_g_fast_write_N(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len);
extern int remote_g_fast_writev_N(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt);
extern GapClient *remote_g_connect_client(GapServer *s, GLock mode);
extern int remote_g_disconnect_client(GapClient *s);
extern void remote_g_close_server(GapServer *s);
extern GapServer *remote_g_open_server(char *database, char *version,
				       int read_only);

#endif /*_GAP_REMOTE_H_*/

