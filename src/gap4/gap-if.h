/*
 * File: gap-if.h
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

#ifndef _GAP_IF_H_
#define _GAP_IF_H_


#include "g-os.h"
#include "g-struct.h"



#define GAP_LOCAL_SERVER 1
#define GAP_REMOTE_SERVER 2



typedef union {
    struct {
	int type;
    } generic;
    struct {
	int type;
	GDB *gdb;
    } local;
    struct {
	int type;
	char *gapserver;
    } remote;
} GapServer;

typedef union {
    struct {
	GapServer *server;
	GLock mode;
    } generic;
    struct {
	GapServer *server;
	GLock mode;
	GClient client;
    } local;
    struct {
	GapServer *server;
	GLock mode;
	int fd;
    } remote;
} GapClient;



extern void gap_set_if_vectors(int use_local);

extern int (*g_lock_file_N)(GapClient *s, GCardinal file_N);
extern int (*g_unlock_file_N)(GapClient *s, GCardinal file_N);
extern GView (*g_lock_N)(GapClient *s, GCardinal file_N, GCardinal rec, GLock lock);
extern int (*g_upgrade)(GapClient *s, GView v, GLock lock);
extern int (*g_unlock)(GapClient *s, GView v);
extern int (*g_abandon)(GapClient *s, GView v);
extern int (*g_flush)(GapClient *s, GView v);
extern int (*g_read)(GapClient *s, GView v, void *buf, GCardinal len);
extern int (*g_readv)(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt);
extern int (*g_write)(GapClient *s, GView v, void *buf, GCardinal len);
extern int (*g_writev)(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt);
extern int (*g_remove)(GapClient *s, GView v);
extern int (*g_view_info)(GapClient *s, GView v, GViewInfo *info);
extern int (*g_rec_info)(GapClient *s, GCardinal file_N, GCardinal rec, GRecInfo *info);
extern int (*g_header_info)(GapClient *s, GCardinal file_N, GHeaderInfo *info);
extern int (*g_fast_read_N)(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len);
extern int (*g_fast_readv_N)(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt);
extern int (*g_fast_write_N)(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len);
extern int (*g_fast_writev_N)(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt);
extern GapClient * (*g_connect_client)(GapServer *s, GLock mode);
extern int (*g_disconnect_client)(GapClient *s);
extern void (*g_close_server)(GapServer *s);
extern GapServer * (*g_open_server)(char *database, char *version,
				    int read_only);



#endif /*_GAP_IF_H_*/

