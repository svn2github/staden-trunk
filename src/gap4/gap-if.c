/*
 * File: gap-if.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: initialise gap-if function pointers
 *
 * Created: 11 February 1993
 * Updated:
 *
 */

#include "gap-if.h"

#include "gap-local.h"
#include "gap-remote.h"
#include "gap-init.h"




int (*g_lock_file_N)(GapClient *s, GCardinal file_N) = NULL;
int (*g_unlock_file_N)(GapClient *s, GCardinal file_N) = NULL;
GView (*g_lock_N)(GapClient *s, GCardinal file_N, GCardinal rec, GLock lock) = NULL;
int (*g_upgrade)(GapClient *s, GView v, GLock lock) = NULL;
int (*g_unlock)(GapClient *s, GView v) = NULL;
int (*g_abandon)(GapClient *s, GView v) = NULL;
int (*g_flush)(GapClient *s, GView v) = NULL;
int (*g_read)(GapClient *s, GView v, void *buf, GCardinal len) = NULL;
int (*g_readv)(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt) = NULL;
int (*g_write)(GapClient *s, GView v, void *buf, GCardinal len) = NULL;
int (*g_writev)(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt) = NULL;
int (*g_remove)(GapClient *s, GView v) = NULL;
int (*g_view_info)(GapClient *s, GView v, GViewInfo *info) = NULL;
int (*g_rec_info)(GapClient *s, GCardinal file_N, GCardinal rec, GRecInfo *info) = NULL;
int (*g_header_info)(GapClient *s, GCardinal file_N, GHeaderInfo *info) = NULL;
int (*g_fast_read_N)(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len) = NULL;
int (*g_fast_readv_N)(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt) = NULL;
int (*g_fast_write_N)(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len) = NULL;
int (*g_fast_writev_N)(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt) = NULL;
GapClient * (*g_connect_client)(GapServer *s, GLock mode) = NULL;
int (*g_disconnect_client)(GapClient *s) = NULL;
void (*g_close_server)(GapServer *s) = NULL;
GapServer * (*g_open_server)(char *database, char *version, int read_only) = NULL;




void gap_set_if_vectors(int use_local)
{
    if (use_local) {
	g_lock_file_N = local_g_lock_file_N;
	g_unlock_file_N = local_g_unlock_file_N;
	g_lock_N = local_g_lock_N;
	g_upgrade = local_g_upgrade;
	g_unlock = local_g_unlock;
	g_abandon = local_g_abandon;
	g_flush = local_g_flush;
	g_read = local_g_read;
	g_readv = local_g_readv;
	g_write = local_g_write;
	g_writev = local_g_writev;
	g_remove = local_g_remove;
	g_view_info = local_g_view_info;
	g_rec_info = local_g_rec_info;
	g_header_info = local_g_header_info;
	g_fast_read_N = local_g_fast_read_N;
	g_fast_readv_N = local_g_fast_readv_N;
	g_fast_write_N = local_g_fast_write_N;
	g_fast_writev_N = local_g_fast_writev_N;
	g_connect_client = local_g_connect_client;
	g_disconnect_client = local_g_disconnect_client;
	g_close_server = local_g_close_server;
	g_open_server = local_g_open_server;
    } else {
	g_lock_file_N = remote_g_lock_file_N;
	g_unlock_file_N = remote_g_unlock_file_N;
	g_lock_N = remote_g_lock_N;
	g_upgrade = remote_g_upgrade;
	g_unlock = remote_g_unlock;
	g_abandon = remote_g_abandon;
	g_flush = remote_g_flush;
	g_read = remote_g_read;
	g_readv = remote_g_readv;
	g_write = remote_g_write;
	g_remove = remote_g_remove;
	g_view_info = remote_g_view_info;
	g_rec_info = remote_g_rec_info;
	g_header_info = remote_g_header_info;
	g_fast_read_N = remote_g_fast_read_N;
	g_fast_readv_N = remote_g_fast_readv_N;
	g_fast_write_N = remote_g_fast_write_N;
	g_fast_writev_N = remote_g_fast_writev_N;
	g_connect_client = remote_g_connect_client;
	g_disconnect_client = remote_g_disconnect_client;
	g_close_server = remote_g_close_server;
	g_open_server = remote_g_open_server;
    }

}
    



