/*
 * File: gap-local.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: interface routines for local server
 *
 * Created: 5 February 1993
 * Updated:
 *
 */

#include <stdio.h> /* IMPORT: NULL */
/*#include <malloc.h>*/


#include "g-request.h"
#include "g-connect.h"
#include "g-db.h"

#include "gap-local.h"
#include "gap-dbstruct.h"
#include "xalloc.h"





int local_g_lock_file_N(GapClient *s, GCardinal file_N)
{
    return g_lock_file_N_(s->local.server->local.gdb, s->local.client, (GFileN)file_N);
}



int local_g_unlock_file_N(GapClient *s, GCardinal file_N)
{
    return g_unlock_file_N_(s->local.server->local.gdb, s->local.client, (GFileN)file_N);
}


 
GView local_g_lock_N(GapClient *s, GCardinal file_N, GCardinal rec, GLock lock)
{
    return g_lock_N_(s->local.server->local.gdb, s->local.client, (GFileN)file_N, rec, lock);
}



/* ARGSUSED */
int local_g_upgrade(GapClient *s, GView v, GLock lock)
{
    return g_upgrade_(s->local.server->local.gdb, s->local.client, v, lock);
}



/* ARGSUSED */
int local_g_unlock(GapClient *s, GView v)
{
    return g_unlock_(s->local.server->local.gdb, s->local.client, v);
}



/* ARGSUSED */
int local_g_abandon(GapClient *s, GView v)
{
    return g_abandon_(s->local.server->local.gdb, s->local.client, v);
}


/* ARGSUSED */
int local_g_flush(GapClient *s, GView v)
{
    return g_flush_(s->local.server->local.gdb, s->local.client, v);
}



/* ARGSUSED */
int local_g_read(GapClient *s, GView v, void *buf, GCardinal len)
{
    return g_read_(s->local.server->local.gdb, s->local.client, v, buf, len );
}



/* ARGSUSED */
int local_g_readv(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt)
{
    return g_readv_(s->local.server->local.gdb, s->local.client, v, vec, vcnt);
}


/* ARGSUSED */
int local_g_write(GapClient *s, GView v, void *buf, GCardinal len)
{
    return g_write_(s->local.server->local.gdb, s->local.client, v, buf, len);
}

/* ARGSUSED */
int local_g_writev(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt)
{
    return g_writev_(s->local.server->local.gdb, s->local.client, v, vec, vcnt);
}

/* ARGSUSED */
int local_g_remove(GapClient *s, GView v)
{
    return g_remove_(s->local.server->local.gdb, s->local.client, v);
}


/* ARGSUSED */
int local_g_view_info(GapClient *s, GView v, GViewInfo *info)
{
    return g_view_info_(s->local.server->local.gdb, s->local.client, v, info);
}

/* ARGSUSED */
int local_g_rec_info(GapClient *s, GCardinal file_N, GCardinal rec, GRecInfo *info)
{
    return g_rec_info_(s->local.server->local.gdb, s->local.client, (GFileN)file_N, rec, info);
}



/* ARGSUSED */
int local_g_header_info(GapClient *s, GCardinal file_N, GHeaderInfo *info)
{
    return g_header_info_(s->local.server->local.gdb, s->local.client, (GFileN)file_N, info);
}



/* ARGSUSED */
int local_g_fast_read_N(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len)
{
    return g_fast_read_N_(s->local.server->local.gdb, s->local.client, (GFileN)file_N, rec, buf, len);
}


/* ARGSUSED */
int local_g_fast_readv_N(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt)
{
    return g_fast_readv_N_(s->local.server->local.gdb, s->local.client, (GFileN)file_N, rec, vec, vcnt);
}



/* ARGSUSED */
int local_g_fast_write_N(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len)
{
    return g_fast_write_N_(s->local.server->local.gdb, s->local.client, (GFileN)file_N, rec, buf, len);
}



/* ARGSUSED */
int local_g_fast_writev_N(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt)
{
    return g_fast_writev_N_(s->local.server->local.gdb, s->local.client, (GFileN)file_N, rec, vec, vcnt);
}







/* ARGSUSED */
GapClient *local_g_connect_client(GapServer *s, GLock mode)
/*
 * Connect client to server
 */
{
    GapClient *c;
    GLock mode_granted;
    static GClient local_client = 0;
    if ( (c = (GapClient *)xmalloc(sizeof(GapClient))) != NULL) {
	local_client++;
	c->local.server = s;
	c->local.client = g_connect_client_(s->local.gdb, local_client, mode, &mode_granted);
	c->local.mode = mode_granted;
    }
    return c;
}



int local_g_disconnect_client(GapClient *s)
/*
 *
 */
{
    int status;
    status = g_disconnect_client_(s->local.server->local.gdb, s->local.client);
    xfree(s);
    return status;
}



GapServer *local_g_open_server(char *database, char *version, int read_only)
{
    GapServer *s;
    char files[GAP_FILES][1024];
    char *fileps[GAP_FILES];
    int i;

    if ( (s=(GapServer*)xmalloc(sizeof(GapServer))) != NULL) {

	for (i=0;i<GAP_FILES;i++) {
	    (void) gap_construct_file(database,file_list[i],version,files[i]);
	    fileps[i] = files[i];
	}
	s->local.type = GAP_LOCAL_SERVER;
	s->local.gdb = g_open_database_(fileps,GAP_FILES, read_only);

	if (s->local.gdb == NULL) {
	    xfree(s);
	    s = NULL;
	}
    }

    return s;
}


void local_g_close_server(GapServer *s)
{
    g_shutdown_database_(s->local.gdb);
    xfree(s);
}

