/*
 * File: gap-remote.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: interface routines for remote server
 *
 * Created: 5 February 1993
 * Updated:
 *
 */

#include <stdio.h> /* IMPORT: NULL */
/*#include <malloc.h>*/

#include "g-request.h"
#include "g-connect.h"
#include "gap-dbstruct.h"
#include "g-db.h"
#include "g-error.h"

#include "gap-remote.h"



int remote_g_lock_file_N(GapClient *s, GCardinal file_N)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



int remote_g_unlock_file_N(GapClient *s, GCardinal file_N)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}


 
GView remote_g_lock_N(GapClient *s, GCardinal file_N, GCardinal rec, GLock lock)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_upgrade(GapClient *s, GView v, GLock lock)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_unlock(GapClient *s, GView v)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_abandon(GapClient *s, GView v)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}




/* ARGSUSED */
int remote_g_flush(GapClient *s, GView v)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_read(GapClient *s, GView v, void *buf, GCardinal len)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_readv(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}


/* ARGSUSED */
int remote_g_write(GapClient *s, GView v, void *buf, GCardinal len)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}

/* ARGSUSED */
int remote_g_writev(GapClient *s, GView v, GIOVec *vec, GCardinal vcnt)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}


/* ARGSUSED */
int remote_g_remove(GapClient *s, GView v)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}


/* ARGSUSED */
int remote_g_view_info(GapClient *s, GView v, GViewInfo *info)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}

/* ARGSUSED */
int remote_g_rec_info(GapClient *s, GCardinal file_N, GCardinal rec, GRecInfo *info)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_header_info(GapClient *s, GCardinal file_N, GHeaderInfo *info)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}






/* ARGSUSED */
int remote_g_fast_read_N(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_fast_readv_N(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}




/* ARGSUSED */
int remote_g_fast_write_N(GapClient *s, GCardinal file_N, GCardinal rec, void *buf, GCardinal len)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



/* ARGSUSED */
int remote_g_fast_writev_N(GapClient *s, GCardinal file_N, GCardinal rec, GIOVec *vec, GCardinal vcnt)
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}







/* ARGSUSED */
GapClient *remote_g_connect_client(GapServer *s, GLock mode)
/*
 * Connect client to server
 */
{
    (void)gerr_set(GERR_NOT_IMPLEMENTED);
    return NULL;
}


int remote_g_disconnect_client(GapClient *s)
/*
 *
 */
{
    return gerr_set(GERR_NOT_IMPLEMENTED);
}



GapServer *remote_g_open_server(char *database, char *version, int read_only)
{
    (void)gerr_set(GERR_NOT_IMPLEMENTED);
    return NULL;
}


void remote_g_close_server(GapServer *s)
{
    (void)gerr_set(GERR_NOT_IMPLEMENTED);
}

