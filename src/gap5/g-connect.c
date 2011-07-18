/*
 * File: g-connect.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: client based server connect/disconnect
 *
 * Created: prior to 18 September 1992
 * Updated:
 *
 */

#include "g-connect.h"

#include "array.h"

#include "g-struct.h"
#include "g-os.h" /* IMPORT: G_LOCK_RW */
#include "g-defs.h" /* IMPORT: G_LOCK_RW */
#include "g-db.h" /* IMPORT: g_client_shutdown */
#include "g-error.h"

GClient g_connect_client_(GDB *gdb, int clientID, GLock mode_requested, GLock *mode_granted)
/*
 * Connect client to server
 */
{
    GCardinal i;

    /* can we accept annother client? */
    if (gdb->ConnectedClients == gdb->Nclient) {
	(void)gerr_set(GERR_MAX_CLIENTS);
	return -1;
    }

    /* check to see that this client isn't already connected */
    for (i=0;i<gdb->Nclient;i++)
	if (arr(Client,gdb->client,i).id==clientID) {
	    (void)gerr_set(GERR_ALREADY_CONNECTED);
	    return -1;
	}

    /* we can accept this client - find a free slot */
    for (i=0;arr(Client,gdb->client,i).id!=-1 && i<gdb->Nclient;i++) ;
    if (i==gdb->Nclient) {
	/* this should not happen */
	(void)gerr_set(GERR_MAX_CLIENTS);
	return -1;
    }

    arr(Client,gdb->client,i).id = clientID;
    /* check mode */
    arr(Client,gdb->client,i).max_lock = mode_requested; /* YUK! should we trust client? */
    *mode_granted = mode_requested;

    gdb->ConnectedClients++;

    return i;
}


int g_disconnect_client_(GDB *gdb, GClient c)

/*
 *
 */
{
    /*
     * Maybe the following should only be called if there
     * are views still open. We could keep a count of the number
     * of open views.
     */
    (void) g_client_shutdown(gdb,c);

    return 0;
}


