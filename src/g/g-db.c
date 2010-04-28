/*
 * File: g-db.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: database opening and shutdown (server)
 *
 * Created: prior to 18-Sep-1992
 * Updated:
 *
 */

#include <staden_config.h>

#include <stdio.h> /* IMPORT: NULL */
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
/*#include <malloc.h>*/ /* IMPORT: malloc */

#include "array.h"

#include "g-db.h"
#include "g-defs.h"
#include "g-os.h"
#include "g-error.h"
#include "g-files.h" /* IMPORT: g_open_file */
#include "g-io.h" /* IMPORT: set_low_level_vector */
#include "g-request.h" /* IMPORT: g_abandon_ */






GDB *g_open_database_(char *fns[], GCardinal Nfns, int read_only)
/*
 * Open a database with a given file name
 */
{
    GDB *gdb;
    GCardinal i;

    /* check arguments */
    if (fns==NULL) {
	(void) gerr_set(GERR_INVALID_ARGUMENTS);
	return NULL;
    }

    if (NULL == (gdb = g_new_gdb()))
	return NULL;

    /*
     * initialise data structs for clients
     */
    gdb->Nclient = G_MAX_CLIENTS;
    if ( (gdb->client = ArrayCreate(sizeof(Client),gdb->Nclient)) == NULL ) {
	g_free_gdb(gdb);
	(void)gerr_set(GERR_OUT_OF_MEMORY);
	return NULL;
    }
    (void)ArrayRef(gdb->client,gdb->Nclient-1);
    for (i=0;i<gdb->Nclient;i++) arr(Client,gdb->client,i).id = -1;


    /*
     * Open the file (now only 1 file allowed, we ignore the rest)
     */
    gdb->gfile = g_open_file(fns[0], read_only);
    if (NULL == gdb->gfile) {
	g_free_gdb(gdb);
	/* g_open_file sets gerrnum */
	return NULL;
    }

    /*
     * allocate views - assume using Gap4, which locks all records.
     */
    gdb->Nview = gdb->gfile->header.num_records;
    if ( (gdb->view = ArrayCreate(sizeof(View),gdb->Nview)) == NULL ) {
	g_free_gdb(gdb);
	(void)gerr_set(GERR_OUT_OF_MEMORY);
	return NULL;
    }
    (void)ArrayRef(gdb->view,gdb->Nview-1);

    /* initialise views */
    for(i=0;i<gdb->Nview;i++) {
	arr(View,gdb->view,i).next = (i-1);
	arr(View,gdb->view,i).flags = G_VIEW_NEW;
    }
    gdb->free_view = gdb->Nview-1;
    

    return gdb;
}



void g_shutdown_database_(GDB *gdb)
/*
 * shut down a database
 */
{
    GFile *g;

    if (gdb==NULL) return;

    /* Save the freetree, as long as we opened in read-write mode */
    if (g = gdb->gfile) {
#ifndef _WIN32
	/* LOW LEVEL IO HERE */
	int mode = fcntl(g->fdaux, F_GETFL, 0);
	if ((mode & O_ACCMODE) & O_RDWR) {
#endif
	    int recsize;

	    /* LOW LEVEL IO HERE */
	    recsize = (g->header.format == G_32BIT)
		? sizeof(AuxIndex32)
		: sizeof(AuxIndex);
	    lseek(g->fdaux, sizeof(AuxHeader) +
		  g->header.num_records * recsize, SEEK_SET);
#ifdef CACHE_FREETREE
	    /* Save a cached copy of the freetree for fast startup next time */
	    if (g->header.format == G_32BIT)
		freetree_save_int4(g->freetree, g->fdaux, g->header.last_time);
	    else
		freetree_save_int8(g->freetree, g->fdaux, g->header.last_time);
#endif
#ifndef _WIN32
	}
#endif
    }

    g_free_gdb(gdb);
}





int g_client_shutdown(GDB *gdb, GClient c)
/*
 * disengage client from a database
 */
{
    GCardinal i;
    int err;
    if (gdb==NULL) return gerr_set(GERR_INVALID_ARGUMENTS);

    /*
     * Remove all view for each file of database from this client
     */
    for(i=0;i<gdb->Nview;i++) {
	if (arr(View,gdb->view,i).flags && !(arr(View,gdb->view,i).flags & G_VIEW_FREE) ) {
	    if (arr(View,gdb->view,i).client == c) {
		/* abandon for this client */
		(void) g_abandon_(gdb,c,i);
	    }
	}
    }


    /*
     * Remove all file locks for this client
     */
    err = g_remove_client(gdb->gfile, c);

    /* disentangle client c from gdb*/
    arr(Client,gdb->client,c).id = -1;
    gdb->ConnectedClients--;

    return err;
}


void panic_shutdown_(char *file, int line)
/*
 * When something fatal happens we need to shut down database in a rather
 * crude fashion.
 */
{
    fprintf(stderr,"** \n");
    fprintf(stderr,"** Panic in file %s at line %d\n",file,line);
    fprintf(stderr,"** A fatal error has occurred - shutting down immediately\n");
    fprintf(stderr,"** \n");

    /*
     * maybe we should dump the contents of the database here - so that we
     * can resume later?
     */

    exit(1);
}
