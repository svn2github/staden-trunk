/*
 * File: g-struct.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: allocation and deallocation of g-struct.h structures
 *
 * Created: 14 October 1992
 * Updated:
 *
 */


#include <stdio.h> /* IMPORT: NULL */
#include <unistd.h>
#include <errno.h>
/*#include <malloc.h>*/

#include "array.h"
#include "freetree.h"
#include "g-files.h" /* IMPORT: g_close_file */
#include "g-struct.h"
#include "g-io.h"
#include "xalloc.h"



static void g_destroy_index(GFile *gfile)
/*
 * destroy gfile index
 */
{
    if (gfile!=NULL && gfile->idx!=NULL) {
	ArrayDestroy(gfile->idx);
	gfile->idx=NULL;
    }
}




/*
 * GFile
 */

GFile *g_new_gfile(int bitsize)
/*
 * create and initialise a new gfile structure
 */
{
    GFile *gfile;
    int endian = 1;

    if (NULL == ( gfile = (GFile *)xmalloc(sizeof(GFile))))
	return NULL;

    gfile->fname = NULL;
    gfile->fd = gfile->fdaux = -1;
    gfile->freetree = NULL;
    gfile->Nidx = 0; gfile->idx = NULL;
    gfile->flock_client = -1; /* An invalid client number */
    gfile->flock_status = G_FLOCK_NONE;
    gfile->check_header = 1;
    if ( *(char *)&endian ) {
	gfile->low_level_vector = (bitsize == G_64BIT)
	    ? low_level_vectors_swapped64
	    : low_level_vectors_swapped32;
	gfile->swapped = 1;
    } else {
	gfile->low_level_vector = (bitsize == G_64BIT)
	    ? low_level_vectors64
	    : low_level_vectors32;
	gfile->swapped = 0;
    }

    return gfile;
}



void g_free_gfile(GFile *gfile)
/*
 * free gfile structure
 */
{
    if (gfile != NULL) {
	if (gfile->fname != NULL) xfree(gfile->fname);
	/* LOW LEVEL IO HERE */
	errno = 0;
	if (gfile->fd != -1) close(gfile->fd);
	/* LOW LEVEL IO HERE */
	if (gfile->fdaux != -1) close(gfile->fdaux);
	if (gfile->idx != NULL) g_destroy_index(gfile);
	if (gfile->freetree != NULL) freetree_destroy(gfile->freetree);
	xfree(gfile);
    }
}



GDB *g_new_gdb()
/*
 * create and initialise a new gdb structure
 */
{
    GDB *gdb;
    if ( ( gdb = (GDB *)xmalloc(sizeof(GDB)) ) != NULL )	{
        gdb->gfile = NULL;
        gdb->client = NULL;
        gdb->Nclient = gdb->ConnectedClients = 0;
        gdb->view = NULL;
        gdb->Nview = 0;
	gdb->free_view = -1;
    }
    return gdb;
}


void g_free_gdb(GDB *gdb)
/*
 * free gdb structure
 */
{
    if (gdb != NULL) {
	if (gdb->gfile != NULL) {
	    g_close_file(gdb->gfile);
	    gdb->gfile = NULL;
	}
	if (gdb->client != NULL) { ArrayDestroy(gdb->client); gdb->client = NULL; }
	if (gdb->view != NULL) { ArrayDestroy(gdb->view); gdb->view = NULL; }
	xfree(gdb);
    }
}





/*
 * Views
 */

GView g_new_view(GDB *gdb)
/*
 * allocate a new view
 */
{
    GView i;

    i = gdb->free_view;
    if (i==-1) {
	(void)ArrayRef(gdb->view,gdb->Nview);
	i = gdb->Nview++;
    } else {
	/* adjust freetree */
	gdb->free_view = arr(View,gdb->view,i).next;
    }

    arr(View,gdb->view,i).next = -1;
    arr(View,gdb->view,i).flags = G_VIEW_NEW;
    arr(View,gdb->view,i).lcache.rec = -1;
    return i;

}


void g_free_view(GDB *gdb, GView view)
/*
 * free a new view
 */
{
    if (gdb==NULL || view==-1) return;

    if (view < 0 || view >= gdb->Nview ||
	(arr(View,gdb->view,view).flags & G_VIEW_FREE)) return;

    arr(View,gdb->view,view).flags = G_VIEW_FREE;
    arr(View,gdb->view,view).next = gdb->free_view;
    gdb->free_view = view;

}



