/*
 * File: gap-init.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: open and close the gap server interface
 *
 * Created: 9 February 1993
 * Updated:
 *
 */


#include <stdio.h>
#include <stdlib.h> /* IMPORT: getenv */
#include <string.h>
/*#include <malloc.h>*/

#include "gap-init.h"

#include "gap-if.h"
#include "gap-io.h"

static int LOCAL = -1;
static char *gap_server;



void gap_init(void)
{

#define GAP_SERVER "GAP_SERVER"

    /* determine local mode */
    if (LOCAL==-1) {
	gap_server = (char *)getenv(GAP_SERVER);
	LOCAL = gap_server==NULL || strlen(gap_server)==0;
	gap_set_if_vectors(LOCAL);

	gap_io_init();
    }

}


int gap_server_is_local(void)
{
    gap_init();

    return LOCAL;
}


GapServer *gap_open_server(char *database, char *version, int read_only)
/*
 * Open a gap database server
 */
{
    gap_init();

    return g_open_server(database,version,read_only);
}



void gap_shutdown_server(GapServer *s)
{
    if (s!=NULL) {
	g_close_server(s);
    }
}
