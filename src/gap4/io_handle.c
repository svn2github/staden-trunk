#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>

#include "IO.h"
#include "misc.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "xalloc.h"
#include "os.h" /* 6/1/99 johnt - for sysconf() */
#include "text_output.h"


/*************************************************************
 * FORTRAN handle stuff
 *************************************************************/

static int NHandles = 0;
static GapIO **Handles = NULL;

/*
 * Called from the SIGTERM handler; sent out during a system shutdown.
 *
 * This attempts to shutdown all plots. When this fails, it quits anyway.
 * The database is flushed, then closed, and finally the process exits.
 */
static void shutdown_handles(int sig) {
    int i, cnum;
    reg_quit rq;
    GapIO *io;

    for (i = 0; i < NHandles; i++) {
	if (!Handles[i])
	    continue;

	log_file(NULL, "Received SIGTERM");

	/* An open database, so quit the displays (if possible) */
	io = Handles[i];
	for (cnum = 1; cnum <= NumContigs(io); cnum++) {
	    rq.job = REG_QUIT;
	    rq.lock = REG_LOCK_WRITE;
	    contig_notify(io, cnum, (reg_data *)&rq);
	    if (!(rq.lock & REG_LOCK_WRITE)) {
		log_file(NULL, "A display has write lock - data may be lost");
	    }
	}

	/* close (which also flushes) */
	close_db(io);
    }

    /* Finally, exit gracefully */
    exit(0);
}

static int initialise_handle(void)
{
    int i;
    static int done_init_handle = 0;
    
    if (done_init_handle)
	return 0;

    done_init_handle++;

    /* get max number of open files */
    if (-1 == (NHandles = sysconf(_SC_OPEN_MAX)))
	return 1; 

    NHandles /=  GAP_FILES;
    if (!NHandles) return 1;

    if ( (Handles = (GapIO **)xmalloc(sizeof(GapIO *) * NHandles)) == NULL) {
	NHandles = 0;
	return 1;
    }
    
    for (i=0; i<NHandles; i++) Handles[i] = NULL;

    /* Add signal handler to shutdown when SIGTERM is received */
    signal(SIGTERM, shutdown_handles);

    return 0;
}


int 
get_free_handle(GapIO *io)
/*
 * find a free entry in the Exp array
 * returns -1 if there is none
 */
{
    int i;
    
    (void) initialise_handle();
    
    if (!NHandles) return -1; /* no slots! */
    for (i=0; i<NHandles && Handles[i]!=NULL; i++) ;

    if (i == NHandles) {
	return -1;

    } else {

	Handles[i] = io;
	return i+1;
    }
    /* return (i==NHandles)?-1:i; */
}


void remove_handle(f_int *handle) {
    Handles[*handle - 1] = NULL;
}


static int check_handle(f_int *handle)
{
    return (handle == NULL ||
	    (int) (*handle) <= 0 ||
	    (int) (*handle) > NHandles);
}


GapIO *io_handle(f_int *HANDLE)
{

    if ( check_handle(HANDLE) ) {
	return NULL;
    }
    return (GapIO *) Handles[(int)(*HANDLE)-1];
}

f_int *handle_io(GapIO *io) {
    static f_int handle;
    
    for (handle = 0; handle < NHandles; handle++)
	if (Handles[handle] == io) {
	    handle++;
	    return &handle;
	}

    return NULL;
}
