#include "misc.h"
#include "IO.h"
#include "fort.h"
#include "list_proc.h"

int 
disassemble_readings(int handle,                                       /* in */
		     char *list,                                       /* in */
                     int iall,                                         /* in */
                     int iopt)                                         /* in */
{
    GapIO *io;
    f_int ngels, nconts;
    char *gel;
    f_int iok, fiopt = iopt, fiall = iall;
    f_int *array;
    int max_gel;

    if ( (io = io_handle(&handle)) == NULL){	
	return -1;
    }

    max_gel = find_max_gel_len(io, 0, 0)+1;

    ngels = NumReadings(io);
    nconts = NumContigs(io);

    if (-1 == set_active_list(list)) {
	return -1;
    }

    if ((array = (f_int *)xmalloc(max_gel * sizeof(f_int)))==NULL){
	return(-1);
    }
    if ((gel = (char *)xmalloc(max_gel * sizeof(char)))==NULL){
	return(-1);
    }

    remgbc_(&io_relpos(io,1), &io_length(io,1), &io_lnbr(io,1),
	    &io_rnbr(io,1), &ngels, &nconts, &io_dbsize(io), gel,
            &max_gel, &handle, &iok, array, &fiall, &fiopt,
            (f_int)max_gel);

    flush2t(io);

    xfree(array);
    xfree(gel);

    return 0;
            
} /* end disassemble_readings */

