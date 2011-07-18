/*
 * File: g-db.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for g-db.h
 *
 * Created: prior to 18-Sep-1992
 * Updated:
 *
 */

#ifndef _G_DB_H_
#define _G_DB_H_

#include "g-struct.h"

extern GDB *g_open_database_(char *fns[], GCardinal Nfns, int read_only);
/*
 * Open a database with a given file name
 */



extern void g_shutdown_database_(GDB *gdb);
/*
 * shut down a database
 */



#define panic_shutdown() panic_shutdown_(__FILE__,__LINE__)
/*
 * When something fatal happens we need to shut down database in a rather
 * crude fashion.
 */


extern void panic_shutdown_(char *file, int line);
/*
 * When something fatal happens we need to shut down database in a rather
 * crude fashion.
 */


extern int g_client_shutdown(GDB *gdb, GClient c);
/*
 * disengage client from a database
 */

#endif /*_G_DB_H_*/
