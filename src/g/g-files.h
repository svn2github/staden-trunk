/*
 * File: g-files.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created: prior to 18 September 1992
 * Updated:
 *
 */

#ifndef _G_FILES_H_
#define _G_FILES_H_

#include "g-struct.h"


extern GFile *g_open_file(char *fn, int read_only);
/*
 * Open a file and its associated index
 */


extern void g_close_file(GFile *g);
/*
 * Close a file and its associated index
 */

extern int g_write_aux_index(GFile *gfile, GCardinal rec);
/*
 * Read a record from the index of the aux file
 */


extern int g_write_aux_header(GFile *gfile);
/*
 * Write the header of the aux file
 */


extern char *g_filename(GFile *gfile);
/*
 * Returns the file name of an open file
 */


extern int g_remove_client(GFile *gfile, GClient client);
/*
 * Remove all locks in this file for this client
 */

extern int g_check_header(GFile *gfile);
/*
 * Checks whether the on-disk copy matches the in-memory copy of the
 * header.
 */

#endif /*_G_FILES_H_*/
