/*
 * File: gap-init.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for opening and closeing the gap server interface
 *
 * Created: 9 February 1993
 * Updated:
 *
 */

#ifndef _GAP_INIT_H_
#define _GAP_INIT_H_

#include "gap-if.h"

extern int gap_server_is_local(void);
extern GapServer *gap_open_server(char *database, char *version,int read_only);
extern void gap_shutdown_server(GapServer *s);

#endif /*_GAP_INIT_H_*/
