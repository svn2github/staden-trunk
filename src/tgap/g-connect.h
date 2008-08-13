/*
 * File: g-connect.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for client based server connect/disconnect
 *
 * Created: prior to 18 September 1992
 * Updated:
 *
 */

#ifndef _G_CONNECT_H_
#define _G_CONNECT_H_

#include "g-struct.h"

extern GClient g_connect_client_(GDB *gdb, int client, GLock mode_requested, GLock *mode_granted);
/*
 * Connect client to server
 */

extern int g_disconnect_client_(GDB *gdb, GClient c);
/*
 *
 */


#endif /*_G_CONNECT_H_*/
