/*
 * File: g-defs.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: gap server defaults
 *
 * Created: prior to 18 September 1992
 * Updated:
 *
 */

#ifndef _G_DEFS_H_
#define _G_DEFS_H_
/* aux file suffix */
#define G_AUX_SUFFIX ".aux"

/* default file permissions */
#define G_DEF_PERMS 0666

/* maximum number of clients */
#define G_MAX_CLIENTS 8


/*
 * define lock modes
 */

#define G_LOCK_NONE (GLock) 0
#define G_LOCK_RO (GLock) 1
#define G_LOCK_RW (GLock) 2
#define G_LOCK_EX (GLock) 3




#endif /*_G_DEFS_H_*/
