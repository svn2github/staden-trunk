/*
 * File: gap-create.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for gap-create.c
 *
 * Created: 21 January 1993
 * Updated:
 *
 */

#ifndef _GAP_CREATE_H_
#define _GAP_CREATE_H_


/*extern int gap_create_db(char *project,char *version);*/

/*extern int gap_init_db(char *project,char *version, int read_only);*/

extern int gap_new_db(char *project,char *version, int read_only);

extern int cpdb(char *base, char *from, char *to);

void set_db_bitsize(int bitsize);

extern int maxdb;

#endif /*_GAP_CREATE_H_*/


    
