/*
 * File:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */
#ifndef _GAP_IO_H_
#define _GAP_IO_H_

#include <stdio.h>

#include "g-error.h"


/*************************************************************
 * Low level IO routines
 *************************************************************/


extern int (* GAP_READ) (GapClient *s, GView v, void *buf, int len, GCardinal type_check, int size);
extern int (* GAP_WRITE) (GapClient *s, GView v, void *buf, int len, GCardinal type_check, int size);


extern void gap_io_init(void);






#endif /*_GAP_IO_H_*/
