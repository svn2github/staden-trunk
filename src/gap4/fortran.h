/*
 * File: fortran.h
 * Version:
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
#ifndef _FORTRAN_H_
#define _FORTRAN_H_

#include "fort.h"

/* This defines the interfaces between fortran and c */

/*
 * COMMON /DEVILS/ IDEVR,IDBSIZ,RELPG
 */

/*
 * Map onto fortran common block holding RELPG array, IDBSIZ and IDEVR
 */



/*
 * Snatch Fortran Common Block
 * Don't EVER EVER let Rodger see this!!
 */
extern struct {
    int_f handle;
    int_f idbsiz;
} devils_;

#endif /* _FORTRAN_H_ */

