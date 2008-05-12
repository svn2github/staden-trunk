/*
 * File: gap-dbstruct.h
 *
 * Author: Staden Package team
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: structure of gap database files
 *
 * Created: 28 October 1992
 * Updated:
 *
 */

#ifndef _GAP_DBSTRUCT_H_
#define _GAP_DBSTRUCT_H_

#include <sys/types.h>

#include "g-os.h"




#define GAP_DB_VERSION 3
/*
 * history
 * V0. Reset to year 0, 11 November 1993
 * V1. 26/04/94  Tags are now stored relative to the 'real sequence' start
 *               rather than the 'used sequence' start.
 * V2. 24/02/95  Primer type now holds values 0 to 4 (instead of 0-3). Version
 *               update is simply to allow for removing dead code in future
 *               if desired.
 * V3. 28/09/98  Notes are added. These are unpositional tags for attaching
 *		 information to entire readings, contigs, or the database.
 */

#define GAP_READ_LEN	        30000
#define GAP_FILES 1
#define DATABASE_FILE		0
#define GAP_DATABASE_FILE	0


/*
 * Macros
 */
/* GDatabase.data_class */
#define GAP_DNA            0
#define GAP_PROTEIN        1
/* GReadings.sense */
#define GAP_SENSE_ORIGINAL 0
#define GAP_SENSE_REVERSE  1
/* GReadings.strand */
#define GAP_STRAND_FORWARD 0
#define GAP_STRAND_REVERSE 1
/* GReadings.primer */
#define GAP_PRIMER_UNKNOWN 0
#define GAP_PRIMER_FORWARD 1
#define GAP_PRIMER_REVERSE 2
#define GAP_PRIMER_CUSTFOR 3
#define GAP_PRIMER_CUSTREV 4
/* GReadings.chemistry */
/*	Bit 0 is 1 for terminator, 0 for primer */
#define GAP_CHEM_TERMINATOR	(1<<0)
/*	Bits 1 to 4 inclusive are the type (any one of, not bit pattern) */
#define GAP_CHEM_TYPE_MASK	(15<<1)
#define GAP_CHEM_TYPE_UNKNOWN	(0<<1)
#define GAP_CHEM_TYPE_ABI_RHOD	(1<<1)
#define GAP_CHEM_TYPE_ABI_DRHOD	(2<<1)
#define GAP_CHEM_TYPE_BIGDYE2	(3<<1)
#define GAP_CHEM_TYPE_ET	(4<<1)
#define GAP_CHEM_TYPE_LICOR	(5<<1)
#define GAP_CHEM_TYPE_MB_ET	(6<<1)
#define GAP_CHEM_TYPE_BIGDYE1	(7<<1)
#define GAP_CHEM_TYPE_BIGDYE3	(8<<1)
/* GVectors.level */
#define GAP_LEVEL_UNKNOWN  0
#define GAP_LEVEL_CLONE    1
#define GAP_LEVEL_SUBCLONE 2

/*
0	unknown		primer
1	unknown		terminator
2	ABI rhodamine	primer
3	ABI rhodamine 	terminator
4	ABI dRhodamine 	primer
5	ABI dRhodamine 	terminator
6	ABI BigDye v2 	primer
7	ABI BigDye v2 	terminator
8	Energy transfer	primer
9	Energy transfer	terminator
10	Licor 		primer
11	Licor 		terminator
12	MegaBACE ET 	primer
13	MegaBACE ET 	terminator
14	ABI BigDye v1 	primer
15	ABI BigDye v1 	terminator
16	ABI BigDye v3 	primer
17	ABI BigDye v3 	terminator
*/




extern char *file_list[];

extern size_t block_sizes[];

extern GCardinal max_recs[];


#endif /*_GAP_DBSTRUCT_H_*/
