#ifndef _GAP_GLOBALS_H_
#define _GAP_GLOBALS_H_

#include <tcl.h>
#include "os.h"
#include "tcl_utils.h"


int init_globals(Tcl_Interp *interp);


/* Gap4 global variable accessor functions */
int    gap4_global_get_quality_cutoff( void );
void   gap4_global_set_quality_cutoff( int qc );
float  gap4_global_get_consensus_cutoff( void );
void   gap4_global_set_consensus_cutoff( float cc );
int    gap4_global_get_consensus_mode( void );
void   gap4_global_set_consensus_mode( int cm );
int    gap4_global_get_maxseq( void );
void   gap4_global_set_maxseq( int ms );
int    gap4_global_get_gopenval( void );
void   gap4_global_set_gopenval( int ov );
int    gap4_global_get_gextendval( void );
void   gap4_global_set_gextendval( int ev );
double gap4_global_get_template_size_tolerance( void );
void   gap4_global_set_template_size_tolerance( double ts );





/* C Variable 			  Tcl variable */
/*---------------------------------------------*/



/*-------------------------------------*/
/* Globals exported to other libraries */
/*-------------------------------------*/

/* On windows, these must only be accessed via the accessor functions above */

extern int   maxseq;			/* "maxseq" */
extern float consensus_cutoff;		/* "consensus_cutoff" */
extern int   consensus_mode;		/* "consensus_mode" */
extern int   quality_cutoff;		/* "quality_cutoff" */
extern int   gopenval;			/* "align_open_val" */
extern int   gextendval;		/* "align_extend_val" */



/*---------------------------------------------*/
/* Globals not yet exported to other libraries */
/*---------------------------------------------*/

extern int chem_as_double;		/* "chem_as_double" */
extern Tcl_Obj *gap_defs;
extern int maxdb;			/* "maxdb" - see gap-create.c,
						   write-only */
extern int ignore_checkdb;		/* "ignore_checkdb" */
extern int consensus_iub;		/* "consensus_iub" -  see qual.c */
					/* "gap_auto_flush" - static to gap-tcl.c */
					/* "gap_fatal_errors" - see gap-error.c */
extern int exec_notes;			/* "exec_notes" - see notes.c */
extern int rawdata_note;		/* "rawdata_note" - see notes.c */

extern double template_size_tolerance;  /* "template_size_tolerance" */
extern int    min_vector_len;		/* Minimum length of SVEC tag */
extern int    template_check_flags;	/* TEMP_OFLAG_* for template.c */

#define CONSENSUS_MODE_FREQ	  0
#define CONSENSUS_MODE_WEIGHTED	  1 	/* as FREQ, determined by qual_cut >= 0 */
#define CONSENSUS_MODE_CONFIDENCE 2



#endif /* _GAP_GLOBALS_H_ */
