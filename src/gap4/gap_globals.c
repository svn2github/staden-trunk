/*
 * Initialises global variables for use in gap.
 *
 * and gap-error.c (gap_fatal_errors)
 * See also gap-tcl.c (auto_flush)
 */

#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <tclXkeylist.h>

#include "read_matrix.h"
#include "gap_globals.h"
#include "misc.h"
#include "dna_utils.h"
#include "align.h"
#include "consen.h"
#include "fort.h"
#include "align_lib.h"
#include "hash_lib.h"
#include "notedb.h"
#include "genetic_code.h"

/* nucleotide weight matrix for alignment routines */
static char *nt_order = "ACGTURYMWSKDHVB-*";
static int **nt_matrix;

/* Globals */
float consensus_cutoff = 0;
int quality_cutoff = -1;
int chem_as_double = 0;
/* char *gap_defs = NULL; */
Tcl_Obj *gap_defs = NULL;
static Tcl_Obj *defs_name;
int maxseq = 100000;
int idm = 5;
int ignore_checkdb = 0;
int consensus_mode = CONSENSUS_MODE_CONFIDENCE;
int exec_notes = 0;
int rawdata_note = 0;
int gopenval = 8;
int gextendval = 1;
int min_vector_len = 0;
int template_check_flags = 0;
double template_size_tolerance = 1;

static char *gap_defs_trace(ClientData cd, Tcl_Interp *interp,
			    char *n1, char *n2, int flags);


/* Gap4, global variable accessor functions */
int   gap4_global_get_quality_cutoff( void ) 		{ return quality_cutoff; }
void  gap4_global_set_quality_cutoff( int qc )		{ quality_cutoff=qc; }
float gap4_global_get_consensus_cutoff( void )		{ return consensus_cutoff; }
void  gap4_global_set_consensus_cutoff( float cc )	{ consensus_cutoff=cc; }
int   gap4_global_get_consensus_mode( void )		{ return consensus_mode; }
void  gap4_global_set_consensus_mode( int cm )		{ consensus_mode=cm; }
int   gap4_global_get_maxseq( void )                    { return maxseq; }
void  gap4_global_set_maxseq( int ms )                  { maxseq=ms; }
int   gap4_global_get_gopenval( void )			{ return gopenval; }
void  gap4_global_set_gopenval( int ov )		{ gopenval=ov; }
int   gap4_global_get_gextendval( void )		{ return gextendval; }
void  gap4_global_set_gextendval( int ev )		{ gextendval=ev; }
double gap4_global_get_template_size_tolerance ( void ) {
    return template_size_tolerance;
}
void   gap4_global_set_template_size_tolerance ( double ts ) {
    template_size_tolerance=ts;
}


/* TraceVar proc for consensus_cutoff */
static char *change_consensus_cutoff(ClientData clientData, Tcl_Interp *interp,
				     char *name1, char *name2, int flags) {
    char *v = Tcl_GetVar2(interp, name1, name2, TCL_GLOBAL_ONLY);
    if (v) consensus_cutoff = atof(v);

    return NULL;
}


static void init_tcl_notes(Tcl_Interp *interp) {
    int i;
    char buf[1024];

    readInNoteDB();	/* Parse and load NOTEDB */

    sprintf(buf, "%d", note_db_count);
    Tcl_SetVar2(interp, "NoteDB", "num_notes", buf, TCL_GLOBAL_ONLY);

    for (i = 0; i < note_db_count; i++) {
	sprintf(buf, "%d,type", i);
	Tcl_SetVar2(interp, "NoteDB", buf, note_db[i].type,
		    TCL_GLOBAL_ONLY);

	sprintf(buf, "%d,id", i);
	Tcl_SetVar2(interp, "NoteDB", buf, note_db[i].search_id,
		    TCL_GLOBAL_ONLY);

	sprintf(buf, "%d,dt", i);
	Tcl_SetVar2(interp, "NoteDB", buf, note_db[i].default_text,
		    TCL_GLOBAL_ONLY);
    }

    return;
}

/* Main global setup function */
int init_globals(Tcl_Interp *interp) {
    static int done_init = 0;
    extern int gap_fatal_errors;
    char *env;

    if (done_init)
	return 0;
    else
	done_init++;

    /* lookup tables */

    set_char_set(1);    /* 1 == DNA */
    set_dna_lookup(); 	/* general lookup and complementing */
    set_iubc_lookup();	/* iubc codes for restriction enzymes */
    set_hash8_lookupn();	/* used by word8 hashing */
    set_mask_lookup();  /* used to mask/mark consensus */
    init_genetic_code();
    inits_();		/* fortran stuff */
    initlu_(&idm);	/* fortran stuff */

    /* Init Tcl note database */
    init_tcl_notes(interp);

    if (NULL == (env = getenv("STADTABL")))
	verror(ERR_FATAL, "init_globals",
	       "STADTABL environment variable is not set.");
    else {
	char buf[1024];

	sprintf(buf, "%s/align_lib_nuc_matrix", env);
	nt_matrix = create_matrix(buf, nt_order);
	if (nt_matrix)
	    init_W128(nt_matrix, nt_order, 0);
	else
	    verror(ERR_FATAL, "init_globals",
		   "%s: file not found", buf);
    }

    /*
     * gap_defs (a Tcl_Obj pointer)
     *
     * We keep this up to date by creating a write trace on the object and
     * doing an ObjGetVar2 when it changes. This way the object is always
     * valid.
     * Firstly we have to create gap_defs though as initially it doesn't
     * exist.
     */
    {
	Tcl_Obj *val;

	defs_name = Tcl_NewStringObj("gap_defs", -1); /* global */

	val = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
	if (NULL == val)
	    val = Tcl_NewStringObj("", -1);

	gap_defs = Tcl_ObjSetVar2(interp, defs_name, NULL, val,
				  TCL_GLOBAL_ONLY);
	Tcl_TraceVar(interp, "gap_defs", TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,
		     gap_defs_trace, NULL);
    }

    /* consensus_cutoff */
    Tcl_TraceVar(interp, "consensus_cutoff", TCL_TRACE_WRITES|TCL_GLOBAL_ONLY,
		 change_consensus_cutoff, (ClientData)NULL);


    /* quality_cutoff */
    Tcl_LinkVar(interp, "quality_cutoff", (char *)&quality_cutoff,
		TCL_LINK_INT);

    /* chem_as_double */
    Tcl_LinkVar(interp, "chem_as_double", (char *)&chem_as_double,
		TCL_LINK_INT);


    /* gap_fatal_errors */
    Tcl_LinkVar(interp, "gap_fatal_errors", (char *)&gap_fatal_errors,
		TCL_LINK_BOOLEAN);


    /* maxseq */
    Tcl_LinkVar(interp, "maxseq", (char *)&maxseq,
		TCL_LINK_INT);

    /* maxdb */
    Tcl_LinkVar(interp, "maxdb", (char *)&maxdb,
		TCL_LINK_INT);


    /* ignore_checkdb */
    Tcl_LinkVar(interp, "ignore_checkdb", (char *)&ignore_checkdb,
		TCL_LINK_INT);

    /* consensus_mode */
    Tcl_LinkVar(interp, "consensus_mode", (char *)&consensus_mode,
		TCL_LINK_INT);

    /* consensus_iub */
    Tcl_LinkVar(interp, "consensus_iub", (char *)&consensus_iub,
		TCL_LINK_INT);

    /* exec_notes */
    Tcl_LinkVar(interp, "exec_notes", (char *)&exec_notes,
		TCL_LINK_INT);

    /* rawdata_note */
    Tcl_LinkVar(interp, "rawdata_note", (char *)&rawdata_note,
		TCL_LINK_INT);

    /* align_open_cost */
    Tcl_LinkVar(interp, "align_open_cost", (char *)&gopenval,
		TCL_LINK_INT);

    /* align_extend_cost */
    Tcl_LinkVar(interp, "align_extend_cost", (char *)&gextendval,
		TCL_LINK_INT);

    /* template_size_tolerance */
    Tcl_LinkVar(interp, "template_size_tolerance", 
		(char *)&template_size_tolerance,
		TCL_LINK_DOUBLE);

    /* min_vector_len */
    Tcl_LinkVar(interp, "min_vector_len", (char *)&min_vector_len,
		TCL_LINK_INT);

    /* template_check_flags */
    Tcl_LinkVar(interp, "template_check_flags", (char *)&template_check_flags,
		TCL_LINK_INT);


    return TCL_OK;
}

static char *gap_defs_trace(ClientData cd, Tcl_Interp *interp,
			    char *n1, char *n2, int flags) {
    gap_defs = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
    return NULL;
}
