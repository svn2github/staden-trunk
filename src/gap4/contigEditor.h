#ifndef _CEDIT_H_
#define _CEDIT_H_

#include <tk.h>

#include "fort.h"
#include "edStructs.h"

extern int edit_contig(Tcl_Interp *interp, GapIO *io, int cnum, int llino,
		       int pos, float con_cut, int qual_cut,
		       int reveal_cutoffs, char *sets);
extern int join_contig(Tcl_Interp *interp, GapIO *io, int cnum[2],
		       int llino[2], int pos[2], float con_cut, int qual_cut);

#define EDITMODE 1
#define JOINMODE 2

extern int inJoinMode(EdStruct *xx);
extern int editorLocked(EdStruct *xx);
extern int editorLockedPos(EdStruct *xx[2], int force);
extern void db_callback_tk(void *xxv, int type, int seq, int pos,
			   void *pointer);

int Ced_Init(Tcl_Interp *interp);

#endif
