#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>
#include <tk.h>

#include "feature_colour.h"
#include "misc.h"
#include "parse_db.h"

#define TAGDB "FEATCOLDB"
feat_col_db *fcol_db = NULL;
int fcol_db_count = 0;

static void fcoldb_parse(char *fn) {
    pf_spec a[] = {
	{"id",	     PF_STR, offsetof(feat_col_db, search_id)},
	{"fg",	     PF_STR, offsetof(feat_col_db, fg_colour)},
	{"bg",	     PF_STR, offsetof(feat_col_db, bg_colour)},
	{"gf",	     PF_STR, offsetof(feat_col_db, gf_colour)},
	{"gb",	     PF_STR, offsetof(feat_col_db, gb_colour)},
	{"dt",	     PF_STR, offsetof(feat_col_db, default_text)},
	{NULL,	     0,	      0}
    };

    fcol_db = (feat_col_db *)parse_file(fn, a, fcol_db, &fcol_db_count,
					 sizeof(*fcol_db), NULL);
}

static void tidyUpFcolDBFields(int fcol) {

    int len;
    
    if (fcol_db[fcol].search_id == NULL)
	fcol_db[fcol].search_id = fcol_db[fcol].type;

    len =  strlen(fcol_db[fcol].search_id);
    if (len < 4)
	strncpy(fcol_db[fcol].id,"    ",4);
    else
	len = 4;
    strncpy(fcol_db[fcol].id,fcol_db[fcol].search_id,len);
    
    if (fcol_db[fcol].gf_colour == NULL && fcol_db[fcol].fg_colour!=NULL)
	fcol_db[fcol].gf_colour = strdup(fcol_db[fcol].fg_colour);

    if (fcol_db[fcol].fg_colour == NULL && fcol_db[fcol].gf_colour!=NULL)
	fcol_db[fcol].fg_colour = strdup(fcol_db[fcol].gf_colour);

    if (fcol_db[fcol].gb_colour == NULL && fcol_db[fcol].bg_colour!=NULL)
	fcol_db[fcol].gb_colour = strdup(fcol_db[fcol].bg_colour);

    if (fcol_db[fcol].bg_colour == NULL && fcol_db[fcol].gb_colour!=NULL)
	fcol_db[fcol].bg_colour = strdup(fcol_db[fcol].gb_colour);
}

void readInFcolDB(void)
{
    char *path, *p;
    char tmp_path[2000];
    int i;
    int gotfile=0;

    /* 22/1/99 johnt - fallback to to use STADTABL is GFCOLDB isn't defined */
    if (NULL == (path = (char *)getenv(TAGDB))){
        if(getenv("STADTABL")) {
	   strcpy(tmp_path,getenv("STADTABL"));
           strcat(tmp_path,"/");
	   strcat(tmp_path,TAGDB);
           path = tmp_path;
	} else {
	    path = TAGDB;
	}
    }

    do {
	p = strrchr(path, PATHSEP); /* 22/1/99 johnt - path seperator now macro to allow WINNT support */
	if (p) {
	    *p = 0;
	    p++;
	} else
	    p = path;
	    
	if (file_exists(p)) {
	    fcoldb_parse(p);
            gotfile++;
	}
    } while (p != path);
    for (i = 0; i < fcol_db_count; i++) {
	tidyUpFcolDBFields(i);
    }
    if( !gotfile )
	verror(ERR_WARN, "Tag DB", "No Files found - please check GTAGDB environment variable.");

}

void get_fcol_types (void) {
    static int done = 0;

    if (!done) {
	readInFcolDB();
	done = 1;
    }
}
int tcl_get_fcol_array(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    Tcl_DString fcols;
    int i;

    get_fcol_types();

    Tcl_DStringInit(&fcols);
    for (i = 0; i < fcol_db_count; i++) {
	Tcl_DStringStartSublist(&fcols);
	Tcl_DStringAppendElement(&fcols, fcol_db[i].type);
	Tcl_DStringAppendElement(&fcols, fcol_db[i].bg_colour);
	Tcl_DStringEndSublist(&fcols);
    }
    Tcl_DStringResult(interp, &fcols);

    return TCL_OK;
}

int tcl_get_key_array(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    Tcl_DString fkey;
    int i;

    if (!fcol_db) {
	get_fcol_types();
    }
    Tcl_DStringInit(&fkey);
    for (i = 0; i < fcol_db_count; i++) {
	Tcl_DStringStartSublist(&fkey);
	Tcl_DStringAppendElement(&fkey, fcol_db[i].type);
	Tcl_DStringEndSublist(&fkey);
    }
    Tcl_DStringResult(interp, &fkey);

    return TCL_OK;
}
