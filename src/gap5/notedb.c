#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parse_db.h"
#include "notedb.h"
#include "misc.h"

note_db_struct *note_db = NULL;
int note_db_count = 0;

/*
 * ---------------------------------------------------------------------------
 * NOTEDB parse specifications.
 * ---------------------------------------------------------------------------
 */

static void notedb_parse(char *fn) {
    pf_spec a[] = {
	{"id",	     PF_STR, offsetof(note_db_struct, search_id)},
	{"dt",	     PF_STR, offsetof(note_db_struct, default_text)},
        {"fg",       PF_STR, offsetof(note_db_struct, fg_colour)},
        {"bg",       PF_STR, offsetof(note_db_struct, bg_colour)},
        {"gf",       PF_STR, offsetof(note_db_struct, gf_colour)},
        {"gb",       PF_STR, offsetof(note_db_struct, gb_colour)},
	{NULL,	     0,	      0}
    };

    note_db = (note_db_struct *)parse_file(fn, a, note_db, &note_db_count,
					   sizeof(*note_db), NULL);
}


static void tidyUpNoteDBFields(int note)
{
    int len;
    
    if (note_db[note].search_id == NULL)
	note_db[note].search_id = note_db[note].type;

    len =  strlen(note_db[note].search_id);
    if (len < 4)
	strncpy(note_db[note].id,"    ",4);
    else
	len = 4;
    strncpy(note_db[note].id,note_db[note].search_id,len);

    if (note_db[note].gf_colour == NULL && note_db[note].fg_colour!=NULL)
        note_db[note].gf_colour = strdup(note_db[note].fg_colour);

    if (note_db[note].fg_colour == NULL && note_db[note].gf_colour!=NULL)
        note_db[note].fg_colour = strdup(note_db[note].gf_colour);

    if (note_db[note].gb_colour == NULL && note_db[note].bg_colour!=NULL)
        note_db[note].gb_colour = strdup(note_db[note].bg_colour);

    if (note_db[note].bg_colour == NULL && note_db[note].gb_colour!=NULL)
        note_db[note].bg_colour = strdup(note_db[note].gb_colour);
}

/*
 * ---------------------------------------------------------------------------
 * The main external function.
 *
 * We break NOTEDB into segments separated by colons. We read backwards so that
 * the first (left) entries have top priority. All files read are merged.
 * ---------------------------------------------------------------------------
 */

void readInNoteDB(void)
{
    char *path, *p;
    char tmp_path[2000];
    int i;
    int gotfile=0;

    if (NULL == (path = getenv("NOTEDB"))) {
	if (getenv("STADTABL")) {
	    strcpy(tmp_path, getenv("STADTABL"));
	    strcat(tmp_path, "/NOTEDB");
	    path = tmp_path;
	} else {
	    path = "NOTEDB";
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
	    notedb_parse(p);
	    gotfile++;
	}
    } while (p != path);

    for (i = 0; i < note_db_count; i++) {
	tidyUpNoteDBFields(i);
    }

    /* 22/1/99 johnt - give a warning if no files found - as may cause errors later in program */
    if( !gotfile )
	verror(ERR_WARN, "Note DB", "No Files found - please check NOTEDB environment variable.");

}

/*
 * Converts a note type string (eg "REFS") to an index into the database.
 * Returns -1 if none found.
 */
int note_id2index(char *id)
{
    int i;

    if (id==NULL)
	return -1;

    for (i=0; i<note_db_count; i++) {
        if (strncmp(id,note_db[i].id,4)==0)
            return i;
    }

    return -1;
}
