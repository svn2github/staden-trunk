#include <stdio.h>
#include <string.h>

#include "parse_db.h"
#include "tagdb.h"
#include "misc.h"

tag_db_struct *tag_db = NULL;
int tag_db_count = 0;

/*
 * ---------------------------------------------------------------------------
 * TAGDB parse specifications.
 * ---------------------------------------------------------------------------
 */

static void tagdb_parse(char *fn) {
    pf_spec a[] = {
	{"id",	     PF_STR, offsetof(tag_db_struct, search_id)},
	{"fg",	     PF_STR, offsetof(tag_db_struct, fg_colour)},
	{"bg",	     PF_STR, offsetof(tag_db_struct, bg_colour)},
	{"gf",	     PF_STR, offsetof(tag_db_struct, gf_colour)},
	{"gb",	     PF_STR, offsetof(tag_db_struct, gb_colour)},
	{"dt",	     PF_STR, offsetof(tag_db_struct, default_text)},
	{NULL,	     0,	      0}
    };

    tag_db = (tag_db_struct *)parse_file(fn, a, tag_db, &tag_db_count,
					 sizeof(*tag_db), NULL);
}


static void tidyUpTagDBFields(int tag)
{
    int len;
    
    if (tag_db[tag].search_id == NULL)
	tag_db[tag].search_id = tag_db[tag].type;

    len =  strlen(tag_db[tag].search_id);
    if (len < 4)
	strncpy(tag_db[tag].id,"    ",4);
    else
	len = 4;
    strncpy(tag_db[tag].id,tag_db[tag].search_id,len);
    
    if (tag_db[tag].gf_colour == NULL && tag_db[tag].fg_colour!=NULL)
	tag_db[tag].gf_colour = strdup(tag_db[tag].fg_colour);

    if (tag_db[tag].fg_colour == NULL && tag_db[tag].gf_colour!=NULL)
	tag_db[tag].fg_colour = strdup(tag_db[tag].gf_colour);

    if (tag_db[tag].gb_colour == NULL && tag_db[tag].bg_colour!=NULL)
	tag_db[tag].gb_colour = strdup(tag_db[tag].bg_colour);

    if (tag_db[tag].bg_colour == NULL && tag_db[tag].gb_colour!=NULL)
	tag_db[tag].bg_colour = strdup(tag_db[tag].gb_colour);


}

/*
 * ---------------------------------------------------------------------------
 * The main external function.
 *
 * We break TAGDB into segments separated by colons. We read backwards so that
 * the first (left) entries have top priority. All files read are merged.
 * ---------------------------------------------------------------------------
 */
#define TAGDB "GTAGDB"

void readInTagDB(void)
{
    char *path, *p;
    char tmp_path[2000];
    int i;
    int gotfile=0;

    /* 22/1/99 johnt - fallback to to use STADTABL is GTAGDB isn't defined */
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
	    tagdb_parse(p);
            gotfile++;
	}
    } while (p != path);

    for (i = 0; i < tag_db_count; i++) {
	tidyUpTagDBFields(i);
    }

    /* 22/1/99 johnt - give a warning if no files found - as may cause errors later in program */
    if( !gotfile )
	verror(ERR_WARN, "Tag DB", "No Files found - please check GTAGDB environment variable.");

}

