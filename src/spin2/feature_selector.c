#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <tcl.h>
#include <tclXkeylist.h>

#include "parse_db.h"
#include "feature_selector.h"
#include "xalloc.h"
#include "tcl_utils.h"

feature_db_struct *feature_db = NULL;
int feature_db_count = 0;

static void featuredb_parse(char *fn) {

    pf_spec a[] = {
	{"bg",	     PF_STR, offsetof(feature_db_struct, fill)},
	{"bc",	     PF_STR, offsetof(feature_db_struct, outline)},
	{"bw",	     PF_INT, offsetof(feature_db_struct, width)},
	{"sh",	     PF_STR, offsetof(feature_db_struct, shape)},
	{NULL,	     0,	      0}
    };

    feature_db = (feature_db_struct *)parse_file(fn, a, feature_db, 
						 &feature_db_count,
						 sizeof(*feature_db), NULL);
}


void read_feature_file(char *filename)
{
    feature_db = NULL;
    feature_db_count = 0;

    featuredb_parse(filename);
}


/*
 * convert from tcl feature keyed list to a feature_db_struct
 * NB function w() can be found in tk_utils/tcl_utils.c and converts a 
 * non-writable string into a writeable one
 */
#define getli(s, x) \
do { \
    if (TCL_OK == TclX_KeyedListGet(interp, list, w(#x), &ptr)) { \
        (void)Tcl_GetIntFromObj(interp, ptr, &(s).x); \
    } \
} while (0);

#define getls(s, x) \
do { \
    if (TCL_OK == TclX_KeyedListGet(interp, list, w(#x), &ptr)) { \
        (s).x = Tcl_GetStringFromObj(ptr, NULL); \
    } \
} while (0);

/*
 * convert a tcl keyed list into a feature db
 */
feature_db_struct* klist_feature(Tcl_Interp *interp,
				 char **items, 
				 int num_items)
{
    feature_db_struct *feat_db;
    int i;
    Tcl_Obj *ptr;
    Tcl_Obj *list;

    if (NULL == (feat_db = (feature_db_struct *)xmalloc(num_items * 
						   sizeof(feature_db_struct))))
	return NULL;
    
    for (i = 0; i < num_items; i++) {
	if (NULL == (feat_db[i].type = (char *)xmalloc(16 * sizeof(char))))
	    return NULL;
	if (NULL == (feat_db[i].fill = (char *)xmalloc(20 * sizeof(char))))
	    return NULL;
	if (NULL == (feat_db[i].outline = (char *)xmalloc(20 * sizeof(char))))
	    return NULL;
	if (NULL == (feat_db[i].shape = (char *)xmalloc(20 * sizeof(char))))
	    return NULL;

	list = Tcl_NewStringObj(items[i], -1);
	getls(feat_db[i], type);
	getli(feat_db[i], index);
	getls(feat_db[i], fill);
	getls(feat_db[i], outline);
	getli(feat_db[i], width);
	getls(feat_db[i], shape);
    }
    return feat_db;

}

int save_feature_file(Tcl_Interp *interp,
		      FILE *fp,
		      char **items,
		      int num_items)
{
    int i;
    feature_db_struct *feat_db;

    if (NULL == (feat_db = (feature_db_struct *)xmalloc(num_items * 
						   sizeof(feature_db_struct))))
	return -1;

    feat_db = klist_feature(interp, items, num_items);

    for (i = 0; i < num_items; i++) {

	fprintf(fp, "\"%s\": \tbg=\"%s\": \tbc=\"%s\": \tbw=\"%d\": \tsh=\"%s\"\n", feat_db[i].type, feat_db[i].fill, feat_db[i].outline,
		feat_db[i].width, feat_db[i].shape);
    }

    fclose(fp);
    return 0;
}
