#ifndef _FEATURE_SELECTOR_H_
#define _FEATURE_SELECTOR_H_

#include <X11/Intrinsic.h>

typedef struct {
    /* Must be first item */
    char *type;

    /* values taken from TAGDB file */
    int index; 
    char *fill;
    char *outline;
    int width;
    char *shape;

    /* values derived from above */
    Pixel fg_pixel;
    Pixel bg_pixel;
    Pixel gf_pixel;
    Pixel gb_pixel;
} feature_db_struct;


extern feature_db_struct *feature_db;
extern int feature_db_count;

void read_feature_file(char *filename);

int save_feature_file(Tcl_Interp *interp, FILE *fp, char **items, int num_items);
feature_db_struct* klist_feature(Tcl_Interp *interp, char **items, int num_items);

#endif /* _FEATURE_SELECTOR_H_ */
