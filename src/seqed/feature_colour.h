#if !defined(FEAT_COLOUR_H)
#define FEAT_COLOUR_H

#include <X11/Intrinsic.h>

typedef struct {
    /* Must be first item */
    char *type;
    /* values taken from FEATCOLDB file */
    char *search_id;
    char *fg_colour;
    char *bg_colour;
    char *gf_colour;
    char *gb_colour;
    char *default_text;

    /* values derived from above */
    Pixel fg_pixel;
    Pixel bg_pixel;
    Pixel gf_pixel;
    Pixel gb_pixel;
    char id[4];
} feat_col_db;


extern feat_col_db *fcol_db;
extern int fcol_db_count;

void readInFColDB(void);

/*static void fcoldb_parse(char *fn);
  static void tidyUpFcolDBFields(int tag);*/
void readInFColDB(void);
void get_fcol_types (void);
int tcl_get_fcol_array(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv);
int tcl_get_key_array(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv);
#endif

