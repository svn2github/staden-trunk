#ifndef _TAGDB_H_
#define _TAGDB_H_

#include "intrinsic_type.h"

typedef struct {
    /* Must be first item */
    char *type;

    /* values taken from TAGDB file */
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
} tag_db_struct;


extern tag_db_struct *tag_db;
extern int tag_db_count;

void readInTagDB(void);

#endif /* _TAGDB_H_ */
