#ifndef _NOTEDB_H_
#define _NOTEDB_H_

#include "intrinsic_type.h"

typedef struct {
    /* Must be first item */
    char *type;

    /* values taken from NOTEDB file */
    char *search_id;
    char *default_text;
    char *fg_colour;
    char *bg_colour;
    char *gf_colour;
    char *gb_colour;

    /* values derived from above */
    char id[4];
    Pixel fg_pixel;
    Pixel bg_pixel;
    Pixel gf_pixel;
    Pixel gb_pixel;
} note_db_struct;


extern note_db_struct *note_db;
extern int note_db_count;

void readInNoteDB(void);

/*
 * Converts a note type string (eg "REFS") to an index into the database.
 * Returns -1 if none found.
 */
int note_id2index(char *id);

#endif /* _NOTEDB_H_ */
