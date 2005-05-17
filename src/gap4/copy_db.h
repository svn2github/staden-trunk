#ifndef _COPY_DB_H_
#define _COPY_DB_H_

#include "IO.h"

/*
 * Copies database pointed to by 'iof' to the database pointed to by 'iot'.
 *
 * If 'verbose' is set then a running summary of its actions is displayed.
 * If 'errs' is set then we stop on errors, otherwise we'll try to ignore
 * them as best as we can.
 * If 'notags' is defined then no tags are copied from iof. This automatically
 * sets tag lists in readings and contigs to be zero.
 */
int copy_database(GapIO *iof, GapIO *iot, int verbose, int errs, int notags);

/*
 * Interface to copy_database. This is called when we already have a database
 * open and wish to copy it to a new database.
 */
int copy_database_from(GapIO *io, char *name, char *version);


#endif
