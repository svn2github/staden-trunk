/*
 * cah.h - caf to gap5 conversion for tg_index
 *
 * Andrew Whitwham, August 2010
 * Wellcome Trust Sanger Institute
 *
 */

#ifndef _CAF_H_
#define _CAF_H_

#include <hache_table.h>
#include <tg_index.h>

int parse_caf(GapIO *io, char *fn, tg_args *a);

#endif
