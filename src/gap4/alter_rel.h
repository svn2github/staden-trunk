#ifndef _ALTER_REL_H_
#define _ALTER_REL_H_

#include <tcl.h>
#include "IO.h"

void delete_contig(GapIO *io, int contig);
void shift_readings(GapIO *io, int gel, int distance);
int reset_contig_order(GapIO *io);
int annotation_address(GapIO *io, int init, int search, int *type, int *num);

#endif
