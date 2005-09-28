#ifndef _gapDB_h
#define _gapDB_h

#include "list.h"

extern void gap_open_for_read(List *l);
extern void gap_close(List *l);
extern void gap_open_for_write(List *l);
extern void gap_write_header(List *l);
extern void gap_write_gel_data(List *l);
extern void gap_write_contig_data(List *l);
extern List *gap_read_header(void);
extern List *gap_read_gel_data(void);
extern List *gap_read_contig_data(void);

#endif /* _gapDB_h */

