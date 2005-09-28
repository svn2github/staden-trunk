#ifndef _flat_sd_h
#define _flat_sd_h

#include "list.h"

extern void flat_sd_open_for_write(List *l);
extern void flat_sd_close(List *l);
extern void flat_sd_write_contig_data(List *l);
extern void flat_sd_write_gel_data(List *l);
extern void flat_sd_write_header(List *l);
extern void flat_sd_open_for_read(List *l);
extern List *flat_sd_read_header(void);
extern List *flat_sd_read_contig_data(void);
extern List *flat_sd_read_gel_data(void);

#endif /* _flat_sd_h */
