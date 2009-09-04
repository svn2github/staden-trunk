#ifndef _bapDB_h
#define _bapDB_h

#include "list.h"

extern void xdap_late_open_for_read(List *l);
extern void xdap_late_close(List *l);
extern void xdap_late_open_for_write(List *l);
extern void xdap_late_write_header(List *l);
extern void xdap_late_write_gel_data(List *l);
extern void xdap_late_write_contig_data(List *l);
extern List *xdap_late_read_header(void);
extern List *xdap_late_read_gel_data(void);
extern List *xdap_late_read_contig_data(void);

#endif /* _bapDB_h */

