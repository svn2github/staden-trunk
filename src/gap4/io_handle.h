#ifndef _IO_HANDLE_H
#define _IO_HANDLE_H

int get_free_handle(GapIO *io);

void remove_handle(f_int *handle);

GapIO *io_handle(f_int *HANDLE);

f_int *handle_io(GapIO *io);

#endif /* _IO_HANDLE_H */
