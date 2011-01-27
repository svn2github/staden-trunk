#ifndef _WIN_FUNCS_H_
#define _WIN_FUNCS_H_

#include <unistd.h>
#include <fcntl.h>

#include "os.h"

size_t pread(int fd, void *buf, size_t count, off_t offset);
size_t pwrite(int fd, const void *buf, size_t count, off_t offset);

#endif /* _WIN_FUNCS_H_ */
