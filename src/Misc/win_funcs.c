/*
 * Implementations of various POSIX functions that are missing when using
 * mingw on windows.
 */
#include "win_funcs.h"

size_t pread(int fd, void *buf, size_t count, off_t offset) {
  int ret;
  off_t curr = lseek(fd, 0, SEEK_CUR);

  if (lseek(fd, offset, SEEK_SET) == -1)
    return -1;
  ret = read(fd, buf, count);
  lseek(fd, curr, SEEK_SET);

  return ret;
}

size_t pwrite(int fd, const void *buf, size_t count, off_t offset) {
  int ret;
  off_t curr = lseek(fd, 0, SEEK_CUR);

  if (lseek(fd, offset, 0) == -1)
    return -1;
  ret = write(fd, buf, count);
  lseek(fd, curr, SEEK_SET);

  return ret;
}
