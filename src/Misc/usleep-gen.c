#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>

/*
 * This sleeps for `usecond' microseconds. This implementation uses select to
 * wait for input on no file descriptors with `useconds' timeout. However it
 * does not wake early if it gets hit by a signal.
 */
int usleep(unsigned int useconds) {
    struct timeval tv;

    tv.tv_sec  = useconds / 1000000;
    tv.tv_usec = useconds % 1000000;

    if (-1 == select(0, NULL, NULL, NULL, &tv))
	return -1;
    else
	return 0;
}

typedef int int_f;

/*
 * FORTRAN interface
 */
void usleep_(int_f *useconds) {
    (void)usleep((unsigned int)*useconds);
}
