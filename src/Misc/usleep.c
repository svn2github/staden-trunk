#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#ifdef _MSC_VER
#include <windows.h>       /* For Sleep() */
#else
#include <sys/time.h>      /* Unavailable on windows */
#endif


/*
 * This sleeps for `usecond' microseconds. This implementation uses select to
 * wait for input on no file descriptors with `useconds' timeout. However it
 * does not wake early if it gets hit by a signal.
 */
int myusleep(unsigned int useconds) {
#ifdef _MSC_VER
    unsigned int mseconds = useconds/1000;
    if(!mseconds)
	mseconds=1;
    Sleep(mseconds);
    return 0;
#else
#ifdef __MINGW32__
    /* I do not know how to do this under mingw yet */
#else
    struct timeval tv;

    tv.tv_sec  = useconds / 1000000;
    tv.tv_usec = useconds % 1000000;

    if (-1 == select(0, NULL, NULL, NULL, &tv))
	return -1;
    else
	return 0;
#endif
#endif
}

typedef int int_f;

/*
 * FORTRAN interface
 */
void myuslep_(int_f *useconds) {
    (void)myusleep((unsigned int)*useconds);
}
