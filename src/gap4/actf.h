#ifndef _ACTF_H
#define _ACTF_H

extern int actfptr;

int_f actf_(int_f *JOB_p,
	   char *FILNAM_p,
	   int_f *FILNAMLEN_p,
	   char *COPYNUM_p,
	   int_f *KBOUT,
	   int_fl  FILNAM_l,
	   int_fl  COPYNUM_l);

/*
 * actf_lock(int mode, char *file, char *version)
 *    Creates the "file.version.BUSY" file.
 *    read_only is 0 for read/write and 1 for read.
 *    Returns 0 for success or >0 for error (see above)
 */
int actf_lock(int read_only, char *file, char *version, int new);

/*
 * actf_unlock(int readonly, char *file, char *version)
 *    Removes the busy file.
 *    Returns 0 for success or >0 for error (see above)
 */
int actf_unlock(int read_only, char *file, char *version);

#endif
