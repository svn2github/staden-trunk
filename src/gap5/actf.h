#ifndef _ACTF_H
#define _ACTF_H

extern int actfptr;

/*
 * Creates the "file.BUSY" file.
 * read_only is 0 for read/write and 1 for read.
 *
 * Returns 0 for success or >0 for error (see above)
 */
int actf_lock(int read_only, char *file, int new);

/*
 * Removes the busy file.
 *
 * Returns 0 for success or >0 for error (see above)
 */
int actf_unlock(int read_only, char *file);

#endif
