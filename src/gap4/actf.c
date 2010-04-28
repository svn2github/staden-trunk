#include <staden_config.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/*
 * Windows defines ERROR somewhere, which IO.h also does. This is just to avoid
 * a harmless compilation warning.
 */
#ifdef ERROR
#undef ERROR
#endif

#include "fort.h"
#include "xalloc.h"
#include "misc.h"
#include "actf.h"
#include "text_output.h"

/* static error messages for use in actf_() */
static char *actferrlist[] = {
    "unknown error",
    "Error! - unknown JOB to ACTF()",
    "Error creating busy file",
    "Error deleting busy file",
    "Sorry, database busy",
    "Misc. error",
    "Database not found"
    };

/*
 * Display an error.
 */
static void actferr(int errnum) {
    verror(ERR_WARN, "lock-database", "%s\n", actferrlist[errnum-1]);
    return;
}

typedef struct {
    char *pathname;
    char *db_name;
    int fd;
} lock_file_t;
static lock_file_t *lock_files = NULL;
static int numa_lock_files = 0;
static int numu_lock_files = 0;


/*
 *-----------------------------------------------------------------------------
 * C interface to the BUSY file code.
 *-----------------------------------------------------------------------------
 */

/*
 * actf_lock(int mode, char *file, char *version)
 *    Creates the "file.version.BUSY" file.
 *    read_only is 0 for read/write and 1 for read.
 *    Returns 0 for success or >0 for error (see above)
 */
/* 15/1/99 johnt - gethostname returns SOCKET_ERROR on WINNT */
#ifdef _WIN32
#include <winsock.h>
#else
#define SOCKET_ERROR -1
#endif

#if !defined(NOLOCKF)
int test_if_locked(char *fname) {
    int locked = 0;
    int fd;

    if ((fd = open(fname, O_RDONLY, 0)) != -1) {
	int i;
	locked = 1;
	/*
	 * If we locked this file ourselves then lockf() will succeed,
	 * so test for this case also. We do this simply by looking through
	 * the list of locked files to see if it is one we already know about.
	 */
	for (i = 0; i < numu_lock_files; i++) {
	    if (strcmp(lock_files[i].pathname, fname) == 0) {
		break;
	    }
	}

	/*
	 * If i == numu_lock_files then we are not using this db
	 * ourselves. Use lockf() to check if another process has
	 * this db open.
	 */
	if (i == numu_lock_files && lockf(fd, F_TEST, 0) == 0) {
	    locked = 0;
	}
	close(fd);
    }

    return locked;
}
#endif

int actf_lock(int read_only, char *file, char *version,int new) {
    char fname[2048];
    struct stat statbuf;
    int fd;
    char dir[1025];
    char db_name[1025];
    char *cp;
    char db_path[2048];
    char aux_path[2048];
    char hostname[1024];
    int locked;

    /*
     * This may appear useless as we're using dir for creating the BUSY file,
     * which could be done relatively anyway.
     * However Gap4 has an option to change directory, so we need to store
     * the full pathname to enable us to remove the BUSY file when closing
     * the database.
     *
     * NB this code doesn't work on Windows, which has pathnames like
     * "C:/fish/xyzzy.0" (and "C:xyzzy.0", which is a relative pathname).
     * However the windows dialogues will have returned a full pathname anyway
     * so this is only really an issue for starting Gap4 on the command line.
     */

#ifdef _WIN32
    /* 6/1/99 johnt - always seem to have a full path, and stat should
     * work fine on a relative path anyway, so don't do full path extension
     * as it is problematic under windows
     */
    *dir=0;
#else
    if (*file != '/') {
	if (NULL == getcwd(dir, 1024)) {
	    *dir = 0;
	} else {
	    strcat(dir, "/");
	}
    } else {
	*dir = 0;
    }
#endif

    if (cp = strrchr(file, '/')) {
	sprintf(db_name, "%s.%s", cp+1, version);
    } else {
	sprintf(db_name, "%s.%s", file, version);
    }
    sprintf(db_path,  "%s.%s",      file, version);
    sprintf(aux_path, "%s.%s.aux",  file, version);
    sprintf(fname,    "%s%s.%s.BUSY", dir, file, version);


    /* Check for existance of lock on the BUSY file (from older gap4s) */
    locked = 0;
    if (stat(fname, &statbuf) != -1) {
#if !defined(NOLOCKF)
	locked = test_if_locked(fname);
	if (!locked) {
	    vmessage("WARNING! Database has lock file, "
		     "but is no longer in use.\n");
	    log_file(NULL, "Overriding lock file");
	    if (!read_only)
		vmessage("WARNING! Taking ownership of lock.\n");
	}
#else
	locked = 1;
#endif
    }

#if 0
    /*
     * WARNING: The following code (along with code in g/g-files.c) checks
     * for a lock on the database itself. HOWEVER this fails in some cases
     * over NFS claiming that a database is locked when it is not (due to
     * server crashes/resets). Perhaps this is bugs in rpc.lockd/statd, but
     * regardless of where the blame lies, for now it causes more problems
     * than it solves.
     */

    /* Check for the existance of a lock on the DB file itself */
    if (!locked) {
	locked = test_if_locked(db_path);
    }
#endif

    if (locked) {
	if (read_only) {
	    vmessage("WARNING! Database is currently in use\n");
	    return 0;
	} else {
	    /* lock already exists */
	    actferr(5);
	    return 5;
	}
    }

    /* Nothing more to do in read only mode */
    if (read_only)
	return 0;

    if (numu_lock_files >= numa_lock_files) {
	numa_lock_files += 10;
	lock_files = (lock_file_t *)xrealloc(lock_files, numa_lock_files *
					     sizeof(lock_file_t));
	if (lock_files == NULL) {
	    actferr(6);
	    return 6;
	}
    }

    /* Sanity check - done the database actually exist? */
    if (!new) {
        if (stat(db_path, &statbuf) == -1 ||
            stat(aux_path, &statbuf)== -1) {
            actferr(7);
            return 7;
        }
    }

    /* Create the BUSY file */
    if ((fd = open(fname, O_CREAT | O_RDWR | O_TRUNC, 0666)) == -1) {
	actferr(3);
	return 3;
    }
    /* Lock it */
#if !defined(NOLOCKF)
    lockf(fd, F_LOCK, 0);
#endif

#if defined(__MINGW32__)
   strcpy(hostname, "unknown");
#else
    /* Write the hostname and process id to it */
    if (SOCKET_ERROR == gethostname(hostname, 1023))
        strcpy(hostname, "unknown");
    hostname[1023] = 0;
#endif
    sprintf(db_path, "%s %d\n", hostname, (int)getpid());
    write(fd, db_path, strlen(db_path));

    lock_files[numu_lock_files].pathname = strdup(fname);
    lock_files[numu_lock_files].db_name = strdup(db_name);
    lock_files[numu_lock_files].fd = fd;
    numu_lock_files++;

    return 0;
}

/*
 * actf_unlock(int readonly, char *file, char *version)
 *    Removes the busy file.
 *    Returns 0 for success or >0 for error (see above)
 */
int actf_unlock(int read_only, char *file, char *version) {
    int i;
    char db_name[1024], *cp;

    /* Nothing to do for read only mode */
    if (read_only)
	return 0;

    /* do the unlocking */
    if (cp = strrchr(file, '/')) {
	sprintf(db_name, "%s.%s", cp+1, version);
    } else {
	sprintf(db_name, "%s.%s", file, version);
    }
    for (i = 0; i < numu_lock_files; i++) {
	if (0 == strcmp(db_name, lock_files[i].db_name))
	    break;
    }

    if (i == numu_lock_files) {
	actferr(4);
	return 4;
    }

    close(lock_files[i].fd); /* Also clears lockf() */

    if (unlink(lock_files[i].pathname) == -1) {
	actferr(4);
	return 4;
    }

    free(lock_files[i].pathname);
    free(lock_files[i].db_name);

    memcpy(&lock_files[i].pathname, &lock_files[i+1].pathname,
	   (numu_lock_files - i - 1) * sizeof(*lock_files));
    numu_lock_files--;

    return 0;
}

