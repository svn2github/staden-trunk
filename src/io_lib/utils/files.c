/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#include "misc.h"

#include <sys/types.h>
#include <sys/stat.h>
/* Alliant's Concentrix <sys/stat.h> is hugely deficient */
/* Define things we require in this program              */
/* Methinks S_IFMT and S_IFDIR aren't defined in POSIX   */
#ifndef S_ISDIR
#define S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
#endif /*!S_ISDIR*/
#ifndef S_ISREG
#define S_ISREG(m)      (((m)&S_IFMT) == S_IFREG)
#endif /*!S_ISREG*/

int is_directory(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISDIR(buf.st_mode);
}

int is_file(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISREG(buf.st_mode);
}

int file_exists(char * fn)
{
    struct stat buf;
    return ( stat(fn,&buf) == 0);
}

int compressed_file_exists(char *fname)
{
    struct stat buf;
    char fn[2048];

    if (stat(fname, &buf) == 0) return 1;

    sprintf(fn, "%s.gz", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.bz", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.bz2", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.Z", fname);
    if (stat(fn, &buf) == 0) return 1;

    sprintf(fn, "%s.z", fname);
    if (stat(fn, &buf) == 0) return 1;

    return 0;
}

int file_size(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) != 0) return 0;
    return buf.st_size;
}

/*
 * ---------------------------------------------------------------------------
 * File of filename management
 */

FILE *open_fofn(char *files) {
  return fopen(files, "r");
}

char *read_fofn(FILE *fp) {
  char line[256];
  static char name[256];
  
  while (fgets(line, 254, fp)) {
    if (1 == sscanf(line, "%s", name))
      return name;
  }

  return NULL;
}

void close_fofn(FILE *fp) {
  fclose(fp);
}
