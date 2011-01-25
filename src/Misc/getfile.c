#include <staden_config.h>

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "staden_config.h"
#ifdef HAVE_PWD_H
#    include <pwd.h>
#endif

#include "misc.h"
#include "getfile.h"

int expandpath( char *namein, char *nameout);
int expandname ( char *nameout, char *namein );
int my_access(char *filename, char mode);
/* FILE *my_fopen ( char *filename, char *mode); */
/* void errout( char *message); */

int my_access(char *filename, char mode)
    
    /*
     * mode
     * r read access
     * w write access
     *
     * return: NULL == fail
     * 1 == successful read only access
     * 2 == successful write access and file exists
     * 3 == successful write access and file does not exist
     */
    
{
    if ( mode == 'r' ) {
	if ( access (filename, R_OK) == 0 && access (filename, X_OK) )
	    return 1;
	else
	    return 0;
    }
    if ( mode == 'w' ) {
	if ( access (filename, F_OK) == 0) {
	    
	    /*    file exists */
	    
	    if ( access (filename, W_OK) == 0 && access (filename, X_OK) )
		return 2;
	}
	else {
	    
	    /* file does not exist, can i write in its directory? */
	    
	    if ( access (filename, X_OK) )
		return 3;
	}
    }
    return 0;
}

/* function to expand path names:
 * 
 * if name starts with the following symbols
 * 
 * /            unchanged
 * ~/           "home directory name"
 * ~fred        "home directory path for fred"
 * ~fred/etc    "home directory path for fred/etc"
 *
 * for all other cases look in the environment variables for the strings
 * between the /'s and expand them if found. Environment can include $ but
 * it is not necessary
 *
 * Arguments
 *     namein	- unexpanded path
 *     nameout  - expanded path (size FILENAME_MAX+1)
 *
 * Returns:
 *    1 == success, 0 == failure
 */
int expandpath(char *namein, char *nameout) {
    char tempbuf[FILENAME_MAX+1], tempb[FILENAME_MAX+1];
    char *tokptr, *strptr = tempbuf, *slash;
    char *nameoutp = nameout;
    int len = FILENAME_MAX, tmp, i;

    if ( namein == NULL ) return 0;

/*    if (strlen(namein) > FILENAME_MAX) return 0; */
    strncpy ( tempbuf, namein, FILENAME_MAX);
    *nameout = '\0';
    
    /* nasty special case: namein starts with  "/" and it
       gets lost so stick it in nameout */
    
    for (i = 0; namein[i] == '/'; i++) {
	*nameoutp++ = '/';
	len--;
    }
    while ( len > 0 && ( tokptr = strtok ( strptr, "/" )) ) {
	if ( expandname ( tempb, tokptr ) == 0 ) return 0;
	strncpy( nameoutp, tempb, len);
	tmp = strlen(tempb);
	nameoutp += tmp;
	len -= tmp + 1;
	if (len > 1) {
	    strcpy(nameoutp++, "/");
	}
	
	strptr = NULL;
    }

    if (NULL != (slash = strrchr ( nameout, '/' )))
	*slash = '\0';

    return 1;
}


int expandname ( char *nameout, char *namein )
    
    /* function to look for namein in environment and expand if found
     * can deal with $NAME, NAME, ~, ~fred
     *
     * errors:
     *
     * if ~fred not in password file
     * if $NAME not in environment
     */
    
{  
    nameout[FILENAME_MAX] = '\0';
    if ( namein[0] == '~' ) {
	if ( strlen ( namein ) == 1 ) {
/* 6/1/99 - johnt - modified to get HOME directory under Windows NT */
#ifdef _WIN32
	    int i;
	    strncpy(nameout,getenv("HOMEDRIVE"),FILENAME_MAX);
	    strncat(nameout,getenv("HOMEPATH"),FILENAME_MAX-strlen(nameout));
	    for(i=0;nameout[i]!='\0';i++)
		if( nameout[i] == '\\')
		    nameout[i]='/';
#else
	    strncpy ( nameout, getenv("HOME"), FILENAME_MAX);
#endif
	    return 1;
	}
#ifndef _WIN32
	else {
	    struct passwd *pwentry;
	    if ( ( pwentry = getpwnam ( &namein[1] )) == NULL ) return 0;
	    strncpy ( nameout, pwentry->pw_dir, FILENAME_MAX);
	    return 1;
	}
#endif
    }

    if ( namein[0] == '$' ) {
	if ( getenv ( &namein[1] ) ) {
	    strncpy ( nameout, getenv(&namein[1]), FILENAME_MAX);
	    return 1;
	}
	else
	    return 0;
    }

#if defined(NO_DOLLARS_FOR_ENV)
    if ( getenv ( &namein[0] ) ) {
	strncpy ( nameout, getenv(namein), FILENAME_MAX);
	return 1;
    } else {
	strncpy ( nameout, namein, FILENAME_MAX);
	return 1;
    }
#else
    strncpy ( nameout, namein, FILENAME_MAX);
    return 1;
#endif
}


FILE *my_fopen ( char *filename, char *mode)
    
    /* handles all opens for c files */
    
{
    return fopen(filename, mode);
}


int my_mkdir (char *filename, int update)
    
    /* handles directory creation for all c files */
    
{
    mode_t mode = 0777;
    return mkdir (filename, mode );
}

/* display error message */
/*
void errout( char *message)
{
    
    fprintf (stderr, message);
}
*/
