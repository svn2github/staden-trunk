#ifndef _GETFILE_H_
#define _GETFILE_H_

#define MAXTRIES 5

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
extern int my_access(char *filename, char mode);


/* function to expand path names:
 * 
 * if name starts with the following symbols
 * 
 * /            unchanged
 * ~/           "home directory name"
 * ~fred        "home directory path for fred"
 * ~fred/etc    "home directory path for fred/etc"
 * for all other cases look in the environment variables for the strings
 * between the /'s and expand them if found. Environment can include $ but
 * it is not necessary
 * 
 * return 1 for success, 0 for failure
 */
extern int expandpath(char *namein, char *nameout);


/* function to look for namein in environment and expand if found
 * can deal with $NAME, NAME, ~, ~fred
 *
 * errors:
 *
 * if ~fred not in password file
 * if $NAME not in environment
 */
extern int expandname(char *nameout, char *namein);


FILE *my_fopen ( char *filename, char *mode);

#endif
