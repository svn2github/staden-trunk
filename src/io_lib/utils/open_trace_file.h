#ifndef _OPEN_TRACE_FILE_H_
#define _OPEN_TRACE_FILE_H_

/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in RAWDATA (which is a colon separated list).
 *
 * Returns a FILE pointer when found.
 *           NULL otherwise.
 */
FILE *open_trace_file(char *file, char *relative_to);


/*
 * Searches for file in the tar pointed to by tarname. If it finds it, it
 * copies it out and returns a file pointer to the temporary file,
 * otherwise we return NULL.
 *
 * If 'tarname'.index exists we will use this as a fast lookup method,
 * otherwise we just do a sequential search through the tar.
 *
 * Offset specifies a starting search position. Set this to zero if you want
 * to search through the entire tar file, otherwise set it to the byte offset
 * into the file of the tar header block for the desired file to extract.
 * (Note that the tar index file overrides this value.)
 *
 * Returns FILE pointer if found
 *         NULL if not.
 */
FILE *find_file_tar(char *file, char *tarname, size_t offset);

#endif /* _OPEN_TRACE_FILE_H_ */
