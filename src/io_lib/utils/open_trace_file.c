#if !(defined(_MSC_VER) || defined(__MINGW32__))
#define TRACE_ARCHIVE
#define USE_WGET
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <limits.h>
#ifdef TRACE_ARCHIVE
#  include <sys/socket.h>
#  include <netinet/in.h>
#  include <netdb.h>
#  include <sys/time.h>
#  include <errno.h>
#endif
#ifdef USE_WGET
#  include <sys/wait.h>
#endif
#ifndef PATH_MAX
#  define PATH_MAX 1024
#endif

#include "open_trace_file.h"
#include "misc.h"
#include "tar_format.h"
#include "compress.h"
#include "hash_table.h"

/*
 * Supported compression extensions. See the magics array in compress.c for
 * the full structure.
 */
static char *magics[] = {"", ".bz", ".gz", ".Z", ".z", ".bz2", ".sz"};

/*
 * Initially produce a new search path where all "::"s are replaced with
 * a single ":". This is because on windows we need to include colons in
 * the search path, but colon is also our component separator.
 *
 * We explicitly add a "./" to the start of the search path
 *
 * Returns: A new search path with items separated by nul chars. Two nul
 *          chars in a row represent the end of the tokenised path.
 * Returns NULL for a failure.
 *
 * The returned data has been malloced. It is up to the caller to free this
 * memory.
 */
static char *tokenise_search_path(char *searchpath) {
    char *newsearch;
    unsigned int i, j;
    size_t len;

    if (!searchpath)
	searchpath="";

    newsearch = (char *)malloc((len = strlen(searchpath))+5);
    if (!newsearch)
	return NULL;

    newsearch[0] = '.';
    newsearch[1] = '/';
    newsearch[2] = '\0';

    for (i = 0, j = 3; i < len; i++) {
	if (i < len-1 && searchpath[i] == ':' && searchpath[i+1] == ':') {
	    newsearch[j++] = ':';
	    i++;
	    continue;
	}

	if (searchpath[i] == ':') {
	    /* Skip blank path components */
	    if (j && newsearch[j-1] != 0)
		newsearch[j++] = 0;
	} else {
	    newsearch[j++] = searchpath[i];
	}
    }

    newsearch[j++] = 0;
    newsearch[j++] = 0;
    
    return newsearch;
}

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
FILE *find_file_tar(char *file, char *tarname, size_t offset) {
    int num_magics = sizeof(magics) / sizeof(*magics);
    char path[PATH_MAX+101];
    FILE *fp;
    tar_block blk;
    int size;
    int name_len = strlen(file);

    /* Maximum name length for a tar file */
    if (name_len > 100)
	return NULL;

    /* Search the .index file */
    sprintf(path, "%s.index", tarname);
    if (file_exists(path)) {
	FILE *fpind = fopen(path, "r");
	char *cp;
	int tmp_off;
	int found = 0;
	
	if (fpind) {
	    while (fgets(path, PATH_MAX+100, fpind)) {
		if (cp = strchr(path, '\n'))
		    *cp = 0;
		tmp_off = strtol(path, &cp, 10);
		while (isspace(*cp))
		    cp++;
		if (strncmp(cp, file, name_len) == 0) {
		    int i;
		    for (i = 0; i < num_magics; i++) {
			if (strcmp(&cp[name_len], magics[i]) == 0) {
			    offset = tmp_off;
			    found = 1;
			    break;
			}
		    }
		    if (found)
			break;
		}
	    }
	    fclose(fpind);

	    /* Not in index */
	    if (!found)
		return NULL;
	}
    }

    if (NULL == (fp = fopen(tarname, "rb")))
	return NULL;

    /*
     * Search through the tar file (starting from index position) looking
     * for our filename. If there was no index then we start from position 0.
     */
    fseek(fp, offset, SEEK_SET);
    while(fread(&blk, sizeof(blk), 1, fp) == 1) {
	if (!blk.header.name[0])
	    break;

	/* start with the same name... */
	if (strncmp(blk.header.name, file, name_len) == 0) {
	    int len;
	    char data[8192];
	    FILE *fpout;
	    char *fname;
	    int i;

	    /* ... but does it end with a known compression extension? */
	    for (i = 0; i < num_magics; i++) {
		if (strcmp(&blk.header.name[name_len], magics[i]) == 0) {
		    break;
		}
	    }
	    /* ... apparently not? continue then */
	    if (i == num_magics)
		continue;

	    /* Found it - copy out the data to a temporary file */
	    fname = tempnam(NULL, NULL);
	    if (NULL == (fpout = fopen(fname, "wb+"))) {
		remove(fname);
		free(fname);
		fclose(fp);
		return NULL;
	    }
	    remove(fname);
	    free(fname);

	    size = strtol(blk.header.size, NULL, 8);
	    while ((len = fread(data, 1, size > 8192 ? 8192 : size, fp)) > 0) {
		fwrite(data, 1, len, fpout);
		size -= len;
	    } 
	    
	    fclose(fp);
	    fseek(fpout, 0, SEEK_SET);
	    return fpout;
	}

	size = strtol(blk.header.size, NULL, 8);
	fseek(fp, TBLOCK*((size+TBLOCK-1)/TBLOCK), SEEK_CUR);
    }

    fclose(fp);
    return NULL;
}

/*
 * Reads a hash file to look for a filename. The hash file contains the
 * (relative) pathname for the file it is an index for along with the
 * positions and sizes of each file contained within it. The file format
 * of the archive itself is irrelevant provided that the data is not
 * internally compressed in some manner specific to that archive.
 *
 * Return FILE pointer if found
 *        NULL if not
 */
FILE *find_file_hash(char *file, char *hashfile) {
    uint64_t pos;
    uint32_t size;
    FILE *fpout;
    int found;
    char *fname;
    static HashFile *hf = NULL;
    static char hf_name[1024];

    /* Cache an open HashFile for fast accesing */
    if (strcmp(hashfile, hf_name) != 0) {
	if (hf)
	    HashFileClose(hf);
	hf = HashFileOpen(hashfile);

	if (!hf)
	    return NULL;
	strcpy(hf_name, hashfile);
    }

    /* Search */
    found = HashFileQuery(hf, (uint8_t *)file, strlen(file), &pos, &size);
    if (-1 == found)
	return NULL;

    /* Found, so copy the contents out and reopen - yuk */
    fname = tempnam(NULL, NULL);
    if (NULL == (fpout = fopen(fname, "wb+"))) {
	remove(fname);
	free(fname);
	return NULL;
    }
    remove(fname);
    free(fname);

    fseek(hf->afp, pos, SEEK_SET);
    do {
	char buf[8192];
	int sz = size >= 8192 ? 8192 : size;
	int got = fread(buf, 1, sz, hf->afp);
	fwrite(buf, 1, sz, fpout);
	size -= got;
    } while (size);

    fseek(fpout, 0, SEEK_SET);

    return fpout;
}

#ifdef TRACE_ARCHIVE
/*
 * Searches for file in the ensembl trace archive pointed to by arcname.
 * If it finds it, it copies it out and returns a file pointer to the
 * temporary file, otherwise we return NULL.
 *
 * Arcname has the form address:port, eg "titan/22100"
 *
 * Returns FILE pointer if found
 *         NULL if not.
 */
#define RDBUFSZ 8192
FILE *find_file_archive(char *file, char *arcname) {
    char server[1024], *cp;
    int port;
    struct hostent *host;
    struct sockaddr_in saddr;
    int s = 0;
    char msg[1024];
    ssize_t msg_len;
    char buf[RDBUFSZ];
    char *fname;
    FILE *fpout;
    int block_count;

    /* Split arc name into server and port */
    if (!(cp = strchr(arcname, '/')))
	return NULL;
    strncpy(server, arcname, 1023);
    server[MIN(1023,cp-arcname)] = 0;
    port = atoi(cp+1);

    /* Make and connect socket */
    if (NULL == (host = gethostbyname(server))) {
	perror("gethostbyname()");
	return NULL;
    }
    saddr.sin_port = htons(port);
    saddr.sin_family = host->h_addrtype;
    memcpy(&saddr.sin_addr,host->h_addr_list[0], host->h_length);
    if ((s = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) == -1) {
	perror("socket()");
	return NULL;
    }
    if (connect(s, (struct sockaddr *)&saddr, sizeof(saddr)) == -1) {
	perror("connect()");
	return NULL;
    }

    /* The minimal message to send down is "--scf tracename" */
    sprintf(msg, "--scf %.*s\n", 1000, file);
    msg_len = strlen(msg);
    if (send(s, msg, msg_len, 0) != msg_len) {
	/*
	 * partial request sent, but requests are short so if this
	 * happens it's unlikely we'll cure it by sending multiple
	 * fragments.
	 */
	/* close(s); */
	return NULL;
    }

    /*
     * Create a temporary file, open it, and unlink it so that on a crash
     * or close disk space is freed.
     */
    fname = tempnam(NULL, NULL);
    if (NULL == (fpout = fopen(fname, "wb+"))) {
	remove(fname);
	free(fname);
	fclose(fpout);
	close(s);
	return NULL;
    }
    remove(fname);
    free(fname);

    /*
     * Read the data back, in multiple blocks if necessary and write it
     * to our temporary file. We use a blocking read with a low timeout to
     * prevent locking up the application indefinitely.
     */
    {
	struct timeval tv = {0, 10000};
	setsockopt(s, SOL_SOCKET, SO_RCVTIMEO, (char *)&tv, sizeof(tv));
    }
    errno = 0;
    block_count = 200;
    while ((msg_len = read(s, buf, RDBUFSZ)) > 0 ||
	   (errno == EWOULDBLOCK && --block_count)) {
	errno = 0;
	if (msg_len > 0)
	    fwrite(buf, 1, msg_len, fpout);
    }
    close(s);

    if (!block_count) {
	fclose(fpout);
	return NULL;
    }

    rewind(fpout);

    return fpout;
}
#endif

#ifdef USE_WGET
FILE *find_file_url(char *file, char *url) {
    char buf[8192], *cp;
    FILE *fp;
    int pid;
    int maxlen = 8190 - strlen(file);
    char *fname = tempnam(NULL, NULL);
    int status;

    /* Expand %s for the trace name */
    for (cp = buf; *url && cp - buf < maxlen; url++) {
	if (*url == '%' && *(url+1) == 's') {
	    url++;
	    cp += strlen(strcpy(cp, file));
	} else {
	    *cp++ = *url;
	}
    }
    *cp++ = 0;

    /* Execute wget */
    if ((pid = fork())) {
	waitpid(pid, &status, 0);
    } else {
	execlp("wget", "wget", "-q", "-O", fname, buf, NULL);
    }

    /* Return a filepointer to the result (if it exists) */
    fp = !status ? fopen(fname, "rb+") : NULL;
    remove(fname);
    free(fname);

    return fp;
}
#endif

/*
 * Searches for file in the directory 'dirname'. If it finds it, it opens
 * it. This also searches for compressed versions of the file in dirname
 * too.
 *
 * Returns FILE pointer if found
 *         NULL if not
 */
static FILE *find_file_dir(char *file, char *dirname) {
    char path[PATH_MAX+1], path2[PATH_MAX+1];
    size_t len = strlen(dirname);
    int num_magics = sizeof(magics) / sizeof(*magics);
    int i;

    if (dirname[len-1] == '/')
	len--;

    /* Special case for "./" */
    if (len==1 && *dirname == '.')
	sprintf(path, "%s", file);
    else 
	sprintf(path, "%.*s/%s", (int)len, dirname, file);

    /* Is it lurking under another, compressed, name? */
    for (i = 0; i < num_magics; i++) {
	sprintf(path2, "%s%s", path, magics[i]);
	if (file_exists(path2)) {
	    return fopen_compressed(path2, NULL);
	}
    }

    return NULL;
}

/*
 * ------------------------------------------------------------------------
 * Public functions below.
 */

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
FILE *open_trace_file(char *file, char *relative_to) {
    char *newsearch;
    char *ele;
    FILE *fp;

    /* Look in the same location as the incoming 'relative_to' filename */
    if (relative_to) {
	char *cp;
	char relative_path[PATH_MAX+1];
	strcpy(relative_path, relative_to);
	if (cp = strrchr(relative_path, '/'))
	    *cp = 0;
	if (fp = find_file_dir(file, relative_path))
	    return fp;
    }

    /* Not found it yet? use RAWDATA then */
    if (NULL == (newsearch = tokenise_search_path(getenv("RAWDATA"))))
	return NULL;
    
    /*
     * Step through the search path testing out each component.
     * We now look through each path element treating some prefixes as
     * special, otherwise we treat the element as a directory.
     */
    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
	if (0 == strncmp(ele, "TAR=", 4)) {
	    if (fp = find_file_tar(file, ele+4, 0)) {
		free(newsearch);
		return fp;
	    }

	} else if (0 == strncmp(ele, "HASH=", 5)) {
	    if (fp = find_file_hash(file, ele+5)) {
		free(newsearch);
		return fp;
	    }
#ifdef TRACE_ARCHIVE
	} else if (0 == strncmp(ele, "ARC=", 4)) {
	    if (fp = find_file_archive(file, ele+4)) {
		free(newsearch);
		return fp;
	    }
#endif
#ifdef USE_WGET
	} else if (0 == strncmp(ele, "URL=", 4)) {
	    if (fp = find_file_url(file, ele+4)) {
		free(newsearch);
		return fp;
	    }
#endif
	} else {
	    if (fp = find_file_dir(file, ele)) {
		free(newsearch);
		return fp;
	    }
	}
    }

    free(newsearch);

    return NULL;
}
