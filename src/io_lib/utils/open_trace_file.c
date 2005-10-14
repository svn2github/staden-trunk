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
#include "os.h"
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
#include "sff.h"

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

    for (i = 0, j = 0; i < len; i++) {
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

    if (j)
	newsearch[j++] = 0;
    newsearch[j++] = '.';
    newsearch[j++] = '/';
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
 * Returns mFILE pointer if found
 *         NULL if not.
 */
static mFILE *find_file_tar(char *file, char *tarname, size_t offset) {
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

	size = strtol(blk.header.size, NULL, 8);

	/* start with the same name... */
	if (strncmp(blk.header.name, file, name_len) == 0) {
	    char *data;
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

	    /* Found it - copy out the data to an mFILE */
	    if (NULL == (data = (char *)malloc(size)))
		return NULL;
	    if (size != fread(data, 1, size, fp)) {
		free(data);
		return NULL;
	    }
	    return mfcreate(data, size);
	}

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
 * Return mFILE pointer if found
 *        NULL if not
 */
static mFILE *find_file_hash(char *file, char *hashfile) {
    size_t size;
    static HashFile *hf = NULL;
    static char hf_name[1024];
    char *data;

    /* Cache an open HashFile for fast accesing */
    if (strcmp(hashfile, hf_name) != 0) {
	if (hf)
	    HashFileDestroy(hf);
	hf = HashFileOpen(hashfile);

	if (!hf)
	    return NULL;
	strcpy(hf_name, hashfile);
    }

    /* Search */
    if (NULL == (data = HashFileExtract(hf, file, &size)))
	return NULL;

    /* Found, so copy the contents to a fake FILE pointer */
    return mfcreate(data, size);
}

#ifdef TRACE_ARCHIVE
/*
 * Searches for file in the ensembl trace archive pointed to by arcname.
 * If it finds it, it copies it out and returns a file pointer to the
 * temporary file, otherwise we return NULL.
 *
 * Arcname has the form address:port, eg "titan/22100"
 *
 * Returns mFILE pointer if found
 *         NULL if not.
 */
#define RDBUFSZ 8192
static mFILE *find_file_archive(char *file, char *arcname) {
    char server[1024], *cp;
    int port;
    struct hostent *host;
    struct sockaddr_in saddr;
    int s = 0;
    char msg[1024];
    ssize_t msg_len;
    char buf[RDBUFSZ];
    mFILE *fpout;
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
     * Create a fake FILE (mFILE) and write to it.
     */
    fpout = mfcreate(NULL, 0);

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
	    mfwrite(buf, 1, msg_len, fpout);
    }
    close(s);

    if (!block_count) {
	mfclose(fpout);
	return NULL;
    }

    mrewind(fpout);

    return fpout;
}
#endif

#ifdef USE_WGET
static mFILE *find_file_url(char *file, char *url) {
    char buf[8192], *cp;
    mFILE *fp;
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
    fp = !status ? mfopen(fname, "rb+") : NULL;
    remove(fname);
    free(fname);

    return fp;
}
#endif

/*
 * Takes an SFF file in 'data' and edits the header to ensure
 * that it has no index listed and only claims to contain a single entry.
 * This isn't strictly necessary for the sff/sff.c reading code, but it is
 * the 'Right Thing' to do.
 *
 * Returns an mFILE on success or NULL on failure.
 */
static mFILE *sff_single(char *data, size_t size) {
    *(uint64_t *)(data+8)  = be_int8(0); /* index offset */
    *(uint32_t *)(data+16) = be_int4(0); /* index size */
    *(uint32_t *)(data+20) = be_int4(1); /* number of reads */

    return mfcreate(data, size);
}

/*
 * This returns an mFILE containing an SFF entry.
 *
 * This does the minimal decoding necessary to skip through the SFF
 * container to find an entry. In this respect it is a semi-duplication
 * of sff/sff.[ch], but implemented for efficiency.
 *
 * Having found an entry it packs the common header, the read specific
 * header and the read data into a single block of memory and returns this
 * as an mFILE. In essence it produces a single-read SFF archive. This
 * is then decoded by the normal sff parsing code representing a small
 * amount of redundancy, but one which is swamped by the I/O time.
 */
static mFILE *find_file_sff(char *entry, char *sff) {
    FILE *fp;
    char chdr[65536], rhdr[65536]; /* generous, but worst case */
    uint32_t nkey, nflows, chdrlen, rhdrlen, dlen, index_length, magic;
    uint64_t index_offset, file_pos;
    uint32_t nreads, i;
    size_t entry_len = strlen(entry);
    int bytes_per_flow = 2;
    char *fake_file;

    /* Read the common header */
    if (NULL == (fp = fopen(sff, "rb")))
	return NULL;
    if (31 != fread(chdr, 1, 31, fp))
	return NULL;

    /* Check magic & vers: TODO */
    magic = be_int4(*(uint32_t *)chdr);
    if (magic != SFF_MAGIC)
	return NULL;
    if (memcmp(chdr+4, SFF_VERSION, 4) != 0)
	return NULL;

    /* If we have an index, use it, otherwise search linearly */
    index_offset = be_int8(*(uint64_t *)(chdr+8));
    index_length = be_int4(*(uint32_t *)(chdr+16));
    if (index_length != 0) {
	char index_format[4];
	long orig_pos = ftell(fp);
	fseek(fp, index_offset, SEEK_SET);
	fread(index_format, 1, 4, fp);

	if (memcmp(index_format, ".hsh", 4) == 0) {
	    /* HASH index */
	    HashFile *hf;
	    char *data;
	    size_t size;

	    fseek(fp, -4, SEEK_CUR);
	    hf = HashFileFopen(fp);
	    data = HashFileExtract(hf, entry, &size);
	    HashFileDestroy(hf); /* closes fp */

	    return data ? sff_single(data, size) : NULL;

	} else {
	    /* Unknown index: revert back to a slow linear scan */
	    fseek(fp, orig_pos, SEEK_SET);
	}
    }

    nreads  = be_int4(*(uint32_t *)(chdr+20));
    chdrlen = be_int2(*(uint16_t *)(chdr+24));
    nkey    = be_int2(*(uint16_t *)(chdr+26));
    nflows  = be_int2(*(uint16_t *)(chdr+28));

    /* Read the remainder of the header */
    if (chdrlen-31 != fread(chdr+31, 1, chdrlen-31, fp))
	return NULL;

    file_pos = chdrlen;

    /* Loop until we find the correct entry */
    for (i = 0; i < nreads; i++) {
	uint16_t name_len;
	uint32_t nbases;

	/* Index could be between common header and first read - skip */
	if (file_pos == index_offset) {
	    fseek(fp, index_length, SEEK_CUR);
	    file_pos += index_length;
	}

	/* Read 16 bytes to get name length */
	if (16 != fread(rhdr, 1, 16, fp))
	    return NULL;
	rhdrlen   = be_int2(*(uint16_t *)rhdr);
	name_len = be_int2(*(uint16_t *)(rhdr+2));
	nbases   = be_int4(*(uint32_t *)(rhdr+4));

	/* Read the rest of the header */
	if (rhdrlen-16 != fread(rhdr+16, 1, rhdrlen-16, fp))
	    return NULL;

	file_pos += rhdrlen;

	dlen = (nflows * bytes_per_flow + nbases * 3 + 7) & ~7;

	if (name_len == entry_len  && 0 == memcmp(rhdr+16, entry, entry_len))
	    break;

	/* This is not the read you are looking for... */
	fseek(fp, dlen, SEEK_CUR);
    }

    if (i == nreads) {
	/* Not found */
	return NULL;
    }

    /*
     * Although we've decoded some bits already, we take the more modular
     * approach of packing the sections together and passing the entire
     * data structure off as a single-read SFF file to be decoded fully
     * by the sff reading code.
     */
    if (NULL == (fake_file = (char *)xmalloc(chdrlen + rhdrlen + dlen)))
	return NULL;

    memcpy(fake_file, chdr, chdrlen);
    memcpy(fake_file+chdrlen, rhdr, rhdrlen);
    if (dlen != fread(fake_file+chdrlen+rhdrlen, 1, dlen, fp)) {
	xfree(fake_file);
	return NULL;
    }

    /* Convert to an mFILE and return */
    fclose(fp);
    return sff_single(fake_file, chdrlen+rhdrlen+dlen);
}

/*
 * Searches for file in the directory 'dirname'. If it finds it, it opens
 * it. This also searches for compressed versions of the file in dirname
 * too.
 *
 * Returns mFILE pointer if found
 *         NULL if not
 */
static mFILE *find_file_dir(char *file, char *dirname) {
    char path[PATH_MAX+1], path2[PATH_MAX+1];
    size_t len = strlen(dirname);
    int num_magics = sizeof(magics) / sizeof(*magics);
    int i;
    char *cp;

    if (dirname[len-1] == '/')
	len--;

    /* Special case for "./" */
    if (len==1 && *dirname == '.')
	sprintf(path, "%s", file);
    else 
	sprintf(path, "%.*s/%s", (int)len, dirname, file);

    /*
     * Given a pathname /a/b/c if a/b is a file and not a directory then
     * we'd get an ENOTDIR error. Instead we assume that a/b is an archive
     * and we attempt to work out what type by reading the first and last
     * bits of the file.
     */
    if (cp = strrchr(file, '/')) {
	strcpy(path2, path); /* path contains / too as it's from file */
	*strrchr(path2, '/') = 0;

	if (is_file(path2)) {
	    /* Open the archive to test for magic numbers */
	    char magic[6];
	    FILE *fp;
	    enum archive_type_t {
		NONE, HASH, TAR, SFF
	    } type = NONE;

	    if (NULL == (fp = fopen(path2, "rb")))
		return NULL;
	    memcpy(magic, "\0\0\0\0\0\0", 4);
	    fread(magic, 1, 4, fp);

	    /* .hsh or .sff at start */
	    if (memcmp(magic, ".hsh", 4) == 0)
		type = HASH;
	    else if (memcmp(magic, ".sff", 4) == 0)
		type = SFF;

	    /* Or .hsh at the end */
	    if (NONE == type) {
		fseek(fp, -12, SEEK_END);
		fread(magic, 1, 4, fp);
		if (memcmp(magic, ".hsh", 4) == 0)
		    type = HASH;
	    }

	    /* or ustar 257 bytes in to indicate un-hashed tar */
	    if (NONE == type) {
		fseek(fp, 257, SEEK_SET);
		fread(magic, 1, 6, fp);
		if (memcmp(magic, "ustar\0", 6) == 0)
		    type = TAR;
	    }
	    fclose(fp);

	    switch (type) {
	    case HASH:
		return find_file_hash(cp+1, path2);
	    case TAR:
		return find_file_tar(cp+1, path2, 0);
	    case SFF:
		return find_file_sff(cp+1, path2);
	    case NONE:
		break;
	    }

	    return NULL;
	}
    }

    /* Is it lurking under another, compressed, name? */
    for (i = 0; i < num_magics; i++) {
	sprintf(path2, "%s%s", path, magics[i]);
	if (file_exists(path2)) {
	    return fopen_compressed(path2, NULL);
	    /* return mfopen(path2, "rb"); */
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
 * Returns a mFILE pointer when found.
 *           NULL otherwise.
 */
mFILE *open_trace_mfile(char *file, char *relative_to) {
    char *newsearch;
    char *ele;
    mFILE *fp;

    /* Use RAWDATA first */
    if (NULL == (newsearch = tokenise_search_path(getenv("RAWDATA"))))
	return NULL;
    
    /*
     * Step through the search path testing out each component.
     * We now look through each path element treating some prefixes as
     * special, otherwise we treat the element as a directory.
     */
    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
	int i;
	char *suffix[6] = {"", ".gz", ".bz2", ".sz", ".Z", ".bz2"};
	for (i = 0; i < 6; i++) {
	    char file2[1024];
	    sprintf(file2, "%s%s", file, suffix[i]);

	    if (0 == strncmp(ele, "TAR=", 4)) {
		if (fp = find_file_tar(file2, ele+4, 0)) {
		    free(newsearch);
		    return fp;
		}

	    } else if (0 == strncmp(ele, "HASH=", 5)) {
		if (fp = find_file_hash(file2, ele+5)) {
		    free(newsearch);
		    return fp;
		}
#ifdef TRACE_ARCHIVE
	    } else if (0 == strncmp(ele, "ARC=", 4)) {
		if (fp = find_file_archive(file2, ele+4)) {
		    free(newsearch);
		    return fp;
		}
#endif
#ifdef USE_WGET
	    } else if (0 == strncmp(ele, "URL=", 4)) {
		if (fp = find_file_url(file2, ele+4)) {
		    free(newsearch);
		    return fp;
		}
#endif
	    } else if (0 == strncmp(ele, "SFF=", 4)) {
		if (fp = find_file_sff(file2, ele+4)) {
		    free(newsearch);
		    return fp;
		}
	    } else {
		if (fp = find_file_dir(file2, ele)) {
		    free(newsearch);
		    return fp;
		}
	    }
	}
    }

    free(newsearch);

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

    return NULL;
}


FILE *open_trace_file(char *file, char *relative_to) {
    mFILE *mf = open_trace_mfile(file, relative_to);
    FILE *fp;
    char *fname;

    if (!mf)
	return NULL;

    if (mf->fp)
	return mf->fp;

    /* Otherwise create a temporary file and write the trace out */
    /* Use tempnam() to force the use of TMP environment variable on Windows */
    if (NULL == (fname=tempnam(NULL, NULL)))
	return 0;
    if (NULL == (fp = fopen(fname, "wb+"))){
	remove(fname);
	free(fname);
	return 0;
    }
    remove(fname);
    free(fname);

    /* Copy the data */
    fwrite(mf->data, 1, mf->size, fp);
    rewind(fp);
    mfclose(mf);

    return fp;
}
