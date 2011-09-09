/* ----------------------------------------------------------------------
 * General purpose fasta and fastq reading code.
 *
 * We provide next_fasta and next_fastq iterators on an opened FILE*.
 * Deallocate memory by passing in NULL to next_fast[aq] function.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "tg_gio.h"
#include "tg_index_common.h"
#include "zfio.h"
#include "fasta.h"

typedef struct {
    char *name;
    char *seq;
    char *qual;
    int  max_name_len;
    size_t  max_seq_len;
    size_t  max_qual_len;
    size_t  seq_len;
} fastq_entry_t;


/*
 * Reads a new fasta entry and returns it. The data in this struct is valid
 * until the next call to this function. To force memory free pass in fp
 * as NULL (or ignore it as it'll be reused on the next file we load anyway).
 * Note the struct returned is a fastq one, but with a NULL quality string.
 *
 * Returns NULL on failure or EOF (distinguish via feof(fp)).
 */
#define BLK_SIZE 8192
fastq_entry_t *fasta_next(zfp *fp) {
    static fastq_entry_t *e = NULL;
    char line[BLK_SIZE], *cp;
    size_t l;

    /* Memory free */
    if (!fp) {
	if (e) {
	    if (e->name)
		free(e->name);
	    if (e->seq)
		free(e->seq);
	    free(e);
	}
	e = NULL;

	return NULL;
    }

    if (!e) {
	e = calloc(1, sizeof(*e));
	if (!e)
	    return NULL;
    }

    if (e->qual) {
	free(e->qual);
	e->qual = NULL;
	e->max_qual_len = 0;
    }

    /* Read name */
    cp = e->name;
    do {
	if (NULL == zfgets(line, BLK_SIZE, fp))
	    return NULL;
	l = strlen(line);
	if (*line != '>') {
	    verror(ERR_WARN, "fasta_next",
		   "Sequence does not appear to be in fasta format");
	    return NULL;
	}
	while (e->max_name_len < cp-e->name + l) {
	    ptrdiff_t diff = cp - e->name;
	    e->max_name_len = e->max_name_len ? 2*e->max_name_len : 1024;
	    e->name = realloc(e->name, e->max_name_len);
	    if (!e->name)
		return NULL;
	    cp = e->name + diff;
	}
	strcpy(cp, line + (cp == e->name ? 1 : 0));
	cp += l - (cp == e->name ? 1 : 0);
	if (cp[-1] == '\n')
	    cp[-1] = 0;
    } while (line[l-1] != '\n');

    /* Truncate name at first whitespace */
    cp = e->name;
    while (*cp && !isspace(*cp)) {
	cp++;
    }
    *cp = 0;

    /* Read sequence */
    cp = e->seq;
    while (!zfeof(fp) && zfpeek(fp) != '>') {
	int more;
	do {
	    char *cp2;

	    if (NULL == zfgets(line, BLK_SIZE, fp)) {
		/* Assumed last entry */
		e->seq_len = cp-e->seq;
		*cp++ = 0;
		return e;
	    }

	    l = strlen(line);
	    if (line[l-1] == '\n') {
		l--;
		more = 0;
	    } else {
		more = 1;
	    }

	    while (e->max_seq_len < cp-e->seq + l + 1) {
		ptrdiff_t diff = cp - e->seq;
		e->max_seq_len = e->max_seq_len ? 2*e->max_seq_len : 1024;
		e->seq = realloc(e->seq, e->max_seq_len);
		if (!e->seq)
		    return NULL;
		cp = e->seq + diff;
	    }

	    cp2 = line;
	    while (l-- > 0) {
		if (!isspace(*cp2))
		    *cp++ = *cp2;
		cp2++;
	    }
	} while (more);
    }

    if (cp) {
	e->seq_len = cp-e->seq;
	*cp++ = 0;
    } else { 
	return NULL;
    }

    return e;
}

/*
 * Reads a new fasta entry and returns it. The data in this struct is valid
 * until the next call to this function. To force memory free pass in fp
 * as NULL (or ignore it as it'll be reused on the next file we load anyway).
 *
 * Returns NULL on failure or EOF (distinguish via feof(fp)).
 */
fastq_entry_t *fastq_next(zfp *fp) {
    static fastq_entry_t *e = NULL;
    char line[BLK_SIZE], *cp;
    size_t l;
    int nqual;

    /* Memory free */
    if (!fp) {
	if (e) {
	    if (e->name)
		free(e->name);
	    if (e->seq)
		free(e->seq);
	    if (e->qual)
		free(e->qual);
	    free(e);
	}
	e = NULL;

	return NULL;
    }

    if (!e) {
	e = calloc(1, sizeof(*e));
	if (!e)
	    return NULL;
    }

    /* Read name */
    cp = e->name;
    do {
	do {
	    if (NULL == zfgets(line, BLK_SIZE, fp))
		return NULL;
	/* blank line detection */
	} while (*line == '\n' && cp == e->name);

	if (*line != '@' && cp == e->name) {
	    fprintf(stderr, "Error: sequence name does not start with @\n"
		    "Previous quality line too long?\n");
	    return NULL;
	}

	l = strlen(line);

	while (e->max_name_len < cp-e->name + l) {
	    ptrdiff_t diff = cp - e->name;
	    e->max_name_len = e->max_name_len ? 2*e->max_name_len : 1024;
	    e->name = realloc(e->name, e->max_name_len);
	    if (!e->name)
		return NULL;
	    cp = e->name + diff;
	}
	strcpy(cp, line + (cp == e->name ? 1 : 0));
	cp += l - (cp == e->name ? 1 : 0);
	if (cp[-1] == '\n')
	    cp[-1] = 0;
    } while (line[l-1] != '\n');

    /* Truncate name at first whitespace */
    cp = e->name;
    while (*cp && !isspace(*cp)) {
	cp++;
    }
    *cp = 0;

    /* Read sequence */
    cp = e->seq;
    while (!zfeof(fp) && zfpeek(fp) != '+') {
	int more;
	do {
	    char *cp2;

	    if (NULL == zfgets(line, BLK_SIZE, fp)) {
		/* Assumed last entry */
		if (*cp) {
		    e->seq_len = cp-e->seq;
		    *cp++ = 0;
		} else {
		    e->seq_len = 0;
		}
		return e;
	    }

	    l = strlen(line);
	    if (line[l-1] == '\n') {
		l--;
		more = 0;
	    } else {
		more = 1;
	    }

	    while (e->max_seq_len < cp-e->seq + l + 1) {
		ptrdiff_t diff = cp - e->seq;
		e->max_seq_len = e->max_seq_len ? 2*e->max_seq_len : 1024;
		e->seq = realloc(e->seq, e->max_seq_len);
		if (!e->seq)
		    return NULL;
		cp = e->seq + diff;
	    }

	    cp2 = line;
	    while (l-- > 0) {
		if (!isspace(*cp2))
		    *cp++ = *cp2;
		cp2++;
	    }
	} while (more);
    }
    if (cp) {
	e->seq_len = cp-e->seq;
	*cp++ = 0;
    } else { 
	e->seq_len = 0;
    }

    /* + line: skip */
    if (NULL == zfgets(line, BLK_SIZE, fp) || *line != '+')
	/* eof */
	return NULL;


    /* Read quality, no more than e->seq_len chars */
    cp = e->qual;
    nqual = e->seq_len;
    while (nqual > 0 && !zfeof(fp)) {
	int more;
	do {
	    char *cp2;

	    if (NULL == zfgets(line, BLK_SIZE, fp)) {
		/* Assumed last entry */
		if (*cp)
		    *cp++ = 0;
		return e;
	    }

	    l = strlen(line);
	    if (line[l-1] == '\n') {
		l--;
		more = 0;
	    } else {
		more = 1;
	    }

	    while (e->max_qual_len < cp-e->qual + l + 1) {
		ptrdiff_t diff = cp - e->qual;
		e->max_qual_len = e->max_qual_len ? 2*e->max_qual_len : 1024;
		e->qual = realloc(e->qual, e->max_qual_len);
		if (!e->qual)
		    return NULL;
		cp = e->qual + diff;
	    }

	    cp2 = line;
	    while (l-- > 0) {
		if (!isspace(*cp2)) {
		    *cp++ = *cp2;
		    nqual--;
		}
		cp2++;
	    }
	} while (more);
    }
    if (cp)
	*cp = 0;

    if (nqual != 0) {
	fprintf(stderr, "Error: differing number of sequence and quality "
		"characters for sequence '%s'\n", e->name);
	return NULL;
    }

    return e;
}


/* ----------------------------------------------------------------------
 * tg_index interface below.
 */
int parse_fasta_or_fastq(GapIO *io, char *fn, tg_args *a, int format) {
    int nseqs = 0, ret = 0;
    zfp *fp;
    struct stat sb;
    fastq_entry_t *ent;
    fastq_entry_t *(*next_seq)(zfp *fp);
    contig_t *c = NULL;

    printf("Loading %s...\n", fn);
    if (-1 == stat(fn, &sb) ||
	NULL == (fp = zfopen(fn, "r"))) {
	perror(fn);
	return -1;
    }

    next_seq = (format == 'a') ? fasta_next : fastq_next;

    /* Fetch sequences */
    while ((ent = next_seq(fp))) {
	seq_t seq;
	static int dummy_qual_len;
	static char *dummy_qual = NULL;

	// printf("@%s\n%s\n+\n%s\n", ent->name, ent->seq, ent->qual);
	// printf("%d\tSeq %s len %d / %d\n",
	// nseqs, ent->name, (int)strlen(ent->seq), ent->seq_len);

	if (ent->seq_len <= 0) {
	    verror(ERR_WARN, "parse_fasta_or_fastq",
		   "Sequence named '%s' appears to be blank", ent->name);
	    continue;
	}

	/* Create 1 read contig */
	create_new_contig(io, &c, ent->name, 0);
	
	seq.rec         = 0;
	seq.parent_rec  = 0;
	seq.parent_type = 0;

	seq.pos      = 1;
	seq.flags    = 0;
	seq.seq_tech = STECH_UNKNOWN;
	seq.format   = SEQ_FORMAT_CNF1;
	seq.left     = 1;
	seq.right    = ent->seq_len;

	seq.name_len = strlen(ent->name);
	seq.name     = strdup(ent->name);

	seq.seq      = ent->seq;
	seq.len      = ent->seq_len;

	if (dummy_qual_len < ent->seq_len) {
	    dummy_qual_len = ent->seq_len;
	    dummy_qual = realloc(dummy_qual, dummy_qual_len);
	    if (!dummy_qual)
		return -1;
	}
	
	seq.conf = dummy_qual;

	if (ent->qual) {
	    int i;
	    for (i = 0; i < ent->seq_len; i++) {
		int q = ent->qual[i] - '!';
		if (q < 0)
		    q = 0;
		if (q > 100)
		    q = 100;
		seq.conf[i] = q;
	    }
	} else {
	    memset(dummy_qual, 0, dummy_qual_len);
	}

	seq.trace_name     = NULL;
	seq.trace_name_len = 0;
	seq.alignment      = NULL;
	seq.alignment_len  = 0;
	seq.sam_aux        = NULL;
	seq.aux_len        = 0;
	seq.anno           = NULL;

	save_range_sequence(io,
			    &seq,
			    0,    /* mapping qual */
			    NULL, /* pair array */
			    0,    /* is_pair */
			    NULL, /* template name */
			    c,    /* contig */
			    a,    /* args */
			    GRANGE_FLAG_TYPE_SINGLE,
			    NULL, /* library */
			    NULL  /* bin return rec */
			    );
			 
	if ((++nseqs & 0xff) == 0) {
	    int perc = 0;
	    off_t pos = zftello(fp);

	    perc = 100.0 * pos / sb.st_size;
	    printf("\r%d%c", perc, (nseqs & 0xfff) ? '%' : '*');
	    fflush(stdout);

	    if ((nseqs & 0xfff) == 0)
		cache_flush(io);
	}
    }

    if (!zfeof(fp)) {
	fprintf(stderr, "Terminated early - malformed data file?\n");
	ret = -1;
    }

    printf("Loaded %d sequences\n", nseqs);

    zfclose(fp);
    cache_flush(io);

    return ret;
}
