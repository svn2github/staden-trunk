#include <staden_config.h>

#include <tcl.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#include <tg_gio.h>
#include "gap_cli_arg.h"
#include "import_gff.h"
#include "consensus.h"

/* Maximum GFF line length */
#define MAX_GFF_LINE 8192

/* Maximum number of attributes per entry */
#define MAX_GFF_ATTRIB 100

typedef struct {
    char *key;
    char *val;
} key_val;

typedef struct {
    char   *seqid;
    char   *source;
    char   *type;
    int     start;
    int     end;
    double  score;
    int     strand; /* '+', '-', '.' or '?' */
    int     phase;
    int     n_attrib;
    key_val attrib[MAX_GFF_ATTRIB];
} gff_entry;


/* ------------------------------------------------------------------------
 * Internal functions
 */

/*
 * This decodes %xx hex components of str and writes back to str.
 * (The string modified is guaranteed to never be longer than the hex encoded
 * form.)
 */
static void decode_hex_insitu(char *str) {
    static int hex[256];
    static int hex_init = 0;
    char *out = str;

    if (!str)
	return;
    

    /* Initialise lookup tables */
    if (!hex_init) {
	int i;
	memset(hex, 0, 256*sizeof(*hex));
	for (i = 0; i <= 9; i++) {
	    hex['0'+i] = i;
	}
	for (i = 0; i <= 5; i++) {
	    hex['a'+i] = 10+i;
	    hex['A'+i] = 10+i;
	}

	hex_init = 1;
    }


    /* Decode */
    while (*str) {
	if (*str == '%') {
	    if (!str[1]) {
		fprintf(stderr,"Truncated %% code in unescape_hex_string()\n");
		*out++ = 0;
		return;
	    }
	    *out++ = (hex[str[1]]<<4) | hex[str[2]];
	    str += 3;
	} else {
	    *out++ = *str++;
	}
    }
    *out++ = 0;

    return;
}

/*
 * Parses a (modifiable) line of text and fills out a gff_entry struct
 * pointing to the contents of this line. For efficiency reasons the
 * gff_struct is passed in, but may be supplied as NULL in which case
 * a fresh one will be malloced.
 *
 * The line supplied is modified by this code, both for replacing entry
 * delimiters with nuls in order to get C strings and to replace %hex
 * encodings with their literal characters. (Therefore %00 is not
 * supported).
 *
 * Returns the gff_struct* on success (usually the pointer passed in)
 *         NULL on failure
 */
static gff_entry *parse_gff_entry(char *line, gff_entry *gff) {
    char *cp, *tmp;
    enum state {SEQID, SOURCE, TYPE, START, END, SCORE, STRAND,
                PHASE, KEY, VAL};

    if (!line)
	return NULL;

    if (!gff)
	gff = malloc(sizeof(*gff));

    cp = line;

    /* SeqID */
    gff->seqid = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    decode_hex_insitu(gff->seqid);

    /* Source */
    gff->source = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    decode_hex_insitu(gff->seqid);

    /* Type */
    gff->type = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    decode_hex_insitu(gff->seqid);

    /* Start */
    tmp = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    gff->start = atoi(tmp);

    /* End */
    tmp = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    gff->end = atoi(tmp);

    /* Score */
    tmp = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    gff->score = atof(tmp);

    /* Strand */    
    tmp = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    if (*tmp == '+')
	gff->strand = '+';
    else if (*tmp == '-')
	gff->strand = '-';
    else if (*tmp == '.')
	gff->strand = '.';
    else
	gff->strand = '?';

    /* Phase */
    tmp = cp;
    while (*cp && *cp != '\t') cp++;
    if (*cp != '\t')
	return NULL;
    *cp++ = 0;
    gff->phase = atoi(tmp);
    
    /* Attribs */
    gff->n_attrib = 0;
    do {
	tmp = cp;
	while (*cp && *cp != '\t' && *cp != '=' && *cp != ';' && *cp != '\n')
	    cp++;
	gff->attrib[gff->n_attrib].key = tmp;
	gff->attrib[gff->n_attrib].val = "";
	if (!*cp || *cp == '\n') {
	    gff->n_attrib++;
	    *cp = 0;
	    break;
	}
	*cp++ = 0;
	decode_hex_insitu(tmp);

	tmp=cp;
	while (*cp && *cp != '\t' && *cp != ';' && *cp != '\n') cp++;
	gff->attrib[gff->n_attrib].val = tmp;
	gff->n_attrib++;
	if (!*cp || *cp == '\n') {
	    *cp = 0;
	    decode_hex_insitu(tmp);
	    break;
	}
	*cp++ = 0;
	decode_hex_insitu(tmp);
    } while (*cp != '\n');
    
    return gff;
}

/*
 * Looks for an attribute matching a specific key and returns the value.
 *
 * Return value if found
 *        NULL if not found or upon error.
 */
static char *gff_find_attrib(gff_entry *gff, char *key) {
    int i;
    for (i = 0; i < gff->n_attrib; i++) {
	if (strcmp(gff->attrib[i].key, key) == 0)
	    return gff->attrib[i].val;
    }

    return NULL;
}


static int   *cached_map   = NULL;
static tg_rec cached_crec  = -1;
static int    cached_start = 0;    /* Bounds in cached_map */
static int    cached_end   = 0;

/*
 * Adds a gff entry as a tag.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int gff_add_tag(GapIO *io, gff_entry *gff, int padded,
		       int plus_as_space) {
    tg_rec rec;
    int rec_type;
    char *type, *txt;
    range_t r;
    bin_index_t *bin;
    anno_ele_t *e;
    contig_t *c;
    char type_a[5];
    tg_rec seq_bin;

    r.flags = GRANGE_FLAG_ISANNO;
    r.start = gff->start;
    r.end   = gff->end;

    /* Get tag TYPE */
    if (!(type = gff_find_attrib(gff, "type"))) {
	char *col = gff_find_attrib(gff, "colour");
	if (!col) col = gff_find_attrib(gff, "color");

	strcpy(type_a, "GF00");
	type = type_a;
	if (col) {
	    int c = atoi(col);
	    type[2] = (c / 10) + '0';
	    type[3] = (c % 10) + '0';
	}
    }
    r.mqual = str2type(type);

    /* Get annotation */
    txt = gff_find_attrib(gff, "Note");
    if (!txt)
	txt = gff_find_attrib(gff, "note"); /* but be pragmatic */

    if (txt && plus_as_space) {
	char *cp;
	for (cp = txt; *cp; cp++) {
	    if (*cp == '+')
		*cp = ' ';
	}
    }

    /* Find seqid rec */
    if ((rec = contig_index_query(io, gff->seqid)) >= 0) {
	c = cache_search(io, GT_Contig, rec);
	rec_type = GT_Contig;

	if (!padded) {
	    /*
	     * Compute mapping table for this contig - not so efficient, so
	     * we cache it. Ideally this should be part of the bin
	     * structure itself and permanently stored.
	     */
	    if (cached_crec != c->rec) {
		int i, np;
		char *con = xmalloc(c->end - c->start + 2);

		if (cached_map)
		    xfree(cached_map);
		cached_map = xmalloc((c->end - c->start + 2) * sizeof(int));
		if (!con || !cached_map)
		    return -1;

		calculate_consensus_simple(io, c->rec, c->start, c->end,
					   con, NULL);

		cached_crec  = c->rec;
		cached_start = c->start;
		for (np = 0, i = c->start; i <= c->end; i++) {
		    cached_map[i - c->start - np] = i;
		    if (con[i - c->start] == '*')
			np++;
		}
		cached_end = c->end - np;
		free(con);
	    }

	    /* Update r.start / r.end via unpadded to padded mapping table */
	    if (r.start >= cached_start && r.start <= cached_end)
		r.start = cached_map[r.start - c->start];

	    if (r.end >= cached_start && r.end <= cached_end)
		r.end = cached_map[r.end - c->start];
	    else
		r.end = cached_map[cached_end - c->start] + r.end - cached_end;

	    if (r.end < r.start)
		r.end = r.start;
	}

	seq_bin = 0;

    } else if ((rec = sequence_index_query(io, gff->seqid)) >= 0) {
	int s_start, s_end, s_orient;
	tg_rec s_contig;

	rec_type = GT_Seq;
	r.flags |= GRANGE_FLAG_TAG_SEQ;

	/* Unpadded coord in seq, just count pads in it */
	if (!padded) {
	    seq_t *sorig = cache_search(io, GT_Seq, rec), *s = sorig;
	    int len = s->len < 0 ? -s->len : s->len;
	    int i, npads;

	    /* Complement seq if needed */
	    if (sequence_get_orient(io, rec)) {
		s = dup_seq(s);
		complement_seq_t(s);
	    }

	    /* Compensate for pads in sequence */
	    for (npads = i = 0; i < len && i <= r.start + npads; i++) {
		if (s->seq[i] == '*')
		    npads++;
	    }
	    if (i > r.start + npads) i--;
	    r.start = i;

	    for (; i < len && i <= r.end + npads; i++) {
		if (s->seq[i] == '*')
		    npads++;
	    }
	    if (i > r.end + npads) i--;
	    r.end = i;

	    if (s != sorig)
		free(s);
	}

	/* Adjust start coord to absolute contig positions */
	sequence_get_position2(io, rec, &s_contig, &s_start, &s_end, &s_orient,
			       &seq_bin, NULL, NULL);
	r.start += s_start;
	r.end   += s_start;
	c = cache_search(io, GT_Contig, s_contig);

    } else {
	fprintf(stderr, "Unknown seqid '%s'\n", gff->seqid);
	return -1;
    }

    r.pair_rec = rec;
    r.rec = anno_ele_new(io, 0, rec_type, rec, 0, str2type(type),
			 gff->strand, txt);

    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, r.rec);
    e = cache_rw(io, e);
	
    if (seq_bin)
	bin = bin_add_to_range(io, &c, seq_bin, &r, NULL, NULL, 0);
    else
	bin = bin_add_range(io, &c, &r, NULL, NULL, 0);
    e->bin = bin->rec;

    return 0;
}

/* ------------------------------------------------------------------------
 * External C API
 */

/*
 * Imports a GFF file and produces gap5 annotations.
 * If 'padded' is true then the GFF coordinates represented padded
 * positions in Gap5. Otherwise they are unpadded and we need to
 * correct for this.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int import_gff(GapIO *io, char *fn, int padded, int plus_as_space) {
    FILE *fp;
    char line[MAX_GFF_LINE];
    gff_entry gff;
    int nentry = 0;

    if (NULL == (fp = fopen(fn, "r"))) {
	perror(fn);
	return -1;
    }

    while (fgets(line, MAX_GFF_LINE, fp)) {
	if (line[0] == '#')
	    continue;

	if (!parse_gff_entry(line, &gff)) {
	    verror(ERR_WARN, "parse_gff", "Malformed gff_line");
	    continue;
	}
	
#if 0
    {
	int i;
	printf("Gff:%s:%s:%s:%d:%d:%f:%d:%d:",
	       gff.seqid,  gff.source, gff.type,
	       gff.start,  gff.end,    gff.score,
	       gff.strand, gff.phase);

	for (i = 0; i < gff.n_attrib; i++) {
	    printf(" :%s:=>:%s:", gff.attrib[i].key, gff.attrib[i].val);
	}
	puts("");
    }
#endif

        gff_add_tag(io, &gff, padded, plus_as_space);
	nentry++;
    }

    fclose(fp);

    gio_debug(io, 1, "Processed %d GFF entries\n", nentry);

    /*
     * Remove old cached map if needed incase data is edited before
     * we next call this function.
     */
    if (cached_map) {
	free(cached_map);
	cached_map = NULL;
    }
    cached_crec = -1;

    return 0;
}

/* ------------------------------------------------------------------------
 * Tcl interfaces
 */
typedef struct {
    GapIO *io;
    char  *infile;
    int    padded;
    int    plus_space;
} ig_arg;

int tcl_import_gff(ClientData clientData, Tcl_Interp *interp,
		   int objc, Tcl_Obj *CONST objv[])
{
    ig_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,     offsetof(ig_arg, io)},
	{"-infile",	ARG_STR, 1, NULL,     offsetof(ig_arg, infile)},
	{"-padded",	ARG_INT, 1, "0",      offsetof(ig_arg, padded)},
	{"-plus_space", ARG_INT, 1, "1",      offsetof(ig_arg, plus_space)},
	{NULL,	        0,	 0, NULL,     0}
    };
    int res;

    vfuncheader("Import GFF");
    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    res = import_gff(args.io, args.infile, args.padded, args.plus_space);
    cache_flush(args.io);

    return res == 0 ? TCL_OK : -1;
}

