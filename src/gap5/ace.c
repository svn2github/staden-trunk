#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

#include "tg_gio.h"
#include "ace.h"
#include "dna_utils.h"
#include "tg_index_common.h"
#include "zfio.h"

/*
 * Code for reading the new-ACE format as described at:
 * http://bcr.musc.edu/manuals/CONSED.txt
 */

#define ACE_AS 1
#define ACE_CO 2
#define ACE_BQ 3
#define ACE_AF 4
#define ACE_BS 5 
#define ACE_RD 6 
#define ACE_QA 7 
#define ACE_DS 8 
#define ACE_RT 9 
#define ACE_CT 10
#define ACE_WA 11
#define ACE_WR 12

#define MAX_NAME 128
#define MAX_LINE_LEN 1024

/* Unknown mapping quality */
#define ACE_MQUAL 255

/* AS <number of contigs> <total number of reads in ace file> */
typedef struct {
    int type;
    int ncontigs;
    int nreads;
} ace_as_t;

/* CO <contig name> <# of bases> <# of reads in contig>
   <# of base segments in contig> <U or C> */
typedef struct {
    int type;
    char cname[MAX_NAME];
    int nbases;
    int nreads;
    int nseg;
    int dir; /* 'U'=>0 or 'C'=>1 */
    char *seq;
} ace_co_t;

/* BQ - basequalities */
typedef struct {
    int type;
    int nqual;
    unsigned char *qual;
} ace_bq_t;

/* AF <read name> <C or U> <padded start consensus position> */
typedef struct {
    int type;
    char rname[MAX_NAME];
    int dir; /* 'U'=>0 or 'C'=>1 */
    int start;
} ace_af_t;

/* BS <padded start consensus position> <padded end consensus position>
   <read name> */
typedef struct {
    int type;
    int start;
    int end;
    char rname[MAX_NAME];
} ace_bs_t;

/* RD <read name> <# of padded bases> <# of whole read info items>
   <# of read tags> */
typedef struct {
    int type;
    char rname[MAX_NAME];
    int nbases;
    int ninfo;
    int ntags;
    char *seq;
} ace_rd_t;

/* QA <qual clipping start> <qual clipping end> <align clipping start>
   <align clipping end> */
typedef struct {
    int type;
    int qstart;
    int qend;
    int astart;
    int aend;
} ace_qa_t;

/* DS CHROMAT_FILE: <name of chromat file> PHD_FILE: <name of phd file>
 * TIME: <date/time of the phd file> CHEM: <prim, term, unknown, etc>
 * DYE: <usually ET, big, etc> TEMPLATE: <template name>
 * DIRECTION: <fwd or rev>
 */
typedef struct {
    int type;
    char chromat[MAX_NAME];
    char phd[MAX_NAME];
    time_t time;
    int chem;
    int dye;
    char tname[MAX_NAME];
    int dir;
} ace_ds_t;

/* RT{
 * data
 * }
 * Similarly WR{}, CT{} and WA{}.
 */
typedef struct {
    int type;
    char *text;
} ace_rt_t;

typedef struct {
    int type;
    char *text;
} ace_wr_t;

typedef struct {
    int type;
    char *text;
} ace_ct_t;

typedef struct {
    int type;
    char *text;
} ace_wa_t;

typedef union {
    int type;
    ace_as_t as;
    ace_co_t co;
    ace_bq_t bq;
    ace_af_t af;
    ace_rd_t rd;
    ace_qa_t qa;
    ace_ds_t ds;
    ace_rt_t rt;
    ace_wr_t wr;
    ace_ct_t ct;
    ace_wa_t wa;
} ace_item_t;

/*
 * Loads an ACE format item and fills out the ace_item_t struct.
 * The caller should check ferror and feof to distinguish the exit
 * conditions.
 *
 * Returns a pointer to a static ace_item_t struct. Do not free it.
 *         NULL on failure or EOF.
 */
ace_item_t *next_ace_item(zfp *fp) {
    static ace_item_t ai = {0};
    char line[MAX_LINE_LEN];
    static int last_depadded_len = 0, i;

    /* Free previous item */
    switch (ai.type) {
    case ACE_CO:
	if (ai.co.seq)
	    free(ai.co.seq);
	break;

    case ACE_BQ:
	if (ai.bq.qual)
	    free(ai.bq.qual);
	break;

    case ACE_RD:
	if (ai.rd.seq)
	    free(ai.rd.seq);
	break;

    case ACE_RT:
	if (ai.rt.text)
	    free(ai.rt.text);
	break;

    case ACE_WR:
	if (ai.wr.text)
	    free(ai.wr.text);
	break;

    case ACE_CT:
	if (ai.ct.text)
	    free(ai.ct.text);
	break;

    case ACE_WA:
	if (ai.wa.text)
	    free(ai.wa.text);
	break;
    }

    /* Read in next item */
    do {
	if (NULL == zfgets(line, MAX_LINE_LEN, fp)) {
	    return NULL;
	}
    } while (line[0] == '\n');

    /* Decode it */
    if (strncmp(line, "AS", 2) == 0) {
	ai.type = ACE_AS;
	if (2 != sscanf(line+3, "%d %d", &ai.as.ncontigs, &ai.as.nreads))
	    return NULL;

    } else if (strncmp(line, "CO", 2) == 0) {
	char dir;
	int pos = 0;

	ai.type = ACE_CO;
	if (5 != sscanf(line+3, "%s %d %d %d %c", ai.co.cname, &ai.co.nbases,
			&ai.co.nreads, &ai.co.nseg, &dir))
	    return NULL;
	ai.co.dir = dir == 'U' ? 0 : 1;
	ai.co.seq = (char *)malloc(ai.co.nbases+1);
	last_depadded_len = 0;
	while (NULL != zfgets(line, MAX_LINE_LEN, fp)) {
	    size_t l = strlen(line);
	    if (line[l-1] == '\n')
		l--;
	    if (l == 0)
		break;

	    if (pos+l > ai.co.nbases) {
		fprintf(stderr, "Sequence in CO line is longer than "
			"declared length\n");
		return NULL;
	    }
	    for (i = 0; i < l; i++) {
		if (line[i] != '*')
		    last_depadded_len++;
	    }
	    strncpy(&ai.co.seq[pos], line, l);
	    pos += l;
	}
	ai.co.seq[pos] = 0;
	if (pos != ai.co.nbases) {
	    fprintf(stderr, "Sequence in CO line does not match "
		    "declared length\n");
	    return NULL;
	}

    } else if (strncmp(line, "BQ", 2) == 0) {
	int pos = 0;

	if (ai.type == ACE_CO) {
	    ai.bq.nqual = last_depadded_len;
	} else {
	    fprintf(stderr, "BQ line does not appear immediately after a "
		    "CO line\n");
	    return NULL;
	}

	ai.type = ACE_BQ;
	ai.bq.qual = (unsigned char *)malloc(ai.bq.nqual);
	while (NULL != zfgets(line, MAX_LINE_LEN, fp)) {
	    char *cp1 = line, *cp2;
	    long l;

	    if (*line == '\n') {
		if (pos != ai.bq.nqual) {
		    fprintf(stderr, "Quality in BQ line does not match the "
			    "declared length\n");
		    return NULL;
		}
		break;
	    }

	    for (;;cp1 = cp2) {
		l = strtol(cp1, &cp2, 10);
		if (cp2 == cp1)
		    break;
		
		ai.bq.qual[pos++] = (int)l;

		if (pos > ai.bq.nqual) {
		    fprintf(stderr, "Quality in BQ line is longer than "
			    "declared length\n");
		    return NULL;
		}
	    }
	}

    } else if (strncmp(line, "AF", 2) == 0) {
	char dir;
	char fmt[256];
	ai.type = ACE_AF;
	sprintf(fmt, "%%%ds %%c %%d", MAX_NAME);
	if (3 != sscanf(line+3, fmt, ai.af.rname, &dir, &ai.af.start))
	    return NULL;
	ai.af.dir = dir == 'U' ? 0 : 1;
	
    } else if (strncmp(line, "BS", 2) == 0) {
	ai.type = ACE_BS;
	/* Unimplemented */

    } else if (strncmp(line, "RD", 2) == 0) {
	char fmt[256];
	int pos = 0;

	ai.type = ACE_RD;
	sprintf(fmt, "%%%ds %%d %%d %%d", MAX_NAME);
	if (4 != sscanf(line+3, fmt, ai.rd.rname,
			&ai.rd.nbases, &ai.rd.ninfo, &ai.rd.ntags))
	    return NULL;

	ai.rd.seq = (char *)malloc(ai.rd.nbases+1);
	while (NULL != zfgets(line, MAX_LINE_LEN, fp)) {
	    size_t l = strlen(line);

	    if (line[l-1] == '\n')
		l--;
	    if (l == 0)
		break;

	    if (pos+l > ai.rd.nbases) {
		fprintf(stderr, "Sequence in RD line is longer than "
			"declared length\n");
		return NULL;
	    }
	    strncpy(&ai.rd.seq[pos], line, l);
	    pos += l;
	}

	ai.rd.seq[pos] = 0;
	if (pos != ai.rd.nbases) {
	    fprintf(stderr, "Sequence in RD line does not match "
		    "declared length\n");
	    return NULL;
	}

    } else if (strncmp(line, "QA", 2) == 0) {
	ai.type = ACE_QA;
	if (4 != sscanf(line+3, "%d %d %d %d",
			&ai.qa.qstart, &ai.qa.qend,
			&ai.qa.astart, &ai.qa.aend))
	    return NULL;

    } else if (strncmp(line, "DS", 2) == 0) {
	char *cp1, *cp2;

	ai.type = ACE_DS;

	if ((cp1 = strstr(line, "CHROMAT_FILE:"))) {
	    for (cp1 += 13; *cp1 == ' '; cp1++);
	    cp2 = strchr(cp1, ' ');
	    strncpy(ai.ds.chromat, cp1, cp2-cp1);
	    ai.ds.chromat[cp2-cp1] = 0;
	} else {
	    ai.ds.chromat[0] = 0;
	}

	if ((cp1 = strstr(line, "PHD_FILE:"))) {
	    for (cp1 += 9; *cp1 == ' '; cp1++);
	    cp2 = strchr(cp1, ' ');
	    strncpy(ai.ds.phd, cp1, cp2-cp1);
	    ai.ds.phd[cp2-cp1] = 0;
	} else {
	    ai.ds.phd[0] = 0;
	}

	if ((cp1 = strstr(line, "TEMPLATE:"))) {
	    for (cp1 += 9; *cp1 == ' '; cp1++);
	    cp2 = strchr(cp1, ' ');
	    strncpy(ai.ds.tname, cp1, cp2-cp1);
	    ai.ds.tname[cp2-cp1] = 0;
	} else {
	    ai.ds.tname[0] = 0;
	}

	/* Other parsing unimplemented at the moment */
	ai.ds.time = 0;
	ai.ds.chem = 0;
	ai.ds.dye = 0;
	ai.ds.dir = 0;

    } else if (strncmp(line, "CT{", 3) == 0) {
	ai.type = ACE_CT;
	goto RT_code;

    } else if (strncmp(line, "WA{", 3) == 0) {
	ai.type = ACE_WA;
	goto RT_code;

    } else if (strncmp(line, "WR{", 3) == 0) {
	ai.type = ACE_WR;
	goto RT_code;

    } else if (strncmp(line, "RT{", 3) == 0) {
	size_t allocated, used;
	ai.type = ACE_RT;

	/* Shared between CT, WA and RT as the structures are compatible */
    RT_code:
	ai.rt.text = NULL;
	allocated = used = 0;

	/* Consume up to the next "}" */
	while (NULL != zfgets(line, MAX_LINE_LEN, fp)) {
	    size_t l;
	    if (line[0] == '}' && line[1] == '\n')
		break;

	    l = strlen(line);
	    while (used + l >= allocated) {
		allocated += 8192;
		ai.rt.text = (char *)realloc(ai.rt.text, allocated);
	    }
	    strcpy(&ai.rt.text[used], line);
	    used += l;
	}

    } else {
	ai.type = 0;
	fprintf(stderr, "Unknown ACE line: %s\n", line);
	return NULL;
    }
    
    return &ai;
}

typedef struct {
    char name[MAX_NAME];
    int dir;
    int pos;
} af_line;

/*
 * Parses a new ACE format file passed in.
 *
 * Returns 0 on success
 *	  -1 on error
 */
int parse_ace(GapIO *io, char *ace_fn, tg_args *a) {
    ace_item_t *ai;
    zfp *fp;
    af_line *af = NULL;
    int af_count = 0, seq_count = -1, nseqs = 0, nseqs_tot = 0, ncontigs = 0;
    contig_t *c = NULL;
    tg_pair_t *pair = NULL;
    seq_t seq;
    char tname[1024];
    
    set_dna_lookup(); /* initialise complement table */

    if (NULL == (fp = zfopen(ace_fn, "r")))
	return -1;

    if (a->pair_reads) {
	pair = create_pair(a->pair_queue);
    }

    while (ai = next_ace_item(fp)) {
	switch (ai->type) {
	case ACE_AS:
	    nseqs_tot = ai->as.nreads;
	    break;

	case ACE_CO:
	    /* New contig => flush out last read in old contig */
	    if (seq_count >= 0 && !(seq.left == -1 && seq.right == -1)) {
		save_range_sequence(io, &seq, ACE_MQUAL, pair,
				    (pair && *tname), tname, c, a,
				    GRANGE_FLAG_TYPE_SINGLE,
				    NULL, NULL);
		nseqs++;

		if ((nseqs & 0x3fff) == 0) {
		    printf("\r%5.2f%%", (100.0 * nseqs) / nseqs_tot);
		    //HacheTableStats(io->cache, stdout);
		    fflush(stdout);
		    cache_flush(io);
		}
		seq.left = -1;
		seq.right = -1;
	    }

	    create_new_contig(io, &c, ai->co.cname, a->merge_contigs);

	    if (af)
		free(af);
	    af = (af_line *)calloc(ai->co.nreads, sizeof(*af));
	    seq_count = -1, af_count = 0;
	    if (NULL == af)
		return 1;

	    ncontigs++;

	    //fprintf(stderr, "Processing contig %s\n", ai->co.cname);
	    break;

	case ACE_AF:
	    assert(af);
	    /* Accumulate assembled-from lines */
	    strcpy(af[af_count].name, ai->af.rname);
	    af[af_count].dir = ai->af.dir;
	    af[af_count].pos = ai->af.start;
	    af_count++;
	    break;

	case ACE_RD:
	    if (seq_count++ >= 0 && !(seq.left == -1 && seq.right == -1)) {
		save_range_sequence(io, &seq, ACE_MQUAL, pair,
				    (pair && *tname), tname, c, a,
				    GRANGE_FLAG_TYPE_SINGLE,
				    NULL, NULL);
		nseqs++;

		if ((nseqs & 0x3fff) == 0) {
		    printf("\r%5.2f%%", (100.0 * nseqs) / nseqs_tot);
		    //HacheTableStats(io->cache, stdout);
		    fflush(stdout);
		    cache_flush(io);
		}
	    }

	    /* Add readings, assumed in same order as AF lines */
	    assert(af);
	    if (strcmp(ai->rd.rname, af[seq_count].name)) {
		fprintf(stderr, "AF lines and RD lines not in same order\n");
		fprintf(stderr, "%s\n%s\n",
			ai->rd.rname,
			af[seq_count].name);
		return 1;
	    }

	    /* Reverse compliment if neeed */
	    if (af[seq_count].dir != 0)
		complement_seq(ai->rd.seq, ai->rd.nbases);

	    /* Fill out a seq_t struct */
	    memset(&seq, 0, sizeof(seq));
	    seq.rec = 0;
	    seq.pos = af[seq_count].pos;
	    seq.len = af[seq_count].dir == 0 ? ai->rd.nbases : -ai->rd.nbases;
	    seq.mapping_qual = ACE_MQUAL;
	    seq.left = 1;
	    seq.right = ai->rd.nbases;
	    seq.flags = af[seq_count].dir == 0 ? 0 : SEQ_COMPLEMENTED;
	    seq.name_len = strlen(ai->rd.rname);
	    seq.name = (char *)malloc(seq.name_len+1+2*ai->rd.nbases);
	    strcpy(seq.name, ai->rd.rname);
	    seq.seq = seq.name + seq.name_len + 1;
	    memcpy(seq.seq, ai->rd.seq, ai->rd.nbases);
	    seq.trace_name = NULL;
	    seq.alignment = NULL;
	    seq.conf = seq.seq + ai->rd.nbases;
	    memset(seq.conf, 4, ai->rd.nbases);

	    /* NB: we don't write it out yet until we see DS the line */
	    break;

	case ACE_QA:
	    if (ai->qa.aend == -1 && ai->qa.astart == -1) {
		/* Unaligned read */
		seq.left = seq.right = -1;
	    } else {
		if (seq.flags & SEQ_COMPLEMENTED) {
		    seq.left  = ABS(seq.len) - ai->qa.aend   + 1;
		    seq.right = ABS(seq.len) - ai->qa.astart + 1;
		} else {
		    seq.left  = ai->qa.astart;
		    seq.right = ai->qa.aend;
		}
	    }
	    break;

	case ACE_DS:
	    strncpy(tname, ai->ds.tname, 1024);
	    seq.template_name_len = *tname ? strlen(tname) : seq.name_len;
	    break;
	}
    }

    if (seq_count++ >= 0 && !(seq.left == -1 && seq.right == -1)) {
	save_range_sequence(io, &seq, ACE_MQUAL, pair,
			    (pair && *tname), tname, c, a,
			    GRANGE_FLAG_TYPE_SINGLE,
			    NULL, NULL);
	nseqs++;

	if ((nseqs & 0x3fff) == 0) {
	    printf("\r%5.2f%%", (100.0 * nseqs) / nseqs_tot);
	    //HacheTableStats(io->cache, stdout);
	    fflush(stdout);
	    cache_flush(io);
	}
    }

    puts("");

    if (pair && !a->fast_mode) {    
	finish_pairs(io, pair);
    }
    
    cache_flush(io);
    zfclose(fp);

    if (pair) delete_pair(pair);

    if (af)
	free(af);

    if (c)
	cache_decr(io, c);

    return 0;
}
