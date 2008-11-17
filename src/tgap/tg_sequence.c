#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "tg_gio.h"

#ifndef ABS
#    define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

/*
 * Given a seq_t struct this allocates a new sequence in the database
 * and copies the contents of 's' into it. If 's' is NULL it simply allocates
 * the new sequence and does nothing with it.
 *
 * Returns the record number on success
 *        -1 on failure
 */
int sequence_new_from(GapIO *io, seq_t *s) {
    int rec;

    rec = io->iface->seq.create(io->dbh, s);
    return rec;
}

/*
 * Sets the sequence position
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_position(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->pos = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence length
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_len(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->len = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence left clip
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_left(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->left = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence right clip
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_right(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->right = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence mapping quality
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_mapping_qual(GapIO *io, seq_t **s, uint8_t value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->mapping_qual = value;
    *s = n;

    return 0;
}

/*
 * Sets the record number of the sequence paired with this (0 => none)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_other_end(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->other_end = value;
    *s = n;

    return 0;
}

/*
 * Sets the parent record type (GT_Template, GT_Ligation, etc)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_parent_type(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->parent_type = value;
    *s = n;

    return 0;
}

/*
 * Sets the parent record number.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_parent_rec(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->parent_rec = value;
    *s = n;

    return 0;
}

/*
 * Sets the flags. See SEQ_* defines in tg_struct.h
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_flags(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->flags = value;
    *s = n;

    return 0;
}

/*
 * Sets the sequence technology. See STECH_* defines in tg_struct.h
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_seq_tech(GapIO *io, seq_t **s, int value) {
    seq_t *n;
    if (!(n = cache_rw(io, *s)))
	return -1;

    n->seq_tech = value;
    *s = n;

    return 0;
}

/*
 * Returns the size needed to store confidence values in this sequence.
 * Ie 1 or 4 per base.
 */
#define sequence_conf_size(s) ((s)->format == SEQ_FORMAT_CNF4 ? 4 : 1)

/*
 * Sets a sequence name.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_name(GapIO *io, seq_t **s, char *name) {
    size_t extra_len;
    seq_t *n;
    char *tmp, *cp;

    if (!(n = cache_rw(io, *s)))
	return -1;

    extra_len = strlen(name)+1 + n->trace_name_len+1 +
	n->alignment_len+1 + ABS(n->len)*(1+sequence_conf_size(n));
    n = cache_item_resize(n, sizeof(*n) + extra_len);
    if (NULL == n)
	return -1;

    n->name = (char *)&n->data;
    n->name_len = strlen(name);
    n->trace_name = n->name + n->name_len + 1;
    n->alignment = n->trace_name + n->trace_name_len + 1;
    n->seq = n->alignment + n->alignment_len + 1;
    n->conf = n->seq + ABS(n->len);

    /* Shift and insert name */
    cp = tmp = malloc(extra_len);
    strcpy(cp, name);
    cp += n->name_len+1;
    strcpy(cp, n->trace_name);
    cp += n->trace_name_len;
    strcpy(cp, n->alignment);
    cp += n->alignment_len;
    memcpy(cp, n->seq, n->len);
    memcpy(cp, n->conf, n->len * sequence_conf_size(n));
    memcpy(&n->data, tmp, extra_len);
    free(tmp);

    *s = n;
    return 0;
}

/*
 * Sets a trace name.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_trace_name(GapIO *io, seq_t **s, char *trace_name) {
    size_t extra_len;
    seq_t *n;
    char *tmp, *cp;

    if (!(n = cache_rw(io, *s)))
	return -1;

    if (!trace_name || 0 == strcmp(n->name, trace_name))
	trace_name = "";

    extra_len = n->name_len+1 + strlen(trace_name)+1 +
	n->alignment_len+1 + ABS(n->len)*(1+sequence_conf_size(n));
    n = cache_item_resize(n, extra_len);
    if (NULL == n)
	return -1;

    n->name = (char *)&n->data;
    n->trace_name = n->name + n->name_len + 1;
    n->trace_name_len = strlen(trace_name);
    n->alignment = n->trace_name + n->trace_name_len + 1;
    n->seq = n->alignment + n->alignment_len + 1;
    n->conf = n->seq + ABS(n->len);

    /* Shift and insert name */
    cp = tmp = malloc(extra_len);
    strcpy(cp, n->name);
    cp += n->name_len+1;
    strcpy(cp, trace_name);
    cp += n->trace_name_len;
    strcpy(cp, n->alignment);
    cp += n->alignment_len;
    memcpy(cp, n->seq, n->len);
    memcpy(cp, n->conf, n->len * sequence_conf_size(n));
    memcpy(&n->data, tmp, extra_len);
    free(tmp);

    *s = n;
    return 0;
}

int sequence_set_seq (GapIO *io, seq_t **s, char *seq) {return -1;}
int sequence_set_conf(GapIO *io, seq_t **s, char *conf) {return -1;}


/*
 * Insert into position 'pos' of a sequence, where pos starts from 0
 * (before first base).
 */
int sequence_insert_base(GapIO *io, seq_t **s, int pos,
			 char base, char conf) {
    seq_t *n;
    int p2, alen;

    if (!(n = cache_rw(io, *s)))
	return -1;

    n = cache_item_resize(n, sizeof(*n) + 
			  n->name_len+1 +
			  n->trace_name_len+1 +
			  n->alignment_len+1 +
			  (ABS(n->len)+1)*(1+sequence_conf_size(n)));
    if (NULL == n)
	return -1;

    n->name = (char *)&n->data;
    n->trace_name = n->name + n->name_len + 1;
    n->alignment = n->trace_name + n->trace_name_len + 1;
    n->seq = n->alignment + n->alignment_len + 1;
    n->conf = n->seq + ABS(n->len);

    p2 = n->len < 0
	? -n->len - pos
	: pos;
    alen = ABS(n->len);

    /* Shift */
    memmove(&n->seq[p2+1], &n->seq[p2], alen-p2+alen);
    n->conf++;
    if (n->format == SEQ_FORMAT_CNF4) {
	memmove(&n->conf[(p2+1)*4], &n->conf[p2*4], (alen-p2)*4);
    } else {
	memmove(&n->conf[p2+1], &n->conf[p2], alen-p2);
    }

    if (n->len < 0)
	n->len--;
    else
	n->len++;

    /* Set */
    n->seq [p2] = base;
    if (n->format == SEQ_FORMAT_CNF4) {
	double remainder = -4.34294482*log(2+3*pow(10, conf/10.0));
	switch (base) {
	case 'A': case 'a':
	    n->conf[p2*4+0] = conf;
	    n->conf[p2*4+1] = remainder;
	    n->conf[p2*4+2] = remainder;
	    n->conf[p2*4+3] = remainder;
	    break;
	case 'C': case 'c':
	    n->conf[p2*4+0] = remainder;
	    n->conf[p2*4+1] = conf;
	    n->conf[p2*4+2] = remainder;
	    n->conf[p2*4+3] = remainder;
	    break;
	case 'G': case 'g':
	    n->conf[p2*4+0] = remainder;
	    n->conf[p2*4+1] = remainder;
	    n->conf[p2*4+2] = conf;
	    n->conf[p2*4+3] = remainder;
	    break;
	case 'T': case 't':
	    n->conf[p2*4+0] = remainder;
	    n->conf[p2*4+1] = remainder;
	    n->conf[p2*4+2] = remainder;
	    n->conf[p2*4+3] = conf;
	    break;
	default:
	    n->conf[p2*4+0] = -5;
	    n->conf[p2*4+1] = -5;
	    n->conf[p2*4+2] = -5;
	    n->conf[p2*4+3] = -5;
	    break;
	}
    } else {
	n->conf[p2] = conf;
    }

    if (p2 < n->left)
	n->left++;

    if (p2 < n->right)
	n->right++;

    return 0;
}

/*
 * Delete position 'pos' of a sequence, where pos starts from 0
 */
int sequence_delete_base(GapIO *io, seq_t **s, int pos) {
    seq_t *n;
    int p2, alen;

    if (!(n = cache_rw(io, *s)))
	return -1;

    if (n->len < 0)
	n->len++;
    else
	n->len--;

    if (p2 < n->left)
	n->left--;

    if (p2 < n->right)
	n->right--;

    p2 = n->len < 0
	? -n->len - pos
	: pos;
    alen = ABS(n->len);

    if (p2 >= alen)
	return 0;

    /* Shift */
    memmove(&n->seq[p2], &n->seq[p2+1], alen-p2+alen-1);
    n->conf--;
    if (n->format == SEQ_FORMAT_CNF4) {
	memmove(&n->conf[p2*4], &n->conf[(p2+1)*4], (alen-1-p2)*4);
    } else {
	memmove(&n->conf[p2], &n->conf[p2+1], alen-1-p2);
    }

    return 0;
}

/* ------------------------------------------------------------------------ 
 * Trivial one-off sequence query functions
 */
int seq_pos(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_pos(&s);
}

int seq_len(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_len(&s);
}

int seq_left(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_left(&s);
}

int seq_right(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_right(&s);
}

int seq_mapping_qual(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_mapping_qual(&s);
}

char *seq_name(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_name(&s);
}

char *seq_seq(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_seq(&s);
}

char *seq_conf(GapIO *io, int rec) {
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
    return sequence_get_conf(&s);
}

/* ------------------------------------------------------------------------ */
/* Sequence manipulation - ripped out of Staden Package's seq_utils.c */

/*
 * Reverses and complements a piece of DNA
 */
static int complementary_base[256];
void complement_seq_conf(char *seq, char *conf, int seq_len, int nconf) {
    int i, middle, j;
    char temp, t[4];
    static int init = 0;

    if (!init) {
	for (i = 0; i < 256; i++)
	    complementary_base[i] = i;

	complementary_base['a'] = 't';
	complementary_base['c'] = 'g';
	complementary_base['g'] = 'c';
	complementary_base['t'] = 'a';
	complementary_base['u'] = 'a';
	complementary_base['A'] = 'T';
	complementary_base['C'] = 'G';
	complementary_base['G'] = 'C';
	complementary_base['T'] = 'A';
	complementary_base['U'] = 'A';

	complementary_base['n'] = 'n';
	complementary_base['-'] = '-';
	complementary_base['b'] = 'v';
	complementary_base['d'] = 'h';
	complementary_base['h'] = 'd';
	complementary_base['k'] = 'm';
	complementary_base['m'] = 'k';
	complementary_base['r'] = 'y';
	complementary_base['s'] = 's';
	complementary_base['v'] = 'b';
	complementary_base['w'] = 'w';
	complementary_base['y'] = 'r';

	complementary_base['B'] = 'V';
	complementary_base['D'] = 'H';
	complementary_base['H'] = 'D';
	complementary_base['K'] = 'M';
	complementary_base['M'] = 'K';
	complementary_base['R'] = 'Y';
	complementary_base['S'] = 'S';
	complementary_base['V'] = 'B';
	complementary_base['W'] = 'W';
	complementary_base['Y'] = 'R';
	init = 1;
    }

    middle = seq_len/2;
    if (nconf == 1) {
	for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
	    temp = complementary_base [ (unsigned char) seq[i] ];
	    seq[i] = complementary_base [ (unsigned char) seq[j] ];
	    seq[j] = temp;
	    temp = conf[i];
	    conf[i] = conf[j];
	    conf[j] = temp;
	}
    } else if (nconf == 4) {
	for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
	    temp = complementary_base [ (unsigned char) seq[i] ];
	    seq[i] = complementary_base [ (unsigned char) seq[j] ];
	    seq[j] = temp;
	    t[0] = conf[i*4+0];
	    t[1] = conf[i*4+1];
	    t[2] = conf[i*4+2];
	    t[3] = conf[i*4+3];
	    conf[i*4+0] = conf[j*4+3];
	    conf[i*4+1] = conf[j*4+2];
	    conf[i*4+2] = conf[j*4+1];
	    conf[i*4+3] = conf[j*4+0];
	    conf[j*4+0] = t[3];
	    conf[j*4+1] = t[2];
	    conf[j*4+2] = t[1];
	    conf[j*4+3] = t[0];
	}
    } else {
	fprintf(stderr, "Unsupported number of confidence values per base\n");
    }

    if ( seq_len % 2 )
      seq[middle] = complementary_base [ (unsigned char) seq[middle] ];
}

seq_t *dup_seq(seq_t *s) {
    size_t len = sizeof(seq_t) - sizeof(char *) + s->name_len+1 +
	s->trace_name_len+1 + s->alignment_len+1 +
	(ABS(s->len)+1)*(1+sequence_conf_size(s));
    seq_t *d = (seq_t *)calloc(1, len);

    memcpy(d, s, len);
    d->name = (char *)&d->data;
    d->trace_name = d->name + d->name_len + 1;
    d->alignment = d->trace_name + d->trace_name_len + 1;
    d->seq = d->alignment + d->alignment_len + 1;
    d->conf = d->seq + ABS(d->len);

    return d;
}

void complement_seq_t(seq_t *s) {
    int tmp, alen;

    alen = ABS(s->len);
    complement_seq_conf(s->seq, s->conf, alen,
			s->format == SEQ_FORMAT_CNF4 ? 4 : 1);
    s->len *= -1;

    tmp = s->left;
    s->left  = alen - (s->right-1);
    s->right = alen - (tmp-1);
}

GRec sequence_index_query(GapIO *io, char *name) {
    return io->iface->seq.index_query(io->dbh, name);
}

int sequence_index_update(GapIO *io, char *name, int name_len, GRec rec) {
    char n2[1024];
    GRec r;
    sprintf(n2, "%.*s", name_len, name);

    r = io->iface->seq.index_add(io->dbh, n2, rec);
    if (r == -1)
	return -1;

    if (r != io->db->seq_name_index) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = r;
    }

    return 0;
}

/*
 * Finds the contig number and position of a sequence record number.
 */
int sequence_get_position(GapIO *io, GRec snum, int *contig, int *pos) {
    bin_index_t *bin;
    int bnum, i;
    int offset = 0, found = 0;
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, snum);

    /* Find the position of this sequence within the bin */
    bnum = s->bin;
    bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
        range_t *r = arrp(range_t, bin->rng, i);
	if (r->rec == snum) {
	    found = 1;
	    offset = r->start;
	    break;
	}
    }

    if (!found)
	return -1;

    /* FIXME: Complemented coordinates need more fixes here */
    if (bin->flags & BIN_COMPLEMENTED) {
	offset = bin->size-1 - offset;
    }

    offset += bin->pos;

    /* Find the position of this bin relative to the contig itself */
    while (bin->parent_type == GT_Bin) {
	bnum = bin->parent;
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	offset += bin->pos;
    }
    
    assert(bin->parent_type == GT_Contig);
    *contig = bin->parent;

    *pos = offset;

    return 0;
}

/*
 * Given the record number for a sequence this returns the record
 * number for the contig containing it.
 */
int sequence_get_contig(GapIO *io, GRec snum) {
    bin_index_t *bin;
    int bnum;
    seq_t *s = (seq_t *)cache_search(io, GT_Seq, snum);

    /* Bubble up bins until we hit the root */
    for (bnum = s->bin; bnum; bnum = bin->parent) {
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	if (bin->parent_type != GT_Bin)
	    break;
    }

    assert(bin->parent_type == GT_Contig);
    return bin->parent;
}

/*
 * ---------------------------------------------------------------------------
 * Base editing functions
 */

#define MAX4(ip)         \
((ip)[0] > (ip)[1]       \
 ? ((ip)[0] > (ip)[2]    \
    ? ((ip)[0] > (ip)[3] \
       ? (ip)[0]         \
       : (ip)[3])        \
    : ((ip)[2] > (ip)[3] \
       ? (ip)[2]         \
       : (ip)[3]))       \
 : ((ip)[1] > (ip)[2]    \
    ? ((ip)[1] > (ip)[3] \
       ? (ip)[1]         \
       : (ip)[3])        \
    : ((ip)[2] > (ip)[3] \
       ? (ip)[2]         \
       : (ip)[3])))


/*
 * Operates on position 'pos' in the displayed orientation rather than the
 * stored orientation.
 */
int sequence_get_base(GapIO *io, seq_t **s, int pos, char *base, int *conf) {
    seq_t *n = *s;

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    if (n->len > 0) {
	if (base)
	    *base = n->seq[pos];
	if (conf) {
	    if (n->format != SEQ_FORMAT_CNF4) {
		*conf = n->conf[pos];
	    } else {
		*conf = MAX4(&n->conf[pos*4]);
	    }
	}
    } else {
	if (base)
	    *base = complementary_base[(int)(n->seq[-n->len-1 - pos])];
	if (conf) {
	    if (n->format != SEQ_FORMAT_CNF4) {
		*conf = n->conf[-n->len-1 - pos];
	    } else {
		*conf = MAX4(&n->conf[(-n->len-1 - pos)*4]);
	    }
	}
    }

    return 0;
}

static double logodds2log[256];
static double *lo2l = logodds2log+128;

static double logodds2remainder[256];
static double *lo2r = logodds2remainder+128;

static unsigned char logodds2phred[256];
static unsigned char *lo2ph = logodds2phred+128;

static int lookup_init = 0;

/*
 * A recalibration table for solexa log-odds scores.
 * Derived by looking at a lane of PhiX control data, but not rigorously
 * tested yet.
 */
static int recalibrate[81] = {
    /*-40 */ -35.5156,
    /*-39 */ -30.2016,
    /*-38 */ -29.808,
    /*-37 */ -30.0427,
    /*-36 */ -30.3126,
    /*-35 */ -28.956,
    /*-34 */ -28.3991,
    /*-33 */ -29.4153,
    /*-32 */ -28.9522,
    /*-31 */ -28.9222,
    /*-30 */ -27.8603,
    /*-29 */ -28.0694,
    /*-28 */ -27.0335,
    /*-27 */ -26.5597,
    /*-26 */ -27.4361,
    /*-25 */ -26.33,
    /*-24 */ -25.8657,
    /*-23 */ -24.935,
    /*-22 */ -24.391,
    /*-21 */ -24.114,
    /*-20 */ -23.4288,
    /*-19 */ -22.4242,
    /*-18 */ -22.0634,
    /*-17 */ -21.4384,
    /*-16 */ -20.5132,
    /*-15 */ -19.6556,
    /*-14 */ -19.2086,
    /*-13 */ -17.9905,
    /*-12 */ -16.8977,
    /*-11 */ -15.7831,
    /*-10 */ -15.1582,
    /* -9 */ -13.7783,
    /* -8 */ -12.4399,
    /* -7 */ -11.1721,
    /* -6 */ -10.09,
    /* -5 */ -6.78853,
    /* -4 */ -6.6224,
    /* -3 */ -4.41534,
    /* -2 */ -2.83268,
    /* -1 */ -1.3449,
    /*  0 */ 1.0271,
    /*  1 */ 3.4759,
    /*  2 */ 4.6632,
    /*  3 */ 6.0888,
    /*  4 */ 7.4930,
    /*  5 */ 8.7698,
    /*  6 */ 10.208,
    /*  7 */ 11.143,
    /*  8 */ 12.224,
    /*  9 */ 13.530,
    /* 10 */ 14.604,
    /* 11 */ 15.313,
    /* 12 */ 16.249,
    /* 13 */ 17.304,
    /* 14 */ 18.278,
    /* 15 */ 18.681,
    /* 16 */ 19.367,
    /* 17 */ 20.487,
    /* 18 */ 20.53,
    /* 19 */ 21.029,
    /* 20 */ 21.694,
    /* 21 */ 22.679,
    /* 22 */ 22.689,
    /* 23 */ 22.79,
    /* 24 */ 23.809,
    /* 25 */ 24.230,
    /* 26 */ 25.431,
    /* 27 */ 24.551,
    /* 28 */ 25.035,
    /* 29 */ 25.411,
    /* 30 */ 25.155,
    /* 31 */ 25.84,
    /* 32 */ 26.471,
    /* 33 */ 25.897,
    /* 34 */ 25.9,
    /* 35 */ 26.689,
    /* 36 */ 26.846,
    /* 37 */ 27.288,
    /* 38 */ 26.707,
    /* 39 */ 26.671,
    /* 40 */ 30.590,
};

/*
 * Operates on position 'pos' in the displayed orientation rather than the
 * stored orientation.
 *
 * As per sequence_get_base but conf is an array of 4 confidence values
 * returned as log(probabilities).
 * If only one is present it is taken to be a phred score and we compute
 * the remainder using (1-p)/3. Otherwise we assume we store 4 log-odds
 * scores instead.
 */
int sequence_get_base4(GapIO *io, seq_t **s, int pos,
		       char *base, double *conf) {
    seq_t *n = *s;

    if (!lookup_init) {
	int i;
	lookup_init = 1;

	/* Log odds value to log(P) */
	for (i = -128; i < 128; i++) {
	    double p = 1 / (1 + pow(10, -i / 10.0));
	    lo2l[i] = log(p);
	    lo2r[i] = log((1-p)/3);
	    lo2ph[i] = 10*log(1+pow(10, i/10.0))/log(10)+0.4999;
	}
    }

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    if (base) {
	if (n->len > 0)
	    *base = n->seq[pos];
	else
	    *base = complementary_base[(int)(n->seq[-n->len-1 - pos])];
    }

    if (conf) {
	if (n->len < 0)
	    pos = -n->len-1 - pos;

	if (n->format != SEQ_FORMAT_CNF4) {
	    switch(n->seq[pos]) {
	    case 'A': case 'a':
		conf[0] = lo2l[n->conf[pos]];
		conf[1] = lo2r[n->conf[pos]];
		conf[2] = lo2r[n->conf[pos]];
		conf[3] = lo2r[n->conf[pos]];
		break;
	    case 'C': case 'c':
		conf[0] = lo2r[n->conf[pos]];
		conf[1] = lo2l[n->conf[pos]];
		conf[2] = lo2r[n->conf[pos]];
		conf[3] = lo2r[n->conf[pos]];
		break;
	    case 'G': case 'g':	
		conf[0] = lo2r[n->conf[pos]];
		conf[1] = lo2r[n->conf[pos]];
		conf[2] = lo2l[n->conf[pos]];
		conf[3] = lo2r[n->conf[pos]];
		break;
	    case 'T': case 't':
		conf[0] = lo2r[n->conf[pos]];
		conf[1] = lo2r[n->conf[pos]];
		conf[2] = lo2r[n->conf[pos]];
		conf[3] = lo2l[n->conf[pos]];
		break;
	    default:
		conf[0] = 0;
		conf[1] = 0;
		conf[2] = 0;
		conf[3] = 0;
		break;
	    }
	} else {
	    int i;
	    for (i = 0; i < 4; i++) {
		//double p = pow(10, n->conf[pos*4+i]/10.0);
		//conf[n->len < 0 ? 4-i : i] = p/(1+p);
		double p = n->conf[pos*4+i];
		//p = recalibrate[(int)(p+40+.49)] / 10.0;
		p /= 10.0;
		conf[i] = p*log(10) - log(1+pow(10, p));
	    }
	}

	if (n->len < 0) {
	    double tmp;
	    tmp = conf[0]; conf[0] = conf[3]; conf[3] = tmp;
	    tmp = conf[1]; conf[1] = conf[2]; conf[2] = tmp;
	}

	richard_munge_conf(conf);
	//rob_munge_conf(conf);
    }

    return 0;
}

/*
 * Simulates Richard's 4-confidence value system.
 * This stores the top value as-is. The 3rd and 4th highest listed as a
 * fraction of the 2nd highest and the 2nd highest value is implicit (1 minus
 * the rest).
 */
void richard_munge_conf(double *lp) {
    double p1 = -100, p2 = -100, p3 = -100, p4 = -100, p34;
    int i1 = 0, i2 = 0, i3 = 0, i4 = 0;
    int i;
    double r3, r4;

    /* Sort */
    for (i = 0; i < 4; i++) {
	if (p1 < lp[i]) {
	    p4 = p3;
	    i4 = i3;
	    p3 = p2;
	    i3 = i2;
	    p2 = p1;
	    i2 = i1;
	    p1 = lp[i];
	    i1 = i;
	} else if (p2 < lp[i]) {
	    p4 = p3;
	    i4 = i3;
	    p3 = p2;
	    i3 = i2;
	    p2 = lp[i];
	    i2 = i;
	} else if (p3 < lp[i]) {
	    p4 = p3;
	    i4 = i3;
	    p3 = lp[i];
	    i3 = i;
	} else if (p4 < lp[i]) {
	    p4 = lp[i];
	    i4 = i;
	}
    }

    /* Rescale */
    p1 = exp(p1);
    p2 = exp(p2);
    p3 = exp(p3);
    p4 = exp(p4);

    /* Simulate truncation using 3 bits for ratios */
    r3 = -10*log(p3/p2)/log(10);
    if      (r3 <  1.5) r3 = 0;
    else if (r3 <  4.5) r3 = 3;
    else if (r3 <  8.0) r3 = 6;
    else if (r3 < 12.5) r3 = 10;
    else if (r3 < 17.5) r3 = 15;
    else if (r3 < 25.5) r3 = 20;
    else if (r3 < 40.0) r3 = 30;
    else r3 = 50;
    r3 = pow(10,r3/-10.0);

    r4 = -10*log(p4/p2)/log(10);
    if      (r4 <  1.5) r4 = 0;
    else if (r4 <  4.5) r4 = 3;
    else if (r4 <  8.0) r4 = 6;
    else if (r4 < 12.5) r4 = 10;
    else if (r4 < 17.5) r4 = 15;
    else if (r4 < 25.5) r4 = 20;
    else if (r4 < 40.0) r4 = 30;
    else r4 = 50;
    r4 = pow(10,r4/-10.0);

    /*
     * From ratio r3 and r4 recompute p2, p3, and p4 using:
     *    1-p1 = p2 + r3*p2 + r4*p3
     *
     * => p2(1+r3+r4) = 1-p1
     * => p2 = (1-p1)/(1+r3+r4)
     */
    p2 = (1-p1)/(1+r3+r4);
    p3 = r3*p2;
    p4 = r4*p2;

    /*
    printf("%5.1f %5.1f %5.1f %5.1f => ",
	   10 * log(exp(lp[0]) / (1-exp(lp[0]))) / log(10),
	   10 * log(exp(lp[1]) / (1-exp(lp[1]))) / log(10),
	   10 * log(exp(lp[2]) / (1-exp(lp[2]))) / log(10),
	   10 * log(exp(lp[3]) / (1-exp(lp[3]))) / log(10));
    */

    lp[i1] = log(p1);
    lp[i2] = log(p2);
    lp[i3] = log(p3);
    lp[i4] = log(p4);

    /*
    printf("%5.1f %5.1f %5.1f %5.1f\n",
	   10 * log(exp(lp[0]) / (1-exp(lp[0]))) / log(10),
	   10 * log(exp(lp[1]) / (1-exp(lp[1]))) / log(10),
	   10 * log(exp(lp[2]) / (1-exp(lp[2]))) / log(10),
	   10 * log(exp(lp[3]) / (1-exp(lp[3]))) / log(10));
    */
}

/*
 * Simulates Rob's 4-confidence value system.
 * This basically has top two values and the remaining two set to the average.
 */
void rob_munge_conf(double *lp) {
    double p1 = -100, p2 = -100, p3 = -100, p4 = -100, p34;
    double i1 = 0, i2 = 0;
    int i;

    for (i = 0; i < 4; i++) {
	if (p1 < lp[i]) {
	    p4 = p3;
	    p3 = p2;
	    p2 = p1;
	    i2 = i1;
	    p1 = lp[i];
	    i1 = i;
	} else if (p2 < lp[i]) {
	    p4 = p3;
	    p3 = p2;
	    p2 = lp[i];
	    i2 = i;
	} else if (p3 < lp[i]) {
	    p4 = p3;
	    p3 = lp[i];
	} else if (p4 < lp[i]) {
	    p4 = lp[i];
	}
    }

    p34 = log((1-((exp(p1)+exp(p2))/(exp(p1)+exp(p2)+exp(p3)+exp(p4))))/2.0);

    //printf("%f %f %f %f => ", lp[0], lp[1], lp[2], lp[3]);

    for (i = 0; i < 4; i++) {
	if (i == i1)
	    lp[i] = p1;
	else if (i == i2)
	    lp[i] = p2;
	else
	    lp[i] = p34;
    }

    //printf("%f %f %f %f\n", lp[0], lp[1], lp[2], lp[3]);
}

int sequence_replace_base(GapIO *io, seq_t **s, int pos, char base, int conf) {
    seq_t *n;

    if (!(n = cache_rw(io, *s)))
	return -1;

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    if (n->format != SEQ_FORMAT_CNF4) {
	if (n->len > 0) {
	    n->seq[pos] = base;
	    n->conf[pos] = conf;
	} else {
	    n->seq[-n->len-1 - pos] = complementary_base[(unsigned char)base];
	    n->conf[-n->len-1 - pos] = conf;
	}
    } else {
	double remainder = -4.34294482*log(2+3*pow(10, conf/10.0));
	int p2 = n->len > 0 ? pos : -n->len-1 -pos;

	n->seq[p2] = n->len > 0
	    ? base
	    : complementary_base[(unsigned char)base];

	switch(base) {
	case 'A': case 'a':
	    n->conf[p2*4+0] = conf;
	    n->conf[p2*4+1] = remainder;
	    n->conf[p2*4+2] = remainder;
	    n->conf[p2*4+3] = remainder;
	    break;
	case 'C': case 'c':
	    n->conf[p2*4+0] = remainder;
	    n->conf[p2*4+1] = conf;
	    n->conf[p2*4+2] = remainder;
	    n->conf[p2*4+3] = remainder;
	    break;
	case 'G': case 'g':
	    n->conf[p2*4+0] = remainder;
	    n->conf[p2*4+1] = remainder;
	    n->conf[p2*4+2] = conf;
	    n->conf[p2*4+3] = remainder;
	    break;
	case 'T': case 't':
	    n->conf[p2*4+0] = remainder;
	    n->conf[p2*4+1] = remainder;
	    n->conf[p2*4+2] = remainder;
	    n->conf[p2*4+3] = conf;
	    break;
	default:
	    n->conf[p2*4+0] = -5;
	    n->conf[p2*4+1] = -5;
	    n->conf[p2*4+2] = -5;
	    n->conf[p2*4+3] = -5;
	    break;
	}
    }

    *s = n;

    return 0;
}
