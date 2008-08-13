#include <stdio.h>
#include <string.h>
#include <assert.h>

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
 * Sets a sequence name.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_set_name(GapIO *io, seq_t **s, char *name) {
    seq_t *n;

    if (!(n = cache_rw(io, *s)))
	return -1;

    n = cache_item_resize(n, sizeof(*n) + strlen(name)+1 + 2*ABS(n->len));
    if (NULL == n)
	return -1;

    n->name = (char *)&n->data;
    n->seq = n->name + n->name_len + 1;
    n->conf = n->seq + ABS(n->len);

    /* Shift and insert name */
    /* FIXME: TO DO */

    *s = n;
    return 0;
}

int sequence_set_seq (GapIO *io, seq_t **s, char *seq) {return 0;}
int sequence_set_conf(GapIO *io, seq_t **s, char *conf) {return 0;}


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

    n = cache_item_resize(n, sizeof(*n) + 1+n->name_len
			  + 2*(ABS(n->len)+1));
    if (NULL == n)
	return -1;

    n->name = (char *)&n->data;
    n->seq = n->name + n->name_len + 1;
    n->conf = n->seq + ABS(n->len);

    p2 = n->len < 0
	? -n->len - pos
	: pos;
    alen = ABS(n->len);

    /* Shift */
    memmove(&n->seq[p2+1], &n->seq[p2], alen-p2+alen);
    n->conf++;
    memmove(&n->conf[p2+1], &n->conf[p2], alen-p2);

    if (n->len < 0)
	n->len--;
    else
	n->len++;

    /* Set */
    n->seq [p2] = base;
    n->conf[p2] = conf;

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
    memmove(&n->conf[p2], &n->conf[p2+1], alen-1-p2);

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
void complement_seq_conf(char *seq, char *conf, int seq_len) {
    int i, middle, j;
    char temp;
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
    for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
	temp = complementary_base [ (unsigned char) seq[i] ];
	seq[i] = complementary_base [ (unsigned char) seq[j] ];
	seq[j] = temp;
	temp = conf[i];
	conf[i] = conf[j];
	conf[j] = temp;
    }

    if ( seq_len % 2 )
      seq[middle] = complementary_base [ (unsigned char) seq[middle] ];
}

seq_t *dup_seq(seq_t *s) {
    size_t len = sizeof(seq_t) - sizeof(char *) + s->name_len+1 + 2*ABS(s->len);
    seq_t *d = (seq_t *)calloc(1, len);

    memcpy(d, s, len);
    d->name = (char *)&d->data;
    d->seq = d->name + d->name_len + 1;
    d->conf = d->seq + ABS(d->len);

    return d;
}

void complement_seq_t(seq_t *s) {
    int tmp, alen;

    alen = ABS(s->len);
    complement_seq_conf(s->seq, s->conf, alen);
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

/*
 * Operates on position 'pos' in the displayed orientation rather than the
 * stored orientation.
 */
int sequence_get_base(GapIO *io, seq_t **s, int pos, char *base, int *conf) {
    seq_t *n = *s;

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    if (n->len > 0) {
	*base = n->seq[pos];
	*conf = n->conf[pos];
    } else {
	*base = complementary_base[(int)(n->seq[-n->len-1 - pos])];
	*conf = n->conf[-n->len-1 - pos];
    }

    return 0;
}

int sequence_replace_base(GapIO *io, seq_t **s, int pos, char base, int conf) {
    seq_t *n;

    if (!(n = cache_rw(io, *s)))
	return -1;

    if (pos < 0 || pos >= ABS(n->len))
	return -1;

    if (n->len > 0) {
	n->seq[pos] = base;
	n->conf[pos] = conf;
    } else {
	n->seq[-n->len-1 - pos] = complementary_base[(unsigned char)base];
	n->conf[-n->len-1 - pos] = conf;
    }

    *s = n;

    return 0;
}
