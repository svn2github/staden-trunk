#ifndef _TG_STRUCT_H_
#define _TG_STRUCT_H_

#include <g.h> /* needed only for some typedefs */
#include <array.h>
#include <inttypes.h>

/* FIXME: fixed sized buffers are so last decade! */
/* Suitable for solexa data, eg "slxa_0025_7_0021_3265" */
#define MAX_NAME_LEN 25
#define MAX_SEQ_LEN 100000

/* ----------------------------------------------------------------------
 * Structures in a more database 'on disk' representation. This is basically
 * as it's stored.
 */
 
#define GT_Generic       0
#define GT_Generic       0
#define GT_RecArray      3
#define GT_Bin           5
#define GT_Range         6
#define GT_BTree         7
#define GT_Database     16
#define GT_Contig       17
#define GT_Seq          18
#define GT_Library      19
#define GT_Track        20
#define GT_AnnoEle      21
#define GT_Anno         22
#define GT_SeqBlock     23
#define GT_AnnoEleBlock 24
#define GT_SeqCons      25

typedef struct {
    GCardinal pos;
    GCardinal size;
    GCardinal start;
    GCardinal end;
    GCardinal parent_type;   /* GT_Bin or GT_Contig */ // FIXME: to add
    GCardinal parent;   /* recno */                    // FIXME: to add
    GCardinal child[2]; /* recno */
    GCardinal range;    /* recno */
    GCardinal id;
    GCardinal flags;
    GCardinal track;
    GCardinal nseqs;
    GCardinal rng_free; /* forms a linked list of items in rng that are free */
} GBin;

typedef struct {
    GCardinal start;
    GCardinal end;
    GCardinal rec; /* recno */
    GCardinal mqual; /* mapping quality */
    GCardinal pair_rec; /* paired end data */
    GCardinal flags; /* see below */
    /* Move library here? */
} GRange; /* An element of the bin->rng record */


#define MQUAL_UNKNOWN   255

#define GRANGE_FLAG_ISMASK     (7<<7) /* Sequence, tag, cons, ref, etc */
#define GRANGE_FLAG_ISSEQ      (0<<7)
#define GRANGE_FLAG_ISANNO     (1<<7)
#define GRANGE_FLAG_ISCONS     (2<<7)
#define GRANGE_FLAG_ISREF      (3<<7)
#define GRANGE_FLAG_ISUMSEQ    (4<<7) /* unmapped sequence */
#define GRANGE_FLAG_ISANY      (7<<7) /* Any */

/* For annotation ranges: */
#define GRANGE_FLAG_TAG_SEQ    (1<<1) /* 0=>contig, 1=>sequence */
#define GRANGE_FLAG_COMPOUND   (1<<2) /* true anno. has multiple components */

/* For sequence ranges: */
#define GRANGE_FLAG_TYPE_MASK  (3<<0)
#  define GRANGE_FLAG_TYPE_SINGLE  0
#  define GRANGE_FLAG_TYPE_PAIRED  1
#  define GRANGE_FLAG_TYPE_COMPLEX 2  /* > 2 or 2x forward, etc */
#define GRANGE_FLAG_END_MASK   (1<<2) /* only applicable if not TYPE_COMPLEX */
#  define GRANGE_FLAG_END_FWD  (0<<2)
#  define GRANGE_FLAG_END_REV  (1<<2)
#define GRANGE_FLAG_PEND_MASK  (1<<6) /* as _END_MASK, but pair data */
#  define GRANGE_FLAG_PEND_FWD (0<<6)
#  define GRANGE_FLAG_PEND_REV (1<<6)
#define GRANGE_FLAG_CONTIG     (1<<3) /* pair held within the same contig */
#define GRANGE_FLAG_COMP1      (1<<4) /* true if complemented */
#define GRANGE_FLAG_COMP2      (1<<5) /* true if complemented */

#define GRANGE_FLAG_UNUSED     (1<<10) /* range has been deleted */

typedef struct {
    GCardinal type;
    GCardinal flags;
    GCardinal rec;
} GBinTrack; /* An element of the bin->track record */

typedef struct {
    GCardinal type;  /* duplicated in bin */
    GCardinal flags; /* duplicated in bin */
    GCardinal item_size;
    GCardinal nitems;
    /* + item_size * nitems bytes of actual track data */
} GTrack_Header; /* The track itself referenced by bin->track[]->rec */ 

typedef struct {
    GCardinal start;
    GCardinal end;
    GCardinal bin;
    /* + 1 byte  name length 'L' */
    /* + L bytes name            */
} GContig_header;


/*
 * A global version number used to indicate the maximum version number of
 * any contents in this database.
 *
 * This can be used to indicate to a gap5 version whether it is capable
 * of completely handling this database or not. Possibly we could also use
 * it to force a newer gap5 release to keep writing data in an older
 * backwards compatible manner when editing an old DB.
 */
#define DB_VERSION 1 /* 1.2.6 */

/* Copied from gap4 - mostly unused at present */
typedef struct { 
    GCardinal version;

    /* Arrays */
    GCardinal Ncontigs;		/* N.elements in array */
    GCardinal contig_order;	/* rec. of array of contig rec. nos */

    GCardinal Nlibraries;       /* N.elements in array */
    GCardinal library;          /* rec. of array o library rec. nos */

    /* Indices */
    GCardinal seq_name_index;	/* rec of type GT_Index */
    GCardinal contig_name_index;/* rec of type GT_Index */
} GDatabase; 

/* ----------------------------------------------------------------------
 * In-memory versions of the above structures where we deem it more
 * efficient to encode in a different format.
 *
 * (Also see bin_index_t from binning.h)
 *
 * Note that this object is one single block of memory *at least* as large
 * as sizeof(seq_t). All the pointers in here (names, seq, conf) except for
 * data are tacked onto the end of the struct. Quite literaly speaking
 * the name starts at the &seq.data. The packing order is described as
 * below. Also see seq_decode() in tg_iface_g.c.
 *
 * field	len
 * ---------------------------------------------
 * name         name_len    (at address &seq.data)
 * nul          1
 * trace_name   trace_name_len
 * nul		1
 * alignment    alignmen_len
 * nul          1
 * seq		ABS(len)
 * conf		ABS(len)    (iff format != SEQ_FORMAT_CNF4)
 * conf		4*ABS(len)  (iff format == SEQ_FORMAT_CNF4, in order ACGT,ACGT)
 */
struct seq_block;
typedef struct {
    signed int  pos;  /* left end, regardless of direction */
    signed int len;   /* +ve or -ve indicates direction */
    int bin;
    int bin_index;    /* index to bin->rng array */
    int left, right;  /* clip left/right coordinates */
    int parent_rec;   /* template record or seq record if type == GT_Seq */
    int parent_type;  /* GT_Seq, GT_Template, GT_Ligation, etc */
    int rec;          /* recno of this seq_t */
    unsigned int seq_tech:3;
    unsigned int flags:3;
    unsigned int format:2;
    uint8_t mapping_qual; /* REMOVE? In GRange already. Same for parent_rec */
    int name_len;
    int trace_name_len;
    int alignment_len;
    int aux_len;
    Array anno;       /* Annotations; FIXME */
    char *name;       /* nul terminated name */
    char *trace_name; /* trace name; blank => same as name */
    char *alignment;  /* alignment; blank => obvious guess from pads */
    char *seq;        /* sequence in ASCII format */
    char *conf;       /* 1 or 4 values per base depending on flags */
    char *sam_aux;    /* Auxillary records */
    struct seq_block *block; /* seq_block_t pointer and index into it */
    int idx;

    char *data;       /* packed memory struct; names/al/seq/conf are here */
} seq_t;

/* Maximum size of a block, actual size maybe less if long sequences */
#define SEQ_BLOCK_BITS 10
#define SEQ_BLOCK_SZ (1<<SEQ_BLOCK_BITS)
typedef struct seq_block {
    int    est_size;
    int    rec[SEQ_BLOCK_SZ];
    seq_t *seq[SEQ_BLOCK_SZ];
} seq_block_t;


/* Sequencing technologies for seq_t.seq_tech */
#define STECH_UNKNOWN 0
#define STECH_SANGER  1
#define STECH_SOLEXA  2
#define STECH_SOLID   3
#define STECH_454     4

/* Sequence flags for seq_t.flags */
#define SEQ_COMPLEMENTED (1<<0)
#define SEQ_CONF_PHRED   (1<<1) /* Confidence values in phred-scale?
				   False => log-odds */
#define SEQ_END_MASK     (1<<2)
#define SEQ_END_FWD      (0<<2)
#define SEQ_END_REV      (1<<2)

#define SEQ_FORMAT_MAQ   0      /* 2-bit base, 6-bit conf */
#define SEQ_FORMAT_CNF1  1      /* 8-bit base, 1 confidence values */
#define SEQ_FORMAT_CNF4  2      /* 8-bit base, 4 confidence values */
#define SEQ_FORMAT_CNF6  3      /* 8-bit base, 6 conf (ins/del too) */

typedef struct {
    int rec;
    signed int start, end;
    unsigned int bin;
    int anno_rec; /* FIXME */
    char *name;
    char *data;
} contig_t;


typedef struct index {
    int rec;
    int pos;
    int size;
    int start_used;
    int end_used;
    int parent_type;
    int parent;   /* -1 implies none */
    int child[2]; /* -ve implies saved already with record -child[?] */
    Array rng;    /* NULL => not loaded */
    int rng_rec;
    int bin_id;
    int flags;
    Array track;  /* array of bin_track_t objects */
    int track_rec;
    int nseqs;
    int rng_free; /* forms a linked list of items in rng that are free */
} bin_index_t;

/* Bit flags for bin_index_t.flags */
#define BIN_COMPLEMENTED  (1<<0)
#define BIN_BIN_UPDATED   (1<<1)
#define BIN_RANGE_UPDATED (1<<2)
#define BIN_TRACK_UPDATED (1<<3)
#define BIN_CONS_CACHED   (1<<4) /* a cached consensus is stored here */
#define BIN_CONS_VALID    (1<<5) /* ... and is up to date */

/*
 * We may also wish to hold in range:
 *     other_end (rec)
 *     template (complex case).
 *
 * Eg consider:
 *        
 *           | VISIBLE PORTION VISIBLE PORTION |    
 * A   >>>>--|---------------------------------|--<<<<
 * B   >>>>--|-------------------<<<<          |
 * C         | >>>>----------------------<<<<  |
 * D         |         >>>>--------------------|--<<<<
 *
 * Template A is invisible as it has no end visible within our range query.
 * We can address this partially by looking a certain distance either
 * side of the visible portion, but most cases of type A will be
 * chimeras and very far apart.
 *
 * Template B and D have one end only visible. We may still need to load the
 * sequence struct to figure out where it is.
 *
 * Template C, hopefully the dominant case, will have both ends returned as
 * range structs from a display query. We therefore can reindex the array
 * in memory by record number and identify both ends and the extents of the
 * line to draw. Hence we do not need to fetch the sequence structs directly.
 *
 * The more complex cases involving more than two reads requires a template
 * record and proper two-way linkage. This is only needed for capillary
 * style sequencing.
 */
typedef struct {
    int start;
    int end;
    int rec;
    int mqual; /* Mapping qual */
    int comp;  /* complemented y/n */
    int pair_rec;
    int pair_start;
    int pair_end;
    int pair_mqual;
    int flags;
    int y;     /* nominal display position, not stored on disc */
    int pair_ind; /* -1 if not found, or index into array of rangec_t */
    /* Move library here? */
} rangec_t;

typedef struct {
    int start;
    int end;
    int rec; /* or alternatively an index if range_t is free */
#if 0
    int type;
    union {
	struct {
	    int mqual;
	    int pair_rec;
	    int flags;
	} seq;
	struct {
	    int obj_type;
	    int obj_rec;
	    int anno_rec;
	} anno;
	struct {
	    int a, b, c;
	} generic;
    } u;
#else
    int mqual;
    int pair_rec;
    int flags;
    /* Move library here? */
#endif
} range_t;

/* Decoded from GTrack_header above */
typedef struct {
    int type;      /* GC %, read-depth, etc. See TRACK_* macros */
    int flag;      /* bit 0 => updated or not, 1 => can free */
    int rec;       /* record number */
    int bin_size;  /* size of corresponding bin, in bases */
    int item_size; /* size of each element in data[] */
    int nitems;    /* Number of items in data[] */
    Array data;    /* cached copy of data for this track record */
} track_t;

typedef struct {
    int type;
    int flags;
    int rec;
    track_t *track;
} bin_track_t;

/* Track types */
#define TRACK_UNKNOWN    0
#define TRACK_ALL        0
#define TRACK_READ_DEPTH 1
#define TRACK_CONS_ARR   2

/* Track flag masks */
#define TRACK_FLAG_VALID  (1<<0)
#define TRACK_FLAG_FREEME (1<<1)

/*
 * Annotations:
 * These now have a two-way relationship. A sequence refers to an annotation
 * struct. Likewise an annotation record consists of multiple fragments
 * possibly spanning multiple reads. This means we can annotate pairs of reads
 * that may link in some way with a single shared annotation, or disjoint 
 * regions within a single sequence (eg gene structures).
 */
struct anno_ele_block;
typedef struct {
    int tag_type;  /* short 4-byte tag type code */
    char *comment; /* Possibly blank, but a per-region comment too */
    int rec;       /* record */
    int bin;       /* bin containing this element */
    int obj_type;  /* type of object referred to by record */
    int obj_rec;   /* record number of object the tag is attached to */
    int anno_rec;  /* link to anno_t record below */
    struct anno_ele_block *block;
    int idx; 
    char *data;    /* location of packed comment */
} anno_ele_t;

#define ANNO_ELE_BLOCK_BITS 10
#define ANNO_ELE_BLOCK_SZ (1<<ANNO_ELE_BLOCK_BITS)
typedef struct anno_ele_block {
    int        est_size;
    int        rec[ANNO_ELE_BLOCK_SZ];
    anno_ele_t *ae[ANNO_ELE_BLOCK_SZ];
} anno_ele_block_t;

typedef struct {
    char *key;        /* the tag type and text */
    char *value;
    int rec;          /* rec of this anno */
    Array ele;        /* recs of anno_ele_t */
} anno_t;

/*
 * Libraries:
 * These hold data on how a set of sequences were produced, describing the
 * processes involved in producing our single or paired-end reads.
 *
 * It's not a strict physical library description as it's possible to take
 * one library and run it on multiple sequencing technologies, but rather
 * we describe a single run on a single machine as our level of granularity.
 *
 * The orientation and direction of our forward and reverse reads will
 * vary by library construction type. Standard puc19 capillary libraries
 * have fwd/rev reads in opposite orientations pointing towards one another,
 * similarly for solexa short insert libraries.
 * Long insert libraries on solexa have the fwd and rev reads pointing away
 * from one another at long insert size (with a smaller sub-set pointing
 * towards one another with small insert size).
 * While 454 long insert libraries have the fwd/rev reading pointing in the
 * same orientation. (Reverse first, then forward?)
 *
 * Hence we describe here our expected orientation, for purposes of working
 * out if a pair for a library is consistent or not.
 * rev_dir and rev_orient are expressed as relative to the fwd_dir and
 * fwd_orient.
 *
 * We measure insert sizes from 0   to +255 in steps of 1 (256 values)
 *                              256 to  511 in steps of 2 (128 values)
 *                              512 to 1024 in steps of 4 (128 valyes)
 * etc, up to a maximum value.
 */
#define LIB_BINS 1792

#define LIB_T_INWARD  0 /* Reads point towards one another */
#define LIB_T_OUTWARD 1 /* Reads point outwards from one another */
#define LIB_T_SAME    2 /* Reads are in the same orientation */

typedef struct {
    int rec;             /* DB record */
    int insert_size[3];  /* Mean insert size */
    double sd[3];        /* standard deviation of insert size */
    int machine;         /* Type of machine, see STECH_* defines above */
    int lib_type;        /* Primary LIB_T_ type expected */
    
    /* A distribution summary, in 1s initially, and then 2s, 4s, 8s, etc */
    int size_hist[3][LIB_BINS+1];

    char *name;
    char *data;
} library_t;

typedef struct {
	double yz;   /* t->yzoom / 200.0 */
	int ymin;
	int ymax;
	double wy0;
	double wy1;
} drawing_t;


#endif /* _TG_STRUCT_H_ */
