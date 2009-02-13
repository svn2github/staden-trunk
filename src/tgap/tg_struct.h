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

#define GT_RecArray     3
#define GT_Bin          5
#define GT_Range        6
#define GT_BTree	7
#define GT_Database    16
#define GT_Contig      17
#define GT_Seq         18
#define GT_DNASource   19
#define GT_Track       20

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
    GCardinal Nanno;    /* followed by Nanno record numbers */
} GBin;

typedef struct {
    GCardinal start;
    GCardinal end;
    GCardinal rec; /* recno */
    GCardinal pair_rec; /* paired end data */
    GCardinal flags; /* see below */
} GRange; /* An element of the bin->rng record */

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

/* Copied from gap4 - mostly unused at present */
typedef struct { 
    GCardinal version;

    /* Arrays */
    GCardinal Ncontigs;		/* N.elements in array */
    GCardinal contig_order;	/* rec. of array of contig rec. nos */

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
typedef struct {
    signed int  pos;  /* left end, regardless of direction */
    signed int len;   /* +ve or -ve indicates direction */
    int bin;
    int bin_index;    /* index to bin->rng array */
    int left, right;  /* clip left/right coordinates */
    int parent_rec;   /* template record */
    int parent_type;  /* GT_Template, GT_Ligation, etc */
    int rec;          /* recno of this seq_t */
    unsigned int seq_tech:3;
    unsigned int flags:3;
    unsigned int format:2;
    uint8_t mapping_qual;
    int name_len;
    int trace_name_len;
    int alignment_len;
    Array anno;       /* Annotations */
    char *name;       /* nul terminated name */
    char *trace_name; /* trace name; blank => same as name */
    char *alignment;  /* alignment; blank => obvious guess from pads */
    char *seq;        /* sequence in ASCII format */
    char *conf;       /* 1 or 4 values per base depending on flags */
    char *data;       /* packed memory struct; names/al/seq/conf are here */
} seq_t;

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
    int anno_rec;
    char *name;
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
    Array rng; /* NULL => not loaded */
    int rng_rec;
    int bin_id;
    int flags;
    Array track;    /* array of GTrack objects */
    int track_rec;
    Array anno;     /* array of GAnno objects */
} bin_index_t;

/* Bit flags for bin_index_t.flags */
#define BIN_COMPLEMENTED  (1<<0)
#define BIN_BIN_UPDATED   (1<<1)
#define BIN_RANGE_UPDATED (1<<2)
#define BIN_TRACK_UPDATED (1<<3)
#define BIN_ANNO_UPDATED  (1<<4)

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
    int comp; /* complemented y/n */
    int pair_rec;
    int pair_start;
    int pair_end;
    int flags;
    int y;    /* nominal display position, not stored on disc */
} rangec_t;

typedef struct {
    int start;
    int end;
    int rec;
    int pair_rec;
    int flags;
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

/* Track types */
#define TRACK_UNKNOWN    0
#define TRACK_ALL        0
#define TRACK_READ_DEPTH 1

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
typedef struct {
    char *comment; /* Possibly blank, but a per-region comment too */
    int rec;       /* record */
    int type;      /* type of object referred to by record */
    int start;
    int end;
} anno_ele_t;

typedef struct {
    char *key;        /* the tag type and text */
    char *value;
    int rec;          /* rec of this anno */
    int nele;
    anno_ele_t *ele;  /* ele = &single if nele == 1 */
    anno_ele_t single;
} anno_t;

#endif /* _TG_STRUCT_H_ */
