/*
 * File: edStructs.h
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: Contig editor structures
 *
 * Created:
 * Updated:
 *
 */

#ifndef _EDSTRUCTS_H_
#define _EDSTRUCTS_H_

#include "xalloc.h"

#include "os.h"
#include "fort.h"
#include "IO.h"
#include "tkEditor.h"
#include "tkEdNames.h"
#include "io-reg.h"
#include "primlib.h"
#include "template.h"

/* New (crude, but simple) consensus caching system */
#define CACHE_CONSENSUS_2

/*
 * constant definitions
 */
#define COMMENT_LENGTH 40
#define MAX_DISPLAY_WIDTH 300
#define DEFAULT_DISPLAY_WIDTH 80
#define DB_GELNOLEN 7
#define NAMELEN (DB_NAMELEN + DB_GELNOLEN + 1)
#define MAXEDSTATES 100
#define BASES 6
#define SEQ_LENGTH_INC 20 /* seq. length inc. when reallocing */

/* Max number of displays of this data */
#define MAX_DISP_PROCS 10

/*
 * We implement a circular stack to store undo commands.
 * 
 * MAX_SAVE_EDITS defines the maximum number of recallable commands
 */
#define MAX_SAVE_EDITS 100

/*
 * Maximum number of different status lines at once, and their types.
 */
#define MAX_STATUS_LINES 8

/* These values are indexes into the xx->status[] array */
#define EDITOR_SL_STRAND	0
#define EDITOR_SL_FRAME1p	1
#define EDITOR_SL_FRAME2p	2
#define EDITOR_SL_FRAME3p	3
#define EDITOR_SL_FRAME1c	4
#define EDITOR_SL_FRAME2c	5
#define EDITOR_SL_FRAME3c	6
#define EDITOR_SL_AUTOTRANSLATE	7 /* FIXME */

/* These are the ways to sort the sequences in the editor */
typedef enum {
    POSITION = 0,
    TEMPLATE = 1,
    STRAND   = 2,
    CLONE    = 3,
    ALPHA    = 4,
    NUMERIC  = 5
} editor_sort_t;

typedef f_int tag_id;
typedef f_int comment_id;

struct _EdLink;

enum States {StateDown=0,StateUp};



/*
 ** element in tag list
 */
typedef union {
    char c[4];
    f_int i;
} tag_types;




typedef struct _tagRecord{
    f_int position;		/* position in sequence */
    f_int length;		/* length of tag */
    tag_types type;
    comment_id comment;		/* index to comment */
    tag_id next;		/* link to next in structure */
    f_int sense;		/* sense: 0=original fwd, 1=original rev */
} tagRecord;




typedef struct _tagstruct{
    /*
     ** Data from the original file
     */
    tagRecord tagrec;
    /*
     ** Data for database management
     */
    tag_id original_tag_id;
    char *newcomment;
    int newcommentlen;
    long flags;
    struct _tagstruct *next;
} tagStruct,*tagptr;



/*
 ** type definitions
 */
typedef struct {
    int relPos;			/* relative position in gel */
    int length;			/* length of editable region of sequence */
    int number;			/* record number of sequence in file */
    int complemented;		/* sequence not in original sense */
    char *name;			/* reading name */
    char *sequence;		/* COMPLETE sequence */
    int flags;
    tagStruct *tagList;		/* list of tags for sequence */
    int1 *gap_conf;		/* confidence values */
    int2 *gap_opos;		/* original position */
    int  alloc_length;		/* length allocated for sequence details */
    int  gap_length;		/* length of complete string */
    int  gap_start;		/* start */
    int  gap_end;		/* end */
    int  template;		/* template number */
} DBStruct;


/*
 * Temporary results used by the "select oligos" function
 */
typedef struct {
    int l,r;			/* region of consensus to analyse */
    char *consensus;		/* the unpadded consensus at this region */
    int *opos;			/* original padded consensus positions */
    primlib_state *state;	/* The primer 'state'. Contains results */
    int *opos_start;		/* Pad-adjusted start position */
    int *opos_end;		/* Pad-adjusted end position */
    int curr_oligo;		/* which results[] element we're using */
    int sense;			/* which strand (0 = top, 1 = bottom) */
    int ave_read_len;		/* average read length */
} select_oligo_t;


typedef struct {
    GapIO *io;

    DBStruct *DB;
    int DB_flags;
    int DB_gelCount;		/* number of gels in contig */
    int DB_contigNum;		/* contig line number */
    int *DBlist;		/* static list on DB_gelCount integers */
    int *DBorder;		/* order of sequences in contig */

    /* Lists of display functions and their 'client data' */
    void (*display_func[MAX_DISP_PROCS])(void *xx, int type,
					 int seq, int pos, void *pointer);
    void *display_data[MAX_DISP_PROCS];
    int next_display;

    /* Undo */
    struct undoNode *undo_lists[MAX_SAVE_EDITS];
    int last_undo;
    int num_undo;
    int discarded_undos;
    int edits_made;
    int since_auto_save;
    int store_undo;		/* 1 => store undo info */
    signed int registration_id;
    int open_undo_count;

    /* Reference sequence */
    int reference_seq;		/* 0 for none. Otherwise a seq id */
    int reference_len;		/* length of ref seq (if cyclic, 0 if not) */
    int reference_offset;	/* the base number for the 1st base of ref */
    
    /* Template information */
    template_c **templates;
} DBInfo;

typedef struct {
    char line[MAX_DISPLAY_WIDTH+1];
    char name[NAMELEN+1];
    XawSheetInk colours[MAX_DISPLAY_WIDTH+1];
} EdStatus;

/* REMEMBER TO: update structure initialisation in edUtils.c */
typedef struct _EdStruct {
    DBInfo *DBi;

    int displayPos ;		/* position of top left corner of display */
    int displayYPos;		/* Y position of display */
    int displayWidth;		/* current screen width */
    int displayHeight;		/* current display height */
    int totalHeight;		/* total height of sequences */
    int cursorPos;		/* cursor position */
    int cursorSeq;		/* sequence on which cursor lies */
    int rulerDisplayed;		/* ruler display switch */
    int consensusDisplayed;	/* consensus display switch */
    int fontWidth;		/* width of contig editor font */
    int fontHeight;		/* height of contig editor font */
    int extent_left;	 	/* Start of end base nums of the displayed */
    int extent_right;		/* contig, including cutoffs if visible */

    int editor_id;		/* A unique id for each editor */

    Editor *ed;			/* Sequence display widget */
    edNames *names;		/* The sequence names sheet widget */

    char displayedConsensus[MAX_DISPLAY_WIDTH+1]; /* consensus */
    float displayedConfidence[MAX_DISPLAY_WIDTH+1]; /* consensus */

    int select_made;		/* a selection has been made */
    int select_seq;		/* sequence selected */
    int select_start_pos;	/* start position in sequence */
    int select_end_pos;		/* end position in sequence */
    tagStruct *select_tag;	/* tag selected */

    int reveal_cutoffs;		/* reveal cutoff switch */
    int showDifferences;	/* highlight disagreements switch */
    int compare_strands;	/* strand consensus comparisons */
    int auto_save;		/* whether or not to Save Contig regularly */
    int insert_mode;		/* Whether we're in insert mode */
    int super_edit;		/* Whether we're allowed super edit */

    int qual_bg[10];		/* 10 'Grey' scale colours for quality */
    int qual_below;		/* Colour for bases below quality cutoff */
    float con_cut;		/* Current consensus cutoff */
    int qual_cut;		/* Current quality cutoff */

    struct _EdLink *link;	/* A link to another EdStruct for joining */
    enum States editorState;
    int editorMode;		/* Could use xx->link != NULL instead */

    select_oligo_t *sel_oli;	/* Pointer to saved data for the
				 * "select oligos" function. NULL normally.
				 */

    int *tag_list;		/* Displayed tag types */
    int display_traces;		/* Whether to display traces for Next.Prob */
    int status[MAX_STATUS_LINES]; /* Which status lines to display */
    EdStatus *status_lines;	/* The status lines themselves */
    int status_depth;		/* Number of elements in status_lines */
    int trans_mode;		/* Codon translation mode, set to 1 or 3 */
    int trace_lock;		/* Whether traces are locked to the editor */
    int show_qual;		/* How to display quality information */
    int show_cons_qual;		/* How to display consensus qual information */

    tagStruct *tmp_tag;		/* Temporary tag used in Select Primers */

    int refresh_flags;		/* Optimisation for redisplaying */
    signed int refresh_seq;	/* Optimisation for redisplaying */
    int refresh_pos;		/* Optimisation for redisplaying */
    int diff_traces;		/* Whether to auto-difference traces */
    signed int compare_trace;	/* Seq number to diff traces with */
    int compare_trace_match;	/* Con trace: use only matching traces? */
    int compare_trace_select;	/* Con trace: ignore selected trace? */
    int diff_bg;		/* Colour for show_differences colour mode */
    int diff_qual;		/* Min quality for differences to be shown */
    int compare_trace_algorithm;/* Trace diff: whether to use ABS() */
    int compare_trace_yscale;	/* Trace diff: do we scale in Y? */
    editor_sort_t group_mode;	/* Sorting for seqs: template, name, strand..*/
    int show_edits;		/* Do we highlight edits? */
    int edit_bg[4];		/* Colours for edit types */
    int tmpl_bg[6];		/* Colours for template statuses */
    cursor_t *cursor;		/* Cursor structure */
    int names_xpos;		/* X position of names display */
    int default_conf_r;		/* Confidence for replace bases, -1 for old */
    int default_conf_n;		/* Confidence for new bases, -1 for averaged */

#ifdef CACHE_CONSENSUS_2
    int valid_consensus;	/* Valid consensus flag */
    int consensus_len;		/* Length of cached arrays */
    char *consensus;		/* Cached consensus */
    float *quality;		/* Cached quality */
#endif

    int unpadded_ruler;		/* Ruler shows unpadded positions? 1 if so */

    int consensus_mode;		/* Current consensus algorithm */
    int lines_per_seq;		/* (Height-1) of mini-traces */
    int diff_trace_size;	/* 0 == auto, N == from cursor +- N bases */
} EdStruct;




typedef struct _EdLink {
    EdStruct *xx[2];		/* A pair of EdStructs */
    int locked;			/* Set if these are locked together */
    int lockOffset;		/* And their relative distances */
    tkSheet *diffs;
} EdLink;

/*
 * FIXME - for now it's easiest doing this in the editor - but we need to
 * move this outside at some stage.
 */
typedef struct _aa_pair_t {
    char AA1; /* A pair of amino acid bases */
    char AA2; /* '\0' if only one element is needed. */
} aa_pair_t;

enum mutation_type {
    no_mutation,
    non_coding,
    silent,
    expressed
};

typedef struct _mutation_t {
    char *tag_text_top;
    char *tag_text_bot;
    enum mutation_type type;
    int seq_top;
    int seq_bot;
    int strands; /* bit 0 == top, bit 1 == bottom */
    int conflict; /* true or false */
    aa_pair_t AA_to;
    char nucleotide_from;
    char nucleotide_to;
    char AA_from;
    char tag_type_top[4];
    char tag_type_bot[4];
} mutation_t;

/* Range of locations mutscan checked; as positions on the reference seq */
typedef struct {
    int ref_start;
    int ref_end;
    int sibling;
    int fwd;
} mutation_cov_t;

#include "undo.h"


#endif /*_EDSTRUCTS_H_*/

