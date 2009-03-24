/*
 * File: gap-dbstruct.h
 *
 * Author: Staden Package team
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: structure of gap database files
 *
 * Created: 28 October 1992
 * Updated:
 *
 */

#ifndef _GAP_DBSTRUCT_H_
#define _GAP_DBSTRUCT_H_

#include <sys/types.h>

#include "g-os.h"




#define GAP_DB_VERSION 3
/*
 * history
 * V0. Reset to year 0, 11 November 1993
 * V1. 26/04/94  Tags are now stored relative to the 'real sequence' start
 *               rather than the 'used sequence' start.
 * V2. 24/02/95  Primer type now holds values 0 to 4 (instead of 0-3). Version
 *               update is simply to allow for removing dead code in future
 *               if desired.
 * V3. 28/09/98  Notes are added. These are unpositional tags for attaching
 *		 information to entire readings, contigs, or the database.
 */

#define GAP_READ_LEN	        30000
#define GAP_FILES 1
#define DATABASE_FILE		0
#define GAP_DATABASE_FILE	0


/*
 * Predefined types
 */
#define GT_Unknown   	0
#define GT_Text   	1
#define GT_Data   	2
#define GT_Array  	3
#define GT_Bitmap 	4

/*
 * Other types
 */
#define GT_Database    	16
#define GT_Contigs	17
#define GT_Readings	18
#define GT_Vectors	19
#define GT_Annotations	20
#define GT_Templates	21
#define GT_Clones	22
#define GT_Notes	23



/*
 * Fixed record numbers
 */
#define GR_Database		0
#define GR_Freerecs		1
#define GR_Contigs		2
#define GR_Readings		3
#define GR_Annotations		4
#define GR_Vectors		5
#define GR_Templates	        6
#define GR_Clones		7
#define GR_Vectors_Default	8
#define GR_Templates_Default	9
#define GR_Clones_Default	10
#define GR_Unknown		11
#define GR_Contig_Order		12


/*
 * Macros
 */
/* GDatabase.data_class */
#define GAP_DNA            0
#define GAP_PROTEIN        1
/* GReadings.sense */
#define GAP_SENSE_ORIGINAL 0
#define GAP_SENSE_REVERSE  1
/* GReadings.strand */
#define GAP_STRAND_FORWARD 0
#define GAP_STRAND_REVERSE 1
/* GReadings.primer */
#define GAP_PRIMER_UNKNOWN 0
#define GAP_PRIMER_FORWARD 1
#define GAP_PRIMER_REVERSE 2
#define GAP_PRIMER_CUSTFOR 3
#define GAP_PRIMER_CUSTREV 4
/* GReadings.chemistry */
/*	Bit 0 is 1 for terminator, 0 for primer */
#define GAP_CHEM_TERMINATOR	(1<<0)
/*	Bits 1 to 4 inclusive are the type (any one of, not bit pattern) */
#define GAP_CHEM_TYPE_MASK	(15<<1)
#define GAP_CHEM_TYPE_UNKNOWN	(0<<1)
#define GAP_CHEM_TYPE_ABI_RHOD	(1<<1)
#define GAP_CHEM_TYPE_ABI_DRHOD	(2<<1)
#define GAP_CHEM_TYPE_BIGDYE2	(3<<1)
#define GAP_CHEM_TYPE_ET	(4<<1)
#define GAP_CHEM_TYPE_LICOR	(5<<1)
#define GAP_CHEM_TYPE_MB_ET	(6<<1)
#define GAP_CHEM_TYPE_BIGDYE1	(7<<1)
#define GAP_CHEM_TYPE_BIGDYE3	(8<<1)
#define GAP_CHEM_TYPE_SOLEXA	(9<<1)
#define GAP_CHEM_TYPE_SOLID	(10<<1)
#define GAP_CHEM_TYPE_454	(11<<1)
#define GAP_CHEM_TYPE_OX_NANO	(12<<1)
/* GVectors.level */
#define GAP_LEVEL_UNKNOWN  0
#define GAP_LEVEL_CLONE    1
#define GAP_LEVEL_SUBCLONE 2

/*
0	unknown		primer
1	unknown		terminator
2	ABI rhodamine	primer
3	ABI rhodamine 	terminator
4	ABI dRhodamine 	primer
5	ABI dRhodamine 	terminator
6	ABI BigDye v2 	primer
7	ABI BigDye v2 	terminator
8	Energy transfer	primer
9	Energy transfer	terminator
10	Licor 		primer
11	Licor 		terminator
12	MegaBACE ET 	primer
13	MegaBACE ET 	terminator
14	ABI BigDye v1 	primer
15	ABI BigDye v1 	terminator
16	ABI BigDye v3 	primer
17	ABI BigDye v3 	terminator
18	-
19	Solexa		terminator
20	AB SOLiD	N/A
21	AB SOLiD	N/A <- use this one
22	454		N/A
23	454		N/A <- use this one
24	Oxford Nanopore	N/A
25	Oxford Nanopore	N/A <- use this one
*/


/*
 * Annotations are singly linked lists, terminated with next == 0
 */
typedef struct { 
    GCardinal type;
    GCardinal position; 
    GCardinal length; 
    GCardinal strand; 
    GCardinal annotation; 
    GCardinal next;
} GAnnotations; 

/*
 * Notes are doubly linked lists. The end terminates with next == 0.
 * The start terminates when prev_type != GT_Notes In that case prev is the
 * contig or reading number (or 0 if prev_type == GT_Database) that this
 * note list comes from.
 */
typedef struct {
    GCardinal type;
    GCardinal ctime_top;	/* Placeholder for 64bit time_t */
    GCardinal ctime;		/* creation date - a time_t. */
    GCardinal mtime_top;	/* Placeholder for 64bit time_t */
    GCardinal mtime;		/* modification date - a time_t. */
    GCardinal annotation;
    GCardinal next;		/* 2 way linked list */
    GCardinal prev;
    GCardinal prev_type;	/* GT_{Notes,Contigs,Readings,Database} */
} GNotes;

typedef struct {
    GCardinal name;		/* vector name */
    GCardinal level;		/* 1=clone, 2=subclone, etc */
} GVectors; 






typedef struct {
    GCardinal name;
    GCardinal trace_name;
    GCardinal trace_type;
    GCardinal left;		/* left neighbour */
    GCardinal right;		/* right neighbour */
    GCardinal position;		/* position in contig */
    GCardinal length;		/* total length of reading */
    GCardinal sense;		/* 0 = original, 1 = reverse */
    GCardinal sequence;
    GCardinal confidence;
    GCardinal orig_positions;
    GCardinal chemistry;	/* see comments above (GAP_CHEM_*) */
    GCardinal annotations;	/* start of annotation list */
    GCardinal sequence_length;	/* clipped length */
    GCardinal start;		/* last base of left cutoff */
    GCardinal end;		/* first base of right cutoff */
#ifdef __cplusplus
    /* template is reserved in c++ */
    GCardinal tmplate;		/* aka subclone */
#else
    GCardinal template;		/* aka subclone */
#endif
    GCardinal strand;		/* 0 = forward, 1 = reverse */
    GCardinal primer;		/* 0 = unknown, 1 = forwards, 2 = reverse, 3 = custom forward, 4 = custom reverse*/
    GCardinal notes;		/* Unpositional annotations */
} GReadings; 







typedef struct { 
    GCardinal left; 
    GCardinal right; 
    GCardinal length;
    GCardinal annotations;	/* start of annotation list */
    GCardinal notes;		/* Unpositional annotations */
} GContigs; 



typedef struct { 

    GCardinal version;
    GCardinal maximum_db_size;	/* MAXDB */
    GCardinal actual_db_size;	/* */
    GCardinal max_gel_len;	/* 4096 */
    GCardinal data_class;	/* 0 = DNA, 1 = protein */

    GCardinal num_contigs;	/* number of contigs used */
    GCardinal num_readings;	/* number of readings used */

    /* Bitmaps */
    GCardinal Nfreerecs;	/* number of words (currently 32bits/word) */
    GCardinal freerecs;

    /* Arrays */
    GCardinal Ncontigs;		/* elements in array */
    GCardinal contigs;		/* records that are GT_Contigs */

    GCardinal Nreadings;	/* elements in array */
    GCardinal readings;		/* records that are GT_Readings */

    GCardinal Nannotations;	/* elements in array */
    GCardinal annotations;	/* records that are GT_Annotations */
    GCardinal free_annotations;	/* head of list of free annotations */

    GCardinal Ntemplates;	/* elements in array */
    GCardinal templates;	/* records that are GT_Templates */

    GCardinal Nclones;		/* elements in array */
    GCardinal clones;		/* records that are GT_Templates */

    GCardinal Nvectors;		/* elements in array */
    GCardinal vectors;		/* records that are  GT_Vectors */

    GCardinal contig_order;	/* Array record number */

    GCardinal Nnotes;		/* elements in array */
    GCardinal notes_a;		/* records that are GT_Notes */
    GCardinal notes;		/* Unpositional annotations */
    GCardinal free_notes;	/* SINGLY linked list of free notes */
} GDatabase; 






typedef struct {
    GCardinal name;
    GCardinal strands; /* Number of strands, 1 or 2 */
    GCardinal vector;
    GCardinal clone;
    GCardinal insert_length_min;
    GCardinal insert_length_max;
} GTemplates;



typedef struct {
    GCardinal name;
    GCardinal vector;
} GClones;


extern char *file_list[];

extern size_t block_sizes[];

extern GCardinal max_recs[];





extern char *gap_construct_file(char *database,char *file, char *version, char *fillbuf);



extern void dumpGContigs(GContigs *r);
extern void dumpGReadings(GReadings *r);
extern void dumpGDatabase(GDatabase *r);


#endif /*_GAP_DBSTRUCT_H_*/
