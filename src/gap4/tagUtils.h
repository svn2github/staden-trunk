/*
 * File: tagUtils.h
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */
#ifndef _tagUtils_h
#define _tagUtils_h

#include "edStructs.h"
#include "gap-dbstruct.h"
#include "IO.h"
#include "tagDefs.h"
#include "tkSheet.h"
#include "xalloc.h"

#include "fortran.h"
#include "fort.h"
#include "array.h"

struct anno_list_t {
    GCardinal anno;
    int type;
    int position;
    int length;
    int strand;
};

typedef struct {
    int number;
    int prev;
    int next;
    int from;
    int from_type;
} anno_ptr;

/*
** The following describe two database files:
**     The tag list file
**     The comment list file
**
** the tag list file consists of at least IDBSIZ records.
**
**     1       : header tag for sequence 1
**     2       : header tag for sequence 2
**     ...
**     IDBSIZ-1: header tag for sequence IDBSIZ-1:
**     IDBSIZ  : descriptor record defining MAX_TAG
**     IDBSIZ+1: supplimentary tags
**     ...
**     MAX_TAG : supplimentary tags
**
** the comment list file consists of at least 1 record.
**
**     1       : descriptor record defining MAX_COM
**     2       : supplimentary comment
**     ...
**     MAX_COM : supplimentary comment
**     
*/

#define TAG_MALLOC(s) (char *)xmalloc(s)
#define TAG_FREE(c)   xfree(c)

/*
**
*/


/* useful macros */

#define str2type(s) ((s)[3] + ((s)[2]<<8) + ((s)[1]<<16) + ((s)[0]<<24))

#define type2str(t,s) \
    ( \
     ((s)[0] = ((t) >> 24) & 0xff), \
     ((s)[1] = ((t) >> 16) & 0xff), \
     ((s)[2] = ((t) >>  8) & 0xff), \
     ((s)[3] = ((t) >>  0) & 0xff), \
     ((s)[4] = '\0'), \
     (s))



/* define external routines */
extern void force_comment(GapIO *io, tagStruct *t);

/* comment interface */
extern tagStruct *readTagList(DBInfo *db, int seq, int ref);
extern void writeTagList(EdStruct *xx, int seq);
extern void destroyTagList(tagStruct *s);
extern void destroyFreeTagList(void);

extern void createAnnotation(EdStruct *xx);
extern void editAnnotation(EdStruct *xx);
extern void deleteAnnotation(EdStruct *xx);

extern void getTagSplodge(EdStruct *xx, int seq, int pos, int width, XawSheetInk *ink);

extern void insertTag(EdStruct *xx, int seq, tagStruct *t);
extern tagStruct *newTag(void);
void freeTag(tagStruct* t);
extern char normaliseBase(EdStruct *xx,int seq,char deletedBase);


extern tagStruct *findTag(EdStruct *xx,int seq,int pos);
extern tagStruct *findTagPos(EdStruct *xx, int seq, int pos);

extern char *get_comment(GapIO *io, comment_id cp);

extern void dump_tags(EdStruct *xx, int seq);
extern comment_id put_comment(GapIO *io, char *c);
extern int read_tag(GapIO *io, tag_id n, tagRecord *t);
extern int saveAnnotation(EdStruct *xx, char *type, char *anno, int strand);
extern int invokeTagEditor(EdStruct *xx, int tagid, int seq, int pos,
			   int length, int sense, char *comment, char *type_id,
			   tagStruct *tag);
extern int write_tag(GapIO *io, tag_id n, tagRecord t);
extern tag_id get_free_tag(GapIO *io);
extern void delete_tag_rec(GapIO *io, tag_id t);
extern void create_tag_for_gel(GapIO *io, int gel, int gellen, char *line,
			       int *cache, int cache_len, int *cache_pos);
extern tagStruct *findPreviousTag(EdStruct *xx, int seq, tagStruct *tag);

Array find_tags(GapIO *io, contig_list_t *contigs, int num_contigs,
		char **type, int num_t);


/*************************************************************
 * low level tag edits
 *************************************************************/


extern int _adjust_position_annotation(DBInfo *db, int seq, tagStruct *tag,
				       int pos, int seq_flags, int tag_flags);
/*
** Set the position of the tag
*/


extern int _adjust_length_annotation(DBInfo *db, int seq, tagStruct *tag,
				     int bases, int seq_flags, int tag_flags);
/*
** Set the position of the tag
*/

extern int _delete_annotation(DBInfo *db, int seq, tagStruct *tag,
			      int seq_flags);
/*
** Remove the edit tag that follows the supplied tag
*/

extern int _insert_annotation(DBInfo *db, int seq, tagStruct *tag,
			      tagStruct *newtag, int seq_flags);
/*
 *
 */

extern tagStruct *_create_annotation(EdStruct *xx, int seq, int pos, int length, char *type, char *comment, tagStruct *tag, int sense, int seq_flags);
/*
** Create a new tag with given attributes after tag
*/

extern int _destroy_annotation(DBInfo *db, EdStruct *xx, int seq, tagStruct *tag, int seq_flags);
/*
** Delete the edit tag that follows the supplied tag
**
** Only ever called in the context of undo _create_annotation
*/

extern int _create_insert_edit(EdStruct *xx, int seq, tagStruct *tag, int pos, char base, int seq_flags);
/*
** Create a new tag after the supplied tag to record the insertion of the base
*/

extern int _create_delete_edit(EdStruct *xx, int seq, tagStruct *tag, int pos, char base, int seq_flags);
/*
** Create a new tag after the supplied tag to record the deletion of the base
*/


extern int _destroy_edit(EdStruct *xx, int seq, tagStruct *tag, int seq_flags);
/*
** Delete the edit tag that follows the supplied tag
**
** Only ever called in the context of undo _create_edit
*/

extern int _delete_edit(EdStruct *xx, int seq, tagStruct *tag, int seq_flags);
/*
** Remove the edit tag that follows the supplied tag
*/

extern int _insert_edit (EdStruct *xx, int seq, tagStruct *tag, tagStruct *newtag, int seq_flags);
/*
** Insert newtag after tag
**
** Only ever called in the context of undo _delete_edit
*/


extern int _adjust_position_edit(EdStruct *xx, int seq, tagStruct *tag, int pos, int seq_flags, int tag_flags);
/*
** Set the position of the tag
*/




/*************************************************************
 * undo level tag edits
 *************************************************************/
extern int U_adjust_position_annotation(EdStruct *xx, int seq, tagStruct *tag, int pos);

extern int U_adjust_length_annotation(EdStruct *xx, int seq, tagStruct *tag, int bases);

extern int U_delete_annotation(EdStruct *xx, int seq, tagStruct *tag);

extern int U_create_annotation(EdStruct *xx, int seq, int pos, int length, char *type, char *comment, tagStruct *tag, int sense);


extern int U_delete_edit(EdStruct *xx, int seq, tagStruct *tag);

extern int U_adjust_position_edit(EdStruct *xx, int seq, tagStruct *tag, int pos);


extern GAnnotations *ctagget(GapIO *io, int gel, char *type);
extern GAnnotations *vtagget(GapIO *io, int gel, int num_t, char **type);

extern void shift_contig_tags(GapIO *io, int contig, int posn, int dist);
extern void merge_contig_tags(GapIO *io, int contig1, int contig2, int off);
extern void complement_contig_tags(GapIO *io, int contig);
extern void split_contig_tags(GapIO *io, int cont1, int cont2, int posl,
			      int posr);
extern void remove_contig_tags(GapIO *io, int contig, int posl, int posr);
extern void remove_gel_tags(GapIO *io, int gel, int posl, int posr);
extern int rmanno(GapIO *io, int anno, int lpos, int rpos);
extern Array anno_list(GapIO *io, int type);

extern int tag2values(char *tag, char *type, int *start, int *end, int *strand,
		      char *comment);
extern int values2tag(char *tag, char *type, int start, int end, int strand,
		      char *comment);
extern tag_id first_tag(GapIO *io, int N);
extern void insert_NEW_tag(GapIO *io, int N, int pos, int length, char *type,
			   char *comment, int sense);
extern void tag_shift_for_insert(GapIO *io, int seq, int pos);
extern void tag_shift_for_delete(GapIO *io, int seq, int pos);
extern void blank_tag_rec(GapIO *io, tag_id t);
extern void update_tag(GapIO *io, int N, tag_id anno);
extern void setUpColourMap(Tcl_Interp *interp, Tk_Window tkwin);
extern void tagDeleteBases(EdStruct *xx, int seq, int cursor_pos,
			   int num_bases);
extern void tagInsertBases(EdStruct *xx, int seq, int pos, int num_bases);

/* Fortran interfaces */
extern f_proc_ret mrgtag_(f_int *HANDLE, f_int *CONT1, f_int *CONT2,
			  f_int *OFF);
extern f_proc_ret comtag_(f_int *HANDLE, f_int *CONT, f_int *LEN);
extern f_proc_ret spltag_(f_int *HANDLE, f_int *CONT1, f_int *CONT2,
			  f_int *POSL, f_int *POSR);
extern f_proc_ret shiftt_(f_int *HANDLE, f_int *cont, f_int *posn, f_int *nc);
extern f_proc_ret rmctag_(f_int *handle, f_int *cont, f_int *lpos,
			  f_int *rpos);
extern f_proc_ret rmgtag_(f_int *handle, f_int *gel, f_int *lpos, f_int *rpos);

void tagfil_(f_int *relpg, f_int *lngthg, f_int *lnbr, f_int *rnbr,
	     f_int *ngels, f_int *nconts, f_int *idbsiz,
	     char *namarc, f_int *idev, f_int *verb, f_implicit l_namarc);
void write_tags(GapIO *io, char *fname, int num_tags,
		int *read1, int *pos1, int *read2, int *pos2, int *length);

/*
 * Removes all annotations listed in anno_ac and anno_av.
 *
 * Returns -1 for error, 0 for success.
 */
int rmanno_list(GapIO *io, int anno_ac, int *anno_av);

#endif  /*_tagUtils_h*/
