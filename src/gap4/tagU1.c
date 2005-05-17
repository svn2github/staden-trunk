/*
 * File: tagU1.c
 * Version:
 *
 * Author:
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

#define TAG_CHECK

#include <tk.h>
#include <stdio.h>

#include "tagDefs.h"    
#include "undo.h"
#include "misc.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "contigEditor.h"
#include "tkSheet.h"
#include "tagdb.h"
#include "notedb.h"
#include "fort.h"
#include "xalloc.h"
#include "active_tags.h"
#include "select.h"
#include "dna_utils.h"

static tagStruct *tagFreeList = NULL;


Pixel ColourNameToPixel(Tcl_Interp *interp, Tk_Window tkwin, char *name)
{
    XColor *cp = Tk_GetColor(interp, tkwin, name);

    if (cp) {
	return cp->pixel;
    } else {
	verror(ERR_WARN, "ColourNameToPixel", "Colourmap is full");
	return 0;
    }
}

void setUpColourMap(Tcl_Interp *interp, Tk_Window tkwin)
{
    static int done = 0;
    int i;

    if (done)
	return;
    else
	done = 1;

    for (i=0;i<tag_db_count;i++) {
        tag_db[i].fg_pixel =  (tag_db[i].fg_colour == NULL) ?
            1 : ColourNameToPixel(interp, tkwin, tag_db[i].fg_colour);
        tag_db[i].bg_pixel =  (tag_db[i].bg_colour == NULL) ?
            0 : ColourNameToPixel(interp, tkwin, tag_db[i].bg_colour);
        tag_db[i].gf_pixel =  (tag_db[i].gf_colour == NULL) ?
            1 : ColourNameToPixel(interp, tkwin, tag_db[i].gf_colour);
        tag_db[i].gb_pixel =  (tag_db[i].gb_colour == NULL) ?
            0 : ColourNameToPixel(interp, tkwin, tag_db[i].gb_colour);
    }

    for (i = 0; i < note_db_count; i++) {
	note_db[i].fg_pixel = note_db[i].fg_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].fg_colour) : 1;
	note_db[i].bg_pixel = note_db[i].bg_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].bg_colour) : 0;
	note_db[i].gf_pixel = note_db[i].gf_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].gf_colour) : 1;
	note_db[i].gb_pixel = note_db[i].gb_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].gb_colour) : 0;
    }
}


/*************************************************************
 *
 *************************************************************/

static void link_tag(DBInfo *db, int seq, tagStruct *tag, tagStruct *newtag)
/*
 * Link a tag into the tag list
 */
{
    tagStruct *next;
    
    if (newtag==NULL) return;
    
    /* determine next */
    if (tag == NULL)
	next = DBgetTags(db, seq);
    else
	next = tag->next;
    
    /* link */
    newtag->next = next;
    if (tag == NULL)
	_DBsetTags(db, seq, newtag);
    else
	tag->next = newtag;
    
}




static tagStruct *delink_tag(DBInfo *db, int seq, tagStruct *tag)
/*
 * Remove the tag following the tag from the tag list
 */
{
    tagStruct *oldtag;
    
    /* The tag lists always start with a blank header tag, so NULL is error */
    if (tag == NULL)
	return NULL;

    oldtag = tag->next;
    
    if (oldtag != NULL) {
	/* delink */
	tag->next = oldtag->next;

	/*
	 * Mark this tag to indicate that the next tag is deleted. This speeds
	 * up the writeTagList function. Note though that this isn't reset on
	 * undo, but as it's only an optimisation it doesn't cause problems,
	 * just slow downs.
	 */
	tag->flags |= TAG_NEXT_DELETED;
    }
    
    return oldtag;
}




tagStruct *findTagPos(EdStruct *xx, int seq, int pos)
/*
 * Return a pointer to the last tag at or before the position
 */
{
    tagStruct *t, *last_t;
    /* read and skip over tag header */
    last_t = t = DBgetTags(DBI(xx),seq);
    if (t != NULL) t = t->next;
    
    while ( t != NULL && t->tagrec.position <= pos ) {
	last_t = t;
	t = t->next;
    }
    
    return last_t;
    
}





/*************************************************************
 *
 *************************************************************/




/* low level */

void force_comment(GapIO *io, tagStruct *t)
/*
 * Force comment to be in memory
 */
{
    if (!(t->flags & TAG_COMMENT_IN_MEMORY)) {
	/*
	 * Read in from database
	 */
	if (t->tagrec.comment) {
	    t->newcomment = get_comment(io, t->tagrec.comment);
	    t->newcommentlen = (int)strlen(t->newcomment);
	} else {
	    t->newcomment = (char *) TAG_MALLOC(1);
	    t->newcomment[0] = '\0';
	    t->newcommentlen = 0;
	}
	t->flags |= TAG_COMMENT_IN_MEMORY;
    }
}


/*
 * Tag internal memory management routines
 */
tagStruct *newTag()
{
    tagStruct *t;
    if (tagFreeList == NULL) {
        t = (tagStruct *) TAG_MALLOC(sizeof(tagStruct));
    } else {
	t = tagFreeList;
	tagFreeList = t->next;
    }
    
    /*
     * Null all the fields
     */
    t->tagrec.position = 0;
    t->tagrec.length = 0;
    t->tagrec.comment = 0;
    t->tagrec.type.i = 0x20202020;
    t->tagrec.next = 0;
    t->tagrec.sense = 0;
    t->original_tag_id = 0;
    t->newcomment = NULL;
    t->newcommentlen = 0;
    t->flags = TAG_UNCHANGED;
    t->next = NULL;

    return t;
}







void freeTag(tagStruct* t)
{
    if (t==NULL) return;
    if (t->newcomment) {
	TAG_FREE(t->newcomment);
	t->newcomment = NULL;
    }
    t->newcommentlen = 0;
    t->next = tagFreeList;
    tagFreeList = t;
}


/*
 * Frees memory used up by the global cache of unused tags (tagFreeList).
 * Not essential - we could just leave this around for reuse by subsequent
 * editors.
 */
void destroyFreeTagList(void) {
    tagStruct *t, *u;

    for (t = tagFreeList; t; t=u) {
	u = t->next;

	if (t->newcomment) {
	    xfree(t->newcomment);
	}
	xfree(t);
    }
    tagFreeList = NULL;
}

void destroyTagList(tagStruct *s)
{
    tagStruct *t,*u;

    t=s;
    while (t!=NULL) {
        u=t->next;
	freeTag(t); /* only removes comment, but keeps t on a freelist */
	t=u;
    }
}

/*
 * Tag creation and modification
 */

void insertTag(EdStruct *xx, int seq, tagStruct *t)
/*
 * insert tag, sorting by position
 */
{
    tagStruct *u, *v;
    
    u = (tagStruct *) DBgetTags(DBI(xx),seq);
    v = NULL;
    while (u != NULL &&
	   (u->tagrec.position <= t->tagrec.position) ) {
	v = u;
	u = u->next;
    }
    t->next = u;
    if (v != NULL) {
	v->next = t;
    } else {
	DBsetTags(xx,seq,t);
    }
    
}












void getTagSplodge(EdStruct *xx, int seq, int pos, int width, XawSheetInk *ink)
/*
 * get the hilighting of a sequence from its `pos' base for `width' bases
 * Bases number from 0?
 */
{
    
    int i;
    tagStruct *t;
    int npos,tpos;
    
    if (!xx->tag_list)
	return; /* happens when GTAGDB is set wrongly */

    if (xx->reveal_cutoffs) {
	int length = DB_Length(xx,seq);
	
        /*blank start*/
        for (i=0; i<width && i<-pos; i++)
	    ink[i].sh=sh_light;
	
        /*copy sequence*/
        for (; i<width && (pos+i)<length; i++)
	    ink[i].sh=sh_default;
	
        /*blank end*/
        for (;i<width;i++)
	    ink[i].sh=sh_light;
	
	
    } else
	for (i=0;i<width;i++)
	    ink[i].sh=sh_default;
    
    pos++;

    pos += DB_Start(xx, seq);
    npos = normalisePos2(xx,seq,pos,width);

    t = (tagStruct *) DBgetTags(DBI(xx),seq);
    /* skip over header */
    if (seq && t != NULL) t = t->next;
    
    for (; t != NULL && (t->tagrec.position < npos+width); t = t->next) {
	if (t->tagrec.position + t->tagrec.length > npos ) {
	    int l,r;
	    int db=idToIndex(t->tagrec.type.c);

	    if (!xx->tag_list[db])
		continue;

	    tpos = normalisePos2(xx, seq, t->tagrec.position,
				 t->tagrec.length);

	    if (tpos < pos)
		l=0;
	    else
		l=tpos-pos;

	    if (tpos + t->tagrec.length > pos+width)
		r=width;
	    else
		r=tpos-pos + t->tagrec.length;

	    if (normaliseSense(xx,seq,t->tagrec.sense)!=1) {
		for (i=l;i<r;i++) {
		    if (tag_db[db].fg_colour!=NULL) {
			ink[i].sh|=sh_fg;
			ink[i].fg=tag_db[db].fg_pixel;
		    }
		    if (tag_db[db].bg_colour!=NULL) {
			ink[i].sh|=sh_bg;
			ink[i].bg=tag_db[db].bg_pixel;
		    }
		}
	    } else {
		for (i=l;i<r;i++) {
		    if (tag_db[db].gf_colour!=NULL) {
			ink[i].sh|=sh_fg;
			ink[i].fg=tag_db[db].gf_pixel;
		    }
		    if (tag_db[db].gb_colour!=NULL) {
			ink[i].sh|=sh_bg;
			ink[i].bg=tag_db[db].gb_pixel;
		    }
		}
	    }

	}
    }
}


char normaliseBase(EdStruct *xx,int seq,char deletedBase)
{
    
    if (DB_Comp(xx,seq) == COMPLEMENTED) {
	return complement_base(deletedBase);
    } else
	return deletedBase;
}





tagStruct *findTag(EdStruct *xx,int seq,int pos)
/*
 * Find the tag (if any) at position `pos' in sequence `seq'
 */
{
    int npos = normalisePos2(xx,seq,pos,1/*character*/);
    
    tagStruct *t, *found = NULL;
    t = (tagStruct *) DBgetTags(DBI(xx),seq);
    while (t != NULL) {
	if (t->tagrec.position <= npos &&
	    t->tagrec.position + t->tagrec.length > npos &&
	    xx->tag_list[idToIndex(t->tagrec.type.c)])
	    found = t;
	t = t->next;
    }
    return found;
}




tagStruct *findAllTags(EdStruct *xx,int seq,int pos)
/*
 * Find the tag (if any) at position `pos' in sequence `seq'
 */
{
    static tagStruct *t;
    static int npos;
    if (xx==NULL) {
	if (t != NULL) t = t->next;
    } else {
	npos = normalisePos2(xx,seq,pos,1/*character*/);
	t = (tagStruct *) DBgetTags(DBI(xx),seq);
    }
    
    while (t != NULL) {
	if (t->tagrec.position <= npos &&
	    t->tagrec.position + t->tagrec.length > npos)
	    return t;
	t = t->next;
    }
    return NULL;
}





void dump_tags(EdStruct *xx, int seq)
{
    tagStruct *t = (tagStruct *) DBgetTags(DBI(xx),seq);
    
    fprintf(stderr,"Tags for %s\n",DBgetName(DBI(xx),seq));
    while (t != NULL) {
	
	fprintf(stderr,"%5d: %5d %3d %.4s %5d %c%c%c%c%c%c%c %5d\n",
		t->original_tag_id,
		t->tagrec.position,
		t->tagrec.length,
		t->tagrec.type.c,
		t->tagrec.comment,
		(t->flags & TAG_POSITION_CHANGED) ?'P':'-',
		(t->flags & TAG_LENGTH_CHANGED)   ?'L':'-',
		(t->flags & TAG_TYPE_CHANGED)     ?'T':'-',
		(t->flags & TAG_COMMENT_CHANGED)  ?'C':'-',
		(t->flags & TAG_INSERTED)         ?'I':'-',
		(t->flags & TAG_COMMENT_IN_MEMORY)?'M':'-',
		(t->flags & TAG_NEXT_DELETED)     ?'D':'-',
		t->tagrec.next
		);
	
	t = t->next;
    }
}




/*************************************************************
 *
 *************************************************************/

#include "undo.h"


/*************************************************************
 *
 * New routines to support new edit routines
 *
 *************************************************************/


/************************************************************
 * Low level for annotation manipulation
 ************************************************************/

int _adjust_position_annotation(DBInfo *db, int seq, tagStruct *tag,
				int pos, int seq_flags, int tag_flags)
/*
 * Set the position of the tag
 */
{
    if (tag==NULL) return 1;
    
    tag->tagrec.position = pos;
    
    _DBsetFlags(db,seq,seq_flags);
    tag->flags = tag_flags;
    
    return 0;
}




int _adjust_length_annotation(DBInfo *db, int seq, tagStruct *tag,
			      int bases, int seq_flags, int tag_flags)
/*
 * Set the position of the tag
 */
{
    if (tag==NULL) return 1;
    
    tag->tagrec.length = bases;
    
    _DBsetFlags(db,seq,seq_flags);
    tag->flags = tag_flags;
    
    return 0;
}



int _delete_annotation(DBInfo *db, int seq, tagStruct *tag, int seq_flags)
/*
 * Remove the edit tag that follows the supplied tag
 */
{
    tagStruct *t;
    t = delink_tag(db, seq, tag);
    /*
     * NB: DO NOT FREE THIS! It is used somewhere. Probably when undo() is
     * called. Not freeing causes a memory leak, but it's small (44 bytes).
     *
     * freeTag(t);
     */
    _DBsetFlags(db,seq,seq_flags);

    return 0;
    
}



int _insert_annotation(DBInfo *db, int seq, tagStruct *tag, tagStruct *newtag,
		       int seq_flags)
/*
 *
 */
{
    link_tag(db,seq,tag,newtag);
    _DBsetFlags(db,seq,seq_flags);
    return 0;
    
}



tagStruct *_create_annotation(EdStruct *xx, int seq, int pos, int length,
			      char *type, char *comment, tagStruct *tag,
			      int sense, int seq_flags)
/*
 * Create a new tag with given attributes after tag
 */
{
    tagStruct *t;
    
    /*
     * Create an new insert tag here
     */
    t = newTag();
    if (t==NULL)
	return NULL;
    
    t->flags = TAG_INSERTED;
    
    t->tagrec.position = pos;
    t->tagrec.length = length;
    strncpy(t->tagrec.type.c,type,4);
    t->newcomment = comment;
    if (t->newcomment != NULL) {
	t->newcommentlen = (int)strlen(comment);
	t->flags=TAG_COMMENT_IN_MEMORY|TAG_COMMENT_CHANGED;
    }
    t->tagrec.sense = sense;
    
    link_tag(DBI(xx),seq,tag,t);
    DBsetFlags(xx,seq,seq_flags);
    if (seq > 0)
	RedisplaySeq(xx, seq);
    else
	xx->refresh_flags |= ED_DISP_CONS;

    _select_tag(xx, seq, t);

    return t;
}

int _destroy_annotation(DBInfo *db, EdStruct *xx, int seq, tagStruct *tag,
			int seq_flags)
/*
 * Delete the edit tag that follows the supplied tag
 *
 * Only ever called in the context of undo _create_annotation or removing
 * the temporary tag created by Select Primer.
 */
{
    tagStruct *t;
    
    if (xx && tag && xx->select_tag == tag->next) {
	_select_tag(xx, seq, NULL);
    }

    t = delink_tag(db,seq,tag);
    freeTag(t);

    _DBsetFlags(db,seq,seq_flags);
    if (seq > 0)
	RedisplaySeq(xx, seq);
    else
	xx->refresh_flags |= ED_DISP_CONS;

    return 0;
}



/*************************************************************
 * Mid level for annotation manipulation (undo information)
 *************************************************************/
int U_adjust_position_annotation(EdStruct *xx, int seq, tagStruct *tag, int pos)
/*
 * YUK! should we send just an adjustment - as with cursor positioning?
 */
{
    int tag_flags, old_tag_flags;
    int seq_flags, old_seq_flags;
    UndoStruct *u;
    int old_pos;
    
    if (tag==NULL) return 1;
    old_pos = tag->tagrec.position;
    old_tag_flags = tag->flags;
    old_seq_flags = DB_Flags(xx,seq);
    
    /* UNDO - _adjust_position_annotation(DBI(xx),tag,pos) */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoAdjustPositionAnnotation;
	u->sequence = seq;
	u->info.adjust_position_annotation.tag = tag;
	u->info.adjust_position_annotation.position = old_pos;
	u->info.adjust_position_annotation.tag_flags = old_tag_flags;
	u->info.adjust_position_annotation.seq_flags = old_seq_flags;
	recordUndo(DBI(xx),u);
    }
    
    tag_flags = old_tag_flags | TAG_POSITION_CHANGED;
    seq_flags = old_seq_flags | DB_FLAG_TAG_MODIFIED;
    return _adjust_position_annotation(DBI(xx),seq,tag,pos,seq_flags,tag_flags);
}


int U_adjust_length_annotation(EdStruct *xx, int seq, tagStruct *tag, int bases)

{
    int tag_flags, old_tag_flags;
    int seq_flags, old_seq_flags;
    UndoStruct *u;
    int old_len;
    
    if (tag==NULL) return 1;
    old_len = tag->tagrec.length;
    old_tag_flags = tag->flags;
    old_seq_flags = DB_Flags(xx,seq);
    
    /* UNDO - _adjust_length_annotation(DBI(xx),tag,old_len) */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoAdjustLengthAnnotation;
	u->sequence = seq;
	u->info.adjust_length_annotation.tag = tag;
	u->info.adjust_length_annotation.length = old_len;
	u->info.adjust_length_annotation.tag_flags = old_tag_flags;
	u->info.adjust_length_annotation.seq_flags = old_seq_flags;
	recordUndo(DBI(xx),u);
    }
    
    tag_flags = old_tag_flags | TAG_POSITION_CHANGED;
    seq_flags = old_seq_flags | DB_FLAG_TAG_MODIFIED;
    return _adjust_length_annotation(DBI(xx),seq,tag,bases,seq_flags,tag_flags);
}


int U_delete_annotation(EdStruct *xx, int seq, tagStruct *tag)
{
    UndoStruct *u;
    tagStruct *t;
    int seq_flags, old_seq_flags;
    
    if (tag==NULL) return 1;
    if ((t = tag->next)==NULL) return 1;
    
    old_seq_flags = DB_Flags(xx,seq);
    /* UNDO - _insert_annotation(xx,...) */
    if ((u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoInsertAnnotation;
	u->sequence = seq;
	u->info.insert_annotation.tag = tag;
	u->info.insert_annotation.new_tag = t;
	u->info.insert_annotation.seq_flags = old_seq_flags;
	
	recordUndo(DBI(xx),u);
    }
    
    seq_flags = old_seq_flags | DB_FLAG_TAG_MODIFIED;
    (void) _delete_annotation(DBI(xx),seq,tag,seq_flags);
    if (seq > 0)
	RedisplaySeq(xx, seq);
    else
	xx->refresh_flags |= ED_DISP_CONS;
    xx->refresh_flags |= ED_DISP_SELECTION;
    
    return 0;
}

int U_create_annotation(EdStruct *xx, int seq, int pos, int length, char *type, char *comment, tagStruct *tag, int sense)
{
    UndoStruct *u;
    int seq_flags, old_seq_flags;
    
    old_seq_flags = DB_Flags(xx,seq);
    
    /* UNDO - _destroy_annotation(xx,seq,tag) */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoDestroyAnnotation;
	u->sequence = seq;
	u->info.destroy_annotation.xx = xx;
	u->info.destroy_annotation.tag = tag;
	u->info.destroy_annotation.seq_flags = old_seq_flags;
	recordUndo(DBI(xx),u);
    }
    
    seq_flags = old_seq_flags | DB_FLAG_TAG_MODIFIED;
    _create_annotation(xx,seq,pos,length,type,comment,tag,sense,seq_flags);

    /*
     * Reset seqeunce cursor, to force CURSOR_NOTIFY and a refresh of the
     * editor tag menus (a hack!).
     */
    U_adjust_cursor(xx, 0);
    
    return 0;
}









/************************************************************
 * Functional level
 ************************************************************/





/************************************************************
 * Interface level
 ************************************************************/


void tagInsertBases(EdStruct *xx, int seq, int pos, int num_bases)
/*
 * Adjust tag structure for insert
 */
{
    int npos;
    tagStruct *last_t, *t, *create_here;
    int this_pos, this_len;
    int last_pos;
    
    /*
     * Normalise everything
     */
    pos += DB_Start(xx, seq);
    npos = normalisePos2(xx,seq,pos,0);

    /*
     * Search along tag list
     * 
     * insert tags at:
     *     npos           : insert nbases[0]
     *     npos+1         : insert nbases[1]
     *     ...
     *     npos+num_bases : insert nbases[num_bases-1]
     */
    
    /* read and skip over tag header */
    last_t = t = DBgetTags(DBI(xx),seq);
    if (t != NULL) t = t->next;
    create_here = NULL;
    last_pos = last_t->tagrec.position;
    
    while (t != NULL) {
	if (last_pos < npos && t->tagrec.position >= npos) {
	    /*
	     * We should create tags here
	     */
	    create_here = last_t;
	}
	/*
	 * Grab this now - t's position may change
	 */
	last_pos = t->tagrec.position;
	/* Adjust annotations */
	this_pos = t->tagrec.position;
	this_len = t->tagrec.length;
	if (this_pos >= npos) {
	    /* need to shift right */
	    this_pos += num_bases;
	    U_adjust_position_annotation(xx, seq, t, this_pos);
	} else if ((this_pos + this_len - 1) >= npos) {
	    /* need to lengthen */
	    this_len += num_bases;
	    U_adjust_length_annotation(xx, seq, t, this_len);
	}
	
	last_t = t;
	t = t->next;
    }
    
    
    
}




void tagDeleteBases(EdStruct *xx, int seq, int cursor_pos, int num_bases)
/*
 * Adjust tag structure for delete
 * deletion of num_bases before pos
 */
{
    int pos;
    int npos;
    tagStruct *last_t, *t;
    int this_pos, this_len;
    int was_pos;
    int end, this_end;
    
    /*
     * Normalise everything
     */
    pos = cursor_pos - num_bases + 1;
    pos += DB_Start(xx, seq);
    npos = normalisePos2(xx,seq,pos,num_bases);
    end = npos + num_bases - 1;
    
    /*
     * From here on:
     *   npos is the position of the first base deleted
     *   num_bases is number deleted
     */
    
    
    /*
     * Second pass - adjust positions of all tags to reflect deletion
     */
    /* read and skip over tag header */
    last_t = t = DBgetTags(DBI(xx),seq);
    if (t != NULL) t = t->next;
    
    while (t != NULL) {
	/*
	 * Adjust tags accordingly
	 */
	was_pos = this_pos = t->tagrec.position;
	this_len = t->tagrec.length;
	this_end = this_pos + this_len - 1;
	/* six cases to consider */
	if (this_end < npos) {
	    /* case 1 - not affected */
	} else if (this_pos < npos && this_end <= end) {
	    /* case 2 - right end of annotaion affected */
	    this_len = npos - this_pos;
	    U_adjust_length_annotation(xx, seq, t, this_len);
	} else if (this_pos < npos /*&& this_end > end*/) {
	    /* case 3 - deletion is in middle of annotation */
	    this_len -= num_bases;
	    U_adjust_length_annotation(xx, seq, t, this_len);
	} else if (this_pos <= end && this_end <= end) {
	    /* case 4 - annotation totally deleted */
	    openUndo(DBI(xx));
	    U_adjust_cursor(xx, 0);
	    U_delete_annotation(xx, seq, last_t);
	    U_adjust_cursor(xx, 0);
	    closeUndo(xx, DBI(xx));
	    /* YUK! */
	    t = last_t;
	} else if (this_pos <= end /*&& this_end > end*/) {
	    /* case 5 - left end affected */
	    this_pos = npos;
	    this_len = this_end - end;
	    U_adjust_length_annotation(xx, seq, t, this_len);
	    U_adjust_position_annotation(xx, seq, t, this_pos);
	} else /*if (this_pos > end)*/ {
	    /* case 6 - annotation falls to right of tag */
	    this_pos -= num_bases;
	    U_adjust_position_annotation(xx, seq, t, this_pos);
	}
	
	last_t = t;
	t = t->next;
	
    }
    
    
    
}











/*************************************************************
 *
 *************************************************************/

char *createAnnotation(EdStruct *xx)
{
    int seq,start,length;
    tagStruct *t;
    static tag_id id = 0;
    
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return NULL;
    }

    if (! getSelection(xx, &seq, &start, &length, &t) || length == 0) {
	/* default selection is current cursor position */
	seq = xx->cursorSeq;
	start = xx->cursorPos + DB_Start(xx, seq);
	length = 1;
	if (start > DB_Length2(xx, seq)) {
	    bell();
	    return NULL;
	}
    }
    
    /* Pick a new tag number; -ve to show that it's unsaved */
    id--;

    /* And invoke the editor */
    return invokeTagEditor(xx, id, seq, start, length, 0 /* sense */,
			   "", "NONE",
			   NULL /* tagStruct for this tag */);
}



tagStruct *findPreviousTag(EdStruct *xx, int seq, tagStruct *tag)
/*
 * Find the tag (if any) at position `pos' in sequence `seq'
 */
{
    tagStruct *t;
    tagStruct *last_t;
    last_t = NULL;
    t = (tagStruct *) DBgetTags(DBI(xx),seq);
    while (t != NULL && t != tag) {
	last_t = t;
	t = t->next;
    }
    return last_t;
}





void deleteAnnotation(EdStruct *xx, tagStruct *t)
/*
 * A rather brutal delete
 */
{
    int seq,start,length;
    tagStruct *tag;
    
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return;
    }

    if (t) {
	seq = xx->cursorSeq;
	_select_tag(xx,seq,t);
    } else {
	if (! getSelection(xx, &seq, &start, &length, &t)) {
	    /* default selection is current cursor position */
	    seq = xx->cursorSeq;
	    start = xx->cursorPos + DB_Start(xx, seq);
	    t = NULL;
	}
	if (t==NULL) {
	    t = findTag(xx,seq,start);
	    _select_tag(xx,seq,t);
	    (void) getSelection(xx, &seq, &start, &length, &t);
	    if (t==NULL) return;
	}
    }

    /* We've now got the tag details, so clear selection */
    edSelectClear(xx);
    
    /* find tag prior to t */
    tag = findPreviousTag(xx,seq,t);
    
    openUndo(DBI(xx));
    U_adjust_cursor(xx, 0);
    U_delete_annotation(xx,seq,tag);
    U_adjust_cursor(xx, 0);
    closeUndo(xx, DBI(xx));
    
    /* Force a redisplay */
    redisplaySequences(xx, 1);
}



/*
 * Returns a list of all tag ids under the current cursor position.
 * The list consists of zero or more curly brace enclosed lists of
 * tagStruct pointer, tag type, position, length.
 *
 * The returned value is
 *    NULL for failure
 *    ndstring_t * for success.
 */
dstring_t *listAnnotation(EdStruct *xx) {
    int seq,pos;
    tagStruct *t;
    int npos;
    dstring_t *ds = dstring_create(NULL);

    seq = xx->cursorSeq;
    pos = xx->cursorPos + DB_Start(xx, seq);
    
    npos = normalisePos2(xx, seq, pos, 1/*character*/);
    
    t = (tagStruct *) DBgetTags(DBI(xx),seq);
    while (t != NULL) {
	if (t->tagrec.position <= npos &&
	    t->tagrec.position + t->tagrec.length > npos &&
	    xx->tag_list[idToIndex(t->tagrec.type.c)]) {
	    dstring_appendf(ds, "{%p %.4s %d %d} ",
			    t,
			    t->tagrec.type.c,
			    t->tagrec.position,
			    t->tagrec.length);
	}
	t = t->next;
    }

    return ds;
}

    
/* Edit annotation 't', or the one under the cursor if 't' is null */
char *editAnnotation(EdStruct *xx, tagStruct *t)
{
    int seq,start,length;

    if (t) {
	seq = xx->cursorSeq;
	_select_tag(xx,seq,t);
    } else {
	if (! getSelection(xx, &seq, &start, &length, &t)) {
	    /* default selection is current cursor position */
	    seq = xx->cursorSeq;
	    start = xx->cursorPos + DB_Start(xx, seq);
	    t = findTag(xx,seq,start);
	    _select_tag(xx,seq,t);
	    (void) getSelection(xx, &seq, &start, &length, &t);
	} else if (t==NULL) {
	    t = findTag(xx,seq,start);
	    _select_tag(xx,seq,t);
	    (void) getSelection(xx, &seq, &start, &length, &t);
	}
    }
    if (t==NULL) return NULL;
    
    /*
     * Find current comment
     */
    force_comment(DBI_io(xx), t);

    return invokeTagEditor(xx, t->original_tag_id, seq, t->tagrec.position,
			   t->tagrec.length,
			   normaliseSense(xx, seq, t->tagrec.sense),
			   t->newcomment, t->tagrec.type.c, t);
}












/************************************************************
 * WRITE TAG LIST
 ************************************************************/

/*
 * Writes the current tag list to disk.
 * For efficiency, we don't update tags that have not changed. However seeing
 * what's changed and what hasn't isn't obvious.
 *
 * We need to scan through both memory and disk copies as the memory data
 * may have tags missing (these need deallocating), and may have extra ones.
 * We can't do a merge sort tactic as the positions may have changed slightly.
 */
void writeTagList(EdStruct *xx, int seq) {
    tagStruct *mem_t, *head_t, *last_t = NULL;
    GapIO *io = DBI_io(xx);
    GAnnotations disk_curr_t, disk_last_t;
    int first_tag_id;
    int curr_id, curr_del;
    int at_header = 1;

    /* determine if any tags have changes */
    if (! (DB_Flags(xx,seq) & DB_FLAG_TAG_MODIFIED) ) return;

    /* Load tag list, skipping first 'header tag' */
    mem_t = head_t = DBgetTags(DBI(xx),seq);
    curr_id = first_tag_id = first_tag(io, DB_Number(xx, seq));

    /* Loop through tags in memory */
    for (; mem_t; last_t = mem_t, mem_t = mem_t->next) {

	if (mem_t->original_tag_id) {
	    at_header = 0;
	    curr_id = mem_t->original_tag_id;
	}

	/*
	 * If TAG_NEXT_DELETED is set, the next tags _may_ have been deleted.
	 * This isn't guaranteed at present due to 'undo' not resetting the
	 * information (this is actually trickier than it sounds).
	 */
	if (mem_t->flags & TAG_NEXT_DELETED && curr_id != 0) {
	    tagStruct *mem_next;
	    int search_id;
	    
	    for (mem_next = mem_t->next; mem_next; mem_next = mem_next->next) {
		if (mem_next->original_tag_id != 0)
		    break;
	    }

	    /*
	     * Found the next tag in memory; search for this on disk, deleting
	     * any tags inbetween.
	     */
	    search_id = mem_next ? mem_next->original_tag_id : 0;
	    if (!at_header) { /* skip for first tag */
		tag_read(io, curr_id, disk_curr_t);
		curr_del = disk_curr_t.next;
	    } else {
		curr_del = curr_id;
		curr_id = search_id;
	    }

	    while (curr_del && curr_del != search_id) {
		int link_to;

		tag_read(io, curr_del, disk_curr_t);
		link_to = disk_curr_t.next;

		delete_tag_rec(io, curr_del);

		mem_t->tagrec.next = link_to;
		if (last_t && curr_del == last_t->tagrec.next)
		    last_t->tagrec.next = link_to;

		if (mem_t->original_tag_id) {
		    tag_read(io, mem_t->original_tag_id, disk_curr_t);
		    disk_curr_t.next = link_to;
		    tag_write(io, mem_t->original_tag_id, disk_curr_t);
		} else if (at_header) {
		    update_tag(io, DB_Number(xx, seq), link_to);
		}
		curr_del = link_to;
	    }
	    if (curr_id)
		tag_read(io, curr_id, disk_curr_t);

	    mem_t->flags &= ~TAG_NEXT_DELETED;
	}

	/* Is it modified? Skip if not */
	if ((mem_t->flags & (TAG_POSITION_CHANGED |
			     TAG_LENGTH_CHANGED |
			     TAG_TYPE_CHANGED |
			     TAG_COMMENT_CHANGED |
			     TAG_INSERTED)) == 0) {
	    continue;
	}

	if (mem_t->original_tag_id) {
	    /* It exists, hence we just update the disk image */
	    tag_read(io, mem_t->original_tag_id, disk_curr_t);
 	} else {
	    /* otherwise we need to create it and relink curr_id */
	    curr_id = mem_t->original_tag_id = get_free_tag(io);
	    at_header = 0;
	    memset(&disk_curr_t, 0, sizeof(disk_curr_t));

	    if (last_t) {
		disk_curr_t.next       = last_t->tagrec.next;
		mem_t->tagrec.next = last_t->tagrec.next;

		if (last_t->original_tag_id) {
		    tag_read(io, last_t->original_tag_id, disk_last_t);
		    disk_last_t.next = mem_t->original_tag_id;
		    tag_write(io, last_t->original_tag_id, disk_last_t);
		} else {
		    update_tag(io, DB_Number(xx, seq), mem_t->original_tag_id);
		}
		last_t->tagrec.next = mem_t->original_tag_id;
	    }
	}

	disk_curr_t.type       = str2type(mem_t->tagrec.type.c);
	disk_curr_t.length     = mem_t->tagrec.length;
	disk_curr_t.position   = mem_t->tagrec.position;
	disk_curr_t.strand     = mem_t->tagrec.sense;
	disk_curr_t.annotation = mem_t->tagrec.comment;

	if (mem_t->flags & TAG_COMMENT_CHANGED) {
	    if (mem_t->newcommentlen > 0) {
		mem_t->tagrec.comment = disk_curr_t.annotation =
		    put_comment(io, mem_t->newcomment);
	    } else {
		if (disk_curr_t.annotation) {
		    deallocate(io, disk_curr_t.annotation);
		    mem_t->tagrec.comment = disk_curr_t.annotation = 0;
		}
	    }
	}

	tag_write(io, mem_t->original_tag_id, disk_curr_t);

	mem_t->flags = 0;
    }

    DB_Flags(xx,seq) &= ~DB_FLAG_TAG_MODIFIED;
}

void writeTagList_old(EdStruct *xx, int seq)
/*
 * Write the tag list
 *
 * PASS 1:
 *	For each tag in memory, create a new tag on disk
 *	Upgrade old comment if necessary
 *
 * CHANGE-OVER:
 *	Change over tag lists
 *
 * PASS 2:
 *	For each tag in memory where the comment has not changed,
 *	delink old comment from original tag record
 *	
 * PASS 3:
 *	Reclaim the original taglist on disk
 *
 */
{
    tagStruct *t;
    tag_id newTagList;
    tag_id oldTagList;
    tag_id thisTag, nextTag, headTag;
    tagRecord thisTagRec;
#ifdef TAG_CHECK
    int chkpos = 0;
    int gellen = seq ? DB_Length2(xx,seq) : DB_Length(xx, seq);
#endif /*TAG_CHECK*/
    GapIO *io = DBI_io(xx);
    
    /* determine if any tags have changes */
    if (! (DB_Flags(xx,seq) & DB_FLAG_TAG_MODIFIED) ) return;
    
    /* */
    t = DBgetTags(DBI(xx),seq);
    if (t==NULL) return;
    
    /******************
     * PASS THE FIRST *
     ******************/
    
    /*
     * Assumption: first tag in list is always a header and must be treated
     * specially
     */
    t = t->next;
    
    /* determine the tag_id for the first tag */
    if (t != NULL)
	nextTag = get_free_tag(io);
    else
	nextTag = 0;
    newTagList = nextTag;
    
    while(t != NULL) {

#ifdef TAG_CHECK
	if (t->tagrec.position <= 0 ||
	    t->tagrec.position + t->tagrec.length > gellen + 1)
	    /*
	     * Tags should never run off the end of the gel
	     */
	    verror(ERR_WARN, "writeTagList",  "INVALID TAG POSITION "
		   "seq=%d (%s) tagpos=%d taglen=%d gellen=%d\n",
		    seq,
		    DBgetName(DBI(xx),seq),
		    t->tagrec.position,
		    t->tagrec.length,
		    gellen);

	if (t->tagrec.position < chkpos)
	    /*
	     * Tag positions should never be out of order
	     */
	    verror(ERR_WARN, "writeTagList",  "TAG OUT OF ORDER "
		   "seq=%d (%s) tagpos=%d taglen=%d\n",
		    seq,
		    DBgetName(DBI(xx),seq),
		    t->tagrec.position,
		    t->tagrec.length );
	chkpos = t->tagrec.position;
#endif /*TAG_CHECK*/



	thisTag = nextTag;
	
	if (t->next != NULL)
	    nextTag = get_free_tag(io);
	else
	    nextTag = 0;
	
	/* set up this tag records */
	thisTagRec.position = t->tagrec.position;
	thisTagRec.length   = t->tagrec.length;
	thisTagRec.type.i   = t->tagrec.type.i;
	thisTagRec.next     = nextTag;
	thisTagRec.sense    = t->tagrec.sense;
	
	/* now the hard part - update comment */
	if (t->flags & TAG_COMMENT_CHANGED) {
	    if (t->newcommentlen>0)
		thisTagRec.comment = put_comment(io, t->newcomment);
	    else
		thisTagRec.comment = 0;
	} else
	    thisTagRec.comment = t->tagrec.comment;
	
	
	write_tag(io, thisTag,thisTagRec);
	
	t = t->next;
    }
    
    
    
    /*******************
     * THE CHANGE-OVER *
     *******************/
    t = DBgetTags(DBI(xx),seq);
    headTag = DB_Number(xx,seq);
    
    /* read original */
    oldTagList = first_tag(io, headTag);
    
    /* write the record back */
    update_tag(io, headTag, newTagList);

    
    
    /*******************
     * PASS THE SECOND *
     *******************/
    
    /*
     * Assumption: first tag in list is always a header and must be treated specially
     */
    t = DBgetTags(DBI(xx),seq);
    t = t->next;
    while(t != NULL) {
	/* now the hard part - update comment */
	if ( !(t->flags & TAG_COMMENT_CHANGED) &&
	    t->tagrec.comment ) {
	    read_tag(io, t->original_tag_id,&thisTagRec);
	    thisTagRec.comment = 0;
	    write_tag(io, t->original_tag_id,thisTagRec);
	}
	t = t->next;
    }
    
    
    /******************
     * PASS THE THIRD *
     ******************/
    
    nextTag = oldTagList;
    while ( nextTag ) {
	thisTag = nextTag;
	(void) read_tag(io,  thisTag , &thisTagRec );
        nextTag = thisTagRec.next;
        delete_tag_rec (io, thisTag);
    }
    
    
}


/*************************************************************
 * READ TAG LIST
 *************************************************************/

tagStruct *readTagList(DBInfo *db, int seq, int ref)
/*
 * Seq is DB_Number()
 */
{
    tagStruct *s,*t,*u;
#ifdef TAG_CHECK
    int chkpos = 0;
    int gellen = seq ? _DB_Length2(db,ref) : _DB_Length(db, ref);
#endif /*TAG_CHECK*/

    s=t=newTag();

    t->tagrec.next = first_tag(_DBI_io(db), seq);

    while (t->tagrec.next) {
        u=newTag();
	read_tag(_DBI_io(db), t->tagrec.next,&u->tagrec);
#ifdef TAG_CHECK
	if (u->tagrec.position <= 0 ||
	    u->tagrec.position + u->tagrec.length > gellen + 1)
	    /*
	     * The tag should never run off the end of the gel
	     */
	    verror(ERR_WARN, "readTagList",  "INVALID TAG POSITION "
		   "seq=%d (%s) tagpos=%d taglen=%d gellen=%d\n",
		    ref,
		    DBgetName(db,ref),
		    u->tagrec.position,
		    u->tagrec.length,
		    gellen);

	if (u->tagrec.position < chkpos)
	    /*
	     * Tag should never be out of ordeer
	     */
	    verror(ERR_WARN, "readTagList",  "TAG OUT OF ORDER "
		   "seq=%d (%s) tagpos=%d taglen=%d\n",
		    ref,
		    DBgetName(db,ref),
		    u->tagrec.position,
		    u->tagrec.length );
	chkpos = u->tagrec.position;
#endif /*TAG_CHECK*/

	u->original_tag_id = t->tagrec.next;
	t->next = u;
	t=u;
    }
    t->next = NULL;
    
    return s;
}



