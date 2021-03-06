#include <tk.h>
#include <math.h>
#include <string.h>

#include "canvas_box.h"
#include "tg_gio.h"
#include "misc.h"
#include "io_utils.h"
#include "cs-object.h"
#include "newgap_cmds.h"
#include "consensus.h"
#include "contig_selector.h"
//#include "contigEditor.h"
#include "gap_globals.h"
#include "tagdb.h"
#include "text_output.h"
//#include "tagUtils.h"
#include "active_tags.h"
//#include "template_display.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "gap4_compat.h"
#include "gap_cli_arg.h"

#include "io_lib/hash_table.h"

/* FIXME - we'll only allow one CSPlot currently */
/* It'll be initialised as NULL anyway due to it's static global nature */
HTablePtr csplot_hash[HASHMODULUS] = {0};

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void cs_callback(GapIO *io, tg_rec contig, void *fdata,
			reg_data *jdata);
int
DoClipping(GapIO *io,                                                  /* in */
	   obj_match *match);                                     /* in, out */

int
find_left_position(GapIO *io,
		   tg_rec *order,
		   double wx);

/*
 * plots the results of a search on the plot window of the contig selector
 */
void
PlotRepeats(GapIO *io,
	    mobj_repeat *repeat) {
    int i;
    char cmd[1024];
    int pos1, pos2, end1, end2;
    int64_t x1, y1, x2, y2;
    /* int max_x = 0; */
    int sense1 = 1;
    int sense2 = 1;
    int inum = 0;
    char *colour = repeat->colour;
    int width = repeat->linewidth;
    char *tag_id = repeat->tagname;
    obj_match new_match;
    int cs_id;
    obj_cs *cs;
    int64_t cur_pos = 0;
    tg_rec *order = ArrayBase(tg_rec, io->contig_order);
    HashTable *cstart;
    HashItem *hi;

    /* Build lookup table for fast contig+pos to abs_pos conversion */
    cstart = HashTableCreate(64, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    for (i = 0; i < NumContigs(io); i++) {
	HashData hd;

	//printf("Rec %"PRIrec" pos %"PRId64"\n", order[i], cur_pos);
	hd.i = cur_pos;
	HashTableAdd(cstart, (char *)&order[i], sizeof(order[i]), hd, 0);
	cur_pos += io_cclength(io, order[i]);
    }

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id);

    for (i = 0; i < repeat->num_match; i++){
	obj_match *match = (obj_match *)&repeat->match[i];
	tg_rec rec;

	/* Check if shown */
	if (match->flags & OBJ_FLAG_HIDDEN)
	    continue;

	/* clip length of match if necessary */
	new_match = *match;
	DoClipping(io, &new_match);
	/*
	 * printf("new pos1 %d pos2 %d length %d\n",
	 * new_match.pos1,  new_match.pos2, new_match.length);
	 *
	 * printf("match pos1 %d pos2 %d length %d \n",
	 * match->pos1, match->pos2, match->length);
	 */
	//pos1 = find_position_in_DB(io, ABS(new_match.c1), new_match.pos1);
	//pos2 = find_position_in_DB(io, ABS(new_match.c2), new_match.pos2);

	rec = ABS(new_match.c1);
	hi = HashTableSearch(cstart, (char *)&rec, sizeof(rec));
	if (!hi) return;
	pos1 = hi->data.i + new_match.pos1;
	end1 = hi->data.i + new_match.end1;

	rec = ABS(new_match.c2);
	hi = HashTableSearch(cstart, (char *)&rec, sizeof(rec));
	if (!hi) return;
	pos2 = hi->data.i + new_match.pos2;
	end2 = hi->data.i + new_match.end2;

	/* convert contig code back to sense ie -ve contig number means
	 * match on opposite strand
	 */
	if (new_match.c1 < 0) {
	    sense1 = -1;
	} else {
	    sense1 = 1;
	}
	if (new_match.c2 < 0) {
	    sense2 = -1;
	} else {
	    sense2 = 1;
	}

	/*
	 * draw matches of same sense (+:+ or -:-) as p1,p2
	 *                                              \
	 *                                         p1+len, p2+len
	 * draw matches of different sense (+:- or -:+) as p1+len, p2
	 *                                                     /
	 *                                                 p1, p2+len
	 */
	x1 = pos1;
	x2 = end1;
	if (sense1 == sense2) {
	    y1 = pos2;
	    y2 = end2;
	} else {
	    y1 = end2;
	    y2 = pos2;

	}

	/* need to plot in top half of screen therefore 'x' contig should be
	 * larger than the corresponding 'y' contig
	 */
	/*
	   printf("R:%d@%d,%d@%d(%d) C:%d@%d,%d@%d(%d)\n",
	   match->pos1, match->c1, match->pos2, match->c2, match->length,
	   new_match.pos1, new_match.c1,
	   new_match.pos2, new_match.c2, new_match.length);

	   printf("tag_id %s \n", tag_id);
	*/
	if (pos1 > pos2){
	    sprintf(cmd,"%s create line %"PRId64" %"PRId64" %"PRId64
		    " %"PRId64" -width %d -capstyle round -fill %s "
		    "-tags {num_%"PRIrec" num_%"PRIrec" %s S}",
		    cs->window, x1, y1, x2, y2, width, colour,
		    ABS(new_match.c1), ABS(new_match.c2), tag_id);
	    
	} else {
	    sprintf(cmd,"%s create line %"PRId64" %"PRId64" %"PRId64
		    " %"PRId64" -width %d -capstyle round "
		    "-tags \"num_%"PRIrec" num_%"PRIrec" %s S\" -fill %s",
		    cs->window, y1, x1, y2, x2, width, ABS(new_match.c1),
		    ABS(new_match.c2), tag_id, colour);
	}
	/* printf("cmd %s \n", cmd); */
	if (TCL_ERROR == Tcl_Eval(GetInterp(), cmd))
	    fprintf(stderr, "%s \n", GetInterpResult());

	inum = atoi(GetInterpResult());
	match->inum = inum;
	HashInsert(csplot_hash, inum, match);
    }

    /* scale new matches */
    scaleSingleCanvas(GetInterp(), cs->world, cs->canvas, cs->window, 'b',
		      tag_id);

    HashTableDestroy(cstart, 0);
}

/* draw the contig lines of the contig selector */
int
display_contigs(Tcl_Interp *interp,                                   /* in */
		GapIO *io,                                            /* in */
		char *win_name,                                       /* in */
		char *colour,                                         /* in */
		int width,                                            /* in */
		int tick_wd,                                          /* in */
		int tick_ht,                                          /* in */
		int64_t offset,                                       /* in */
		char *direction)                                      /* in */
{
    char cmd[1024];
    int i;
    int64_t x1 = 1;
    int64_t x2 = x1;
    int64_t y1 = 1;
    int64_t y2 = y1;

    sprintf(cmd, "%s delete all", win_name);
    Tcl_Eval(interp, cmd);

    /* draw first tick */
    if (strcmp(direction, "horizontal")==0){
	sprintf(cmd, "%s create line %"PRId64" %"PRId64" %"PRId64" %"PRId64
		" -fill %s -width %d -tags sep_1\n",
		win_name, x1, offset-tick_ht, x1, offset+tick_ht,
		colour, tick_wd);
    } else if (strcmp(direction, "vertical")==0){
	sprintf(cmd, "%s create line %"PRId64" %"PRId64" %"PRId64" %"PRId64
		" -fill %s -width %d -tags sep_1\n",
		win_name, offset-tick_ht, y1, offset+tick_ht, y1,
		colour, tick_wd);
    }
    Tcl_Eval(interp, cmd);

#ifdef DEBUG
    printf("num contigs %d \n", NumContigs(io));
    for (i = 0; i < NumContigs(io); i++ ){
	printf("i %d %d\n", i, arr(tg_rec, io->contig_order, i));
    }
#endif

    for (i = 0; i < NumContigs(io); i++){
	if (arr(tg_rec, io->contig_order, i) > 0) {
	    int clen = io_cclength(io, arr(tg_rec, io->contig_order, i));
	    if (strcmp(direction, "horizontal")==0){
		x1 = x2;
		x2 = clen + x2;
		/*
		  printf("i %d num %d length %d x1 %d x2 %d \n",
		  i, arr(tg_rec, io->contig_order, i), clen,
		  x1, x2);
		*/
		/* contig line */
		sprintf(cmd,"%s create line %"PRId64" %"PRId64" %"PRId64
			" %"PRId64" -fill %s -width %d "
			"-tags {contig c_%d num_%"PRIrec" hl_%"PRIrec" S}\n",
			win_name, x1, offset, x2, offset,
			colour, width, i+1,
			arr(tg_rec, io->contig_order, i),
			arr(tg_rec, io->contig_order, i));
	    } else if (strcmp(direction, "vertical")==0){
		y1 = y2;
		y2 = clen + y2;
		sprintf(cmd,"%s create line %"PRId64" %"PRId64" %"PRId64
			" %"PRId64" -fill %s -width %d "
			"-tags {contig c_%d num_%"PRIrec" hl_%"PRIrec" S}\n",
			win_name, offset, y1, offset, y2,
			colour, width, i+1,
			arr(tg_rec, io->contig_order, i),
			arr(tg_rec, io->contig_order, i));
	    }
	    /* printf("cmd %s \n", cmd); */
	    Tcl_Eval(interp, cmd);

	    /* Store canvas item number in an array containing contig no. */
	    {
		char aname[1024], aele[50];
		sprintf(aname, "%s.Cnum", win_name);
		sprintf(aele, "%d", i+1);
		Tcl_SetVar2(interp, aname, aele, Tcl_GetStringResult(interp),
			    TCL_GLOBAL_ONLY);
	    }

	    /* tick at end of line */
	    if (strcmp(direction, "horizontal")==0){
		sprintf(cmd, "%s create line %"PRId64" %"PRId64" %"PRId64
			" %"PRId64" -fill %s -width %d -tags sep_%d\n",
			win_name, x2, offset-tick_ht, x2, offset+tick_ht,
			colour, tick_wd, i+2);
	    } else if (strcmp(direction, "vertical")==0){
		sprintf(cmd, "%s create line %"PRId64" %"PRId64" %"PRId64
			" %"PRId64" -fill %s -width %d -tags sep_%d\n",
			win_name, offset-tick_ht, y2, offset+tick_ht, y2,
			colour, tick_wd, i+2);

	    }
	    /* printf("cmd %s \n", cmd); */
	    Tcl_Eval(interp, cmd);
	}
    }
    return TCL_OK;
}

void
update_contig_order(Tcl_Interp *interp,
		    GapIO *io,
		    int cs_id,
		    contig_list_t *contig_array,
		    int num_contigs,
		    int64_t cx)
{
    tg_rec *order = ArrayBase(tg_rec, io->contig_order);
    obj_cs *cs;
    int i, j;
    double wx, wy;
    int64_t left_position;
    char cmd[1024];
    int64_t orig_pos = 0;
    reg_buffer_start rs;
    reg_buffer_end re;
    reg_order ro;

    cs = result_data(io, cs_id);

    CanvasToWorld(cs->canvas, cx, 0, &wx, &wy);

    /*
     * returns the nth contig to the left of the wx, NOT the contig number.
     * If this is to the left of the first contig, returns 0.
     */
    left_position = find_left_position(io, order, wx);

    for (i = 0; i < NumContigs(io); i++) {
	if (order[i] == contig_array[0].contig) {
	    orig_pos = i+1;
	    break;
	}
    }

    /* convert index on order to index on contig num */
    for (i = 0; i < num_contigs; i++) {

	for (j = 0; j < NumContigs(io); j++) {
	    if (order[j] == contig_array[i].contig)
		break;
	}
	ReOrder(io, order, j, left_position);

	if (j > left_position) {
	    left_position++;
	    orig_pos++;
	}
    }

    ro.job = REG_ORDER;
    ro.pos = left_position;

#ifdef HACK
    /* HACK is there a better way of representing this - only need to
     * replot once
     */
    contig_notify(io, 1, (reg_data *)&ro);
#endif

    /* Notify of the start of the flurry of updates */
    rs.job = REG_BUFFER_START;
    for (i = 0; i < num_contigs; i++) {
	contig_notify(io, contig_array[i].contig, (reg_data *)&rs);
    }

    ro.job = REG_ORDER;
    ro.pos = left_position;

    for (i = 0; i< num_contigs; i++)
	contig_notify(io, contig_array[i].contig, (reg_data *)&ro);

    /* Notify the end of our updates */
    re.job = REG_BUFFER_END;
    for (i = 0; i < num_contigs; i++) {
	contig_notify(io, contig_array[i].contig, (reg_data *)&re);
    }

    /* draw larger separator tick to show where contig was moved from */
    sprintf(cmd, "HighlightSeparator %s %"PRId64, cs->hori, orig_pos);
    Tcl_Eval(interp, cmd);
}

void
ReOrder(GapIO *io,
	tg_rec *order,
	int c_from,
	int c_to)
{
    int length;
    tg_rec tmp;
    tmp = order[c_from];

#ifdef DEBUG
    printf("from %d to %d \n", c_from, c_to);
    printf("num_conts %d \n",  NumContigs(io));
  {
    int i;
    for (i = 0; i < NumContigs(io); i++ ){
	printf("reorder i %d %d\n", i, order[i]);
    }
  }
#endif
    /* if shifting array elements up */
    if (c_from < c_to) {
	c_to = c_to - 1;
	length = ABS(c_from - c_to);
	memmove(&order[c_from], &order[c_from]+1,
		length * sizeof(*order));
    }
    /* if shifting array elements down */
    else {
	length = ABS(c_from - c_to);
	memmove(&order[c_to]+1, &order[c_to], length * sizeof(*order));
    }
    order[c_to] = tmp;
#ifdef DEBUG
  {
    int i;
    for (i = 0; i < NumContigs(io); i++ ){
	printf("reorder i %d %d\n", i, order[i]);
    }
  }
#endif
} /* end ReOrder*/

/*
 * returns the nth contig to the left of the wx, ie the order NOT the contig
 * number.
 */
int
find_left_position(GapIO *io,
		   tg_rec *order,
		   double wx)
{

    int num_contigs;
    tg_rec cur_contig;
    int64_t length, prev_len;
    int nearest_contig;
    int i;

    length = 0;
    prev_len = 0;
    num_contigs = NumContigs(io);

    for (i = 0; i < num_contigs; i++) {

	cur_contig = order[i];
	prev_len = length;
	length += ABS(io_cclength(io, cur_contig));
#ifdef DEBUG
	printf("i %d length %d prev %d curcontig %"PRIrec" wx %f\n",
	       i, length, prev_len, cur_contig, wx);
#endif
	if (wx < length) {
	    if (ABS(wx - prev_len) >= ABS(wx - length)) {
		/* nearest length */
		nearest_contig = i+1;
	    } else {
		nearest_contig = i;
	    }
	    return nearest_contig;
	}
    }
    return num_contigs;
}

/* determines the position of a base in terms of the entire database */
int64_t
find_position_in_DB(GapIO *io,
		    tg_rec c_num,
		    int64_t position)
{
    tg_rec *order = ArrayBase(tg_rec, io->contig_order);
    int i;
    int64_t cur_length = 0;
    int64_t cur_contig;

    for (i = 0; i < NumContigs(io); i++){

	cur_contig = order[i];
	if (c_num == cur_contig) {
#ifdef DEBUG
	       printf("position %d cur_length %d c_num %d cur_contig %d\n",
	       position, cur_length, c_num, cur_contig);
#endif
	    return(cur_length + position);
	}
	/* cur_length += io_cclength(io, cur_contig) + 1; */
	cur_length += io_cclength(io, cur_contig);
    }
    return -1;
}
/*
 * find the position (in bases) of the contig selector cursor local to a contig
 */
double
CSLocalCursor(GapIO *io,
	      double wx)
{
    int i;
    int64_t offset = 0;
    int64_t prev_offset = 0;
    int num_contigs;
    tg_rec *order = ArrayBase(tg_rec, io->contig_order);
    tg_rec cur_contig;

    num_contigs = NumContigs(io);
    /*
     * a couple of fudges: if num_contigs is 1 then wx is still wx
     * if wx < 0 then must be to the left of the 1st contig and wx is still wx
     */
    if ((num_contigs == 1) || (wx < 0)) {
	return wx;
    }

    for (i = 0; i < num_contigs; i++) {
	int cstart, cend;

	cur_contig = order[i];
	prev_offset = offset;
	consensus_valid_range(io, cur_contig, &cstart, &cend);
	//offset += ABS(io_cclength(io, cur_contig));
	offset += cend - cstart + 1;
	if ((wx > prev_offset) && (wx <= offset+1)) {
	    return (wx - prev_offset + cstart);
	}
    }
    /* last contig */
    return (wx - offset);
}

int
DoClipping(GapIO *io,                                                  /* in */
	   obj_match *match)                                      /* in, out */
{
    //contig_t *c;
    int c1_start, c1_end;
    int c2_start, c2_end;

    //c = cache_search(io, GT_Contig, ABS(match->c1));
    //c1_start = c->start;
    //c1_end   = c->end;
    consensus_valid_range(io, ABS(match->c1), &c1_start, &c1_end);

    if (match->pos1 < c1_start) match->pos1 = c1_start;
    if (match->pos1 > c1_end)	match->pos1 = c1_end;
    if (match->end1 < c1_start)	match->end1 = c1_start;
    if (match->end1 > c1_end)	match->end1 = c1_end;

    match->pos1 -= c1_start-1;
    match->end1 -= c1_start-1;

    //c = cache_search(io, GT_Contig, ABS(match->c2));
    //c2_start = c->start;
    //c2_end   = c->end;
    consensus_valid_range(io, ABS(match->c2), &c2_start, &c2_end);

    if (match->pos2 < c2_start)	match->pos2 = c2_start;
    if (match->pos2 > c2_end)	match->pos2 = c2_end;
    if (match->end2 < c2_start)	match->end2 = c2_start;
    if (match->end2 > c2_end)	match->end2 = c2_end;

    match->pos2 -= c2_start-1;
    match->end2 -= c2_start-1;

    match->length = MIN(match->end1 - match->pos1,
			match->end2 - match->pos2) + 1;

    /* Set minimum size of line, just for clarity? */

    return 0;
}

void
DrawCSTags(Tcl_Interp *interp,                                       /* in */
	   int x1,
	   int x2,
	   tg_rec tag_num,
	   int tag_type,
	   int offset,
	   char *win_name,
	   int width,
	   tg_rec contig_num,
	   tg_rec read_num)
{
    char type[100];
    char *colour = tag_db[0].bg_colour;
    char cmd[1024], str[5];
    int k;

    sprintf(type, "{tag %s t_%"PRIrec" num_%"PRIrec" rnum_%"PRIrec"}",
	    type2str(tag_type,str), tag_num, contig_num, read_num);

    /* find tag colour in tag_db */
    for (k = 0; k < tag_db_count; k++){
	if (tag_type == str2type(tag_db[k].id)) {
	    colour = tag_db[k].bg_colour;
	    break;

	}
    }

    sprintf(cmd, "%s create rectangle %d %d %d %d "
	    "-fill %s -tags %s -width %d -outline %s\n",
	    win_name, x1, offset, x2 + 1, offset, colour, type, width, colour);

    Tcl_Eval(interp, cmd);
}


#if 0
/*
 * variation on vtagget. In this case, the tag_num is returned
 */
GAnnotations *
get_tag_num(GapIO *io,                                                 /* in */
	    int gel,                                                   /* in */
	    int num_t,                                                 /* in */
	    char **type,                                               /* in */
	    int *tag_num)                                             /* out */
{
    static GAnnotations a;
    static int anno;
    int arg;
    char str[5];

    if (!gel)
	anno = a.next;
    else
	if (-1 == io_read_annotation(io, gel, &anno))
	    return (GAnnotations *)-1;

    while (anno) {
	tag_read(io, anno, a);
	for (arg = 0; arg < num_t; arg++) {
	    if (idToIndex(type[arg]) == idToIndex(type2str(a.type,str))) {
		*tag_num = anno;
		return &a;
	    }
	}
	anno = a.next;
    }

    return (GAnnotations *)NULL;
}
#endif

int
display_cs_tags(Tcl_Interp *interp,                                   /* in */
		GapIO *io,                                            /* in */
		obj_cs *cs)                                           /* in */
{
    char **tag_types = NULL;
    int num_tags, i, cstart, clen;
    HashTable *ttype;
    

    /* get template display tag list */
    /* HACK - put in registration structure ? */
    if (TCL_ERROR == Tcl_VarEval(interp, "GetDefaultTags ", "CONTIG_SEL.TAGS ", NULL)) {
	printf("ERROR %s\n", Tcl_GetStringResult(interp));
    }


    /* Build a hash for fast lookup of tag type */
    if (SetActiveTags2(Tcl_GetStringResult(interp), &num_tags, &tag_types) == -1) {
	return -1;
    }

    if (num_tags == 0) {
	if (tag_types)
	    Tcl_Free((char *)tag_types);
	return 0;
    }

    ttype = HashTableCreate(64, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    for (i = 0; i < num_tags; i++) {
	int t = str2type(tag_types[i]);
	HashData hd;
	hd.i = 1;
	HashTableAdd(ttype, (char *)&t, sizeof(t), hd, 0);
    }
    if (tag_types)
	Tcl_Free((char *)tag_types);


    /* Use contig_iterator to loop through all tags on all contigs */
    for (cstart = i = 0; i < NumContigs(io); i++, cstart += clen) {
	tg_rec crec;
	rangec_t *r;
	contig_iterator *iter;
	int x1, x2;

	if (arr(tg_rec, io->contig_order, i) <= 0) {
	    clen = 0;
	    continue;
	}

	crec = arr(tg_rec, io->contig_order, i);
	clen = io_cclength(io, crec);
	iter = contig_iter_new_by_type(io, crec, 1,
				       CITER_FIRST | CITER_ISTART,
				       CITER_CSTART, CITER_CEND,
				       GRANGE_FLAG_ISANNO);

	while (NULL != (r = contig_iter_next(io, iter))) {
	    int t;
	    x1 = cstart + r->start;
	    x2 = cstart + r->end;

	    /* Check it is a type we've requested to display */
	    t = r->mqual;
	    if (!HashTableSearch(ttype, (char *)&t, sizeof(t)))
		continue;

	    if (r->flags & GRANGE_FLAG_TAG_SEQ) {
		/* read */
		DrawCSTags(interp, x1, x2, r->rec, r->mqual,
			   cs->tag.offset, cs->hori, cs->tag.width,
			   crec, r->pair_rec);
	    } else {
		/* cons */
		DrawCSTags(interp, x1, x2, r->rec, r->mqual,
			   cs->tag.offset+20, cs->hori, cs->tag.width,
			   crec, 0);
	    }
	}
	contig_iter_del(iter);
    }

    HashTableDestroy(ttype, 0);
    return 0;
}

/*
 * plot horizontal contigs
 */
void
update_contig_selector(Tcl_Interp *interp,
		       GapIO *io,
		       obj_cs *cs)
{
    int win_ht;
    char cmd[1024];

    Tcl_VarEval(interp, "winfo height ", cs->hori, NULL);
    win_ht = atoi(Tcl_GetStringResult(interp));

    display_contigs(interp, io, cs->hori, cs->line_colour, cs->line_width,
		    cs->tick->line_width, cs->tick->ht, win_ht/2, "horizontal");

    cs->world->total->x1 = 1;
    cs->world->total->x2 = CalcTotalContigLen(io);
    cs->world->total->y1 = 1;
    cs->world->total->y2 = CalcTotalContigLen(io);

    if (lengthZoom(cs->zoom) <= 1) {

	memcpy(cs->world->visible, cs->world->total, sizeof(d_box));
	SetCanvasCoords(interp, cs->world->visible->x1, cs->world->visible->y1,
			cs->world->visible->x2,cs->world->visible->y2, cs->canvas);
	/* remove all current zooming info */
	freeZoom(&cs->zoom);

	/* add first zoom */
	pushZoom(&cs->zoom, cs->world->visible);
    }

    display_cs_tags(interp, io, cs);

    scaleSingleCanvas(interp, cs->world, cs->canvas, cs->hori, 'x', "all");

    /* rehighlight selected contigs */
    sprintf(cmd, "ReHighlightContigSelection %s %s",
	    io_obj_as_string(io), cs->hori);
    Tcl_Eval(interp, cmd);
}

/*
 * plot vertical contigs and deal with diagonal line in dot plot
 */
void
update_contig_comparator(Tcl_Interp *interp,
			 GapIO *io,
			 obj_cs *cs)
{
    int win_wd;
    char cmd[1024];

    Tcl_VarEval(interp, "winfo width ", cs->vert, NULL);
    win_wd = atoi(Tcl_GetStringResult(interp));

    display_contigs(interp, io, cs->vert, cs->line_colour, cs->line_width,
		    cs->tick->line_width, cs->tick->ht, win_wd/2, "vertical");

    scaleSingleCanvas(interp, cs->world, cs->canvas, cs->vert, 'y', "all");

    sprintf(cmd, "DisplayDiagonal %s %s %s", cs->frame, cs->window,
	    io_obj_as_string(io));
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	printf("update_contig_comparator: %s\n", Tcl_GetStringResult(interp));

}

int
contig_comparator_reg(Tcl_Interp *interp,
		      GapIO *io,
		      obj_cs *cs,
		      char *csp_win,
		      char *csv_win)
{
    int id;
    id = register_id();

    strcpy(cs->hori, cs->window);
    strcpy(cs->vert, csv_win);
    strcpy(cs->window, csp_win);

    /* create list of windows in the contig selector display */
    deleteWindow(cs->win_list, &cs->num_wins, cs->win_list[0]->window);
    addWindow(cs->win_list, &cs->num_wins, cs->window, 'b', id);
    addWindow(cs->win_list, &cs->num_wins, cs->hori,   'x', id);
    addWindow(cs->win_list, &cs->num_wins, cs->vert,   'y', id);

    /* do the first plot */
    update_contig_comparator(interp, io, cs);

/*
    for (i = 1; i <= NumContigs(io); i++) {
	contig_register(io, i, cs_callback, (void *)cs, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ANNO | REG_GENERIC |
			REG_BUFFER | REG_FLAG_INVIS, REG_TYPE_CONTIGSEL);
    }
*/
    return id;
}

int
contig_selector_reg(Tcl_Interp *interp,
		    GapIO *io,
		    char *frame,
		    char *csh_win,
		    tag_s tag,
		    cursor_s cursor,
		    tick_s *tick)
{
    obj_cs *cs;
    int id;

    if (NULL == (cs = (obj_cs *)xmalloc(sizeof(obj_cs))))
	return 0;

    id = register_id();

    cs->line_width = get_default_int(interp, gap5_defs,
				     "CONTIG_SEL.LINE_WIDTH");

    cs->line_colour = get_default_astring(interp, gap5_defs,
					  "CONTIG_SEL.COLOUR");

/*
    cs->tag = tag_struct(interp, "CONTIG_SEL");
    cs->cursor = cursor_struct(interp, "CONTIG_SEL");
    cs->tick = tick_struct(interp, "CONTIG_SEL");
*/
    cs->tag = tag;
    cs->cursor = cursor;
    cs->tick = tick;
    cs->buffer_count = 0;
    cs->do_update = 0;

    cs->vert[0] = '\0';
    strcpy(cs->frame, frame);
    strcpy(cs->window, csh_win);
    strcpy(cs->hori, cs->window);

    /* create list of windows in the contig selector display */
    if (NULL == (cs->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;
    cs->num_wins = 0;
    addWindow(cs->win_list, &cs->num_wins, cs->window, 'x', id);

    if (NULL == (cs->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (cs->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (cs->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (cs->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    initCanvas(interp, cs->canvas, cs->window);
    createZoom(&cs->zoom);

    /* do the first plot */
    update_contig_selector(interp, io, cs);

    contig_register(io, 0, cs_callback, (void *)cs, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
		    REG_NUMBER_CHANGE | REG_ANNO | REG_GENERIC |
		    REG_BUFFER | REG_ORDER | REG_FLAG_INVIS, REG_TYPE_CONTIGSEL);

    return id;

}

/*
 * Removes the template and reading display (and unplots etc).
 */
static void cs_shutdown(GapIO *io, obj_cs *cs) {
    reg_quit rq;

/*
    for (i = 1; i <= NumContigs(io); i++) {
	contig_deregister(io, i, cs_callback, (void *)cs);
    }
*/
    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;

    type_notify(io, REG_TYPE_FIJ,      (reg_data *)&rq);
    type_notify(io, REG_TYPE_READPAIR, (reg_data *)&rq);
    type_notify(io, REG_TYPE_REPEAT,   (reg_data *)&rq);
    type_notify(io, REG_TYPE_CHECKASS, (reg_data *)&rq);
    type_notify(io, REG_TYPE_OLIGO,    (reg_data *)&rq);

    /*
     * need to deregister AFTER done type_notify requests because they need
     * the cs data structure which is deleted during deregistration
     */
    //    for (i = 1; i <= NumContigs(io); i++) {
    //	contig_deregister(io, i, cs_callback, (void *)cs);
    //    }
    contig_deregister(io, 0, cs_callback, (void *)cs);

    if (TCL_ERROR == Tcl_VarEval(GetInterp(), "DeleteContigSelector ",
				 cs->frame, NULL)) {
	printf("cs_shutdown %s\n", GetInterpResult());
    }

    free_win_list(cs->win_list, cs->num_wins);

    xfree(cs->line_colour);
    xfree(cs->canvas);
    xfree(cs->world->visible);
    xfree(cs->world->total);
    xfree(cs->world);
    if (cs->cursor.colour) free(cs->cursor.colour);
    if (cs->tick->colour) free(cs->tick->colour);
    xfree(cs->tick);
    freeZoom(&cs->zoom);
    xfree(cs);
}

void
cs_callback(GapIO *io, tg_rec contig, void *fdata, reg_data *jdata) {
    char cmd[1024];
    obj_cs *cs = (obj_cs *)fdata;

    switch(jdata->job) {
    case REG_BUFFER_START:
	{
#ifdef DEBUG
	    printf("REG_BUFFER_START count %d \n", cs->buffer_count);
#endif
	    cs->buffer_count++;
	    cs->do_update = REG_BUFFER_START;
	    return;
	}
    case REG_BUFFER_END:
	{
#ifdef DEBUG
	    printf("REG_BUFFER_END count %d \n", cs->buffer_count);
#endif
	    cs->buffer_count--;
	    if (cs->buffer_count <= 0) {
		cs->buffer_count = 0;
		if (cs->do_update & REG_LENGTH) {

		} else if (cs->do_update & REG_ANNO) {
		    Tcl_VarEval(GetInterp(), cs->hori, " delete tag", NULL);
		    display_cs_tags(GetInterp(), io, cs);
		    scaleSingleCanvas(GetInterp(), cs->world, cs->canvas,
				      cs->hori, 'x', "tag");
		} else if (cs->do_update & REG_ORDER) {
		    update_contig_selector(GetInterp(), io, cs);
		    if (cs->vert[0] != '\0') {
			update_contig_comparator(GetInterp(), io, cs);
		    }
		}
		cs->do_update = 0;
	    }
	    return;
	}
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Contig selector");
	    return;
	}

    case REG_GET_OPS:
	{
	    /* jdata->get_ops.ops = "Information\0Configure\0"; */
	    return;
	}
    case REG_ANNO:
	{
#ifdef DEBUG
	    printf("contig selector REG_ANNO\n");
#endif
	    if (!cs->do_update) {
		Tcl_VarEval(GetInterp(), cs->hori, " delete tag", NULL);
		display_cs_tags(GetInterp(), io, cs);
		scaleSingleCanvas(GetInterp(), cs->world, cs->canvas,
				  cs->hori, 'x', "tag");
	    } else {
		cs->do_update |= REG_ANNO;
	    }
	    return;
	}
    case REG_ORDER:
	{

#ifdef DEBUG
	    printf("contig selector REG_ORDER %d\n", cs->buffer_count);
#endif
	    if (!cs->do_update) {

		update_contig_selector(GetInterp(), io, cs);
		if (cs->vert[0] != '\0') {
		    update_contig_comparator(GetInterp(), io, cs);
		}
	    } else {
		cs->do_update |= REG_ORDER;
	    }
	    break;
	}
    case REG_QUIT:
	{
	    cs_shutdown(io, cs);
	    return;
	}

    case REG_GENERIC:
	switch(jdata->generic.task) {

	case TASK_WINDOW_ADD:
	    {
		win *winfo = (win *)jdata->generic.data;

		addWindow(cs->win_list, &cs->num_wins, winfo->window,
			  winfo->scroll, winfo->id);
		break;
	    }
	case TASK_WINDOW_DELETE:
	    {
		char *window = (char *)jdata->generic.data;

		deleteWindow(cs->win_list, &cs->num_wins, window);
		break;
	    }
	case TASK_CANVAS_SCROLLX:
	    {
		char *scroll = (char *)jdata->generic.data;

		canvasScrollX(GetInterp(), cs->window, cs->win_list,
			      cs->num_wins, cs->world->visible, cs->canvas,
			      scroll);
		break;
	    }
	case TASK_CANVAS_SCROLLY:
	    {
		char *scroll = (char *)jdata->generic.data;

		canvasScrollY(GetInterp(), cs->window, cs->win_list,
			      cs->num_wins, cs->world->visible, cs->canvas,
			      scroll);

		break;
	    }
	case TASK_CANVAS_RESIZE:
	    {
		char scroll_args[20];
		/* resize template display window */
		resizeCanvas(GetInterp(), cs->window, cs->win_list,
			     cs->num_wins, cs->world->visible,
			     cs->world->total, cs->canvas);
		sprintf(scroll_args, "scroll 0 units");
		canvasScrollX(GetInterp(), cs->window, cs->win_list,
			      cs->num_wins, cs->world->visible, cs->canvas,
			      scroll_args);

		break;
	    }
	case TASK_CANVAS_ZOOMBACK:
	    {

		if (lengthZoom(cs->zoom) <= 2) {
		    freeZoom(&cs->zoom);
		    pushZoom(&cs->zoom, cs->world->total);
		}

		canvasZoomback(GetInterp(), cs->canvas, cs->window, cs->world,
			       cs->win_list, cs->num_wins, &cs->zoom);

		break;
		}
	case TASK_CANVAS_ZOOM:
	    {
		s_zoom *szoom = (s_zoom *)jdata->generic.data;
		canvasZoom(GetInterp(), cs->canvas, cs->window, cs->world,
			   cs->win_list, cs->num_wins, &cs->zoom, szoom->zoom,
			   szoom->scroll);

		break;
	    }
	case TASK_CANVAS_CURSOR_X:
	    {
		char *label;
		int *cx = (int *)jdata->generic.data;
		double local_pos;
		double wx, wy;

		CanvasToWorld(cs->canvas, *cx, 0, &wx, &wy);

		label = get_default_string(GetInterp(), gap5_defs,
					   "CONTIG_SEL.CURSOR1_X");

		canvasCursorX(GetInterp(), cs->canvas, cs->frame, label,
			      cs->cursor.colour, cs->cursor.width, *cx, wx,
			      cs->win_list, cs->num_wins);

		/* fill in local position of cursor in label box */
		local_pos = CSLocalCursor(io, wx);

		label = get_default_string(GetInterp(), gap5_defs,
					   "CONTIG_SEL.CURSOR2_X");

		sprintf(cmd, "%s%s configure -text %d\n", cs->frame, label,
			(int)local_pos);
		Tcl_Eval(GetInterp(), cmd);

		break;
	    }
	case TASK_CANVAS_CURSOR_Y:
	    {
		char *label;
		int *cy = (int *)jdata->generic.data;
		double local_pos;
		double wx, wy;
		double cx1, cy1;

		CanvasToWorld(cs->canvas, 0, *cy, &wx, &wy);
		WorldToCanvas(cs->canvas, wy, 0, &cx1, &cy1);

		label = get_default_string(GetInterp(), gap5_defs,
					   "CONTIG_SEL.CURSOR1_Y");

		canvasCursorY(GetInterp(), cs->canvas, cs->frame, label,
			      cs->cursor.colour, cs->cursor.width, *cy, wy,
			      cs->win_list, cs->num_wins);

		sprintf(cmd, "DrawCanvasCursorX1 %s %s %.20f %s %d\n",
			cs->frame, cs->hori, cx1, cs->cursor.colour,
			cs->cursor.width);
		if (TCL_ERROR == Tcl_Eval(GetInterp(), cmd))
		    printf("%s\n", GetInterpResult());

		/* fill in local position of cursor in label box */
		local_pos = CSLocalCursor(io, wy);

		label = get_default_string(GetInterp(), gap5_defs,
					   "CONTIG_SEL.CURSOR2_Y");

		sprintf(cmd, "%s%s configure -text %d\n", cs->frame, label,
			(int)local_pos);
		Tcl_Eval(GetInterp(), cmd);

		break;
	    }
	case TASK_CANVAS_CURSOR_DELETE:
	    {
		int i;
		for (i = 0; i < cs->num_wins; i++) {
		    Tcl_VarEval(GetInterp(), cs->win_list[i]->window,
				" delete cursor_x cursor_x1 cursor_y", NULL);
		}
		break;
	    }
#if 0
	case TASK_CS_REDRAW:
	    {
		/* HACK - never used */
		int i, id = register_id();

		contig_deregister(io, 0, cs_callback, fdata);
		contig_register(io, 0, cs_callback, fdata, id,
				REG_REQUIRED |
				REG_DATA_CHANGE |
				REG_OPS |
				REG_NUMBER_CHANGE |
				REG_ANNO |
				REG_GENERIC |
				REG_FLAG_INVIS |
				REG_BUFFER,
				REG_TYPE_CONTIGSEL);
		break;
	    }
#endif
	    break;
	}
	break;
    case REG_JOIN_TO:
    case REG_LENGTH:
    case REG_DELETE:
    case REG_COMPLEMENT:
    case REG_NUMBER_CHANGE:
#ifdef DEBUG
	printf("contig selector REG_REDRAW %d\n", cs->buffer_count);
#endif
	update_contig_selector(GetInterp(), io, cs);
	if (cs->vert[0] != '\0') {
	    update_contig_comparator(GetInterp(), io, cs);
	}
	/* update tcl globals, CurContig, LREG and RREG */
	sprintf(cmd, "ContigParams %s", io_obj_as_string(io));
	Tcl_Eval(GetInterp(), cmd);


#ifdef HACK
	printf("COM %s \n", cs->com);
	if (cs->buffer_count) {
	    cs->do_update = 1;
	} else {
	    Tcl_Eval(cs->interp, cs->com);
	}
#endif
	break;
    }
}


/*
 * Delete all contig comparator displays.
 */
typedef struct {
    GapIO *io;
    int id;
} rd_arg;

int tk_clear_cp(ClientData clientData, Tcl_Interp *interp,
		int objc, Tcl_Obj *CONST objv[]) {
    obj_cs *cs;
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_quit rq;


    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;

    type_notify(args.io, REG_TYPE_FIJ,      (reg_data *)&rq);
    type_notify(args.io, REG_TYPE_READPAIR, (reg_data *)&rq);
    type_notify(args.io, REG_TYPE_REPEAT,   (reg_data *)&rq);
    type_notify(args.io, REG_TYPE_CHECKASS, (reg_data *)&rq);
    type_notify(args.io, REG_TYPE_OLIGO,    (reg_data *)&rq);

    cs = result_data(args.io, args.id);
    strcpy(cs->window, cs->hori);
    cs->vert[0] = '\0';
    return TCL_OK;
}

typedef struct {
    char *result;
    char *colour;
    char *csplot;
    int width;
} conf_arg;

/*
 * Configures a result.
 */
int tk_matchresult_configure(ClientData clientData, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    conf_arg args;
    cli_args a[] = {
	{"-result",  ARG_STR, 1, NULL, offsetof(conf_arg, result)},
	{"-colour",   ARG_STR, 1, "",   offsetof(conf_arg, colour)},
	{"-width",    ARG_INT, 1, "-1", offsetof(conf_arg, width)},
	{"-csplot",   ARG_STR, 1, NULL, offsetof(conf_arg, csplot)},
        {NULL,        0,       0, NULL, 0}
    };
    mobj_repeat *r;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    /*
     * Find the result associated with this tag. 'result' is a C pointer
     * held in an implementation defined manner as a tcl string.
     */
    r = (mobj_repeat *)TclPtr2C(args.result);

    /* Adjust configurations */
    if (*args.colour) {
	strncpy(r->colour, args.colour, COLOUR_LEN-1);
    }

    if (args.width != -1) {
	r->linewidth = args.width;
    }

    return TCL_OK;
}
