#include <tk.h>
#include <math.h>
#include <string.h>

#include "canvas_box.h"
#include "IO.h"
#include "misc.h"
#include "io_utils.h"
#include "io-reg.h"
#include "cs-object.h"
#include "newgap_cmds.h"
#include "contig_selector.h"
#include "complement.h"
#include "contigEditor.h"
#include "gap_globals.h"
#include "tagdb.h"
#include "text_output.h"
#include "tagUtils.h"
#include "active_tags.h"
#include "template_display.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"

/* FIXME - we'll only allow one CSPlot currently */
/* It'll be initialised as NULL anyway due to it's static global nature */
HTablePtr csplot_hash[HASHMODULUS] = {0};

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void cs_callback(GapIO *io, int contig, void *fdata,
			reg_data *jdata);
int
DoClipping(GapIO *io,                                                  /* in */
	   obj_match *match);                                     /* in, out */

int
find_left_position(GapIO *io,
		   GCardinal *order,
		   double wx);

/*
 * plots the results of a search on the plot window of the contig selector
 */
void
PlotRepeats(GapIO *io,
	    mobj_repeat *repeat) {
    int i;
    char cmd[1024];
    int pos1, pos2;
    int x1, y1, x2, y2;
    /* int max_x = 0; */
    int sense1 = 1;
    int sense2 = 1;
    int inum;
    char *colour = repeat->colour;
    int width = repeat->linewidth;
    char *tag_id = repeat->tagname;
    obj_match new_match;
    int cs_id;
    obj_cs *cs;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id, 0);

    for (i = 0; i < repeat->num_match; i++){
	obj_match *match = (obj_match *)&repeat->match[i];

	/* Check if shown */
	if (match->flags & OBJ_FLAG_HIDDEN)
	    continue;

	/* clip length of match if necessary */
	new_match = *match;
	DoClipping(io, &new_match);
	/*
	* printf("new pos1 %d pos2 %d length %d\n",
	* new_match.pos1,  new_match.pos2, new_match.length);

	* printf("match pos1 %d pos2 %d length %d \n",
	* match->pos1, match->pos2, match->length);
	*/
	pos1 = find_position_in_DB(io, abs(new_match.c1), new_match.pos1);
	pos2 = find_position_in_DB(io, abs(new_match.c2), new_match.pos2);

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
	x2 = pos1 + new_match.length;
	if (sense1 == sense2) {
	    y1 = pos2;
	    y2 = pos2 + new_match.length;
	} else {
	    y1 = pos2 + new_match.length;
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
	    sprintf(cmd,"%s create line %d %d %d %d -width %d -capstyle round "
		    "-tags {num_%d num_%d %s S} -fill %s",
		    cs->window, x1, y1, x2, y2, width, abs(new_match.c1),
		    abs(new_match.c2), tag_id, colour);
	} else {
	    sprintf(cmd,"%s create line %d %d %d %d -width %d -capstyle round "
		    "-tags \"num_%d num_%d %s S\" -fill %s",
		    cs->window, y1, x1, y2, x2, width, abs(new_match.c1),
		    abs(new_match.c2), tag_id, colour);
	}
	/* printf("cmd %s \n", cmd); */
	if (TCL_ERROR == Tcl_Eval(GetInterp(), cmd))
	    printf("%s \n", GetInterp()->result);

	inum = atoi(GetInterp()->result);
	match->inum = inum;
	HashInsert(csplot_hash, inum, match);
    }

    /* scale new matches */
    scaleSingleCanvas(GetInterp(), cs->world, cs->canvas, cs->window, 'b',
		      tag_id);
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
		int offset,                                           /* in */
		char *direction)                                      /* in */
{
    char cmd[1024];
    int i;
    int x1 = 1;
    int x2 = x1;
    int y1 = 1;
    int y2 = y1;

    sprintf(cmd, "%s delete all", win_name);
    Tcl_Eval(interp, cmd);

    /* draw first tick */
    if (strcmp(direction, "horizontal")==0){
	sprintf(cmd, "%s create line %d %d %d %d "
		"-fill %s -width %d -tags sep_1\n",
		win_name, x1, offset-tick_ht, x1, offset+tick_ht,
		colour, tick_wd);
    } else if (strcmp(direction, "vertical")==0){
	sprintf(cmd, "%s create line %d %d %d %d "
		"-fill %s -width %d -tags sep_1\n",
		win_name, offset-tick_ht, y1, offset+tick_ht, y1,
		colour, tick_wd);
    }
    /* printf("cmd %s \n", cmd); */
    Tcl_Eval(interp, cmd);

#ifdef DEBUG
    printf("num contigs %d \n", NumContigs(io));
    for (i = 0; i < NumContigs(io); i++ ){
	printf("i %d %d\n", i, arr(GCardinal, io->contig_order, i));
    }
#endif

    for (i = 0; i < NumContigs(io); i++){
	if (arr(GCardinal, io->contig_order, i) > 0) {
	    int clen = io_clength(io, arr(GCardinal, io->contig_order, i));
	    if (strcmp(direction, "horizontal")==0){
		x1 = x2;
		x2 = clen + x2;
		/*
		  printf("i %d num %d length %d x1 %d x2 %d \n",
		  i, arr(GCardinal, io->contig_order, i), clen,
		  x1, x2);
		*/
		/* contig line */
		sprintf(cmd,"%s create line %d %d %d %d "
			"-fill %s -width %d "
			"-tags {contig c_%d num_%d hl_%d S}\n",
			win_name, x1, offset, x2, offset,
			colour, width, i+1,
			arr(GCardinal, io->contig_order, i),
			arr(GCardinal, io->contig_order, i));
	    } else if (strcmp(direction, "vertical")==0){
		y1 = y2;
		y2 = clen + y2;
		sprintf(cmd,"%s create line %d %d %d %d "
			"-fill %s -width %d "
			"-tags {contig c_%d num_%d hl_%d S}\n",
			win_name, offset, y1, offset, y2,
			colour, width, i+1,
			arr(GCardinal, io->contig_order, i),
			arr(GCardinal, io->contig_order, i));
	    }
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
		sprintf(cmd, "%s create line %d %d %d %d "
			"-fill %s -width %d -tags sep_%d\n",
			win_name, x2, offset-tick_ht, x2, offset+tick_ht,
			colour, tick_wd, i+2);
	    } else if (strcmp(direction, "vertical")==0){
		sprintf(cmd, "%s create line %d %d %d %d "
			"-fill %s -width %d -tags sep_%d\n",
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
		    int *contig_array,
		    int num_contigs,
		    int cx)
{
    GCardinal *order = ArrayBase(GCardinal, io->contig_order);
    obj_cs *cs;
    int i, j;
    double wx, wy;
    int left_position;
    char cmd[1024];
    int orig_pos;
    reg_buffer_start rs;
    reg_buffer_end re;
    reg_order ro;

    cs = result_data(io, cs_id, 0);

    CanvasToWorld(cs->canvas, cx, 0, &wx, &wy);

    /*
     * returns the nth contig to the left of the wx, NOT the contig number.
     * If this is to the left of the first contig, returns 0.
     */
    left_position = find_left_position(io, order, wx);

    for (i = 0; i < NumContigs(io); i++) {
	if (order[i] == contig_array[0]) {
	    orig_pos = i+1;
	    break;
	}
    }

    /* convert index on order to index on contig num */
    for (i = 0; i < num_contigs; i++) {

	for (j = 0; j < NumContigs(io); j++) {
	    if (order[j] == contig_array[i])
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
	contig_notify(io, contig_array[i], (reg_data *)&rs);
    }

    ro.job = REG_ORDER;
    ro.pos = left_position;

    for (i = 0; i< num_contigs; i++)
	contig_notify(io, contig_array[i], (reg_data *)&ro);

    /* Notify the end of our updates */
    re.job = REG_BUFFER_END;
    for (i = 0; i < num_contigs; i++) {
	contig_notify(io, contig_array[i], (reg_data *)&re);
    }

    /* draw larger separator tick to show where contig was moved from */
    sprintf(cmd, "HighlightSeparator %s %d", cs->hori, orig_pos);
    Tcl_Eval(interp, cmd);
}

void
ReOrder(GapIO *io,
	GCardinal *order,
	int c_from,
	int c_to)
{
    int length;
    GCardinal tmp;
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
	length = abs(c_from - c_to);
	memmove(&order[c_from], &order[c_from]+1,
		length * sizeof(int));
    }
    /* if shifting array elements down */
    else {
	length = abs(c_from - c_to);
	memmove(&order[c_to]+1, &order[c_to], length * sizeof(int));
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
		   GCardinal *order,
		   double wx)
{

    int num_contigs;
    int cur_contig;
    int length, prev_len;
    int nearest_contig;
    int i;

    length = 0;
    prev_len = 0;
    num_contigs = NumContigs(io);

    for (i = 0; i < num_contigs; i++) {

	cur_contig = order[i];
	prev_len = length;
	length += ABS(io_clength(io, cur_contig));
#ifdef DEBUG
	printf("i %d length %d prev %d curcontig %d wx %f\n", i, length,
	       prev_len, cur_contig, wx);
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
int
find_position_in_DB(GapIO *io,
		    int c_num,
		    int position)
{
    GCardinal *order = ArrayBase(GCardinal, io->contig_order);
    int i;
    int cur_length = 0;
    int cur_contig;

    for (i = 0; i < NumContigs(io); i++){

	cur_contig = order[i];
	if (c_num == cur_contig) {
#ifdef DEBUG
	       printf("position %d cur_length %d c_num %d cur_contig %d\n",
	       position, cur_length, c_num, cur_contig);
#endif
	    return(cur_length + position);
	}
	/* cur_length += io_clength(io, cur_contig) + 1; */
	cur_length += io_clength(io, cur_contig);
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
    int offset = 0;
    int prev_offset = 0;
    int num_contigs;
    GCardinal *order = ArrayBase(GCardinal, io->contig_order);
    int cur_contig;

    num_contigs = NumContigs(io);
    /*
     * a couple of fudges: if num_contigs is 1 then wx is still wx
     * if wx < 0 then must be to the left of the 1st contig and wx is still wx
     */
    if ((num_contigs == 1) || (wx < 0)) {
	return wx;
    }

    for (i = 0; i < num_contigs; i++) {

	cur_contig = order[i];
	prev_offset = offset;
	offset += ABS(io_clength(io, cur_contig));
	if ((wx > prev_offset) && (wx <= offset+1)) {
	    return (wx - prev_offset);
	}
    }
    /* last contig */
    return (wx - offset);
}

int
DoClipping(GapIO *io,                                                  /* in */
	   obj_match *match)                                      /* in, out */

{
    int length[4];
    int min_length = INT_MAX;
    int c1_len, c2_len, i;

    length[0] = match->length;
    length[1] = length[0];
    length[2] = length[0];
    length[3] = length[0];

    if (match->pos1 <= 0) {
	length[0] = match->length + match->pos1 - 1;
	if (length[0] < 1) length[0] = 1;
	match->pos1 = 1;
    }
    if (match->pos2 <= 0) {
	length[1] = match->length + match->pos2 - 1;
	if (length[1] < 1) length[1] = 1;
	match->pos2 = 1;
    }

    /* clip the length of contig 1 if necessary */
    c1_len = io_clength(io, abs(match->c1));
    if (match->pos1 + match->length > c1_len) {
	length[2] = c1_len - match->pos1;
	if (length[2] < 1) length[2] = 1;
    if (match->pos1 > c1_len)
	match->pos1 = c1_len;
    }

    /* clip the length of contig 1 if necessary */
    c2_len = io_clength(io, abs(match->c2));
    if (match->pos2 + match->length > c2_len) {
	length[3] = c2_len - match->pos2;
	if (length[3] < 1) length[3] = 1;
    if (match->pos2 > c2_len)
	match->pos2 = c2_len;
    }

    /* select smallest of clipped lengths */
    for (i = 0; i < 4; i++) {
	if (length[i] < min_length) {
	    min_length = length[i];
	}
    }

    match->length = min_length;

    return 0;
}

void
DrawCSTags(Tcl_Interp *interp,                                       /* in */
	   int x1,
	   int x2,
	   int tag_num,
	   GAnnotations *annotation,
	   int offset,
	   char *win_name,
	   int width,
	   int contig_num,
	   int read_num)
{
    char type[100];
    char *colour = tag_db[0].bg_colour;
    char cmd[1024], str[5];
    int k;                                                       /* counter */

    sprintf(type, "{tag %s t_%d num_%d rnum_%d}",
	    type2str(annotation->type,str), tag_num, contig_num, read_num);

    /* find tag colour in tag_db */
    for (k = 0; k < tag_db_count; k++){

	if (annotation->type == str2type(tag_db[k].id)) {

/*
	    sprintf(type, "{tag %s t_%d num_%d}", tag_db[k].search_id,
		    tag_num, contig_num);
*/
	    colour = tag_db[k].bg_colour;
	    break;

	}  /* end if */

    } /* end for */


    sprintf(cmd, "%s create rectangle %d %d %d %d "
	    "-fill %s -tags %s -width %d -outline %s\n",
	    win_name, x1, offset, x2 + 1, offset, colour, type, width, colour);

    Tcl_Eval(interp, cmd);

    /* printf("cmd %s \n", cmd); */
}

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


int
display_cs_tags(Tcl_Interp *interp,                                   /* in */
		GapIO *io,                                            /* in */
		obj_cs *cs)                                           /* in */
{
    GAnnotations *annotation;                       /* annotation structure */
    char **tag_types = NULL;
    int num_tags;
    int tag_pos;                                         /* tag LH position */
    int whole_reading = 0;                /* if TRUE display entire reading */
    int r_pos, r_len;                     /* reading LH position and length */
    int c_num;
    int r_num;
    GContigs contig;
    GReadings reading;
    int x1, x2;
    int tag_num = 0;

   /* get template display tag list */
    /* HACK - put in registration structure ? */
    if (TCL_ERROR == Tcl_VarEval(interp, "GetDefaultTags ", "CONTIG_SEL.TAGS ", NULL)) {
	printf("ERROR %s\n", interp->result);
    }

    if (SetActiveTags2(interp->result, &num_tags, &tag_types) == -1) {
	return -1;
    }

    if (num_tags == 0) {
	if (tag_types)
	    Tcl_Free((char *)tag_types);
	return 0;
    }

    for (c_num = 1; c_num <= NumContigs(io); c_num++) {

	/* reading tags */
	contig_read(io, c_num, contig);

	for (r_num=contig.left; r_num; r_num=reading.right) {

	    gel_read(io, r_num, reading);
	    /* like vtagget except also returns tag_num */
	    annotation = get_tag_num(io, r_num, num_tags,
				     tag_types, &tag_num);

	    while (annotation && annotation != (GAnnotations *)-1){

		/* if reading has been complemented, find new pos of tag */
		if (reading.sense) {

		    tag_pos = (reading.position - reading.start) +
			(reading.length - annotation->position -
			 annotation->length + 1);
		    tag_pos = find_position_in_DB(io, c_num, tag_pos);

		} else {

		    tag_pos = annotation->position - (reading.start -
						      reading.position +1);
		    tag_pos = find_position_in_DB(io, c_num, tag_pos);

		} /* end if */

	       /*
		   printf("num %d tag_pos %d len %d \n",
		   c_num, tag_pos, annotation->length);
		 */
		SetReadingPosLen(whole_reading, io, r_num, &r_pos, &r_len);
		r_pos = find_position_in_DB(io, c_num, r_pos);
		CalcXCoords(tag_pos, annotation->length, &x1, &x2);

		/* clip tag at cutoff data */
		x1 = MAX(x1, r_pos);
		x2 = MIN(x2, r_pos + r_len - 1);
		if (x2 >= x1) {
		    DrawCSTags(interp, x1, x2, tag_num, annotation,
			       cs->tag.offset,
			       cs->hori, cs->tag.width, c_num, r_num);
		}
		annotation = get_tag_num(io, 0, num_tags, tag_types, &tag_num);
	    } /* end while */
	} /* end for each reading */

	/* consensus tags */
	annotation = get_tag_num(io, -c_num, num_tags, tag_types, &tag_num);

	while (annotation && annotation != (GAnnotations *)-1){

	    tag_pos = annotation->position;
	    tag_pos = find_position_in_DB(io, c_num, tag_pos);

	    CalcXCoords(tag_pos, annotation->length, &x1, &x2);

	    DrawCSTags(interp, x1, x2, tag_num, annotation, cs->tag.offset+20,
		       cs->hori, cs->tag.width, c_num, 0);
	    annotation = get_tag_num(io, 0, num_tags, tag_types, &tag_num);

	} /* end while */

    } /* end for each contig */

    if (tag_types)
	Tcl_Free((char *)tag_types);
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
    win_ht = atoi(interp->result);

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
    sprintf(cmd, "ReHighlightContigSelection %d %s", *handle_io(io), cs->hori);
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
    win_wd = atoi(interp->result);

    display_contigs(interp, io, cs->vert, cs->line_colour, cs->line_width,
		    cs->tick->line_width, cs->tick->ht, win_wd/2, "vertical");

    scaleSingleCanvas(interp, cs->world, cs->canvas, cs->vert, 'y', "all");

    sprintf(cmd, "DisplayDiagonal %s %s %d", cs->frame, cs->window,
	    *handle_io(io));
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	printf("update_contig_comparator: %s\n", interp->result);

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
    int i;

    if (NULL == (cs = (obj_cs *)xmalloc(sizeof(obj_cs))))
	return 0;

    id = register_id();

    cs->line_width = get_default_int(interp, gap_defs,
				     "CONTIG_SEL.LINE_WIDTH");

    cs->line_colour = get_default_astring(interp, gap_defs,
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

    for (i = 1; i <= NumContigs(io); i++) {
	contig_register(io, i, cs_callback, (void *)cs, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ANNO | REG_GENERIC |
			REG_BUFFER | REG_ORDER | REG_FLAG_INVIS, REG_TYPE_CONTIGSEL);
    }
    return id;

}

/*
 * Removes the template and reading display (and unplots etc).
 */
static void cs_shutdown(GapIO *io, obj_cs *cs) {
    int i;
    reg_quit rq;

/*
    for (i = 1; i <= NumContigs(io); i++) {
	contig_deregister(io, i, cs_callback, (void *)cs);
    }
*/
    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;

    type_notify(io, REG_TYPE_FIJ,      (reg_data *)&rq, 1);
    type_notify(io, REG_TYPE_READPAIR, (reg_data *)&rq, 1);
    type_notify(io, REG_TYPE_REPEAT,   (reg_data *)&rq, 1);
    type_notify(io, REG_TYPE_CHECKASS, (reg_data *)&rq, 1);
    type_notify(io, REG_TYPE_OLIGO,    (reg_data *)&rq, 1);

    /*
     * need to deregister AFTER done type_notify requests because they need
     * the cs data structure which is deleted during deregistration
     */
    for (i = 1; i <= NumContigs(io); i++) {
	contig_deregister(io, i, cs_callback, (void *)cs);
    }

    if (TCL_ERROR == Tcl_VarEval(GetInterp(), "DeleteContigSelector ",
				 cs->frame, NULL)) {
	printf("cs_shutdown %s\n", GetInterp()->result);
    }

    free_win_list(cs->win_list, cs->num_wins);

    xfree(cs->line_colour);
    xfree(cs->canvas);
    xfree(cs->world->visible);
    xfree(cs->world->total);
    xfree(cs->world);
    if (cs->cursor.colour) free(cs->cursor.colour);
    if (cs->tick->colour) free(cs->tick->colour);
    freeZoom(&cs->zoom);
    xfree(cs);
}

void
cs_callback(GapIO *io, int contig, void *fdata, reg_data *jdata) {
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

		label = get_default_string(GetInterp(), gap_defs,
					   "CONTIG_SEL.CURSOR1_X");

		canvasCursorX(GetInterp(), cs->canvas, cs->frame, label,
			      cs->cursor.colour, cs->cursor.width, *cx, wx,
			      cs->win_list, cs->num_wins);

		/* fill in local position of cursor in label box */
		local_pos = CSLocalCursor(io, wx);

		label = get_default_string(GetInterp(), gap_defs,
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
		char cmd[1024];
		double cx1, cy1;

		CanvasToWorld(cs->canvas, 0, *cy, &wx, &wy);
		WorldToCanvas(cs->canvas, wy, 0, &cx1, &cy1);

		label = get_default_string(GetInterp(), gap_defs,
					   "CONTIG_SEL.CURSOR1_Y");

		canvasCursorY(GetInterp(), cs->canvas, cs->frame, label,
			      cs->cursor.colour, cs->cursor.width, *cy, wy,
			      cs->win_list, cs->num_wins);

		sprintf(cmd, "DrawCanvasCursorX1 %s %s %.20f %s %d\n",
			cs->frame, cs->hori, cx1, cs->cursor.colour,
			cs->cursor.width);
		if (TCL_ERROR == Tcl_Eval(GetInterp(), cmd))
		    printf("%s\n", GetInterp()->result);

		/* fill in local position of cursor in label box */
		local_pos = CSLocalCursor(io, wy);

		label = get_default_string(GetInterp(), gap_defs,
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
	case TASK_CS_REDRAW:
	    {
		/* HACK - never used */
		int i, id = register_id();

		for (i = 1; i <= NumContigs(io); i++) {
		    contig_deregister(io, i, cs_callback, fdata);
		    contig_register(io, i, cs_callback, fdata, id,
				    REG_REQUIRED |
				    REG_DATA_CHANGE |
				    REG_OPS |
				    REG_NUMBER_CHANGE |
				    REG_ANNO |
				    REG_GENERIC |
				    REG_FLAG_INVIS |
				    REG_BUFFER,
				    REG_TYPE_CONTIGSEL);
		}
		break;
	    }
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
	sprintf(cmd, "ContigParams %d", *handle_io(io));
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
