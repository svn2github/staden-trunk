#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <itk.h>

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "misc.h"
#include "array.h"
#include "tcl_utils.h"
#include "renzyme_box.h"
#include "renzyme_search.h"
#include "graphic_editor.h"
#include "getfile.h"
#include "text_output.h"
#include "renz_utils.h" /*seq_utils */
#include "sequence_formats.h"
#include "parse_feature.h"
#include "dna_utils.h"
#include "spin_globals.h"
#include "editor_reg.h"
#include "nip_results.h"
#include "nip_globals.h"
#include "seq_results.h"
#include "restriction_enzyme_map.h"
#include "canvas_box.h"
#include "renzyme_map_canvas.h"
#include "read_sequence.h"
#include "editor.h"
#include "seq_reg_structs.h"

/* needed for update_arg */

#define NAMEDEL "/"
#define SEQDEL "/"
#define CUT '\''
#define PAD 'N'
#define MAXLINE 1024
#define MAXRSEQ 10
#define WHITESPACE " \t\n"
#define COMMENT ';'
/* above defines come from seq_utils/renz_utils.c */

extern int compare_rmatch(const void *p1, const void *p2);
extern void printWorld(WorldPtr *world);
extern void printCanvas(CanvasPtr *c);

void ed_nip_renz_print_func(void *obj, seq_reg_plot *plot);
void ed_nip_renz_text_func(void *obj);
void static ed_renz_callback(int seq_num, void *obj, editor_reg_data *jdata);

/*
 * find the position of the left hand base of where the enzyme target sequence
 * cuts the dna
 */
int find_matches(RENZYME *rs,                                         /* in */
		 char *sequence,                                      /* in */
		 int sequence_len,                                    /* in */
		 int sequence_type,                                   /* in */
		 int rs_id,                                           /* in */
		 R_MATCH **match,                                     /* out */
		 int *total_matches)                                  /* out */
{
    
    #define SIZE_HASH 256

    int *seq_hash_values;
    int *matches, max_matches = MAXMATCHES;
    int num_matches, k;
    int cnt = 0;
    int array_size = MAXMATCHES;
    int array_inc = MAXMATCHES;
    int last_word[SIZE_HASH];
    int word_count[SIZE_HASH];

    if ( ! (seq_hash_values = (int *) xmalloc ( sizeof(int)*(sequence_len)))) {
	return -2;
    }
    if ( ! (matches = (int *) xmalloc ( sizeof(int)*(max_matches) ))) {
	return -2;
    }
    /* hash sequence */
    hash_dna(sequence, sequence_len, seq_hash_values, last_word, word_count);
    
    if (sequence_type == 2 || sequence_type == 4)
	sequence_type = 1;
    if (sequence_type == 1 || sequence_type == 3)
	sequence_type = 0;

    /* for each recognition sequence find the matches between rec seq and dna */
    dna_search(sequence, sequence_len, rs->rec_seq,
	       strlen (rs->rec_seq), sequence_type,
	       seq_hash_values, last_word, word_count, matches, 
	       max_matches, &num_matches);
    /* store the match results in array 'match' */
    for (k = 0; k < num_matches; k++) {
	(*match)[cnt].enz_name = rs_id;
	/* store the cut position in matches as required for */
	if (matches[k] + rs->cut_pos_1 == sequence_len) {
	    (*match)[cnt].cut_pos1 = sequence_len;	   
	    (*match)[cnt].cut_pos2 = sequence_len;	   
	} else {
	    (*match)[cnt].cut_pos1 = matches[k] + rs->cut_pos_1;
	    if (rs->cut_pos_2 != 0) {
		(*match)[cnt].cut_pos2 = matches[k] 
		    + strlen(rs->rec_seq) - rs->cut_pos_2;
	    } else {
		(*match)[cnt].cut_pos2 = matches[k] 
		    + strlen(rs->rec_seq) - rs->cut_pos_1;
	    }
	}
	cnt++;
	/* cnt is less than MAXMATCHES 
	 * if it is larger, then need to allocate more space to array
	 */
	if (cnt >= array_size) {
	    /* set array_size to be count plus arbitary increment */
	    array_size = cnt + array_inc;
	    if (NULL == ( (*match) = (R_MATCH *)realloc((*match), 
							array_size * 
							sizeof(R_MATCH))))
		return 0;
	    /* clear the new memory */
	    memset(&(*match)[cnt], 0, sizeof(R_MATCH) * array_inc); 
	}
    }
    *total_matches = cnt;
    xfree(seq_hash_values);
    xfree(matches);
    return 1;
}

/*
 * find the position of the left hand base of where the enzymes target sequence
 * cuts the dna
 */
int find_matches_all (RENZYMES *rss,                                         /* in */
		      char *sequence,                                       /* in */
		      int sequence_len,                                     /* in */
		      int sequence_type,                                    /* in */
		      R_MATCH **match,                                      /* out */
		      int *total_matches)                                   /* out */
{

    #define SIZE_HASH 256

    int *seq_hash_values;
    int *matches, max_matches = MAXMATCHES;
    int num_matches, k;
    int cnt = 0;
    int array_size = MAXMATCHES;
    int array_inc = MAXMATCHES;
    int last_word[SIZE_HASH];
    int word_count[SIZE_HASH];
    int i, num_enzymes;
    RENZYME *rs;

    if ( ! (seq_hash_values = (int *) xmalloc ( sizeof(int)*(sequence_len)))) {
	return -2;
    }
    if ( ! (matches = (int *) xmalloc ( sizeof(int)*(max_matches) ))) {
	return -2;
    }
    /* hash sequence */
    hash_dna(sequence, sequence_len, seq_hash_values, last_word, word_count);
    
    if (sequence_type == 2 || sequence_type == 4)
	sequence_type = 1;
    if (sequence_type == 1 || sequence_type == 3)
	sequence_type = 0;
    
    num_enzymes = rss->used;
    /* for each enzyme */
    for (i = 0; i < num_enzymes; i++) {
	rs = rss->renzyme[i];
    /* for each recognition sequence find the matches between rec seq and dna */
	dna_search(sequence, sequence_len, rs->rec_seq,
		   strlen (rs->rec_seq), sequence_type,
		   seq_hash_values, last_word, word_count, matches, 
		   max_matches, &num_matches);
	/* store the match results in array 'match' */

	for (k = 0; k < num_matches; k++) {
	    (*match)[cnt].enz_name = i;
	    /* store the cut position in matches as required for */
	    if (matches[k] + rs->cut_pos_1 == sequence_len) {
		(*match)[cnt].cut_pos1 = sequence_len;	   
		(*match)[cnt].cut_pos2 = sequence_len;	   
	    } else {
		(*match)[cnt].cut_pos1 = matches[k] + rs->cut_pos_1;
		if (rs->cut_pos_2 != 0) {
		    (*match)[cnt].cut_pos2 = matches[k] 
			+ strlen(rs->rec_seq) - rs->cut_pos_2;
		} else {
		    (*match)[cnt].cut_pos2 = matches[k] 
			+ strlen(rs->rec_seq) - rs->cut_pos_1;
		}
	    }
	    cnt++;
	    /* cnt is less than MAXMATCHES 
	     * if it is larger, then need to allocate more space to array
	     */
	    if (cnt >= array_size) {
		/* set array_size to be count plus arbitary increment */
		array_size = cnt + array_inc;
		if (NULL == ( (*match) = (R_MATCH *)realloc((*match), 
							    array_size * 
							    sizeof(R_MATCH))))
		    return 0;
		/* clear the new memory */
		memset(&(*match)[cnt], 0, sizeof(R_MATCH) * array_inc); 
	    }
	}
    }
    *total_matches = cnt;
    xfree(seq_hash_values);
    xfree(matches);
    return 1;
}

RENZYME *renzyme_copy ( RENZYME *r2 ) {

    RENZYME *r;
    char *name, *rs, *pt, *cr, *rst;

    r = init_renzyme ();

    if (NULL == (name = (char *)malloc((strlen (r2->name )+ 1) * sizeof(char))))
	goto err;
    if (NULL == (rs = (char *)malloc((strlen (r2->rec_seq )+ 1) * sizeof(char))))
	goto err;
    if (NULL == (pt = (char *)malloc((strlen (r2->prototype )+ 1) * sizeof(char))))
	goto err;
    if (NULL == (cr = (char *)malloc((strlen (r2->supplier_codes )+ 1) * sizeof(char))))
	goto err;
    if (NULL == (rst = (char *)malloc((strlen (r2->rec_seq_text )+ 1) * sizeof(char))))
	goto err;
 
    strcpy (name, r2->name);
    name[strlen (r2->name)] = 0;
    strcpy (rs, r2->rec_seq);  
    rs[strlen (r2->rec_seq)] = 0;
    strcpy (rst, r2->rec_seq_text);  
    rst[strlen (r2->rec_seq_text)] = 0;
    strcpy (pt, r2->prototype);
    pt[strlen (r2->prototype)] = 0;
    strcpy (cr, r2->supplier_codes);
    cr[strlen (r2->supplier_codes)] = 0;

    r->name = name;
    r->rec_seq = rs;
     r->rec_seq_text = rst;
    r->prototype = pt;
    r->supplier_codes = cr;
    r->cut_pos_1 = r2->cut_pos_1;
    r->cut_pos_2 = r2->cut_pos_2;
    r->av_frag_size = r2->av_frag_size;
    return r;
 err:
    if (name) xfree(name);
    if (rs) xfree(rs);
    if (rst) xfree(rst);
    if (pt) xfree(pt);
    if (cr) xfree(cr);
    return NULL;
}


/*
 * tcl interface to plotting lines on restriction enzyme display
 * stand alone and template displays
 */
void
EdPlotStickMap(Tcl_Interp *interp,
	     char *win_name,
	     int cut_site,
	     int xoffset,
	     int yoffset,
	     int tick_ht,
	     int tick_wd,
	     char *colour,
	     int index,
	     int sequence_len,
	     int seq_id)
{
    char cmd[1024];
    int xcut;

    /* Force cut sites off the end of the contig to be displayed. */
    xcut = MAX(cut_site, 0);
    xcut = MIN(sequence_len, xcut);
    
    /* printf ("xoffset=%d  xcut=%d\n", xoffset, xcut); */
    sprintf(cmd,"%s create line %d %d %d %d -fill %s -width %d "
	    "-tag {S renz re_%d pos_%d seqid_%d}",
	    win_name, xcut+xoffset, yoffset, xcut+xoffset, 
	    yoffset + tick_ht, colour, tick_wd, index, cut_site, seq_id);
    Tcl_Eval(interp, cmd);
}

void ed_plot_renz_matches(Tcl_Interp *interp,
			  char *re_win,
			  char *names_win,
			  int text_offset,
			  char *t_colour,
			  int yoffset,
			  int num_enzymes,
			  RENZYMES *r_enzyme,
			  ruler_s *ruler,
			  int sequence_len,
			  int num_matches,
			  R_MATCH *match,
			  tick_s *tick,
			  char *frame,
			  WorldPtr *world,
			  CanvasPtr *canvas,
			  win **win_list,
			  int num_wins,
			  StackPtr **zoom,
			  int seq_id,
			  SELECTION *sel)
{
    char cmd[1024];
    int item;
    int seq_len;
    int offset;
    int y_offset;
    int t_offset;
    int i;
    int cut_site;
    int *done;

    if (NULL == (done = (int *)xcalloc(num_enzymes, sizeof(int)))) {
	return;
    }

    sprintf(cmd, "%s delete all", re_win);
    Tcl_Eval(interp, cmd);

    sprintf(cmd, "%s delete all", names_win);
    Tcl_Eval(interp, cmd);

    t_offset = text_offset;
    y_offset = yoffset;

    seq_len = ruler->end - ruler->start + 1;

    /*printf ("seq_len=%d\n", seq_len);*/
   
    /* for each selected enzyme, display its name and cut sites */
    for (item = 0; item < num_enzymes; item++) {
        offset = yoffset;
	sprintf(cmd,"%s create text 10 %d -text %s -anchor w -fill %s "
		"-font enzyme_font -tag {S re_%d seqid_%d}",
		names_win, t_offset, r_enzyme->renzyme[item]->name, t_colour, 
		item, seq_id);
	Tcl_Eval(interp, cmd);
	sprintf(cmd, "%s create line %d %d %d %d -tag contig -fill %s",
		re_win, ruler->start, y_offset, ruler->end, y_offset, ruler->colour);
	Tcl_Eval(interp, cmd);
	 
	for (i = 0; i < num_matches; i++) {
	    if (item == match[i].enz_name) {	    
		offset = (item * tick->ht) + yoffset;
                cut_site = ruler->start - 1 + match[i].cut_pos1 - 1;
		EdPlotStickMap(interp, re_win, cut_site, 0, offset, 
			       tick->ht, tick->line_width, tick->colour, item,
			       ruler->end - ruler->start + 1, seq_id);
	    }
	}
	/* draw a selection on the re_win, if it needed */
	if (sel != NULL) {
	if (sel->sel_first != 0 && sel->sel_last != 0) {
	    if (sel->name_first != NULL && sel->name_last != NULL) {
		if (!strcmp (sel->name_first, r_enzyme->renzyme[item]->name)) {
		    sprintf (cmd, "%s create rectangle %d %d %d %d -fill green -outline green -tags \"fragment dummy seqid_%d xL_%d xR_%d name1_%s name2_%s\"", 
			     re_win, sel->sel_first+2, offset+2, sel->sel_last-2, offset+tick->ht-2, 
			     seq_id, sel->sel_first,sel->sel_last,sel->name_first, sel->name_last);
		}
	    } else {
		sprintf (cmd, "%s create rectangle %d %d %d %d -fill green -outline green -tags \"dummy\"", 
			 re_win, sel->sel_first, offset, sel->sel_last, offset+tick->ht); 
	    }
	    Tcl_Eval(interp, cmd);
	}
	}
	y_offset += tick->ht;
	t_offset += tick->ht;
    }
    /* draw final line */
    sprintf(cmd, "%s create line %d %d %d %d -tag contig -fill %s",
	    re_win, ruler->start, y_offset, ruler->end, y_offset, ruler->colour);
    Tcl_Eval(interp, cmd);
    
    /* add back selection rectangles */
    if (TCL_ERROR == Tcl_VarEval(interp, "ReSelectRect ", frame, " ",
				 names_win, NULL))
	verror(ERR_WARN, "plot_renz_matches", "%s\n", interp->result);

    world->total->x1 = ruler->start;
    world->total->x2 = ruler->end;
    world->total->y1 = 1;
    world->total->y2 = y_offset;

    memcpy(world->visible, world->total, sizeof(d_box));
    world->visible->y2 = canvas->height;
   
    /* defined in tk_utils/canvas_box.c */
    SetCanvasCoords(interp, world->visible->x1, world->visible->y1, 
		    world->visible->x2, world->visible->y2, canvas);

    draw_single_ruler(interp, ruler, canvas, ruler->start, ruler->end, 1);

    scaleCanvas(interp, win_list, num_wins, "all", world->visible, canvas);
    scrollRegion(interp, win_list, num_wins, world->total, canvas);
   
    /* remove all current zooming info */
	freeZoom(zoom);
    /* add first zoom */
	pushZoom(zoom, world->visible);
    /* defined in tk_utils/canvas_box.c */   
}

void ed_renz_shutdown(RENZYMES *r_enzyme, 
		   int num_enzymes,
		   R_MATCH *match,
		   CanvasPtr *canvas,
		   WorldPtr *world,
		   StackPtr *zoom)
{
    int i;

    /* Free memory */
    if (r_enzyme) {
	for (i = 0; i < num_enzymes; i++) {
	    xfree(r_enzyme->renzyme[i]->name);
	    xfree(r_enzyme->renzyme[i]->rec_seq);
	    xfree(r_enzyme->renzyme[i]->rec_seq_text);
	    xfree(r_enzyme->renzyme[i]->prototype);
	    xfree(r_enzyme->renzyme[i]->supplier_codes);
	}
	xfree(r_enzyme);
    }

    xfree(match);

    if (canvas)
	xfree(canvas);
    
    if (world->visible)
	xfree(world->visible);

    if (world->total)
	xfree (world->total);
    
    if (world)
	xfree(world);

    freeZoom(&zoom);
}

void editor_renz_shutdown(Tcl_Interp *interp, editor_result *result, int seq_num){

    EDITOR_RECORD *er;
    ed_renz_res *data = result->data;
    out_canvas_e *output = result->output;
    char cmd[1024];
    int ed_id, seq_id;
    
    /* need to deregister sequence */
    editor_deregister(seq_num, ed_renz_callback, (editor_result *)result);
    ed_delete_cursor(seq_num, output->cursor->id, 0);
    
    seq_id = GetEdenId (seq_num);
    ed_id = GetEdIdFromSeqId (seq_id);
    er = GetEditor (ed_id);
    er->graphical = 0;
    if (er->text == 0) {
	delete_editor (ed_id);	
    }
    
    sprintf(cmd, "DeleteREnzPlot %s %s\n", data->frame, data->re_win);
    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	verror(ERR_WARN, "restriction enzymes", "shutdown %s\n", 
	       interp->result);
    }

    ed_renz_shutdown(data->r_enzyme, data->num_enzymes, data->match,
		  data->canvas, data->world, data->zoom);
    free(data->text_colour);
    free(data->tick->colour);
    free(data->cursor.colour);
    free(data->ruler->tick.t.colour);
    free(data->ruler->colour);

    xfree(result->data);
    /* xfree(result->input); */
    xfree(result->output);
    xfree(result);
}

int ed_renz_replot (Tcl_Interp *interp, 
		    out_canvas_e *output, 
		    editor_result *result,
		    ed_renz_res *data) 
{

    R_MATCH *match;
    int total_matches;
    int seq_num, seq_len;
    char *seq;

    seq_num = GetEdenNum (data->seq_id);   
    seq = GetEdenSeq (seq_num);         /* sequence and its length has been changed,*/
    seq_len = GetEdenLength (seq_num);  /* so we inquire every time before plot */    
   
    if (NULL == (match = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
	return -1;
    find_matches_all (data->r_enzyme, seq, seq_len, data->sequence_type, &match, &total_matches );
    /* printf ("total_matches=%d\n", total_matches);  */
    data->match = match;
    data->num_match = total_matches;
    data->sequence_len = seq_len; /* do we need it? we have seq_id, so we can get it from seq_id */
    data->ruler->end = seq_len;
    ed_plot_renz_matches(output->interp, data->re_win, data->names_win,
		      data->text_offset, data->text_colour, 
		      data->yoffset, data->num_enzymes, data->r_enzyme,
		      data->ruler, data->sequence_len,
		      data->num_match, data->match,
		      data->tick, data->frame, data->world,
		      data->canvas, data->win_list, data->num_wins,
		      &data->zoom, data->seq_id, data->sel);    
    if (total_matches == 0) {
	return -1;
    }
    
    return 0;
}

int ed_renz_reg (Tcl_Interp *interp,
		 int seq_id,
		 out_canvas_e *output,
		 char *frame,
		 char *names_win,
		 char *re_win,
		 char *inlist,
		 int num_items,
		 int start,
		 int end,
		 int text_offset,
		 char *text_fill,
		 tick_s *tick,
		 int yoffset,
		 ruler_s *ruler,
		 cursor_s cursor)
{
    editor_result *result;
    RENZYMES *r_enzyme;
    SELECTION *selection;
    EDITOR_RECORD *er;
    ed_renz_res *data;
    int id, ed_id, group_id, member_id;
    int seq_num, seq_type;
    seq_reg_cursor_notify cn;
    int line_width;
    char **sel = NULL;
    int num_sel;
     
    if (NULL == (result = (editor_result *)xmalloc(sizeof(editor_result))))
	return -1;
    if (NULL == (data = (ed_renz_res *)xmalloc(sizeof(ed_renz_res))))
	return -1;
  
    seq_num = GetEdenNum(seq_id); /* get position of the sequence in registration list */

    seq_type = GetEdenType (seq_num);
    result->data = data;
    result->seq_id[HORIZONTAL] = seq_id;
    result->seq_id[VERTICAL] = -1;
 
    id = get_editor_reg_id();
    result->id = id; 
    result->output = (void *)output;
    result->pr_func = ed_nip_renz_print_func;
    result->op_func = ed_renz_callback;
    result->txt_func = ed_nip_renz_text_func;

    /* data structure */
    strcpy(data->re_win, re_win);
    strcpy(data->frame, frame);
    strcpy(data->names_win, names_win);
    
    data->tick = tick;
    data->ruler = ruler;
    data->cursor = cursor; /* cursor_s */
    data->sequence_len = GetEdenLength (seq_num);
    data->sequence_type = seq_type;
   
    /* create list of windows in the restriction enzyme display */
    if (NULL == (data->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;

    data->num_wins = 0;
    addWindow(data->win_list, &data->num_wins, data->re_win, 'b', id);
    addWindow(data->win_list, &data->num_wins, data->ruler->window, 'x', id);
    addWindow(data->win_list, &data->num_wins, data->names_win, 'y', id);
   
    if (NULL == (data->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (data->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (data->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (data->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    initCanvas(interp, data->canvas, data->re_win);
    createZoom(&data->zoom);

    /* create selecing Renzyme name array */
    if (Tcl_SplitList(interp, inlist, &num_sel, &sel) != TCL_OK)
      return -1;
    
    r_enzyme = get_selected_renzyme ( num_sel, sel);   
    data->r_enzyme = r_enzyme;
    data->num_enzymes = num_sel; /* num_items = num_sel ? */;    
    data->tick = tick;
    data->yoffset = yoffset;
    data->text_offset = text_offset;   
    data->text_colour = strdup(text_fill);
    /*data->seq_id = seq_num;*/
    data->seq_id = seq_id;
    data->match = NULL;
    data->num_match = 0;
    selection = init_selection ();
    data->sel = selection;
   
    line_width = get_default_int(interp, nip_defs, w("NIP.CURSOR.LINE_WIDTH"));
   
    /* private=0: share cursor for every plot */
    output->cursor = editor_create_cursor(seq_num, 0, NULL, line_width, 1, HORIZONTAL);
    output->cursor_visible = 0;
    
    ed_id = GetEdIdFromSeqId (seq_id);
    member_id = GetMemberIdFromSeqId (seq_id);
    group_id = GetGroupIdFromSeqId (seq_id);
    er = GetEditor (ed_id);
    er->graphical = 1;
    output->cursor->frame_name = NULL;
   
    output->cursor->posy = member_id;
    set_editor_cursor (ed_id, output->cursor);

    /* move cursor to start position if this is our own cursor */
    if (output->cursor->refs == 1) {
	output->cursor->abspos = start;
	output->cursor->posy = 1;
    }
    
    editor_register(seq_num, ed_renz_callback, (void *)result, SEQ_PLOT_TEMP, id);

    ed_renz_replot (interp, output, result, data);
    /* if above returns -1, may want to
     * to editor_renz_shutdown(interp, result, seq_num); ?
     */
    
    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = output->cursor;
    cn.cursor->job = CURSOR_MOVE;
    cn.selection = selection;
    editor_notify(seq_num, (editor_reg_data *)&cn);
    
    return id;
}

int CreateGraphicEditor (ClientData clientData, Tcl_Interp *interp, int argc, char **argv){

    ed_renz_arg args;
    cursor_s cursor; /* element is width and colour */
    tick_s *tick;
    ruler_s *ruler;
    int id;
    out_canvas_e *output;
    Tcl_DString input_params;
    int seq_num, seq_id;

    cli_args a[] = {
	{"-frame",	 ARG_STR, 1, NULL, offsetof(ed_renz_arg, frame)},
	{"-win_names",	 ARG_STR, 1, NULL, offsetof(ed_renz_arg, win_name)},
	{"-window",	 ARG_STR, 1, NULL, offsetof(ed_renz_arg, plot)},
	{"-win_ruler",	 ARG_STR, 1, NULL, offsetof(ed_renz_arg, win_ruler)},
	{"-enzymes",	 ARG_STR, 1, NULL, offsetof(ed_renz_arg, inlist)},
	{"-num_enzymes", ARG_INT, 1, NULL, offsetof(ed_renz_arg, num_items)},
	{"-text_offset", ARG_INT, 1, NULL, offsetof(ed_renz_arg, text_offset)},
	{"-text_fill",   ARG_STR, 1, NULL, offsetof(ed_renz_arg, text_fill)},
	{"-tick_height", ARG_INT, 1, "-1", offsetof(ed_renz_arg, tick_ht)},
	{"-tick_width",  ARG_INT, 1, "-1", offsetof(ed_renz_arg, tick_wd)},
	{"-tick_fill",   ARG_STR, 1,   "", offsetof(ed_renz_arg, tick_fill)},
	{"-cursor_width",ARG_INT, 1, "-1",  offsetof(ed_renz_arg, cursor_wd)},
	{"-cursor_fill", ARG_STR, 1,  "",  offsetof(ed_renz_arg, cursor_fill)},
	{"-yoffset",	 ARG_INT, 1, NULL, offsetof(ed_renz_arg, yoffset)},
	{"-seq_id",	 ARG_INT, 1, NULL, offsetof(ed_renz_arg, seq_id)},
	{"-start",	 ARG_INT, 1, "1",  offsetof(ed_renz_arg, start)},
	{"-end",	 ARG_INT, 1, "-1", offsetof(ed_renz_arg, end)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (output = (out_canvas_e *)xmalloc(sizeof(out_canvas_e))))
	return TCL_OK;

    set_char_set(DNA);       
    seq_id = args.seq_id;
    
    /* get register num */
    seq_num = GetEdenNum (seq_id);

    if (args.end == -1) {
	args.end = GetEdenLength (seq_num); 
    } 
    vfuncheader("restriction enzyme plot");
    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    /*vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"  
      "enzymes: %s\n", "NAME", args.start, args.end,args.inlist);*/
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    cursor = cursor_struct(interp, tk_utils_defs, "R_ENZ", args.cursor_wd, args.cursor_fill);
    tick = tick_struct(interp, tk_utils_defs, "R_ENZ", args.tick_wd, args.tick_ht, args.tick_fill);
    /*printf ("line_width=%d\n",tick->line_width); */
    ruler = ruler_struct(interp, tk_utils_defs, "R_ENZ", 0);
    ruler->start = args.start;
    ruler->end = args.end;
    strcpy(ruler->window, args.win_ruler);
    output->interp = interp;
   
    id = ed_renz_reg(interp, args.seq_id, output, args.frame, 
		     args.win_name, args.plot, args.inlist, args.num_items,
		     args.start, args.end, args.text_offset, args.text_fill,
		     tick, args.yoffset, ruler, cursor);
   
    vTcl_SetResult(interp, "%d", id);

    return TCL_OK;
}


void
ed_nip_renz_print_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_canvas_e *output = result->output;
    ed_renz_res *data = result->data;

#ifdef DEBUG_START
    printf("START ed_nip_renz_print_func\n");
#endif
    ed_plot_renz_matches(output->interp, data->re_win, data->names_win,
		      data->text_offset, data->text_colour,
		      data->yoffset, data->num_enzymes, data->r_enzyme,
		      data->ruler, data->sequence_len,
		      data->num_match, data->match,
		      data->tick, data->frame, data->world,
		      data->canvas, data->win_list, data->num_wins,
		      &data->zoom, data->seq_id, data->sel);

}

void
ed_nip_renz_text_func(void *obj)
{

}

void
EdExpandRSeq(int match,                                                  /* in */
	     int cut_pos,                                                /* in */
	     char *sequence,                                             /* in */
	     int sequence_len,                                           /* in */
	     int sequence_type,                                          /* in */
	     char *in_seq,                                               /* in */
	     char *out_seq)                                             /* out */
{
    int i, j;
    int cnt = 0;
    int start = 0;
    int length;
    int inseq_len = 0;

    /* Match is the position of the cut site. We want to go back to the first
     * base of the recognition sequence.
     */
    match--; /* starts at offset 1, we want 0 */
    if (cut_pos > 0) {
	int count = cut_pos;
	while (count > 0) {
	    while (--match && match >= 0 && sequence[match] == '*')
		;
	    count--;
	}
    } else {
	match -= cut_pos;
    }

    inseq_len = strlen(in_seq);

    /* only happens if at the end of a circular dna */
    if ((match < 0) && (sequence_type == 1)) { 
	match = sequence_len + match;
    }

    /* padding characters before target string */
    if (cut_pos < 0) {
	length = inseq_len;
	start = cut_pos;
    }
    /* padding characters after target string */
    else if (cut_pos >= inseq_len) {
	length = cut_pos + 1;
    } else {
	length = inseq_len;
    }

    for (j = 0, i = start; i < length; i++) {
	if (i == cut_pos) {
	    out_seq[cnt++] = CUT;
	    if (cut_pos >= inseq_len)
		break;
	}

	if (sequence_type == 0) { /* linear dna */
	    /* at end */
	    while (1) {
		if (match+i+j >= sequence_len || match+i < 0) {
		    out_seq[cnt++] = PAD;
		} else {
		    if (sequence[match+i+j] == '*') {
			j++;
			continue;
		    }
		    out_seq[cnt++] = sequence[match+i+j];
		} 
		break;
	    }

	} else { /* circular dna */
	    while (sequence[(match+i+j+sequence_len)%sequence_len] == '*')
		j++;
	    out_seq[cnt++] = sequence[(match+i+j+sequence_len)%sequence_len];
	}
    }

    /* terminate out_seq */
    out_seq[cnt] = '\0';
}


static int ed_compare_rmatch(const void *p1, const void *p2)
{
    R_MATCH *r1 = (R_MATCH *)p1, *r2 = (R_MATCH *)p2;
    if ((*r1).cut_pos1 < (*r2).cut_pos1) 
	return (-1);
    else if ((*r1).cut_pos1 == (*r2).cut_pos1)
	return (0);
    else
	return (1);
}

static int ed_compareint(const void *p1, const void *p2)
{
    int *i1 = (int *)p1, *i2 = (int *)p2;
    if (*i1 < *i2) 
	return (-1);
    else if (*i1 == *i2)
	return (0);
    else
	return (1);
}

void EdFindFragments(int num_matches,                            /* in */
		     R_MATCH *tmp_match,                         /* in */
		     int sequence_len,                           /* in */
		     int sequence_type,                          /* in */
		     int *fragment)                              /* out */
{
    int i;

    /* if sequence is a circle */
    if (sequence_type == 1) {
	fragment[0] = sequence_len - tmp_match[num_matches-1].cut_pos1 + 
	    tmp_match[0].cut_pos1;
	
	for (i = 1; i < num_matches; i++) {
	    
	    fragment[i] = tmp_match[i].cut_pos1 - tmp_match[i-1].cut_pos1;
	}
    } else {
	fragment[0] = tmp_match[0].cut_pos1 - 1;
	for (i = 1; i < num_matches; i++) {
	    fragment[i] = tmp_match[i].cut_pos1 - tmp_match[i-1].cut_pos1;
	}
	fragment[num_matches] = sequence_len - tmp_match[num_matches-1].cut_pos1 + 1;
    }
}



int EdOrderOnPosition (RENZYMES *r_enzyme,
		       R_MATCH *match,
		       int total_matches,
		       char *sequence,
		       int sequence_len,
		       int sequence_type,
		       int lreg)
{
    int k;
    R_MATCH *tmp_match;
    char r_seq[MAXLINE];
    int *fragment;
    int *lengths;
    char fbuf[1024], lbuf[1024];
    int fragments_printed = 0;

    /* make temporary copy of match data */
    if (NULL == (tmp_match = (R_MATCH *)xmalloc(total_matches * sizeof(R_MATCH))))
	return 0;
    memcpy(tmp_match, match, total_matches * sizeof(R_MATCH));

    /* sort on position over all the enzymes */
    qsort((void *) tmp_match, total_matches, sizeof(R_MATCH), ed_compare_rmatch);
 
    vmessage("%10s%20s%34s%9s%8s\n", "Name", "Sequence", "Position", 
	   "Fragment", "lengths");

    /* 
     * malloc total_matches +1 because for linear sequences, have an extra
     * fragment ie if cut n times, get n+1 fragments
     */
    if (NULL == (fragment = (int *)xmalloc((total_matches+1) * sizeof(int))))
	return 0;
    if (NULL == (lengths = (int *)xmalloc((total_matches+1) * sizeof(int))))
	return 0;

    EdFindFragments(total_matches, tmp_match, sequence_len, sequence_type, 
		  fragment);

    /* linear */
    if (sequence_type == 0) { 
	memcpy(lengths, fragment, (total_matches + 1) *sizeof(int));
	qsort((void *) lengths, (total_matches +1), sizeof(int), ed_compareint);
    }
    /* circular */
    else {
	memcpy(lengths, fragment, (total_matches) *sizeof(int));
	qsort((void *) lengths, (total_matches), sizeof(int), ed_compareint);
    }
    
/*
    memcpy(lengths, fragment, total_matches*sizeof(int));
    qsort((void *) lengths, total_matches, sizeof(int), ed_compareint);
*/
    for (k = 0; k < total_matches; k++) {

	EdExpandRSeq(tmp_match[k].cut_pos1,
		     /*r_enzyme[tmp_match[k].enz_name].cut_site[tmp_match[k].enz_seq],*/
		     r_enzyme->renzyme[tmp_match[k].enz_name]->cut_pos_1,  /*FIXME*/
		     sequence, 
		     sequence_len,
		     sequence_type,
		     /*r_enzyme[tmp_match[k].enz_name].seq[tmp_match[k].enz_seq],*/
		     r_enzyme->renzyme[tmp_match[k].enz_name]->rec_seq,
		     r_seq);
	    
	if (fragment[k] >= 0 && fragment[k] <= sequence_len) {
	    sprintf(fbuf, "%7d", fragment[k]);
	    fragments_printed++;
	} else
	    sprintf(fbuf, "%7s", "-");
	if (lengths[k] >= 0)
	    sprintf(lbuf, "%7d", lengths[k]);
	else
	    sprintf(lbuf, "%7s", "-");
	    
	vmessage("%5d %-15s %-32s%10d%s%s \n",
		 k+1,
		 r_enzyme->renzyme[tmp_match[k].enz_name]->name,
		 r_seq,
		 tmp_match[k].cut_pos1 + lreg - 1,
		 fbuf, lbuf);
    }
    /* print last fragment if linear sequence */
    if (sequence_type == 0) {
	if (fragment[total_matches] >= 0)
	    vmessage("%71d%7d \n",
		     fragment[total_matches],
		     lengths[total_matches]);
	else if (fragments_printed < 2)
	    vmessage("%71d%7d \n",
		     lengths[total_matches],
		     lengths[total_matches]);
	else
	    vmessage("%71s%7d \n",
		     "-",
		     lengths[total_matches]);
    }
    xfree(tmp_match);
    xfree(fragment);
    xfree(lengths);

    return 1;
}


int EdPrintEnzymeByEnzyme(RENZYMES *r_enzyme,                           /* in */
			  R_MATCH *match,                              /* in */
			  int total_matches,                           /* in */
			  int num_enzymes,                             /* in */
			  char *sequence,                              /* in */
			  int sequence_len,                            /* in */
			  int sequence_type,                           /* in */
			  int lreg)                                    /* in */
{
    int i,k;
    int cnt =0;
    char r_seq[MAXLINE];
    R_MATCH *tmp_match;
    int start = 0;
    int num_matches;
    int j = 0;
    int *fragment;
    int *lengths;
    char fbuf[1024], lbuf[1024];
    int fragments_printed = 0;

    if (0 == num_enzymes)
      return 1;

    /* printf("total matches %d\n", total_matches); */
    if (NULL == (tmp_match = (R_MATCH *)xmalloc(total_matches * sizeof(R_MATCH))))
	return 0;
    
    /* printf("num enz %d\n", num_enzymes); */

    for (i = 0; i < num_enzymes; i++) {
	while ((j < total_matches) && (i == match[j].enz_name)) {
	    memcpy(&tmp_match[cnt], &match[j], sizeof(R_MATCH));
	    cnt++;
	    j++;
	}
	if (0 == (num_matches = j - start))
	    continue;

	/* 
	 * malloc total_matches +1 because for linear sequences, have an extra
	 * fragment ie if cut n times, get n+1 fragments
	 */
	if (NULL == (fragment = (int *)xmalloc((num_matches +1)* sizeof(int))))
	    return 0;
	if (NULL == (lengths = (int *)xmalloc((num_matches +1) * sizeof(int))))
	    return 0;
	
	qsort((void *) tmp_match, num_matches, sizeof(R_MATCH), 
	      ed_compare_rmatch);
	
	vmessage("  Matches found= %5d \n", num_matches);
	vmessage("%10s%20s%34s%9s%8s\n", "Name", "Sequence", "Position", 
	       "Fragment", "lengths");

	EdFindFragments(num_matches, tmp_match, sequence_len, sequence_type,
		      fragment);

	/* linear */
	if (sequence_type == 0) { 
	    memcpy(lengths, fragment, (num_matches + 1) *sizeof(int));
	    qsort((void *) lengths, (num_matches +1), sizeof(int), ed_compareint);
	}
	/* circular */
	else {
	    memcpy(lengths, fragment, (num_matches) *sizeof(int));
	    qsort((void *) lengths, (num_matches), sizeof(int), ed_compareint);
	}

	for (k =0; k < num_matches; k++) {
	    /*ExpandRSeq(tmp_match[k].cut_pos,
		       r_enzyme[tmp_match[k].enz_name].cut_site[tmp_match[k].enz_seq],
		       sequence, 
		       sequence_len,
		       sequence_type,
		       r_enzyme[tmp_match[k].enz_name].seq[tmp_match[k].enz_seq],
		       r_seq);*/
	    EdExpandRSeq(tmp_match[k].cut_pos1,
			 r_enzyme->renzyme[tmp_match[k].enz_name]->cut_pos_1,  /*FIXME*/
			 sequence, 
			 sequence_len,
			 sequence_type,
			 r_enzyme->renzyme[tmp_match[k].enz_name]->rec_seq,
			 r_seq);
	    
	    if (fragment[k] > 0 && fragment[k] <= sequence_len) {
		sprintf(fbuf, "%7d", fragment[k]);
		fragments_printed++;
	    } else
		sprintf(fbuf, "%7s", "-");
	    if (lengths[k] > 0)
		sprintf(lbuf, "%7d", lengths[k]);
	    else
		sprintf(lbuf, "%7s", "-");
	    
	    vmessage("%5d %-15s %-32s%10d%s%s \n",
		     k+1,
		     r_enzyme->renzyme[tmp_match[k].enz_name]->name,
		     r_seq,
		     tmp_match[k].cut_pos1 + lreg - 1,
		     fbuf, lbuf);
	}

	/* print last fragment if linear sequence */
	if (sequence_type == 0) {
	    if (fragment[num_matches] > 0)
		vmessage("%71d%7d \n",
			 fragment[num_matches],
			 lengths[num_matches]);
	    else if (fragments_printed < 2)
		vmessage("%71d%7d \n",
			 lengths[num_matches],
			 lengths[num_matches]);
	    else
		vmessage("%71s%7d \n",
			 "-",
			 lengths[num_matches]);
	}
	cnt = 0;
	num_matches = 0;
	start = j;

	xfree(fragment);
	xfree(lengths);
    }

    vmessage("\nZero cutters:\n");
    cnt = j = start = 0;
    for (i = 0; i < num_enzymes; i++) {
	while ((j < total_matches) && (i == match[j].enz_name)) {
	    cnt++;
	    j++;
	}
	if (0 == (num_matches = j - start)) {
	    vmessage("      %s\n", r_enzyme->renzyme[i]->name);
	}
	cnt = 0;
	num_matches = 0;
	start = j;
     }

   xfree(tmp_match);
    return 1;
}


void ed_nip_renz_info(int seq_num,
		      ed_renz_res *data,
		      int start,
		      int print_opt)
{
    int seq_len;
    char *seq;
    int sequence_type;
    
    seq = GetEdenSeq(seq_num);
    seq_len = GetEdenLength(seq_num);
    /*sequence_type = GetSequenceStructure(seq_num);*/
    sequence_type = GetEdenType(seq_num);

    vfuncheader("Restriction enzymes result list");
    
    vmessage("Sequence %s\n", GetEdenName(seq_num));

    vmessage("Number of enzymes = %d\n", data->num_enzymes);
    vmessage("Number of matches = %d\n", data->num_match);

    if (print_opt == 0) {
	start_message();
	EdPrintEnzymeByEnzyme(data->r_enzyme, data->match, data->num_match, 
			      data->num_enzymes, seq, seq_len, 
			      sequence_type, start);
	end_message(data->re_win);
    } else {
	start_message();
	EdOrderOnPosition(data->r_enzyme, data->match, data->num_match,
			seq, seq_len, sequence_type, start);
	end_message(data->re_win);
    }
}

void static ed_renz_callback(int seq_num, void *obj, editor_reg_data *jdata) {

    editor_result *result = (editor_result *) obj;
    ed_renz_res *data = result->data;

    switch(jdata->job) {     
    case SEQ_CHANGED:
	{	   
	    out_canvas_e *output = result->output;
       
	    /* data->sel = jdata->selected.selection; */
	    data->sel = NULL;
	    ed_renz_replot (output->interp, output, result, data);
	    break;
	}
    case SEQ_SELECTED:
	{
	    out_canvas_e *output = result->output;
	   
	    data->sel = jdata->selected.selection;
	    ed_renz_replot (output->interp, output, result, data);
	    /* printf ("graphicl_editor:SEQ_SELECTED\n"); */
	    
	    break;
	}
    case SEQ_CURSOR_NOTIFY:
	{
	    ed_renz_res *r = result->data;
	    out_canvas_e *output = result->output;
	   
	    ed_nip_canvas_cursor_refresh (output->interp, GetEdenId(seq_num),
					  jdata->cursor_moved.cursor, 
					  output->cursor, r->canvas, r->win_list, 
					  r->num_wins, result->id, 
					  &output->cursor_visible, r->world, 1);
	    
	    /*  printf("graphical_editor: SEQ_CURSOR_NOTIFY\n"); */
	    break;
	}
	 /*case SEQ_QUIT:
	 case SEQ_DELETE:*/
    case GRAPHIC_EDITOR_QUIT:
	{
	    out_canvas_e *output = result->output;
	    editor_renz_shutdown(output->interp, result, seq_num);
	    /* printf ("graphic_editor:GRAPHIC_QUIT\n"); */
	    break;
	} 
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
#ifdef DEBUG
	case OUTPUT:
	    jdata->info.result = (void *)output;
	    break;
	case INPUT: 
	    jdata->info.result = (void *)input;
	    break;
	case NPTS: 
	    jdata->info.result = (void *)num_pts;
	    break;
	case INDEX: 
	    jdata->info.result = (void *)id;
	    break;
#endif
	case RESULT:
	    jdata->info.result = (void *)result;
	    break;
	} /* switch SEQ_RESULT_INFO */
	break;

    case SEQ_GENERIC:
	switch (jdata->generic.task) {
	case TASK_NIP_RENZ_INFO: 
	    {
		/* info on all cuts by a specific restriction enzyme */
		ed_renz_res *r;
		int i;
	        R_MATCH *tmp_match;
		int cnt = 0;
		int *enz_name = (int *)jdata->generic.data;
		int seq_len;
		char *seq;
		int sequence_type;
		int seq_num;

		r = result->data;
#ifdef DEBUG
		printf("TASK_RENZ_INFO name %d\n", *enz_name);
#endif
		seq_num = GetEdenNum(result->seq_id[0]);
		
		seq = GetEdenSeq (seq_num);
		seq_len = GetEdenLength (seq_num);   
		sequence_type = GetEdenType (seq_num); 
		/*FIXME */
		if (NULL == (tmp_match = (R_MATCH *)malloc(r->num_match * 
							   sizeof(R_MATCH))))
		    return;
		/* 
		 * make tmp match array containing only those matches for the 
		 * selected enzyme
		 */
		for (i = 0; i < r->num_match; i++) {
		    if (*enz_name == r->match[i].enz_name) {
			tmp_match[cnt++] = r->match[i];
		    }
		}

		/*start_message();
		PrintEnzymeByEnzyme(r->r_enzyme, tmp_match, cnt,
				    r->num_enzymes, seq, seq_len, 
				    sequence_type, r->ruler.start);
		end_message(r->re_win);
		xfree(tmp_match);*/
		break;
	    }
	}
	break; /* SEQ_GENERIC */
	
    } /* switch */
}

int RenzymeCmds_Init(Tcl_Interp *interp) {
 
    Tcl_CreateCommand(interp, "create_graphic_editor", CreateGraphicEditor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
       
    return TCL_OK;
}
