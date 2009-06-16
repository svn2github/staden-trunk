#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "io-reg.h"
#include "canvas_box.h"
#include "gap_canvas_box.h"
#include "template_display.h"
#include "restriction_enzymes.h"
#include "restriction_enzyme_map.h"
#include "qual.h"
#include "gap_globals.h"
#include "stop_codon.h"
#include "misc.h"
#include "tcl_utils.h"
#include "edUtils.h"
#include "tclXkeylist.h"
#include "dna_utils.h"

#define MAXMATCHES 10000
/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void stop_codon_callback(GapIO *io, int contig,
				void *fdata,
				reg_data *jdata);

int
FindStopCodons(int strand,                                            /* in */
	       char *sequence,                                        /* in */
	       int sequence_len,                                      /* in */
	       int sequence_type,                                     /* in */
	       int num_codons,                                        /* in */
	       char **codon,                                          /* in */
	       R_Match **match,                                      /* out */
	       int *total_matches)                                   /* out */
{
#define PLUS 1
#define COMP 2
#define BOTH 3
#define CODON_LEN 3

    int compareint(int *, int *);
    int *seq_hash_values;
    int *matches, max_matches;
    int num_matches;
    int i, k;
    int cnt = 0;
    int array_size = MAXMATCHES;
    int array_inc = MAXMATCHES;
    int last_word[SIZE_HASH];
    int word_count[SIZE_HASH];
    int start_codon = 0;
    int end_codon;
    int *depad_to_pad;
    int unpadded_len;

    end_codon = num_codons - 1;

    if (strand == COMP) {
	start_codon = num_codons;
    }
    if ((strand == COMP) || (strand == BOTH)) {
	end_codon = num_codons * 2 - 1;
    }

    if ( ! (seq_hash_values = (int *) xmalloc (sizeof(int)*(sequence_len)))){
	return -2;
    }

    /* match gets realloced as it goes */
    if (NULL == (*match = (R_Match*)xcalloc(MAXMATCHES, sizeof(R_Match))))
	return 0;
    
    /* matches isn't realloced, so we use worst case. It's free soon though */
    max_matches = sequence_len+1;
    if ( ! (matches = (int *) xmalloc ( sizeof(int)*(max_matches) ))) {
	return -2;
    }

    /*
     * Hashing code requires non-padded sequence. We also need this to
     * determine the correct reading frame and ORF lengths.
     */
    depad_to_pad = (int *)xcalloc(sequence_len+1, sizeof(int));
    unpadded_len = sequence_len;
    depad_seq(sequence, &unpadded_len, depad_to_pad);
    hash_dna(sequence, unpadded_len, seq_hash_values, last_word,
	     word_count);

     for (i = start_codon; i <= end_codon; i++) {
	/* find the matches between rec seq and dna */

	/* error = search_dna (sequence, sequence_len, codon[i], 
			    CODON_LEN, sequence_type, 
			    matches, max_matches, 
			    &num_matches, seq_hash_values); */
	
	dna_search(sequence, unpadded_len, codon[i], CODON_LEN, 
		   sequence_type, seq_hash_values, last_word, word_count,
		   matches, max_matches, &num_matches);
	    
	/* store the match results in array 'match' */
	for (k = 0; k < num_matches; k++) {

	    /* find reading frame */
	    (*match)[cnt].enz_name = matches[k]%CODON_LEN;
	    /* set reading frame '0' to be 3 */
	    if ((*match)[cnt].enz_name == 0) {
		(*match)[cnt].enz_name = CODON_LEN;
	    }
	    (*match)[cnt].cut_pos = depad_to_pad[matches[k]]%sequence_len;
	    /* for complementary strand, have reading frames 4, 5, 6 */
	    if (i >= num_codons) {
		(*match)[cnt].enz_name += CODON_LEN;
		(*match)[cnt].cut_pos = (depad_to_pad[matches[k]]%sequence_len) + 2;
	    }
	    /* reference to sequence of stop codon */
	    (*match)[cnt].enz_seq = i;
	    /* store the cut position in matches as required for 
	     * displaying 
	     */

	    cnt++;
	    /* cnt is less than MAXMATCHES 
	     * if it is larger, then need to allocate more space to array
	     */
	    if (cnt >= array_size) {
		/* set array_size to be count plus arbitary increment */
		array_size = cnt + array_inc;
		
		if (NULL == ( (*match) = (R_Match *)realloc((*match), 
							    array_size * 
							    sizeof(R_Match))))
		    return 0;
		/* clear the new memory */
		
		memset(&(*match)[cnt], 0, sizeof(R_Match) * array_inc); 
	    }
	}
    }

    *total_matches = cnt;
    xfree(seq_hash_values);
    xfree(matches);
    xfree(depad_to_pad);
    return 1;
}

void
display_stop_codons(Tcl_Interp *interp,
		    GapIO *io,
		    obj_codon *s)
{
    char cmd[1024];
    char *rf[6];
    int text_offset = 25;
    int i;
    int y_offset;
    int cut_site;
    char *t_colour;
    int t_offset;
    int start_rf = 0;
    int end_rf = 2;

    sprintf(cmd, "%s delete all", s->window);
    Tcl_Eval(interp, cmd);

    sprintf(cmd, "%s delete all", s->names_win);
    Tcl_Eval(interp, cmd);

    t_colour = get_default_string(interp, gap_defs, "STOP_CODON.TEXT_COLOUR");

    rf[0] = "\"frame 1 +\"";
    rf[1] = "\"frame 2 +\"";
    rf[2] = "\"frame 3 +\"";
    rf[3] = "\"frame 1 -\"";
    rf[4] = "\"frame 2 -\"";
    rf[5] = "\"frame 3 -\"";

    /* set up start and end reading frames */
    if (s->strand == COMP)
	start_rf = 3;
    if ((s->strand == COMP) || (s->strand == BOTH))
	end_rf = 5;

    y_offset = s->yoffset;
    t_offset = text_offset;

    for (i = start_rf; i <= end_rf; i++) {
	if (s->strand == BOTH && i >= 3) {
	    t_offset = text_offset + ((i+1) * s->tick->ht);
	    y_offset = s->yoffset + (i * s->tick->ht);
	} else {
	    t_offset = text_offset + ((i%3) * s->tick->ht);
	    y_offset = s->yoffset + ((i%3) * s->tick->ht);
	}

	sprintf(cmd, "%s create line %d %d %d %d -fill %s -tag contig",
		s->window, s->start, y_offset, s->end, y_offset, 
		s->ruler->colour);
	Tcl_Eval(interp, cmd);
	sprintf(cmd, "%s create text 10 %d -text %s -anchor w -fill %s", 
		s->names_win, t_offset, rf[i], t_colour);
	Tcl_Eval(interp, cmd);
    }

    /* draw last lines */
    y_offset += s->tick->ht;
    sprintf(cmd, "%s create line %d %d %d %d -fill %s -tag contig",
	    s->window, s->start, y_offset, s->end, y_offset, s->ruler->colour);
    Tcl_Eval(interp, cmd);
    if (s->strand == BOTH) {
	y_offset += s->tick->ht;
	sprintf(cmd, "%s create line %d %d %d %d -fill %s -tag contig",
		s->window, s->start, y_offset, s->end, y_offset, 
		s->ruler->colour);
	Tcl_Eval(interp, cmd);
    }

    for (i = 0; i < s->num_match; i++) {
	
	cut_site = s->start - 1 + s->c_match.match[i].cut_pos;

	if (s->strand == BOTH && s->c_match.match[i].enz_name > 3) {
	    y_offset = s->yoffset + ((s->c_match.match[i].enz_name) * s->tick->ht);
	} else {
	    y_offset = s->yoffset + ((s->c_match.match[i].enz_name-1)%3) * s->tick->ht;
	}
	PlotStickMap(interp, s->window, cut_site, cut_site, 0,  y_offset,
		     s->tick->ht, s->tick->line_width, s->tick->colour, 
		     s->c_match.match[i].enz_seq,
		     1, io_clength(io, s->c_match.contig));
    }

    s->world->total->x1 = s->start;
    s->world->total->x2 = s->end;
    s->world->total->y1 = 1;
    s->world->total->y2 = y_offset;

    memcpy(s->world->visible, s->world->total, sizeof(d_box));
    SetCanvasCoords(interp, s->world->visible->x1, s->world->visible->y1, 
		    s->world->visible->x2, s->world->visible->y2, s->canvas);

    draw_single_ruler(interp, s->ruler, s->canvas, s->ruler->start, 
		      s->ruler->end, 1);

    scaleCanvas(interp, s->win_list, s->num_wins, "all", s->world->visible, 
		s->canvas);
    scrollRegion(interp, s->win_list, s->num_wins, s->world->total, s->canvas);

    /* remove all current zooming info */
    freeZoom(&s->zoom);

    /* add first zoom */
    pushZoom(&s->zoom, s->world->visible);
}

int
codon_reg(Tcl_Interp *interp,
	  int strand,
	  GapIO *io,
	  char *frame,
	  char *names_win,
	  char *sc_win,
	  int contig_num,
	  int lreg,
	  int rreg,
	  tick_s *tick,
	  int yoffset,
	  ruler_s *ruler,
	  cursor_s xhair) 
{
    int sequence_type = 0;   /* circle if 1; linear if 0 */
    obj_codon *s;
    int id;
    int num_codons = 3;
    char cmd[1024];


    if (NULL == (s = (obj_codon *)xmalloc(sizeof(obj_codon)))) {
	return 0;
    }
    if (NULL == (s->codon = (char **)xmalloc((num_codons * 2) * 
					     sizeof(char*)))) {
	return 0;
    }

    id = register_id();
    s->id = id;

    strcpy(s->window, sc_win);
    strcpy(s->frame, frame);
    strcpy(s->names_win, names_win);

    s->interp = interp;
    s->xhair = xhair;
    s->ruler = ruler;

    /*
     * create list of windows in the stop codon display.
     */
    if (NULL == (s->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;
    s->num_wins = 0;
    addWindow(s->win_list, &s->num_wins, s->window, 'x', s->id);
    addWindow(s->win_list, &s->num_wins, s->ruler->window, 'x', s->id);
    addWindow(s->win_list, &s->num_wins, s->names_win, 'n', s->id);

    if (NULL == (s->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (s->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (s->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (s->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    initCanvas(interp, s->canvas, s->window);
    createZoom(&s->zoom);

    s->codon[0] = "TAA";
    s->codon[1] = "TGA";
    s->codon[2] = "TAG";
    s->codon[3] = "TTA";
    s->codon[4] = "TCA";
    s->codon[5] = "CTA";

    s->c_match.contig = contig_num;
    s->start = lreg;
    s->end = rreg;
    s->c_match.match = NULL;
    s->num_match = 0;

    s->num_codons = num_codons;
    s->tick = tick;
    s->yoffset = yoffset;
    s->sequence_type = sequence_type;
    s->strand = strand;
    s->editor_avail = type_to_result(io, REG_TYPE_EDITOR, contig_num);

    s->cursor = create_contig_cursor(io, contig_num, 0, id);
    s->cursor_visible = 0;

    sprintf(cmd, "%s.buttons.refresh configure -state %s",
	    s->frame, s->editor_avail ? "normal" : "disabled");
	Tcl_Eval(interp, cmd);

    stop_codon_replot(interp, io, s, NULL);

    contig_register(io, contig_num, stop_codon_callback, (void *)s, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS 
		    | REG_CURSOR_NOTIFY | REG_NUMBER_CHANGE | REG_REGISTERS | 
		    REG_GENERIC, REG_TYPE_STOPCODON);
    canvas_cursor_refresh(s->interp, io, s->c_match.contig,
			  s->cursor, s->cursor,
			  s->canvas, s->win_list, s->num_wins, s->id,
			  0, &s->cursor_visible, s->world, 1);

    return id;
}


/*
 * Remove a result from the codon mapping window
 */
void codon_shutdown(Tcl_Interp *interp,
		    GapIO *io,
		    obj_codon *s)
{
    char cmd[1024];
   
    /*
     * Remove from the registration lists.
     */
    contig_deregister(io, s->c_match.contig, stop_codon_callback, (void *)s);
     
    delete_contig_cursor(io, s->c_match.contig, s->cursor->id, 0);

    /*
     * Pop down configuration window if visible
     */
    sprintf(cmd, "DeleteCodonPlot %s %s\n", s->frame, s->window);
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	printf("%s\n", Tcl_GetStringResult(interp));

    /* Free memory */
    if (s->c_match.match)
	xfree(s->c_match.match);

    if (s->codon)
	xfree(s->codon);

    if (s->canvas)
	xfree(s->canvas);

    if (s->world->visible)
	xfree(s->world->visible);

    if (s->world->total)
	xfree (s->world->total);

    if (s->world)
	xfree(s->world);

    if (s->xhair.colour) free(s->xhair.colour);
    if (s->tick->colour) free(s->tick->colour);
    xfree(s->tick);

    free_win_list(s->win_list, s->num_wins);
    free_ruler_struct(s->ruler);

    freeZoom(&s->zoom);

    xfree(s);
}

/* registration replotting */
int
stop_codon_replot(Tcl_Interp *interp,
		  GapIO *io, 
		  obj_codon *s,
		  char *seq)
{
    char *sequence;
    int sequence_len;
    R_Match *match;
    int total_matches;

    if (seq) {
	sequence_len = strlen(seq);
	sequence = seq;
    } else {
	/* generate a consensus sequence of the contig */
	sequence_len = s->end - s->start + 1;

	if (NULL == (sequence = (char *)malloc((sequence_len + 1) * 
					       sizeof(char ))))
	    return 0;

	calc_consensus(s->c_match.contig, s->start, s->end, CON_SUM, sequence, 
		       NULL, NULL, NULL, consensus_cutoff, quality_cutoff, 
		       database_info, io);
    }

    /* free the old match data */
    if (s->c_match.match)
	xfree(s->c_match.match);

    FindStopCodons(s->strand, sequence, sequence_len, s->sequence_type, 
		   s->num_codons, s->codon, &match, &total_matches);

    s->c_match.match = match;
    s->num_match = total_matches;

    display_stop_codons(interp, io, s);
    
    /* update window size */
    Tcl_VarEval(interp, "update idletasks", NULL);

    if (!seq)
	xfree(sequence);
    return 1;
}

static void
update_obj_codon(GapIO *io, obj_codon *s, int contig) {
    int newl;

    newl = ABS(io_clength(io, s->c_match.contig));

    /* if need to replot, lose original lreg to rreg info */
    s->start = 1;
    s->end = newl;

    /* reset ruler from and to */
    s->ruler->start = 1;
    s->ruler->end = newl;
 
}

/* callback */
static void
stop_codon_callback(GapIO *io, 
		  int contig,
		  void *fdata,
		  reg_data *jdata) 
{
    obj_codon *s = (obj_codon *)fdata;

    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Stop codons");
	    return;
	}
    case REG_JOIN_TO:
#ifdef DEBUG
	printf("stop codon JOIN TO contig %d join %d\n", contig, jdata->join.contig);
#endif
	s->c_match.contig = jdata->join.contig;
	/* NOTE - flow through complement and length */

    case REG_COMPLEMENT:
    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("stop codon COMPLEMENT LENGTH contig %d\n", contig);
#endif
	    update_obj_codon(io, s, contig);
	    stop_codon_replot(s->interp, io, s, NULL);
	    return;
	}
    case REG_QUIT:
    case REG_DELETE:
	{
	    codon_shutdown(s->interp, io, s);
	    return;
	}

    case REG_GET_OPS:
	jdata->get_ops.ops = "Remove\0";
	/*
	   jdata->get_ops.ops = 
	   "Output enzyme by enzyme\0"
	   "Output ordered on position\0"
		"SEPARATOR\0"
		"Remove\0";
		*/
	return;

    case REG_INVOKE_OP:
	switch (jdata->invoke_op.op) {

	case 0: /* Remove */
	    codon_shutdown(s->interp, io, s);
	    break;
	}
	break;
    case REG_PARAMS:
	/* jdata->params.string = s->params; */
	return;

    case REG_CURSOR_NOTIFY:
	canvas_cursor_refresh(s->interp, io, s->c_match.contig,
			      jdata->cursor_notify.cursor, s->cursor,
			      s->canvas, s->win_list, s->num_wins, s->id,
			      0, &s->cursor_visible, s->world, 1);
	return;

    case REG_NUMBER_CHANGE:
#ifdef DEBUG
	printf("stop codon NUMBER_CHANGE %d\n", jdata->number.number);
#endif
	/* FIXME - update tcl saved contig */
	s->c_match.contig = jdata->number.number; 
	return; 

    case REG_REGISTER:
	if (s->editor_avail == 0 &&
	    jdata->c_register.type == REG_TYPE_EDITOR) {
	    char cmd[1024];

	    s->editor_avail = 1;
	    sprintf(cmd, "%s.buttons.refresh configure -state normal",
		    s->frame);
	    Tcl_Eval(s->interp, cmd);
	}
	break;

    case REG_DEREGISTER:
	if (s->editor_avail != 0 &&
	    jdata->c_register.type == REG_TYPE_EDITOR) {
	    s->editor_avail = type_to_result(io, REG_TYPE_EDITOR,
					     s->c_match.contig);
	    if (0 == s->editor_avail) {
		char cmd[1024];

		sprintf(cmd, "%s.buttons.refresh configure -state disabled",
			s->frame);
		Tcl_Eval(s->interp, cmd);
	    }
	}
	break;
    case REG_GENERIC:
	switch (jdata->generic.task) {
	case TASK_CANVAS_SCROLLX: 
	    {
		char *scroll = (char *)jdata->generic.data;
		
		canvasScrollX(s->interp, s->window, s->win_list, s->num_wins,
			      s->world->visible, s->canvas, scroll);
		break;
	    }
	case TASK_CANVAS_ZOOMBACK:
	    {
		canvasZoomback(s->interp, s->canvas, s->window, s->world,
			       s->win_list, s->num_wins, &s->zoom);
		
		/* redraw ruler ticks by redrawing ruler */
		draw_single_ruler(s->interp, s->ruler, s->canvas, 
				  s->ruler->start, s->ruler->end, 1); 
		scaleSingleCanvas(s->interp, s->world, s->canvas, 
				  s->ruler->window, 'x', "all");
		canvas_cursor_refresh(s->interp, io, s->c_match.contig,
				      s->cursor, s->cursor, s->canvas, 
				      s->win_list, s->num_wins, s->id, 0, 
				      &s->cursor_visible, s->world, 0);
		break;
	    }
	case TASK_CANVAS_ZOOM:
	    {
		s_zoom *szoom = (s_zoom *)jdata->generic.data;
		canvasZoom(s->interp, s->canvas, s->window, s->world,
			   s->win_list, s->num_wins, &s->zoom, szoom->zoom,
			   szoom->scroll);
		/* redraw ruler ticks by redrawing ruler */
		draw_single_ruler(s->interp,s->ruler, s->canvas, 
				  s->ruler->start, s->ruler->end, 1);
		
		scaleSingleCanvas(s->interp, s->world, s->canvas, 
				  s->ruler->window, 'x', "all");
		
		canvas_cursor_refresh(s->interp, io, s->c_match.contig,
				      s->cursor, s->cursor, s->canvas, 
				      s->win_list, s->num_wins, s->id, 0, 
				      &s->cursor_visible, s->world, 0);
		break;
	    }
	case TASK_CANVAS_CURSOR_X:
	    {    
		char *label;
		int *cx = (int *)jdata->generic.data;
		double wx, wy;

		CanvasToWorld(s->canvas, *cx, 0, &wx, &wy);
		label = get_default_string(s->interp, gap_defs,
					   "CODON.CURSOR");

		canvasCursorX(s->interp, s->canvas, s->frame,label, 
			      s->xhair.colour, s->xhair.width, *cx, wx, 
			      s->win_list, s->num_wins);

		break;
	    }
	case TASK_CANVAS_RESIZE:
	    {
		char scroll_args[20];
		/* resize stop codon display window */
		resizeCanvas(s->interp, s->window, s->win_list, s->num_wins,
			     s->world->visible, s->world->total, s->canvas);
		sprintf(scroll_args, "scroll 0 units");
		canvasScrollX(s->interp, s->window, s->win_list, s->num_wins,
			      s->world->visible, s->canvas, scroll_args);
	    break;
	    }
	case TASK_CANVAS_WORLD:
	    {
		task_world_t *tw = (task_world_t *)jdata->generic.data;
		int cx = tw->canvasx;
		double wx, wy;
		
		CanvasToWorld(s->canvas, cx, 0, &wx, &wy);
		tw->basex = wx;
		break;
	    }
	}
    }
}


