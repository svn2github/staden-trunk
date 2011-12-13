#include <stdlib.h>
#include <tk.h>
#include <limits.h>

#include "io-reg.h"
#include "template_display.h"
#include "gap-dbstruct.h"
#include "canvas_box.h"
#include "gap_canvas_box.h"
#include "io_utils.h"
#include "tagdb.h"
#include "tagUtils.h"
#include "misc.h"
#include "text_output.h"
#include "template.h"
#include "gap_globals.h"
#include "contig_selector.h"
#include "active_tags.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "vlen.h"
#include "FtoC.h"
#include "ruler_tick.h"
#include "os.h"

/*
 *---------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------
 */
static void template_callback(GapIO *io, int contig, void *fdata,
			      reg_data *jdata);


/*****************************************************************************/
/*                                   CalcXCoords                             */
/*****************************************************************************/
/* calculate the start and end x-coords of a line of initial position pos and
 * length len and convert to pixel coords using factor
 */
void
CalcXCoords(int pos,		                                       /* in */
	    int len,		                                       /* in */
	    int *x1,		                                      /* out */
	    int *x2)		                                      /* out */
{

    *x1 = pos;
    /* *x2 = (pos + len); */
    /* eg reading starts at 1 and is 205 long, goes from position 1 to 205 */
    *x2 = (pos + len - 1);

} /* end CalcXCoords */

/*****************************************************************************/
/*                                CalcYDepth                                 */
/*****************************************************************************/
/* algorithm #2
 * calculate the y depth and the corresponding level each line item will be
 * positioned
 */
void
CalcYDepth(int num,		                                       /* in */
	   PlotRec *PArray,	                                  /* in, out */
	   int max_depth,	                                       /* in */
	   int *num_levels)	                                      /* out */
{
    int i;			                                  /* counter */
    int cnt = 1;		                            /* level counter */
    int found = FALSE;        /* found is TRUE if line been assigned a level */
    int delta_x = 10;		/* distance in pixels between adjacent lines */
    int *level;

    *num_levels = 0;

    /* level array contains the right hand most x2 coord */
    if ((level =  (int *)xmalloc((max_depth+1) * sizeof(int)))==NULL) {
	return;
    }

    /* initialise level array to a large negative number */
    for (i = 1; i <= max_depth; i++){
	level[i] = INT_MIN;
    } /* end for */

    level[cnt] = PArray[0].l.x2;
    PArray[0].l.y1 = cnt;
    PArray[0].l.y2 = cnt;

    /* for each item */
    for (i = 1; i < num; i++) {
	while (!found) {
	    if (level[cnt] > (PArray[i].l.x1) - delta_x) {
		cnt++;
	    } else {
		level[cnt] = PArray[i].l.x2;

		PArray[i].l.y1 = cnt;
		PArray[i].l.y2 = cnt;
		found = TRUE;
	    } /* end if */
	} /* end while */

	if (*num_levels < cnt) {
	    *num_levels = cnt;
	}

	cnt = 1;
	found = FALSE;

    } /* end for */

    if (*num_levels == 0)
	(*num_levels)++;

    xfree(level);
} /* end CalcYDepth */

void CalcYDepthTemplate(int num,	                               /* in */
			DPlotRec **PArray,	                  /* in, out */
			int start_depth,                               /* in */
			int max_depth,	                               /* in */
			int *num_levels)	                      /* out */
{
    int i;			                                  /* counter */
    /* int cnt = 1;*/		                            /* level counter */
    int cnt = start_depth;                                  /* level counter */
    int found = FALSE;        /* found is TRUE if line been assigned a level */
    int delta_x = 10;		/* distance in pixels between adjacent lines */
    int *level;

    *num_levels = 0;

    /* level array contains the right hand most x2 coord */
    if ((level =  (int *)xmalloc((max_depth+1) * sizeof(int)))==NULL) {
	return;
    }

    /* initialise level array to a large negative number */
    for (i = 1; i <= max_depth; i++){
	level[i] = INT_MIN;
    } /* end for */

    level[cnt] = PArray[0]->l.x2;

    PArray[0]->l.y1 = cnt;
    PArray[0]->l.y2 = cnt;

    /* for each item */
    for (i = 1; i < num; i++) {
	while (!found) {
	    if (level[cnt] > (PArray[i]->l.x1) - delta_x) {
		cnt++;
	    } else {
		level[cnt] = PArray[i]->l.x2;
		PArray[i]->l.y1 = cnt;
		PArray[i]->l.y2 = cnt;
		found = TRUE;
	    } /* end if */
	} /* end while */

	if (*num_levels < cnt) {
	    *num_levels = cnt;
	}

	cnt = start_depth;
	found = FALSE;

    } /* end for */

    if (*num_levels == 0)
	*num_levels = cnt;

    xfree(level);
} /* end CalcYDepth */


/*****************************************************************************/
/*                                CalcYCoords                                */
/*****************************************************************************/
/* algorithm #2
 * calculate the y coords for each line segment
*/
void CalcTemplateYCoords(int num,	                               /* in */
			 DPlotRec *PArray,                        /* in, out */
			 int depth,
			 int height)                                   /* in */
{
    int i;			                                  /* counter */
    float yincr;

    yincr = (float)height / (depth + 1);
    /* special case for a single template */
    if (depth == 1) {
	yincr = 20;
    }

    for (i = 0; i < num; i++) {
#ifdef DEBUG
	printf("depth %d yincr %f y1 %f",
	       depth, yincr, PArray[i].l.y1);
#endif
	PArray[i].l.y1 = height - (yincr * PArray[i].l.y1);
	PArray[i].l.y2 = height - (yincr * PArray[i].l.y2);
#ifdef DEBUG
	printf(" y2 %f \n", PArray[i].l.y1);
#endif
    } /* end for */
} /* end CalcYCoords */


/*****************************************************************************/
/*                         CalcReadingYDepth                                 */
/*****************************************************************************/
/* algorithm #2
 * calculate the y depth and the corresponding level each reading will be
 * positioned
 */
void
CalcReadingYDepth(GapIO *io,                                           /* in */
		  int *contig_array,
		  int num_contigs,
		  DPlotRec PArray[],	                          /* in, out */
		  int *num_levels)	                              /* out */
{
    int i, c;			                                  /* counter */
    int cnt = 1;		                            /* level counter */
    int found = FALSE;        /* found is TRUE if line been assigned a level */
    int delta_x = 10;		 /* distance in bases between adjacent lines */
    float *level;		                        /* array of y levels */
    int max_depth;

    max_depth = NumReadings(io);
    *num_levels = 0;

    /* level array contains the right hand most x2 coord */
    if ((level = (float *)xmalloc((max_depth+1) * sizeof(float)))==NULL){
	return;
    }

    /* initialise level array to a large negative number */
    for (i = 0; i <= max_depth; i++){
	level[i] = INT_MIN;
    } /* end for */

    for (c = 0; c < num_contigs; c++) {
	for (i=io_clnbr(io, contig_array[c]); i; i=io_rnbr(io, i)) {

	    if (PArray[i].l.x1 == 0 && PArray[i].l.x2 == 0)
		continue;

	    while (!found) {
		if (level[cnt] > (PArray[i].l.x1) - delta_x) {
		    cnt++;
		} else {
		    level[cnt] = PArray[i].l.x2;
		    PArray[i].l.y1 = cnt;
		    PArray[i].l.y2 = cnt;
		    found = TRUE;
		} /* end if */
	    } /* end while */
	    if (*num_levels < cnt) {
		*num_levels = cnt;
	    }
	    cnt = 1;
	    found = FALSE;
	} /* end for */
    }

    if (*num_levels == 0)
	(*num_levels)++;
    xfree(level);

} /* end CalcReadingYDepth */


/*****************************************************************************/
/*                         CalcReadingYCoords                                */
/*****************************************************************************/
/* algorithm #2
 * calculate the y coords for each line segment
*/
void
CalcReadingYCoords(GapIO *io,                                          /* in */
		   int *contig_array,                                  /* in */
		   int num_contigs,
		   DPlotRec *RArray,	                               /* in */
		   DPlotRec *ReadArray,	                              /* out */
		   int depth,
		   int height,
		   int *num_readings)                                 /* out */
{
    int i, c;
    float yincr;

    *num_readings = 0;
    yincr = (float)height / (depth + 1);
    /* special case for a single reading */
    if (depth == 1) {
	yincr = 20;
    }

    for (c = 0; c < num_contigs; c++) {
	for (i = io_clnbr(io, contig_array[c]); i; i = io_rnbr(io, i)) {

	    if (RArray[i].colour) {
		RArray[i].l.y1 = height - (yincr * RArray[i].l.y1);
		RArray[i].l.y2 = height - (yincr * RArray[i].l.y2);

		ReadArray[(*num_readings)++] = RArray[i];
		RArray[i].type = NULL;
	    }
	} /* end for */
    }
} /* end CalcReadingYCoords */

int
inContigList(int *contig_array,
	     int num_contigs,
	     int contig)
{
    int i;
    for (i = 0; i < num_contigs; i++) {
	if (contig == contig_array[i]) {
	    return 1;
	}
    }
    return 0;
}

int
getContigIndex(int *contig_array,
	       int num_contigs,
	       int contig)
{
    int i;
    for (i = 0; i < num_contigs; i++) {
	if (contig == contig_array[i]) {
	    return i;
	}
    }
    return -1;
}

/* template for left contig has to go to the right
 * AND
 * template for right contig has to go to the left
 */
int
checkTemplateConsistency(template_c *t_l,
			 template_c *t_r)
{
    /* forward reading, guessed end */
    /*   ---->                      */
    /* st------------------ end     */

    /* reverse reading, guessed start */
    /*    ---->                       */
    /* end------------------start     */

    /* forward reading, guessed start */
    /*               <-----           */
    /* st------------------ end       */

    /* reverse reading, guessed end   */
    /*               <------          */
    /* end------------------start     */

    if ((t_l->start <= t_l->end) == (t_r->start <= t_r->end)) {

#ifdef DEBUG
	printf("left and right contig are consistent! %d %d \n", t_l->num, t_r->num);
#endif
	return(1);
    }
    
#ifdef DEBUG
	printf("left and right contig are NOT consistent! %d %d\n", t_l->num, t_r->num);
#endif
    return(0);
}

int
FindSpanningReadPair(template_c *t_l,
		     template_c *t_r)
{

    if ((t_l->flags & TEMP_FLAG_GUESSED_END) && (t_r->flags & TEMP_FLAG_GUESSED_START))
	return 1;

   if ((t_l->flags & TEMP_FLAG_GUESSED_START) && (t_r->flags & TEMP_FLAG_GUESSED_END))
       return 1;

    return 0;
}

int
FindSpanningTemplates(GapIO *io,
		      template_c **tarr,
		      int *contig_array,
		      int num_contigs,
		      template_p *t_pos,
		      int *t_num)
{
    template_c *t;
    int i;
    item_t *it, *new_contig_l;
    int cnt = 0;
    int num_r;
    int tmp;
    int num_readings = NumReadings(io);
    int contig;
    gel_cont_t *gc;

    for (i = 0; i < num_readings; i++) {
	t_pos[i].t = NULL;
	t_pos[i].start = 0;
	t_pos[i].end = 0;
	t_pos[i].consist = -1;
	t_pos[i].num_r = 0;
    }

    check_all_templates(io, tarr);

    /*
     * find all the templates which span contigs in our contig_array
     */
    for (i = 1; i <= Ntemplates(io); i++){
	t = tarr[i];

	/* if template in contig_array && */
	/* if template spans two or more contigs */
	if (t && t->flags & TEMP_FLAG_SPANNING) {

	    /* Find first valid gel_cont in this item list */
	    new_contig_l = head(t->gel_cont);
	    for (it = new_contig_l; it; it = it->next) {
		gel_cont_t *gc = (gel_cont_t *)(it->data);
		if (inContigList(contig_array, num_contigs,
				 gc->contig)) {
		    break;
		}
	    }
	    new_contig_l = it;

	    /*
	     * Loop around the list starting at new_contig_l looking for
	     * all readings in this contig (which are implicitly on this
	     * template).
	     *
	     * We call get_template_positions once and create one
	     * t_pos array entry provided that there are 2 or more
	     * reads on this template from this contig.
	     */
	    num_r = 0;
	    while (new_contig_l) {
		contig = ((gel_cont_t *)(new_contig_l->data))->contig;
		it = new_contig_l;
		new_contig_l = NULL;

		/* Get the template positions & consistency */
		tmp = t->consistency;
		get_template_positions(io, t, contig);
		t->consistency |= tmp;

		/*
		 * If we find a reading on this template for another contig
		 * then we remember that list point ready for the next
		 * loop.
		 *
		 * We also mark the other readings for this contig
		 * so that we avoid processing them again later.
		 * We do this by temporarily negating gc->contig.
		 */
		for (; it; it = it->next) {
		    gel_cont_t *gc = (gel_cont_t *)(it->data);

		    if (gc->contig == contig) {
			gc->contig = -gc->contig;
		    } else {
			if (!new_contig_l &&
			    gc->contig > 0 &&
			    inContigList(contig_array, num_contigs,
					     gc->contig)) {
			    new_contig_l = it;
			}
		    }
		}

		/*
		 * Allocate an entry in the t_pos array and initialise
		 * the bits that we know.
		 * num_r here is redundant and we could remove it if
		 * we fix the code that uses it.
		 */
		if (t_pos[cnt].t == NULL) {
		    if ((t_pos[cnt].t =
			 (template_c *)xmalloc(sizeof(template_c)))==NULL) {
			return -1;
		    }
		}
		memcpy(t_pos[cnt].t, t, sizeof(template_c));


		t_pos[cnt].contig = contig;
		t_pos[cnt].t_num = i;
		t_pos[cnt].start = 0; /* initialise for later */
		t_pos[cnt].end = 0; /* initialise for later */
		cnt++;
		num_r++;
	    } /* while(new_contig_l) */

	    /*
	     * Only keep this template information if there are at
	     * least 2 different reads on it.
	     */
	    if (num_r < 2) {
		cnt--;
	    } else {
		t_pos[cnt-num_r].num_r = num_r;
	    }

	    /*
	     * Now make all contigs in the read/contig pairs
	     * positive again.
	     */
	    for (it = head(t->gel_cont); it; it = it->next) {
		gc = (gel_cont_t *)(it->data);
		if (gc->contig < 0)
		    gc->contig = -gc->contig;
	    }
	}
    }
    *t_num = cnt;
    return 0;
}

/*
 * go through spanning templates from left to right contig
 * decide if the template is consistent across the two contigs ie is
 * going in the correct direction
 * if it is consistent, find the gap between contigs and start and end of
 *     the template
 * if it is not consistent, find the "new" start and end of the template
 */
void
FindSpanningTemplatePositions(GapIO *io,
			      int *contig_array,
			      int num_contigs,
			      template_p *t_pos,
			      int num_templates,
			      template_o *offset)
{
    int i, j;
    int c1_index, c2_index;
    int len1, end1, end2;

    for (i = 0; i < num_templates; i++) {

	/*loop through each reading on template */
	for (j = i+1; j < i+t_pos[i].num_r; j++) {

	    c1_index = getContigIndex(contig_array, num_contigs,
				      t_pos[i].contig);
	    c2_index = getContigIndex(contig_array, num_contigs,
				      t_pos[j].contig);

	    /* check contigs are next to one another */
	    if (abs(c1_index - c2_index) == 1) {
		if (c1_index < c2_index) {

		    t_pos[i].consist = checkTemplateConsistency(t_pos[i].t, t_pos[j].t);
		    t_pos[j].consist = t_pos[i].consist;
#ifdef DEBUG
		    printf("CONSIST %d %d\n", t_pos[i].t_num, t_pos[i].consist);

		    printf("i %d start %d end %d min %d max %d\n",
			   i, t_pos[i].t->start, t_pos[i].t->end,
			   t_pos[i].t->min, t_pos[i].t->max);
		    printf("j %d start %d end %d min %d max %d\n",
			   j, t_pos[j].t->start, t_pos[j].t->end,
			   t_pos[j].t->min, t_pos[j].t->max);
#endif
		    if (t_pos[i].consist == 1) {
			len1 = io_clength(io, t_pos[i].contig);
			end1 = MAX(MAX(t_pos[i].t->start,
					   t_pos[i].t->end), t_pos[i].t->max);
			end2 = MAX(MAX(t_pos[j].t->start,
				       t_pos[j].t->end), t_pos[j].t->max);

			t_pos[j].length = (end1 - len1) - end2;

			if (offset) {
			    offset[c2_index].gap += (end1 - len1) - end2;
			    offset[c2_index].cnt++;
			}
			t_pos[i].start = MIN(MIN(t_pos[i].t->start,
						 t_pos[i].t->end),
					     t_pos[i].t->min);
			t_pos[j].end = end2;
#ifdef DEBUG
			printf("len1 %d end1 %d end2 %d\n", len1, end1, end2);
			printf("start %d end %d len %d offset %d\n",
			       t_pos[i].start, t_pos[j].end, t_pos[j].length, (end1 - len1) - end2);

#endif
		    }
		} else {
#ifdef DEBUG
		    printf("********************GOT HERE *************\n");
#endif
		    t_pos[i].consist = checkTemplateConsistency(t_pos[j].t, t_pos[i].t);
		    t_pos[j].consist = t_pos[i].consist;
		    /*
		     * if the contigs are not consistent with respect to
		     * the direction of the templates, draw a line between
		     * the two readings
		     */

		    if (t_pos[i].consist == 1) {
			len1 = io_clength(io, t_pos[j].contig);
			end1 = MAX(MAX(t_pos[j].t->start,
				       t_pos[j].t->end), t_pos[j].t->max);

			end2 = MAX(MAX(t_pos[i].t->start,
				       t_pos[i].t->end), t_pos[i].t->max);
			t_pos[i].length = (end1 - len1) - end2;

			if (offset) {
			    offset[c1_index].gap += (end1 - len1) - end2;
			    offset[c1_index].cnt++;
			}
			t_pos[j].start = MIN(MIN(t_pos[j].t->start,
						 t_pos[j].t->end),
					     t_pos[j].t->min);
			t_pos[i].end = end2;
		    }
		}
	    } else if (abs(c1_index - c2_index) > 1) {
		/*
		 * templates spanning non-adjacent contigs in contig_array
		 */
#ifdef DEBUG
		printf("i %d j %d min %d max %d \n", i, j, t_pos[i].t->min,
		       t_pos[j].t->max);
#endif
		t_pos[i].start = t_pos[i].t->min;
		t_pos[j].end = t_pos[j].t->max;
	    }
	}
    }
}

void
FindTemplatePositionChanges(GapIO *io,
			    c_offset *contig_offset,
			    template_p *t_pos,
			    int num_templates,
			    template_d *t_changes)
{
    int i, j;
    int start1, start2, end1, end2;
    int consist;

    /* indexed on template number */
    for (i = 1; i <= Ntemplates(io); i++) {
	t_changes[i].diff = 0;
	t_changes[i].start = 0;
	t_changes[i].end = 0;
	t_changes[i].consist = 1;
	t_changes[i].readpair = 0;
    }

    /* find start and end for non consistent templates */
    for (i = 0; i < num_templates; i++) {
	for (j = i+1; j < i+t_pos[i].num_r; j++) {

#ifdef DEBUG
	    printf("num %d i %d j %d num_r %d \n", t_pos[i].t_num, i, j, t_pos[i].num_r);
	    printf("consist %d %d \n", t_pos[i].consist, t_pos[j].consist);
#endif
	    t_changes[t_pos[i].t_num].readpair =
	      FindSpanningReadPair(t_pos[i].t, t_pos[j].t);

	    if ((!t_pos[i].consist && !t_pos[j].consist) ||
		(t_pos[i].consist == -1 && t_pos[j].consist == -1)) {

		if (t_pos[i].contig == t_pos[j].contig)
		    continue;

		t_changes[t_pos[i].t_num].consist = 0;

		/* spanning templates not adjacent to each other */
		if (t_pos[i].consist == -1 && t_pos[j].consist == -1) {
		    consist = checkTemplateConsistency(t_pos[i].t, t_pos[j].t);
		    if (consist) {
			/* non adjacent but consistent readings */
			t_changes[t_pos[i].t_num].consist = 2;
		    } else {
			/* non adjacent and inconsistent readings */
			t_changes[t_pos[i].t_num].consist = 0;
		    }
		}

		start1 = t_pos[i].t->min + contig_offset[t_pos[i].contig].offset;
		start2 = t_pos[j].t->min + contig_offset[t_pos[j].contig].offset;
		end1 = t_pos[i].t->max + contig_offset[t_pos[i].contig].offset;
		end2 = t_pos[j].t->max + contig_offset[t_pos[j].contig].offset;
		t_pos[i].start = t_pos[j].start = 0;
		t_pos[i].end = t_pos[j].end = 0;

#ifdef DEBUG
		printf("s1 %d s2 %d e1 %d e2 %d \n", start1, start2, end1, end2);
		printf("i %d min %d max %d j %d min %d max %d\n", i,
		       t_pos[i].t->min, t_pos[i].t->max, j, t_pos[j].t->min,
		       t_pos[j].t->max);
#endif
		if (start1 < start2)
		    t_pos[i].start = t_pos[i].t->min;
		else
		    t_pos[j].start = t_pos[j].t->min;

		if (end1 > end2)
		    t_pos[i].end = t_pos[i].t->max;
		else
		    t_pos[j].end = t_pos[j].t->max;

	    }
	}
    }
    /*
     * determine the new start and ends of spanning templates
     */
    for (i = 0; i < num_templates; i++) {
	/*
	 * note here that the offset will be different for the start and
	 * end because this depends on the contig the end of the template
	 * is in
	 */
	if (t_pos[i].start) {
	    t_changes[t_pos[i].t_num].start = t_pos[i].start +
		contig_offset[t_pos[i].contig].offset;
	}
	if (t_pos[i].end) {
	    t_changes[t_pos[i].t_num].end = t_pos[i].end +
		contig_offset[t_pos[i].contig].offset;
	}
#ifdef DEBUG
	printf("num %d start %d end %d\n", t_pos[i].t_num,
	       t_changes[t_pos[i].t_num].start,
	       t_changes[t_pos[i].t_num].end);
#endif
    }

}

int FindTemplatePositions(GapIO *io,
			  c_offset *contig_offset,
			  int *contig_array,
			  int num_contigs,
			  template_c **tarr,
			  template_d **t_changes)
{
    int i;
    template_p *t_pos = NULL;
    int num_templates, num_readings;

    if ((t_pos = (template_p *)xmalloc(NumReadings(io) *
				       sizeof(template_p)))==NULL){
	return -1;
    }

    if (NULL == (*t_changes = (template_d *)xmalloc((Ntemplates(io)+1) *
						   sizeof(template_d)))){
	return -1;
    }

    /* indexed on template number */
    for (i = 1; i <= Ntemplates(io); i++) {
	(*t_changes)[i].diff = 0;
	(*t_changes)[i].start = 0;
	(*t_changes)[i].end = 0;
	(*t_changes)[i].consist = 1;
	(*t_changes)[i].readpair = 0;
    }

    FindSpanningTemplates(io, tarr, contig_array, num_contigs, t_pos,
			  &num_templates);

    FindSpanningTemplatePositions(io, contig_array, num_contigs, t_pos,
				  num_templates, NULL);

    FindTemplatePositionChanges(io, contig_offset, t_pos, num_templates,
				*t_changes);

    num_readings = NumReadings(io);
    if (t_pos) {
	for (i = 0; i < num_readings; i++) {
	    if (t_pos[i].t)
		xfree(t_pos[i].t);
	}
	xfree(t_pos);
    }
    return 0;
}


int contigOffsets(GapIO *io,
		  template_c **tarr,
		  c_offset *contig_offset,
		  int *contig_array,
		  int num_contigs,
		  int calc_contig_pos,
		  template_d *t_changes)
{
    int i, j, k;
    template_p *t_pos = NULL;
    template_o *offset = NULL;
    int num_templates;
    int num_readings;

     vfuncgroup(2, "Template display");

    if ((t_pos = (template_p *)xmalloc(NumReadings(io) *
				       sizeof(template_p)))==NULL){
	return -1;
    }

    /*
     * find all the templates which span contigs in our contig_array
     */

    FindSpanningTemplates(io, tarr, contig_array, num_contigs, t_pos,
			  &num_templates);

    if ((offset = (template_o *)xmalloc(num_contigs * sizeof(template_o)))==NULL){
	return -1;
    }
    for (i = 0; i < num_contigs; i++) {
	offset[i].gap = 0;
	offset[i].cnt = 0;
    }

    FindSpanningTemplatePositions(io, contig_array, num_contigs, t_pos,
				  num_templates, offset);

    /*
     * find the average gap between contigs
     */
    contig_offset[contig_array[0]].offset = 0;

    for (i = 1; i < num_contigs; i++) {
	if (!calc_contig_pos)
	    offset[i].gap = 0;

	if (offset[i].gap == 0) {
	    offset[i].average = 0;
	} else {
	    offset[i].average = (double)offset[i].gap / offset[i].cnt;
	}
	contig_offset[contig_array[i]].gap = (int)offset[i].average;

	contig_offset[contig_array[i]].offset =
	    contig_offset[contig_array[i-1]].offset +
		(int)offset[i].average +
		    ABS(io_clength(io, contig_array[i-1]));
#ifdef DEBUG
	printf("i %d contig %d gap %d num_templates %d ave %f offset %d\n",
	       i,
	       contig_array[i],
	       offset[i].gap,
	       offset[i].cnt,
	       offset[i].average,
	       contig_offset[contig_array[i]].offset);
#endif
    }

    /*
     * determine the new start and ends of spanning templates
     */
    FindTemplatePositionChanges(io, contig_offset, t_pos, num_templates,
				t_changes);

    for (i = 1; i < num_contigs; i++) {
	GTemplates t;
	GReadings r;
	item_t *item;
	char name[DB_NAMELEN + 1];
	char name2[DB_NAMELEN + 1];
	int length;

	strcpy(name, get_contig_name(io, ABS(contig_array[i-1])));
	strcpy(name2, get_contig_name(io, ABS(contig_array[i])));
	vmessage("Contig %s(%d) and Contig %s(%d) \n",
		 name, io_clnbr(io, ABS(contig_array[i-1])),
		 name2,io_clnbr(io, ABS(contig_array[i])));

	/* print out the spanning templates */
	for (j = 0; j < num_templates; j++) {

	    /*
	     * only write out info for those templates spanning adjacent
	     * contigs and are "consistent" ie used in the calculation of
	     * distance between them
	     */
	    if (t_pos[j].contig == contig_array[i-1]) {

		for (k = j+1; k < j+t_pos[j].num_r; k++) {
		    if (t_pos[k].contig == contig_array[i] &&
			t_changes[t_pos[j].t_num].consist) {

			template_read(io, t_pos[j].t_num, t);

			TextRead(io, t.name, name, sizeof(name)-1);
			length = t_changes[t_pos[j].t_num].end -
			    t_changes[t_pos[j].t_num].start + 1;

			vmessage("Template %12s(%4d) length %d\n",
				     name, t_pos[j].t_num, length);

			for (item = head(t_pos[j].t->gel_cont); item;
			     item = item->next) {
			    gel_cont_t *gc = (gel_cont_t *)item->data;

			    strcpy(name, io_rname(io, gc->read));
			    gel_read(io, gc->read, r);

			    vmessage("Reading %*s(%+5d%c), pos %6d%+5d, contig %4d\n",
				     DB_NAMELEN, name,
				     gc->read * (r.sense ? -1 : 1),
				     "?FRfr"[PRIMER_TYPE_GUESS(r)],
				     r.position, r.end - r.start - 1,
				     chain_left(io, gc->read));

			}
		    }
		}
	    }
	}

	vmessage("Gap between contigs = %d\n",
		 contig_offset[contig_array[i]].gap);
	vmessage("Offset of contig %s (%d) from the beginning = %d\n\n",
		 name2, 
		 io_clnbr(io, ABS(contig_array[i])),
		 contig_offset[contig_array[i]].offset);
    }

    if (offset)
	xfree(offset);

    num_readings = NumReadings(io);
    if (t_pos) {
	for (i = 0; i < num_readings; i++) {
	    if (t_pos[i].t)
		xfree(t_pos[i].t);
	}
	xfree(t_pos);
    }
    return 0;
}

int
CalcContigOffsets(GapIO *io,
		  c_offset *contig_offset,  /* out */
		  int *contig_array,
		  int num_contigs,
		  int calc_contig_pos,
		  template_c ***tarr,        /* out */
		  template_d **t_changes)   /* out */
{
    int i;

    if (Ntemplates(io) == 0)
	return -1;

    /* Initialise templates */
    if (NULL == (*tarr = init_template_checks(io, num_contigs, contig_array,
					      1)))
	return -1;

    /*
     * This is how we could turn on checking of inter-contig distances
     * between sequences.
     *
     * template_check_set_flags(io, *tarr, TEMP_OFLAG_INTERDIST, 0);
     */

    check_all_templates(io, *tarr);

    if (NULL == (*t_changes = (template_d *)xmalloc((Ntemplates(io)+1) *
						   sizeof(template_d)))){
	return -1;
    }

    /*
     * store changes to template info
     */
    for (i = 1; i <= Ntemplates(io); i++){
	(*t_changes)[i].diff = 0;
	(*t_changes)[i].start = 0;
	(*t_changes)[i].end = 0;
	(*t_changes)[i].consist = 1;
    }

    /* find "proper" offsets between adjacent contigs in the display */
    if (-1 == (contigOffsets(io, *tarr, contig_offset, contig_array,
			     num_contigs, calc_contig_pos, *t_changes))) {
	return -1;
    }
   return 0;
}


void
templatePosition(template_c *t,
		 DPlotRec *TArray,
		 int t_num,
		 int cnt,
		 int c_num,
		 int status,
		 int start,
		 int end,
		 int *min_x1,
		 int *max_x2,
		 DPlotRec **TArray1, /* out */
		 DPlotRec **TArray2, /* out */
		 int *cnt1,          /* out */
		 int *cnt2)          /* out */

{
    char *colour;

    TArray[cnt].num = t_num;

/*    sprintf(TArray[cnt].type, "{template te_%d c_%d S}", t_num, c_num); */

    TArray[cnt].l.x1 = start;
    TArray[cnt].l.x2 = end;
#ifdef DEBUG
    printf("t_num %d status %d\n", t_num, status);
#endif
    if (status & CONTRADICTORY) {
#ifdef DEBUG
	printf("colour contradictory\n");
#endif
	colour =
	    get_default_string(GetInterp(), gap_defs,
				"TEMPLATE.CONTRADICT_COLOUR");
    } else if (status & ONE_READING) {
#ifdef DEBUG
	printf("colour one reading\n");
#endif
	colour =
	    get_default_string(GetInterp(), gap_defs,
				"TEMPLATE.ONE_READING_COLOUR");
    } else if (status & DIFF_CONTIGS) {
#ifdef DEBUG
	printf("colour diff contigs\n");
#endif
	colour =
	    get_default_string(GetInterp(), gap_defs,
				"TEMPLATE.DIFF_CONTIGS_COLOUR");
    } else if (status & FORW_REV_READINGS) {
#ifdef DEBUG
	printf("colour forw rev readings\n");
#endif
	colour =
	    get_default_string(GetInterp(), gap_defs,
				"TEMPLATE.FORW_REV_COLOUR");
    } else if (status & SPAN_CONTIG) {
#ifdef DEBUG
	printf("colour span contig\n");
#endif
	/* spanning contig, consistent */
	colour =
	    get_default_string(GetInterp(), gap_defs,
				"TEMPLATE.SPAN_CONTIG_COLOUR");
    } else {
	/* spanning contig, inconsistent */
#ifdef DEBUG
	printf("colour span contig inconsistent\n");
#endif
	colour =
	    get_default_string(GetInterp(), gap_defs,
				"TEMPLATE.SPAN_CONTIG_INCONS_COLOUR");
    }

    TArray[cnt].colour = colour;
    strcpy(TArray[cnt].arrow, "none");

    /* determine the right hand most x coord */
    if (*max_x2 < TArray[cnt].l.x2) {

	*max_x2 = TArray[cnt].l.x2;

    } /* end if */

    /* determine the left hand most x coord */
    if (*min_x1 > TArray[cnt].l.x1) {

	*min_x1 = TArray[cnt].l.x1;

    } /* end if */

    if (status & CONTRADICTORY || status & SPAN_CONTIG_INCONS) {
        TArray2[(*cnt2)++] = &TArray[cnt];
    } else {
        TArray1[(*cnt1)++] = &TArray[cnt];
    }
}

/*
 * Set status (corresponds to colour in template display).
 * The order here counts as it sets the precedence for which
 * colour should be displayed.
 */
int
getStatus(template_c *t)
{
    int status;

    /* default value eg when two readings at the same end of the contig */
    status = ONE_READING;

    if (head(t->gel_cont)->next == NULL) {
#ifdef DEBUG
	printf("status = ONE_READING %d\n", t->num);
#endif
	status = ONE_READING;
    }
    if (!(t->flags & (TEMP_FLAG_GUESSED_START |
		      TEMP_FLAG_GUESSED_END))) {

#ifdef DEBUG
	printf("status = FORW_REV_READINGS %d\n", t->num);
#endif
	status = FORW_REV_READINGS;
    }

    if (t->flags & TEMP_FLAG_SPANNING){
#ifdef DEBUG
	printf("status = DIFF_CONTIGS %d\n", t->num);
#endif
	status = DIFF_CONTIGS;
    }

    if (t->consistency & (TEMP_CONSIST_DIST |
			  TEMP_CONSIST_PRIMER |
			  TEMP_CONSIST_STRAND)) {
#ifdef DEBUG
	printf("status = CONTRADICTORY %d\n", t->num);
#endif
	status |= CONTRADICTORY;
    }
    return(status);
}

int
comparePlotRec(DPlotRec **r1, DPlotRec **r2)
{
    if ((*r1)->l.x1 < (*r2)->l.x1)
	return (-1);
    else if ((*r1)->l.x1 == (*r2)->l.x1)
	return (0);
    else
	return (1);
}

/*****************************************************************************/
/*                                CalcTemplates                              */
/*****************************************************************************/
/* calculates the x & y coords of each template
   assigns template_num to be the template index
   assigns type to be "template"
   assigns the colour dependent on the status of the template ie
   reading at both ends in contig, one reading in a different contig, just one
   reading

   Returns an allocated array with an element for each reading detailing
   whether it's template is displayed.
*/
int
CalcTemplates(GapIO *io,	                                      /*  in */
	      c_offset *contig_offset,
	      int *contig_array,
	      int num_contigs,	                                      /*  in */
	      template_d *t_changes,                                  /*  in */
	      template_c **tarr,                                      /* out */
	      DPlotRec *TArray,	                                      /* out */
	      DPlotRec **TArray1,                                     /* out */
	      DPlotRec **TArray2,                                     /* out */
	      int *num_templates,                                     /* out */
	      int *min_x1,	                                      /* out */
	      int *max_x2,	                                      /* out */
	      int *max_y2,                                            /* out */
	      int not_single,	                                      /*  in */
	      int read_pairs,
	      int span_read_pairs,
	      int consist_read_pairs,
	      int win_height)	                                      /*  in */
{
    int i;			                                  /* counter */
    int depth, max_depth;   /* depth & max depth of templates in y direction */
    int cnt = 0;	                                 /* template counter */
    template_c *t;
    int status;
    int offset;
    gel_cont_t *gc;
    gel_cont_t *c_num;
    int x1, x2;
    item_t *item;
    int prev_cnum;
    Tcl_DString tmp;
    int len;
    int cnt1 = 0;
    int cnt2 = 0;
    int ntemplates = Ntemplates(io);

    if (ntemplates == 0)
	return -1;

    /* for each template */
    for (i = 1; i <= ntemplates; i++){
	t = tarr[i];

	/* if the template is in the contig */
	if (t) {
	    item_t *it;

	    /* Check if we only want consistent read pairs shown */
	    if (t->consistency && consist_read_pairs)
		continue;

	    /*
	     * user display configurations: display templates with > 1 reading,
	     * display templates with read pairs
	     */
	    status = getStatus(t);

	    
	    /* first check if template has > 1 reading */
	    if ((not_single || read_pairs || span_read_pairs) &&
		(NULL == head(t->gel_cont) || NULL == head(t->gel_cont)->next)) {
		continue;
	    }

#ifdef DEBUG
	    printf("TEMPLATE readpair %d forw %d contra %d diff %d\n", 
		   t_changes[i].readpair, 
		   status & FORW_REV_READINGS, status & CONTRADICTORY,
		   status & DIFF_CONTIGS);
#endif
	    if ((read_pairs || span_read_pairs) && !t_changes[i].readpair &&
		!(status & FORW_REV_READINGS) && !(status & CONTRADICTORY)) {
		continue;
	    }


	    /* display read pairs that span contigs */
	    if (span_read_pairs && !(status & DIFF_CONTIGS)) {
		continue;
	    }

#ifdef DEBUG
	    if (t->flags & TEMP_FLAG_SPANNING) {
		for (it = head(t->gel_cont); it; it=it->next) {
		    gc = (gel_cont_t *)(it->data);

		    get_template_positions(io, t, gc->contig);
		    printf("%3d: %04d-%04d, %04d-%04d, 0x%02x, 0x%x:",
			   t->num, t->start, t->end, t->min, t->max,
			   t->consistency, t->flags);
		    printf(" %02d.%03d\n", gc->contig, gc->read);
		}
	    }
#endif
	    it = head(t->gel_cont);
	    gc = (gel_cont_t *)(it->data);
	    offset = contig_offset[gc->contig].offset;

	    /* templates spanning contigs in contig_array */
	    if (t_changes[i].start && t_changes[i].end) {
		GTemplates te;

#ifdef DEBUG
		printf("templates spanning contigs in contig array %d %d\n", t->num, status);
#endif
		/* check to see if new template length is within range */
		template_read(io, i, te);

		if (t_changes[i].start < t_changes[i].end) {
		    x1 = t_changes[i].start;
		    x2 = t_changes[i].end;
		} else {
		    x1 = t_changes[i].end;
		    x2 = t_changes[i].start;
		}
		if (t_changes[i].consist == 0) {
		    status = SPAN_CONTIG_INCONS;
		} else if (status & CONTRADICTORY) {
		    status = SPAN_CONTIG_INCONS;
		} else {
		    if (x2 - x1 < te.insert_length_min ||
			x2 - x1 > te.insert_length_max) {
			status = SPAN_CONTIG_INCONS;
		    } else {
			status = SPAN_CONTIG;
		    }
		}
		prev_cnum = -1;
		Tcl_DStringInit(&tmp);

		vTcl_DStringAppend(&tmp, "{template te_%d ", i);

		for (item = head(t->gel_cont); item; item=item->next) {
		    c_num = (gel_cont_t *)(item->data);

		    if (c_num->contig != prev_cnum) {
			vTcl_DStringAppend(&tmp, "c_%d ", c_num->contig);
			prev_cnum = c_num->contig;
		    }

		}

		/* strcat(TArray[cnt].type, " S}"); */
		vTcl_DStringAppend(&tmp, " S}");
		len = strlen(Tcl_DStringValue(&tmp));
		if (NULL == (TArray[cnt].type = (char *)xmalloc(len+1))) {
		    verror(ERR_FATAL, "CalcTemplates", "out of memory");
		    return -1;
		}
		sprintf(TArray[cnt].type, "%s", Tcl_DStringValue(&tmp));
		Tcl_DStringFree(&tmp);

		/* templatePosition(t, TArray, i, cnt, gc->contig, status,
		   x1, x2, min_x1, max_x2); */
		templatePosition(t, TArray, i, cnt, gc->contig, status,
				 x1, x2, min_x1, max_x2, TArray1, TArray2,
				 &cnt1, &cnt2);
		cnt++;
	    } else {
		/* 
		 * templates not spanning contigs or spanning contigs which 
		 * are not displayed 
		 */
		if (span_read_pairs) 
		    continue;

#ifdef DEBUG
		printf("templates not spanning contigs or spanning contigs which are not displayed %d %d\n", t->num, status);
#endif

		if (t->flags & TEMP_FLAG_SPANNING) {
		    /* other spanning templates */
		    get_template_positions(io, t, gc->contig);
		}

		x1 = MIN(MIN(t->start, t->end), t->min) + offset;
		x2 = MAX(MAX(t->start, t->end), t->max) + offset;

#ifdef DEBUG
		printf("start %d end %d min %d max %d\n", t->start, t->end, t->min, t->max);

		printf("2 i %d c %d x1 %d x2 %d offset %d consist %d\n",
		       i, gc->contig, x1, x2, offset, t->consistency);
#endif
		if (NULL == (TArray[cnt].type = (char *)xmalloc(40))) {
		    verror(ERR_FATAL, "CalcTemplates", "out of memory");
		    return -1;
		}
		sprintf(TArray[cnt].type, "{template te_%d c_%d S}",
			i, gc->contig);
		templatePosition(t, TArray, i, cnt, gc->contig, status,
				 x1, x2, min_x1, max_x2, TArray1, TArray2,
				 &cnt1, &cnt2);
		cnt++;
	    }
	}

    } /* end for */

    *num_templates = cnt;

    qsort((void *) TArray1, cnt1, sizeof(DPlotRec*),
	  (int (*)(const void *, const void *))comparePlotRec);

    qsort((void *) TArray2, cnt2, sizeof(DPlotRec*),
	  (int (*)(const void *, const void *))comparePlotRec);

    max_depth = Ntemplates(io)+1;

    /* determine actual number of y levels required */
    /* y spacing alg #2: columns of readings of different heights */

    CalcYDepthTemplate(cnt1, TArray1, 1, max_depth, &depth);
    CalcYDepthTemplate(cnt2, TArray2, depth+1, max_depth, &depth);

    CalcTemplateYCoords(*num_templates, TArray, depth, win_height);
    *max_y2 = win_height;

    return 0;
} /* end CalcTemplates */


/*****************************************************************************/
/*                          FindReadingYCoords                               */
/*****************************************************************************/
/* find the y coords of a reading from the template the reading is on */
void
FindReadingYCoords(GapIO *io,	                                       /* in */
		   template_c **tarr,
		   DPlotRec *TArray,                                  /* in */
		   DPlotRec *RArray,                                  /* in */
		   DPlotRec *ReadArray,                               /* out */
		   int *num_readings,                                 /* out */
		   int num_templates)                                  /* in */
{
    int i, t;			                                 /* counters */
    int gel;
    item_t *ip;

    for (i = 0; i < num_templates; i++){

	t = TArray[i].num;

	for (ip = head(tarr[t]->gel_cont); ip; ip = ip->next) {
	    gel = ((gel_cont_t *)ip->data)->read;

	    /*
	     * need this check to guard against templates having readings that
	     * are not in the contigs we are interested in. For these cases,
	     * the reading colour will not have been set so use this as a
	     * simple check
	     */

	    /* printf("read %d t %d y1 %d y2 %d colour %s type %s\n",
		   gel, t, RArray[gel].l.y1, RArray[gel].l.y2,
		   RArray[gel].colour ? RArray[gel].colour : "(nul)",
		   RArray[gel].type ? RArray[gel].type : "(nul)"); */

	    if (RArray[gel].colour) {
		RArray[gel].l.y1 = TArray[i].l.y1;
		RArray[gel].l.y2 = TArray[i].l.y2;

		ReadArray[(*num_readings)++] = RArray[gel];
		RArray[gel].type = NULL;
	    }
	} /* end for */

    } /* end for */
} /* end FindReadingYCoords */


/*****************************************************************************/
/*                                  CalcReadings                             */
/*****************************************************************************/
/* calculates the x coords of each reading
   assigns reading_num to be the reading index
   assigns type to be "reading"
   assigns colour dependent on primer type
*/
void
CalcReadings(GapIO *io,		                                       /* in */
	     int contig_num,	                                       /* in */
	     int offset,                                               /* in */
	     template_d *t_changes,                                   /* in */
	     template_c **tarr,
	     int not_single,
	     int read_pairs,
	     int span_read_pairs,
	     DPlotRec *RArray,	                                      /* out */
	     int *max_x2,	                                      /* out */
	     int *min_x1)	                                      /* out */
{
    int i;			                                  /* counter */
    GReadings reading;		                        /* reading structure */
    char *colour;
    int cnt = 0;
    template_c *t;
    int status;

    *max_x2 = 0;		        /* largest x2 value */
    *min_x1 = INT_MAX;		        /* smallest x1 value */

    for (i = io_clnbr(io, contig_num); i; i = io_rnbr(io, i)) {

	gel_read(io, i, reading);

	/* only plot relevant readings if only displaying a subset */
	t = tarr[reading.template];
	if ((not_single || read_pairs) &&
	    (NULL == head(t->gel_cont) || NULL == head(t->gel_cont)->next)) {
	    continue;

	}
	status = getStatus(t);

#ifdef DEBUG
	printf("READING readpair %d forw %d contra %d diff %d\n", 
	       t_changes[i].readpair, 
	       status & FORW_REV_READINGS, status & CONTRADICTORY,
	       status & DIFF_CONTIGS);
#endif
	if ((read_pairs || span_read_pairs) && !t_changes[reading.template].readpair &&
	    !(status & FORW_REV_READINGS) && !(status & CONTRADICTORY)) {

	    continue;
	}

	/* display read pairs that span contigs */
	if (span_read_pairs && !(status & DIFF_CONTIGS)) {
	  continue;
	}

	{
	    int tx1, tx2;
	    CalcXCoords(reading.position + offset, reading.sequence_length,
			&tx1, &tx2);
	    RArray[i].l.x1 = tx1;
	    RArray[i].l.x2 = tx2;
	}
	cnt++;
	RArray[i].num = i;

	if (NULL == (RArray[i].type = (char *)xmalloc(40))) {
	    verror(ERR_FATAL, "CalcReadings", "out of memory");
	    return;
	}
	sprintf(RArray[i].type, "{reading r_%d num_%d S}", i, contig_num);
	if (reading.sense) {
	    strcpy(RArray[i].arrow, "first");
	} else {
	    strcpy(RArray[i].arrow, "last");
	}
	switch(PRIMER_TYPE(reading)) {
	case GAP_PRIMER_UNKNOWN:
	    colour = get_default_string(GetInterp(), gap_defs,
					 "TEMPLATE.PRIMER_UNKNOWN_COLOUR");
	    RArray[i].colour = colour; /* 4 red */
	    break;

	case GAP_PRIMER_FORWARD:
	    colour = get_default_string(GetInterp(), gap_defs,
					 "TEMPLATE.PRIMER_FORWARD_COLOUR");
	    RArray[i].colour = colour;  /* 5 blue */
	    break;

	case GAP_PRIMER_REVERSE:
	    colour = get_default_string(GetInterp(), gap_defs,
					 "TEMPLATE.PRIMER_REVERSE_COLOUR");
	    RArray[i].colour = colour;  /* 8 grey */
	    break;

	case GAP_PRIMER_CUSTFOR:
	    colour = get_default_string(GetInterp(), gap_defs,
					 "TEMPLATE.PRIMER_CUSTOM_FOR_COLOUR");
	    RArray[i].colour = colour;
	    break;

	case GAP_PRIMER_CUSTREV:
	    colour = get_default_string(GetInterp(), gap_defs,
					 "TEMPLATE.PRIMER_CUSTOM_REV_COLOUR");
	    RArray[i].colour = colour;
	}

/*
	printf("max_x2 %d c2 %d min_x1 %d x1 %d\n", *max_x2, RArray[i].l.x2,
	       *min_x1, RArray[i].l.x1);
*/
	if (*max_x2 < RArray[i].l.x2) {
	    *max_x2 = RArray[i].l.x2;
	} /* end if */

	if (*min_x1 > RArray[i].l.x1) {
	    *min_x1 = RArray[i].l.x1;
	} /* end if */

    } /* end for */

} /* end CalcReadings */


/*****************************************************************************/
/*                          SetReadingPosLen                                 */
/*****************************************************************************/
/* set reading position and length depending on whether whole reading or
   displayed reading has been chosen */
void
SetReadingPosLen(int whole_reading,                                    /* in */
		 GapIO  *io,                                           /* in */
		 int reading_num,                                      /* in */
		 int *r_pos,                                          /* out */
		 int *r_len)                                          /* out */
{
    GReadings reading;

    gel_read(io, reading_num, reading);

    if (whole_reading) {

	*r_pos = reading.position - reading.start;
	*r_len = reading.length;

    } else {

	*r_pos = reading.position;
	*r_len = reading.sequence_length;

    } /* end if */

} /* end SetReadingPosLen */

/*****************************************************************************/
/*                             CalcTotalContigLen                            */
/*****************************************************************************/
/* return the total length of all the contigs in a database */
int64_t
CalcTotalContigLen(GapIO *io)                                          /* in */
{
    int i;
    int64_t total_len = 0;

    for (i = 1; i <= NumContigs(io); i++){
	total_len = total_len + ABS(io_clength(io, i));
    }

    return total_len;

} /* end CalcTotalContigLen */

/*****************************************************************************/
/*                             CalcLongContig                                */
/*****************************************************************************/
/* find the contig number of the longest contig in the database */
int64_t
CalcLongContig(GapIO *io)                                              /* in */
{
    int i;
    int contig_num = 0;
    int64_t len, longest = 0;

    for (i = 1; i <= NumContigs(io); i++){
	len = ABS(io_clength(io, i));
	if (len > longest) {
	    contig_num = i;
	    longest = len;
	}
    }

    return contig_num;

} /* end CalcLongContig */

/*****************************************************************************/
/*                             CalcTotalReadingLen                           */
/*****************************************************************************/
/* return the total length of all the readings in a database */
int64_t
CalcTotalReadingLen (GapIO *io,                                        /* in */
		     int ngels)                                        /* in */
{
    int i;
    int64_t total_len = 0;

    /* for each reading in the database */
    for (i = 1; i <= ngels; i++) {
 	total_len = total_len + ABS(io_length(io, i));
    } /* end for */

    return(total_len);

} /* end CalcTotalReadingLen */


/*****************************************************************************/
/*                             CalcContigNumLen                              */
/*****************************************************************************/
/* find the contig length & number of a contig containing reading_num */
void
CalcNumLenContig(GapIO *io,                                            /* in */
		 int reading_num,                                      /* in */
		 int *contig_num,                                     /* out */
		 int *contig_len)                                     /* out */
{
    int i, j;

    for (i = 1; i <= io->db.num_contigs; i++){
	for (j=io_clnbr(io, i); j; j=io_rnbr(io, j)) {

	    if (j == reading_num) {
		*contig_num = i;
		*contig_len = ABS(io_clength(io, i));
		break;
	    }

	} /* end for */

    } /* end for */

} /* end CalcNumLenContig */


/*
 * Removes the template and reading display (and unplots etc).
 */
static void template_shutdown(GapIO *io, obj_template_disp *t) {
    int i;
    char cmd[1024];

    for (i = 0; i < t->num_contigs; i++) {
	contig_deregister(io, t->contig[i], template_callback, (void *)t);
	delete_contig_cursor(io, t->contig[i], t->cursor[i]->id, 0);
    }
    xfree(t->contig);

    sprintf(cmd, "DeleteTemplateDisplay %s %s %d\n", t->frame, t->t_win,
	    t->id);
    if (TCL_ERROR == Tcl_Eval(t->interp, cmd))
	printf("template_shutdown:%s\n", Tcl_GetStringResult(t->interp));


    if (t->tarr)
	uninit_template_checks(io, t->tarr);

    if (t->contig_offset)
	xfree(t->contig_offset);

    free_win_list(t->win_list, t->num_wins);

    if (t->readings) {
	for (i = 0; i < t->num_readings; i++) {
	    xfree(t->readings[i].type);
	}
	xfree(t->readings);
    }

    if (t->ruler_coord) {
	for (i = 0; i < t->num_contigs; i++) {
	    xfree(t->ruler_coord[i].type);
	}
	xfree(t->ruler_coord);
    }
    xfree(t->canvas);
    xfree(t->world->visible);
    xfree(t->world->total);
    xfree(t->world);
    if (t->ruler->colour) free(t->ruler->colour);
    if (t->ruler->tick.t.colour) free(t->ruler->tick.t.colour);
    if (t->xhair.colour) free(t->xhair.colour);
    xfree(t->ruler);
    freeZoom(&t->zoom);
    xfree(t->cursor);
    xfree(t->cursor_visible);
    xfree(t);
}

static void
template_join(Tcl_Interp *interp,
	      GapIO *io,
	      obj_template_disp *t,
	      int old_contig,
	      int new_contig)
{
    int i, j;
    int length, old_found = 0;
    int remove_dup = -1;
    static char old_contig_name[DB_NAMELEN+1];
    static char new_contig_name[DB_NAMELEN+1];
    char cmd[1024];

    strcpy(old_contig_name,
	   get_read_name(io, io_clnbr(io, old_contig)));
    strcpy(new_contig_name,
	   get_read_name(io, io_clnbr(io, new_contig)));


    /* replace the old contig with the new contig */
    for (i = 0; i < t->num_contigs; i++) {
	if (ABS(t->contig[i]) == old_contig) {
	    old_found = 1;

	    /* check for duplicates and set flag */
	    for (j = 0; j < t->num_contigs; j++) {
		 if (ABS(t->contig[j]) == new_contig) {
		     remove_dup = i;
		     break;
		 }
	     }
	    /* no duplicates, don't set flag */
	    t->contig[i] = new_contig;
	    break;
	}
    }

    if (remove_dup != -1) {
	/* remove any duplicates */
	delete_contig_cursor(io, old_contig, t->cursor[remove_dup]->id, 0);
	canvas_cursor_refresh(interp, io, new_contig,
			      t->cursor[remove_dup], t->cursor[remove_dup],
			      t->canvas, t->win_list, t->num_wins,
			      t->id, t->contig_offset[new_contig].offset,
			      &t->cursor_visible[remove_dup], t->world, 1);

	length = t->num_contigs - remove_dup - 1;
	memmove(&(t->contig)[remove_dup],
		&(t->contig)[remove_dup+1],
		length * sizeof(int));
	memmove(&(t->cursor)[remove_dup],
		&(t->cursor)[remove_dup+1],
		length * sizeof(int));
	memmove(&(t->cursor_visible)[remove_dup],
		&(t->cursor_visible)[remove_dup+1],
		length * sizeof(int));
	t->num_contigs--;
    }

    if (old_found) {
	sprintf(cmd,"%s.menubar.[menu_path {View.Quality Plot}] delete \"contig %s\" ",
		t->frame, old_contig_name);
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    printf("template_join:%s\n", Tcl_GetStringResult(interp));
	sprintf(cmd,"%s.menubar.[menu_path {View.Restriction Enzyme Plot}] delete \"contig %s\" ",
		t->frame, old_contig_name);
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    printf("template_join:%s\n", Tcl_GetStringResult(interp));
    }
}


static void
template_renumber(GapIO *io,
		  obj_template_disp *t,
		  int old_contig,
		  int new_contig)
{
    int i;
    for (i = 0; i < t->num_contigs; i++) {
	if (ABS(t->contig[i]) == old_contig) {
	    t->contig[i] = t->contig[i] > 0 ? new_contig :
		-new_contig;
	    break;
	}
    }
}

/*
 * Cursor updates
 */
static void template_update_cursor(GapIO *io,
				   obj_template_disp *t,
				   int contig,
				   cursor_t *cursor,
				   int show_cursor)
{
    int index, i;

    index = 0;
    for (i = 0; i < t->num_contigs; i++) {
	if (t->contig[i] == contig) {
	    index = i;
	    break;
	}
    }

    canvas_cursor_refresh(t->interp, io, t->contig[index],
			  cursor, t->cursor[index],
			  t->canvas, t->win_list, t->num_wins,
			  t->id, t->contig_offset[t->contig[index]].offset,
			  &t->cursor_visible[index], t->world, show_cursor);
}

void template_update_cursors(GapIO *io,
			     obj_template_disp *t,
			     int show_cursor)
{
    int i;
    for (i = 0; i < t->num_contigs; i++) {
	template_update_cursor(io, t, t->contig[i], t->cursor[i], show_cursor);
    }
}

/*
 * Callback for the template plot
 */
static void template_callback(GapIO *io, int contig, void *fdata,
			      reg_data *jdata) {
    obj_template_disp *t = (obj_template_disp *)fdata;
    char cmd[1024];
    int i;
#ifdef DEBUG
    printf("template_callback job %d contig %d\n", jdata->job, contig);
#endif
    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Template display");
	    return;
	}
    case REG_COMPLEMENT:
	{
	    Tcl_DString c_list;
	    /* relabel with new contig */

#ifdef HACK
	    int right;
	    right = io_crnbr(io, ABS(contig));
 	    printf("RIGHT %d %d\n", right, contig);

	    sprintf(cmd,"%s.menubar.view.opts.quality entryconfigure \"contig %s\" -label \"contig %s\"", t->frame, get_read_name(io, right),
		    get_contig_name(io, ABS(contig)));
	    if (TCL_ERROR == Tcl_Eval(t->interp, cmd))
		printf("update_template_view_menu:%s\n", Tcl_GetStringResult(t->interp));

	    sprintf(cmd,"%s.menubar.view.opts.renz entryconfigure \"contig %s\" -label \"contig %s\"", t->frame, get_read_name(io, right),
		    get_contig_name(io, ABS(contig)));
	    if (TCL_ERROR == Tcl_Eval(t->interp, cmd))
		printf("update_template_view_menu:%s\n", Tcl_GetStringResult(t->interp));
#endif
	    Tcl_DStringInit(&c_list);

	    Tcl_DStringAppendElement(&c_list, "UpdateTemplateContigList");
	    sprintf(cmd, "%d", *handle_io(io));
	    Tcl_DStringAppendElement(&c_list, cmd);
	    Tcl_DStringAppendElement(&c_list, t->frame);
	    Tcl_DStringStartSublist(&c_list);
	    for (i = 0; i < t->num_contigs; i++) {
		Tcl_DStringAppendElement(&c_list, get_contig_name(io,
							  ABS(t->contig[i])));
	    }
	    Tcl_DStringEndSublist(&c_list);

	    if (TCL_ERROR == Tcl_Eval(t->interp, Tcl_DStringValue(&c_list)))
		printf("template complement: %s \n", Tcl_GetStringResult(t->interp));
	    Tcl_DStringFree(&c_list);

	    if (!t->do_update)
		update_template_display(t->interp, io, t, 1);
	    else
		t->do_update |= REG_LENGTH;

	    return;

	}

    case REG_JOIN_TO:
#ifdef DEBUG
	printf("template JOIN_TO contig %d join %d\n",
	       contig, jdata->join.contig);
#endif
	template_join(t->interp, io, t, contig, jdata->join.contig);
	return;
    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("template LENGTH %d\n", contig);
#endif
	    if (!t->do_update)
		update_template_display(t->interp, io, t, 1);
	    else
		t->do_update |= REG_LENGTH;

	    return;
	}
    case REG_BUFFER_START:
	{
#ifdef DEBUG
	    printf("template REG_BUFFER_START\n");
#endif
	    t->buffer_count++;
	    t->do_update = 1; /* initialise do_update */
	    return;
	}

    case REG_BUFFER_END:
	{
#ifdef DEBUG
	    printf("template REG_BUFFER_END count %d \n", t->buffer_count);
#endif
	    t->buffer_count--;

	    if (t->buffer_count <= 0) {

		if (t->do_update & REG_LENGTH) {
		    update_template_display(t->interp, io, t, 1);
		} else if (t->do_update & TASK_CANVAS_REDRAW) {
		    int *recalc = (int *)jdata->generic.data;
		    update_template_display(t->interp, io, t, *recalc);
		} else if (t->do_update & REG_ANNO) {
		    display_reading_tags(t->interp, io, t);
		    scaleSingleCanvas(t->interp, t->world, t->canvas,
				      t->window, 'b', "tag");
		    if (strcmp(t->window, t->ruler->window) != 0) {
			scaleSingleCanvas(t->interp, t->world, t->canvas,
					  t->ruler->window, 'x', "tag");
		    }
		}
		t->buffer_count = 0;
		t->do_update = 0;
	    }
	    return;
	}
    case REG_ANNO:
	{
#ifdef DEBUG
	    printf("template REG_ANNO \n");
#endif
	    if (!t->do_update) {
		display_reading_tags(t->interp, io, t);
		scaleSingleCanvas(t->interp, t->world, t->canvas,
				  t->window, 'b', "tag");
		if (strcmp(t->window, t->ruler->window) != 0) {
		    scaleSingleCanvas(t->interp, t->world, t->canvas,
					  t->ruler->window, 'x', "tag");
		}
	    } else {
		t->do_update |= REG_ANNO;
	    }
	    return;
	}
    case REG_QUIT:
    case REG_DELETE:
	{
	    /* if deleting a contig shutdown the display */
#ifdef DEBUG
	    printf("template QUIT\n");
#endif
	    template_shutdown(io, t);
	    return;
	}
    case REG_GET_OPS:
	jdata->get_ops.ops = "Remove\0";
	return;

    case REG_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0:
	    template_shutdown(io, t);
	}
	return;

    case REG_PARAMS:
	return;

    case REG_CURSOR_NOTIFY:
	template_update_cursor(io, t, contig, jdata->cursor_notify.cursor, 1);
	return;

    case REG_NUMBER_CHANGE:
	/*
	 * only called in io_delete_contig in IO2.c but never executed
	 * because we always quit displays beforehand
	 */
#ifdef DEBUG
	printf("template NUMBER_CHANGE\n");
#endif
	template_renumber(io, t, contig, jdata->number.number);
	return;

    case REG_HIGHLIGHT_READ: {
	char buf[1024];
	int width;
#ifdef DEBUG
	printf("template HIGHLIGHT_READ\n");
#endif
	if (jdata->highlight.seq > 0) { /* can't do the consensus */
	    if (jdata->highlight.val == 1) {
		width = t->line_bold;
	    } else {
		width = t->line_width;
	    }
	    sprintf(buf, "%s itemconfig r_%d -width %d",
			t->window, jdata->highlight.seq, width);

	    if (Tcl_Eval(t->interp, buf) == TCL_ERROR) {
		puts(Tcl_GetStringResult(t->interp));
	    }
	}
	return;
    }
    case REG_GENERIC:
	/* printf("template GENERIC %d\n", jdata->generic.task); */
	switch (jdata->generic.task) {
	case TASK_WINDOW_ADD:
	    {
		win *winfo = (win *)jdata->generic.data;

#ifdef DEBUG
		printf("TASK_WINDOW_ADD %d \n", t->num_wins);
		for (i = 0; i < t->num_wins; i++) {
		    printf("before win %s \n", t->win_list[i]->window);
		}
#endif
		addWindow(t->win_list, &t->num_wins, winfo->window,
			  winfo->scroll, winfo->id);

		/*
		 * if this is the only window, need to set up the canvas struct
		 */
		if (t->num_wins == 1) {
		    strcpy(t->window, t->win_list[0]->window);
		    initCanvas(t->interp, t->canvas, t->win_list[0]->window);
		    SetCanvasCoords(t->interp,
				    t->world->visible->x1,
				    t->world->visible->y1,
				    t->world->visible->x2,
				    t->world->visible->y2, t->canvas);
		}
#ifdef DEBUG
		printf("TASK_WINDOW_ADD %d \n", t->num_wins);
		for (i = 0; i < t->num_wins; i++) {
		    printf("after win %s \n", t->win_list[i]->window);
		}
#endif

		break;
	    }
	case TASK_WINDOW_DELETE:
	    {
		char *window = (char *)jdata->generic.data;
#ifdef DEBUG
		int i;
		printf("TASK_WINDOW_DELETE %d \n", t->num_wins);
		for (i = 0; i < t->num_wins; i++) {
		    printf("before win %s \n", t->win_list[i]->window);
		}
#endif
		deleteWindow(t->win_list, &t->num_wins, window);
		/*
		 * if delete the window in t->window, must set it to be
		 * another window in the list
		 */
		if (strcmp(t->window, window) == 0 && t->num_wins > 0) {
		    strcpy(t->window, t->win_list[0]->window);
		    initCanvas(t->interp, t->canvas, t->win_list[0]->window);
		    SetCanvasCoords(t->interp,
				    t->world->visible->x1,
				    t->world->visible->y1,
				    t->world->visible->x2,
				    t->world->visible->y2, t->canvas);
		}

#ifdef DEBUG
		for (i = 0; i < t->num_wins; i++) {
		    printf("after win %s \n", t->win_list[i]->window);
		}
#endif
		break;
	    }
	case TASK_CANVAS_SCROLLX:
	    {
		char *scroll = (char *)jdata->generic.data;

		canvasScrollX(t->interp, t->window, t->win_list, t->num_wins,
			      t->world->visible, t->canvas, scroll);
		break;
	    }
	case TASK_CANVAS_SCROLLY:
	    {
		char *scroll = (char *)jdata->generic.data;

		canvasScrollY(t->interp, t->window, t->win_list, t->num_wins,
			      t->world->visible, t->canvas, scroll);

		break;
	    }
	case TASK_CANVAS_RESIZE:
	    {
		char scroll_args[20];
		/* resize template display window */

		resizeCanvas(t->interp, t->window, t->win_list, t->num_wins,
			     t->world->visible, t->world->total, t->canvas);
		sprintf(scroll_args, "scroll 0 units");
		canvasScrollX(t->interp, t->window, t->win_list, t->num_wins,
			      t->world->visible, t->canvas, scroll_args);

		break;
	    }
	case TASK_CANVAS_WORLD:
	    {
		task_world_t *tw = (task_world_t *)jdata->generic.data;
		int cx = tw->canvasx;
		double wx, wy;

		CanvasToWorld(t->canvas, cx, 0, &wx, &wy);
		tw->basex = wx - t->contig_offset[tw->cnum].offset;
		break;
	    }
	case TASK_CANVAS_ZOOMBACK:
	    {
		reg_generic gen;

		/*
		 * to allow zoom not to be reset even when add/delete
		 * templates/readings
		 */
		/* listZoom(t->zoom); */
		if (lengthZoom(t->zoom) <= 2) {
		    freeZoom(&t->zoom);
		    pushZoom(&t->zoom, t->world->total);
		}
		/* listZoom(t->zoom); */

		canvasZoomback(t->interp, t->canvas, t->window, t->world,
			   t->win_list, t->num_wins, &t->zoom);


		/* redraw ruler ticks */
		Tcl_VarEval(t->interp, t->ruler->window, " delete tick", NULL);
		gen.job = REG_GENERIC;
		gen.task = TASK_DISPLAY_TICKS;
		result_notify(io, t->id, (reg_data *)&gen, 0);
		break;
	    }
	case TASK_CANVAS_ZOOM:
	    {
		s_zoom *szoom = (s_zoom *)jdata->generic.data;
		reg_generic gen;

		canvasZoom(t->interp, t->canvas, t->window, t->world,
			   t->win_list, t->num_wins, &t->zoom, szoom->zoom,
			   szoom->scroll);

		/* redraw ruler ticks */
		Tcl_VarEval(t->interp, t->ruler->window, " delete tick", NULL);
		gen.job = REG_GENERIC;
		gen.task = TASK_DISPLAY_TICKS;
		result_notify(io, t->id, (reg_data *)&gen, 0);

		break;
	    }
	case TASK_CANVAS_CURSOR_X:
	    {
		char *label;
		int *cx = (int *)jdata->generic.data;
		double local_pos;
		double wx, wy;

		CanvasToWorld(t->canvas, *cx, 0, &wx, &wy);

		label = get_default_string(t->interp, gap_defs,
					   "TEMPLATE.CURSOR1");

		canvasCursorX(t->interp, t->canvas, t->frame, label,
			      t->xhair.colour, t->xhair.width, *cx, wx,
			      t->win_list, t->num_wins);

		/* fill in local position of cursor in label box */
		local_pos = TemplateLocalCursor(t->id, t->contig_offset,
						t->contig,
						t->num_contigs, wx);

		label = get_default_string(t->interp, gap_defs,
					   "TEMPLATE.CURSOR2");
		sprintf(cmd, "%s%s configure -text %d\n", t->frame, label,
			(int)local_pos);
		Tcl_Eval(t->interp, cmd);

		break;
	    }
	case TASK_CANVAS_CURSOR_DELETE:
	    {
		int i;
		for (i = 0; i < t->num_wins; i++) {
		    Tcl_VarEval(t->interp, t->win_list[i]->window,
				" delete cursor_x", NULL);
		}
		break;
	    }
	case TASK_CANVAS_REDRAW:
	    {
		int *recalc = (int *)jdata->generic.data;
		/* despite many requests, only redraw once */

#ifdef DEBUG
		printf("TASK_CANVAS_REDRAW update %d recalc %d\n",
		       t->do_update, *recalc);
#endif
		if (!t->do_update) {
		    update_template_display(t->interp, io, t, *recalc);
		} else {
		    t->do_update |= TASK_CANVAS_REDRAW;
		}
		break;
	    }
	case TASK_DISPLAY_RULER:
	    {
		ruler_s *ruler = (ruler_s *)jdata->generic.data;

#ifdef DEBUG
		printf("TASK_DISPLAY_RULER\n");
#endif
		if (t->ruler_coord) {
		    for (i = 0; i < t->num_contigs; i++) {
			xfree(t->ruler_coord[i].type);
		    }
		    xfree(t->ruler_coord);
		    t->ruler_coord = NULL;
		}
		display_ruler(t->interp, io, t->canvas, t->contig_offset,
			      t->contig, t->num_contigs, t->configs[RULER],
			      t->configs[TICKS],
			      t->ruler, t->frame, &t->ruler_coord);

#ifdef FIXME
		if (NULL == (canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
		    return;
		/*
		 * turn ruler off, zoom, then turn ruler on, ruler is
		 * displayed in the wrong place.
		 */

		canvas->width = t->canvas->width;
		canvas->height = t->canvas->height;
		canvas->ax = t->canvas->ax;
		canvas->ay = t->canvas->ay;
		canvas->bx = t->canvas->bx;
		canvas->by = t->canvas->by;
		canvas->x = t->canvas->x;
		canvas->y = t->canvas->y;
#endif

		scaleSingleCanvas(t->interp, t->world, t->canvas,
				  ruler->window, 'x', "all");

		template_update_cursors(io, t, 0);
		break;
	    }
	case TASK_DISPLAY_TICKS:
	    {
		/* int *ticks = (int *)jdata->generic.data; */
		char cmd[100];

		if (!t->configs[TICKS] || !t->configs[RULER])
		    return;

		for (i = 0; i < t->num_contigs; i++) {
		    /*
		    display_ruler_ticks(t->interp, t->canvas,
					t->contig_offset[t->contig[i]].offset,
					t->ruler_coord[i].l.y1,
					t->ruler, t->ruler_coord[i].l.x1,
					t->ruler_coord[i].l.x2);
					*/
		    display_ruler_ticks(t->interp, t->canvas,
					t->contig_offset[t->contig[i]].offset,
					t->ruler_coord[i].l.y1,
					t->ruler, 1, t->ruler_coord[i].l.x2 -
					t->ruler_coord[i].l.x1 + 1);
		}
		scaleSingleCanvas(t->interp, t->world, t->canvas,
				  t->ruler->window, 'x', "tick");
		sprintf(cmd, "RulerWindowSize %d %s %s ", 1, t->frame,
			t->ruler->window);
		Tcl_Eval(t->interp, cmd);
		template_update_cursors(io, t, 0);

		break;
	    }
	}
    }
}

void
template_config(Tcl_Interp *interp,
		char *frame,
		int *config_array)
{
    char config[1024];
    int i;

    /*
     * initialise config_array
     */
    for (i = 0; i < NUM_CONFIGS; i++) {
	config_array[i] = 0;
    }

    sprintf(config, "config%s.template", frame);
    config_array[TEMPLATES] = atoi(Tcl_GetVar(interp, config, TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[TEMPLATES],
		TCL_LINK_INT);

    sprintf(config, "config%s.reading", frame);
    config_array[READINGS] = atoi(Tcl_GetVar(interp, config, TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[READINGS],
		TCL_LINK_INT);

    sprintf(config, "config%s.multi_template", frame);
    config_array[MULTI_TEMPLATES] = atoi(Tcl_GetVar(interp, config,
						    TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[MULTI_TEMPLATES],
		TCL_LINK_INT);

    sprintf(config, "config%s.read_pairs", frame);
    config_array[READ_PAIRS] = atoi(Tcl_GetVar(interp,config,TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[READ_PAIRS],
		TCL_LINK_INT);

    sprintf(config, "config%s.ruler", frame);
    config_array[RULER] = atoi(Tcl_GetVar(interp, config, TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[RULER],
		TCL_LINK_INT);

    sprintf(config, "config%s.ticks", frame);
    config_array[TICKS] = atoi(Tcl_GetVar(interp, config, TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[TICKS],
		TCL_LINK_INT);

    sprintf(config, "config%s.span_read_pairs", frame);
    config_array[SPAN_READ_PAIRS] = atoi(Tcl_GetVar(interp,config,TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[SPAN_READ_PAIRS],
		TCL_LINK_INT);

    sprintf(config, "config%s.consist_read_pairs", frame);
    config_array[CONSIST_READ_PAIRS] = atoi(Tcl_GetVar(interp,config,TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[CONSIST_READ_PAIRS],
		TCL_LINK_INT);


    sprintf(config, "config%s.calc_contig_pos", frame);
    config_array[CALC_CONTIG_POS] = atoi(Tcl_GetVar(interp,config,TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[CALC_CONTIG_POS],
		TCL_LINK_INT);

}

/*
 * Registers the template display.
 */
int template_reg(Tcl_Interp *interp,
		 GapIO *io,
		 int *contig_array,
		 int num_contigs,
		 char *frame,
		 char *t_win,
		 char *r_win,
		 ruler_s *ruler,
		 cursor_s xhair,
		 int line_width,
		 int line_bold)
{
    obj_template_disp *t;
    int id;
    int i;

    if (NULL == (t = (obj_template_disp *)xmalloc(sizeof(obj_template_disp))))
	return 0;

    if (NULL == (t->contig_offset = (c_offset *)xmalloc((1+io->db.num_contigs)
							* sizeof(c_offset))))
	return -1;
    if (NULL == (t->cursor =
		 (cursor_t **)xmalloc(num_contigs * sizeof(*(t->cursor)))))
	return -1;
    if (NULL == (t->cursor_visible =
		 (int *)xmalloc(num_contigs * sizeof(*(t->cursor_visible)))))
	return -1;

    id = register_id();

    strcpy(t->frame, frame);
    strcpy(t->window, t_win);
    strcpy(t->t_win, t_win);
    t->contig = contig_array;
    t->num_contigs = num_contigs;
    t->id = id;
    t->line_width = line_width;
    t->line_bold = line_bold;
    t->xhair = xhair;
    t->do_update = 0;
    t->buffer_count = 0;
    t->interp = interp;
    t->tarr = NULL;
    t->readings = NULL;
    t->ruler_coord = NULL;

    t->ruler = ruler;
    t->ruler->start = -1;
    t->ruler->end = -1;
    sprintf(t->ruler->window, "%s", r_win);

    for (i = 0; i < num_contigs; i++) {
	t->cursor_visible[i] = 0;
	t->cursor[i] = create_contig_cursor(io, contig_array[i], 0, id);
    }

    /* create list of windows in the template display */
    if (NULL == (t->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;

    t->num_wins = 0;
    addWindow(t->win_list, &t->num_wins, t->window, 'b', t->id);
    addWindow(t->win_list, &t->num_wins, t->ruler->window, 'x', t->id);

    if (NULL == (t->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (t->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (t->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (t->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    initCanvas(interp, t->canvas, t->window);
    createZoom(&t->zoom);
    template_config(interp, t->frame, t->configs);

    /*
     * do the first plot and do this BEFORE registrations because otherwise
     * any contig editors already up will try to tell the as yet uninitialised
     * template display of a cursor_notify request
     */
    update_template_display(interp, io, t, 1);

    /*
     * to allow zoom not to be reset even when add/delete templates/readings
     */
    /* add first zoom */
    pushZoom(&t->zoom, t->world->visible);

    for (i = 0; i < num_contigs; i++) {

	contig_register(io, contig_array[i], template_callback, (void *)t, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_CURSOR_NOTIFY | REG_NUMBER_CHANGE |
			REG_BUFFER_START | REG_BUFFER_END |REG_ANNO |
			REG_JOIN_TO | REG_GENERIC |
			REG_HIGHLIGHT_READ, REG_TYPE_TEMPLATE);
	template_update_cursor(io, t, contig_array[i], t->cursor[i], 1);
    }

    return id;
}

void
plot_lines(Tcl_Interp *interp,
	   PlotRec *array,
	   int num,
	   char *win_name,
	   int line_width)
{
    char *cmd;
    int i;
    int len, alen;

    alen = 1024;
    if (NULL == (cmd = (char *)xmalloc(alen * sizeof(char)))) {
	return;
    }

    for (i = 0; i < num; i++) {

	len = flen("%s create line %d %d %d %d "
		   "-fill {%s} -tags %s -width %d -arrow %s\n",
		   win_name,
		   array[i].l.x1,
		   array[i].l.y1,
		   array[i].l.x2,
		   array[i].l.y2,
		   array[i].colour,
		   array[i].type,
		   line_width,
		   array[i].arrow);

	if (len > alen) {
	    alen = len;
	    if (NULL == (cmd = (char *)xrealloc(cmd, alen * sizeof(char))))
		return;
	}

	sprintf(cmd, "%s create line %d %d %d %d "
		"-fill {%s} -tags %s -width %d -arrow %s\n",
		win_name,
		array[i].l.x1,
		array[i].l.y1,
		array[i].l.x2,
		array[i].l.y2,
		array[i].colour,
		array[i].type,
		line_width,
		array[i].arrow);

	Tcl_Eval(interp, cmd);
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
    }
    xfree(cmd);
}

void
plot_dlines(Tcl_Interp *interp,
	    DPlotRec *array,
	    int num,
	    char *win_name,
	    int line_width)
{
    char *cmd;
    int i;
    int len, alen;

    alen = 1024;
    if (NULL == (cmd = (char *)xmalloc(alen * sizeof(char)))) {
	return;
    }

    for (i = 0; i < num; i++) {

	len = flen("%s create line %f %f %f %f "
		   "-fill {%s} -tags %s -width %d -arrow %s\n",
		   win_name,
		   array[i].l.x1,
		   array[i].l.y1,
		   array[i].l.x2,
		   array[i].l.y2,
		   array[i].colour,
		   array[i].type,
		   line_width,
		   array[i].arrow);

	if (len > alen) {
	    alen = len;
	    if (NULL == (cmd = (char *)xrealloc(cmd, alen * sizeof(char))))
		return;
	}

	sprintf(cmd, "%s create line %f %f %f %f "
		"-fill {%s} -tags %s -width %d -arrow %s\n",
		win_name,
		array[i].l.x1,
		array[i].l.y1,
		array[i].l.x2,
		array[i].l.y2,
		array[i].colour,
		array[i].type,
		line_width,
		array[i].arrow);

	Tcl_Eval(interp, cmd);
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
    }
    xfree(cmd);
}

int
display_templates(Tcl_Interp *interp,
		  GapIO *io,
		  obj_template_disp *t,
		  template_d *t_changes)
{
    int num_templates;		/* number of templates in the contig */
    int t_x1 = INT_MAX;	        /* min template x coord */
    int t_x2 = 0;		/* max template x coord */
    int r_x1 = INT_MAX;	        /* min readings x coord */
    int r_x2 = 0;		/* max reading x coord */
    int t_y2 = 0;
    int r_y2 = 0;
    int i;
    DPlotRec *TArray, **TArray1, **TArray2;
    DPlotRec *RArray;
    DPlotRec *ReadArray;
    int num_readings, check;
    int depth;                         /* depth of templates in y direction */
    char cmd[1024];
    int ntemplates = Ntemplates(io);

    /* first check if template window is displayed and return if it isn't */
    if (strcmp(t->window, t->t_win) != 0)
	return 0;

    /* check that there are templates or readings to be displayed */
    if (!t->configs[TEMPLATES] && !t->configs[READINGS]) {

	sprintf(cmd, "DeleteTemplatePlot %d %d %s %s",
		*handle_io(io), t->id, t->frame, t->t_win);
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    printf("display_templates: %s \n", Tcl_GetStringResult(interp));

	deleteWindow(t->win_list, &t->num_wins, t->window);
	if (t->num_wins > 0) {
	    strcpy(t->window, t->win_list[0]->window);
	} else {
	    *t->window = '\0';
	}

	t->world->total->x1 = t_x1;
	t->world->total->x2 = t_x2;
	t->world->total->y1 = t_y2;
	t->world->total->y2 = t_y2;

	return 0;
    }

    num_readings = NumReadings(io);

    if ((RArray = (DPlotRec *)xcalloc((num_readings+1), sizeof(DPlotRec)))==NULL){
	return -1;
    }
    if ((TArray = (DPlotRec *)xcalloc((ntemplates+1),
				      sizeof(DPlotRec))) ==NULL){
	return -1;
    }
    if ((TArray1 = (DPlotRec **)xcalloc((ntemplates+1),
				     sizeof(DPlotRec*))) == NULL){
	return -1;
    }

    if ((TArray2 = (DPlotRec **)xcalloc((Ntemplates(io)+1),
				     sizeof(DPlotRec*))) ==NULL){
	return -1;
    }
    TArray1[0] = &TArray[0];
    TArray2[0] = &TArray[0];

    /*
     * need to check that colour/type has been set later on so initialise here
     */
    for (i = 0; i < num_readings+1; i++) {
	RArray[i].colour = NULL;
	RArray[i].type = NULL;
    }


    /* display templates */
    if (t->configs[TEMPLATES]){

	check = CalcTemplates(io, t->contig_offset, t->contig, t->num_contigs,
			      t_changes,
			      t->tarr, TArray, TArray1, TArray2,
			      &num_templates, &t_x1, &t_x2, &t_y2,
			      t->configs[MULTI_TEMPLATES],
			      t->configs[READ_PAIRS],
			      t->configs[SPAN_READ_PAIRS],
			      t->configs[CONSIST_READ_PAIRS],
			      t->canvas->height);

	if (-1 == check) {
	    printf("ERROR: in calctemplates \n");
	    return -1;
	}
	plot_dlines(interp, TArray, num_templates, t->window, t->line_width);
    } /* end if */

    /* display readings */
    if (t->configs[READINGS]){

	if ((ReadArray = (DPlotRec *)xcalloc((num_readings+1),
					     sizeof(DPlotRec)))==NULL) {
	    return -1;
	}

	num_readings = 0;
	for (i = 0; i < t->num_contigs; i++) {

	    /* find readings on current contig */
	    CalcReadings(io, t->contig[i],
			 t->contig_offset[t->contig[i]].offset, t_changes,
			 t->tarr, t->configs[MULTI_TEMPLATES],
			 t->configs[READ_PAIRS], t->configs[SPAN_READ_PAIRS],
			 RArray, &r_x2, &r_x1);

	}


	/* if have template defined then find Y coords using template */
	if (t->configs[TEMPLATES]) {

	    FindReadingYCoords(io, t->tarr, TArray, RArray, ReadArray,
			       &num_readings, num_templates);

	} else {

	    /* y spacing alg #2: columns of readings of different heights */
	    CalcReadingYDepth(io, t->contig, t->num_contigs, RArray, &depth);
	    /* r_y2 = depth; */
	    r_y2 = t->canvas->height;
	    CalcReadingYCoords(io, t->contig, t->num_contigs, RArray,
			       ReadArray, depth, t->canvas->height,
			       &num_readings);

	}
	if (t->readings) {
	    for (i = 0; i < t->num_readings; i++) {
		xfree(t->readings[i].type);
	    }
	    xfree(t->readings);
	}
	t->readings = ReadArray;
	t->num_readings = num_readings;
	plot_dlines(interp, ReadArray, num_readings, t->window, t->line_width);

	/* highlight reading list */
	sprintf(cmd, "SelectReadingList %d ", *handle_io(io));
	Tcl_Eval(interp, cmd);

    } /* end if */

    if (t_x1 < r_x1)
	t->world->total->x1 = (double)t_x1;
    else
	t->world->total->x1 = (double)r_x1;

    if (t_x2 > r_x2)
	t->world->total->x2 = (double)t_x2;
    else
	t->world->total->x2 = (double)r_x2;

    t->world->total->y1 = (double)1;
    if (t_y2 > r_y2) {
	t->world->total->y2 = (double)t_y2;
    } else {
	t->world->total->y2 = (double)r_y2;
    }

#ifdef DEBUG
    printf("t_x1 %d r_x1 %d t_x2 %d r_x2 %d t_y2 %d r_y2 %d\n",
	   t_x1, r_x1, t_x2, r_x2, t_y2, r_y2);
#endif
    if (TArray[0].type) {
	for (i = 0; i < num_templates+1; i++) {
	    xfree(TArray[i].type);
	}
    }

    for (i = 0; i < num_readings+1; i++) {
	if (RArray[i].type)
	    xfree(RArray[i].type);
    }

    xfree(TArray);
    xfree(TArray1);
    xfree(TArray2);
    xfree(RArray);
    return 0;
}

void
DrawReadingTags(Tcl_Interp *interp,                                    /* in */
		int x1,
		int y,
		int x2,
		int tag_num,
		GAnnotations *annotation,
		char *win_name,
		int line_width,
		int contig_num)
{
    char type[30];
    char *colour = tag_db[0].bg_colour;
    char cmd[1024], str[5];
    int k;                                                        /* counter */

    sprintf(type, "{tag %s t_%d num_%d}",
	    type2str(annotation->type, str), tag_num, contig_num);

    /* find tag colour in tag_db */
    for (k = 0; k < tag_db_count; k++){

	if (annotation->type == str2type(tag_db[k].id)) {
	    colour = tag_db[k].bg_colour;
	    break;
	}  /* end if */

    } /* end for */

    sprintf(cmd, "%s create rectangle %d %d %d %d "
	    "-fill {%s} -tag %s -width %d -outline {%s}\n",
	    win_name, x1, y, x2, y, colour, type, line_width, colour);

    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	printf("%s\n", Tcl_GetStringResult(interp));
    }
    /* printf("cmd %s \n", cmd); */
}

void
display_consensus_tags(Tcl_Interp *interp,
		       GapIO *io,
		       int num_tags,
		       char **tags,
		       int c_num,
		       int contig_offset,
		       char *win_ruler,
		       int ruler_offset,
		       int tag_line_width)
{
    GAnnotations *annotation;                       /* annotation structure */
    int x1, x2;
    int tag_num, tag_pos;

    /* consensus tags */
    annotation = get_tag_num(io, -c_num, num_tags, tags, &tag_num);

    while (annotation && annotation != (GAnnotations *)-1){

	tag_pos = annotation->position;
	CalcXCoords(tag_pos, annotation->length, &x1, &x2);

	x1 += contig_offset;
	x2 += contig_offset;

	DrawReadingTags(interp, x1, ruler_offset, x2,
			tag_num, annotation, win_ruler,
			tag_line_width, c_num);

	annotation = get_tag_num(io, 0, num_tags, tags, &tag_num);

    } /* end while */
}

/*
 * draw tags on readings and ruler if the readings are currently displayed
 */
int
display_reading_tags(Tcl_Interp *interp,
		     GapIO *io,
		     obj_template_disp *t)
{
    GAnnotations *annotation;                       /* annotation structure */
    int i, c, r;
    GReadings reading;		                        /* reading structure */
    char **tags = NULL;
    int num_tags;
    int tag_num, tag_pos;
    int r_len, r_pos;
    int x1, x2;
    int *r_offsets = NULL;
    int c_index;

    if (!t->configs[READINGS] && !t->configs[RULER])
	return 0;

    /* get template display tag list */
    /* HACK - put in registration structure ? */
    Tcl_VarEval(interp, "GetDefaultTags ", "TEMPLATE.TAGS ", t->frame, NULL);

    if (SetActiveTags2(Tcl_GetStringResult(interp), &num_tags, &tags) == -1) {
	return -1;
    }

    if (num_tags == 0) {
	if (tags)
	    Tcl_Free((char *)tags);
	return 0;
    }

    /*
     * display consensus tags on the ruler
     */
    if (t->configs[RULER]) {
	Tcl_VarEval(interp, t->ruler->window, " delete tag", NULL);
	for (i = 0; i < t->num_contigs; i++) {
	    display_consensus_tags(interp, io, num_tags, tags, t->contig[i],
				   t->contig_offset[t->contig[i]].offset,
				   t->ruler->window,
				   t->ruler->tag.offset+
				   t->ruler_coord[i].l.y1,
				   t->ruler->tag.width);
	}

    }

    /*
     * display reading tags on readings and on the ruler (if displayed)
     */
    if (!t->configs[READINGS]) {
	Tcl_Free((char *)tags);
	return 0;
    }

    if (t->configs[RULER]) {
	/*
	 * ruler offsets in template ruler display. Indexed on contig number
	 */
	if (NULL == (r_offsets = (int *)xmalloc((NumContigs(io)+1) *
						sizeof(int)))) {
	    if (tags)
		Tcl_Free((char *)tags);
	    return -1;
	}

	for (i = 0; i < t->num_contigs; i++) {
	    r_offsets[t->contig[i]] = t->ruler_coord[i].l.y1;
	}
    }

    Tcl_VarEval(interp, t->window, " delete tag", NULL);
    for (i = 0; i < t->num_readings; i++) {
	r = t->readings[i].num;
	c = rnumtocnum(io, r);

	c_index = getContigIndex(t->contig, t->num_contigs, c);

	gel_read(io, r, reading);
	/* get first tag */

	annotation = get_tag_num(io, r, num_tags, tags, &tag_num);

	while (annotation && annotation != (GAnnotations *)-1){

	    /* if reading has been complemented, find new pos of tag */
	    if (reading.sense) {

		tag_pos = (reading.position - reading.start) +
		    (reading.length - annotation->position -
		     annotation->length + 1);

	    } else {

		tag_pos = annotation->position - (reading.start -
						  reading.position +1);

	    } /* end if */

	    r_pos = reading.position;
	    r_len = reading.sequence_length;
	    CalcXCoords(tag_pos, annotation->length, &x1, &x2);

	    /* clip tag at cutoff data */
	    x1 = MAX(x1, r_pos);
	    x2 = MIN(x2, r_pos + r_len - 1);
	    if (x2 >= x1) {
		x1 += t->contig_offset[c].offset;
		x2 += t->contig_offset[c].offset;

		DrawReadingTags(interp, x1, t->readings[i].l.y1, x2, tag_num,
				annotation, t->window, t->ruler->tag.width, c);

		if (t->configs[RULER]) {
		    DrawReadingTags(interp, x1,
				    r_offsets[c] - t->ruler->tag.offset, x2,
				    tag_num, annotation, t->ruler->window,
				    t->ruler->tag.width, c);
		}
	    }

	    /* get next annotation */
	    annotation = get_tag_num(io, 0, num_tags, tags, &tag_num);
	} /* end while */
    } /* end for */

    if (t->configs[RULER])
	xfree(r_offsets);
    Tcl_Free((char *)tags);

    return 0;
}

/*
 * draws the template display from scratch, ie removes all previous zooming
 * information
 */
int
update_template_display(Tcl_Interp *interp,
			GapIO *io,
			obj_template_disp *t,
			int recalculate)
{
    int ruler_length;
    template_d *t_changes = NULL;

    Tcl_VarEval(interp, t->window, " delete template", NULL);
    Tcl_VarEval(interp, t->window, " delete reading", NULL);
    Tcl_VarEval(interp, t->window, " delete tag", NULL);

/*
    printf("total1 %f %f %f %f\n", t->world->total->x1,
	   t->world->total->y1, t->world->total->x2, t->world->total->x2);
*/
    if (recalculate) {
	if (t->tarr)
	    uninit_template_checks(io, t->tarr);
	CalcContigOffsets(io, t->contig_offset, t->contig, t->num_contigs,
			  t->configs[CALC_CONTIG_POS], &t->tarr, &t_changes);
    } else {
	/* only bother finding template positions if displaying templates */
	/* if (t->configs[TEMPLATES]) { */
	/* but also need to know for readings whether to display readpairs */
	    FindTemplatePositions(io, t->contig_offset, t->contig,
				  t->num_contigs, t->tarr, &t_changes);
	/* } */
    }
    if (-1 == (display_templates(interp, io, t, t_changes)))
	return -1;

    ruler_length = t->contig_offset[t->contig[t->num_contigs-1]].offset +
	io_clength(io, t->contig[t->num_contigs-1]);

#ifdef DEBUG
    printf("RULER length %d %f %f win %s\n", ruler_length, t->world->total->x1,
	   t->world->total->x2, t->window);
#endif
    /* must ensure that the world->total is at least ruler length */
    if (t->world->total->x1 > 1) {
	t->world->total->x1 = 1;
    }
    if (t->world->total->x2 < ruler_length) {
	t->world->total->x2 = ruler_length;
    }

 /*
  * memcpy(t->world->visible, t->world->total, sizeof(d_box));
  * SetCanvasCoords(interp, t->world->visible->x1, t->world->visible->y1,
  *		    t->world->visible->x2, t->world->visible->y2, t->canvas);
  */
    /*
     * to allow zoom not to be reset even when add/delete templates/readings
     */
    if (lengthZoom(t->zoom) <= 1) {
	memcpy(t->world->visible, t->world->total, sizeof(d_box));
	SetCanvasCoords(interp, t->world->visible->x1, t->world->visible->y1,
			t->world->visible->x2, t->world->visible->y2,
			t->canvas);
    } else {
/*
	printf("total2 %f %f %f %f\n", t->world->total->x1,
	       t->world->total->y1, t->world->total->x2, t->world->total->x2);
*/
	SetCanvasCoords(interp, t->world->visible->x1, t->world->visible->y1,
		    t->world->visible->x2, t->world->visible->y2, t->canvas);

    }

    if (t->ruler_coord) {
	int i;
	for (i = 0; i < t->num_contigs; i++) {
	    xfree(t->ruler_coord[i].type);
	}
	xfree(t->ruler_coord);
	t->ruler_coord = NULL;
    }
    display_ruler(interp, io, t->canvas, t->contig_offset, t->contig,
		  t->num_contigs, t->configs[RULER], t->configs[TICKS],
		  t->ruler, t->frame, &t->ruler_coord);

    display_reading_tags(interp, io, t);

    if (t->configs[TEMPLATES] || t->configs[READINGS]) {
	scaleSingleCanvas(t->interp, t->world, t->canvas, t->window, 'b', "all");
    }

    if (t->configs[RULER]) {
	scaleSingleCanvas(t->interp, t->world, t->canvas, t->ruler->window, 'x', "all");
    }

    template_update_cursors(io, t, 0);
    if (t_changes)
	xfree(t_changes);

    return 0;
}

/*
 * find the position (in bases) of the template cursor local to a contig
 */
double TemplateLocalCursor(int id,
			   c_offset *contig_offset,
			   int *contig_array,
			   int num_contigs,
			   double wx)
{
    int i;
    int offset = 0;
    int prev_offset = 0;

    /*
     * a couple of fudges: if num_contigs is 1 then wx is still wx
     * if wx < 0 then must be to the left of the 1st contig and wx is still wx
     */
    if ((num_contigs == 1) || (wx < 0)) {
	return wx;
    }

    for (i = 1; i < num_contigs; i++) {

	prev_offset = offset;
	offset = contig_offset[contig_array[i]].offset;

	if ((wx > prev_offset) && (wx <= offset)) {
	    return (wx - prev_offset);
	}
    }
    /* last contig */
    return (wx - offset);
}


int
template_find_left_position(GapIO *io,
			    int *contig_array,
			    int num_contigs,
			    c_offset *contig_offset,
			    double wx)
{
    int cur_contig;
    int length, prev_len;
    int nearest_contig;
    int i;

    length = 0;
    prev_len = 0;

    for (i = 0; i < num_contigs; i++) {

	cur_contig = contig_array[i];
	prev_len = length;
	if (i+1 == num_contigs) {
	    length += ABS(io_clength(io, contig_array[i]));
	} else {
	    length = contig_offset[contig_array[i+1]].offset;
	}
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


void
ReOrderContigs(int *order,
	       cursor_t **cursor,
	       int *vis,
	       int c_from,
	       int c_to)
{
    int tmpo, tmpv;
    int length;
    cursor_t *tmpc;

    tmpo = order[c_from];
    tmpv = vis[c_from];
    tmpc = cursor[c_from];

    /* if shifting array elements up */
    if (c_from < c_to) {
	 if (c_to != 1) {
             c_to = c_to - 1;
	 }
	length = abs(c_from - c_to);
	memmove(&order[c_from],  &order[c_from]+1,  length * sizeof(*order));
	memmove(&vis[c_from],    &vis[c_from]+1,    length * sizeof(*vis));
	memmove(&cursor[c_from], &cursor[c_from]+1, length * sizeof(*cursor));
    }
    /* if shifting array elements down */
    else {
	length = abs(c_from - c_to);
	memmove(&order[c_to]+1,  &order[c_to],  length * sizeof(*order));
	memmove(&vis[c_to]+1,    &vis[c_to],    length * sizeof(*vis));
	memmove(&cursor[c_to]+1, &cursor[c_to], length * sizeof(*cursor));
    }
    order[c_to] = tmpo;
    vis[c_to] = tmpv;
    cursor[c_to] = tmpc;

}

void
update_template_contig_order(Tcl_Interp *interp,
			     GapIO *io,
			     int template_id,
			     int cx,
			     int *contig_list,
			     int num_contigs)
{
    reg_generic gen;
    obj_template_disp *t;
    double wx, wy;
    int left_index;
    int i, j;
    int recalc;
    /* char cmd[1024]; */

    t = result_data(io, template_id, 0);

    CanvasToWorld(t->canvas, cx, 0, &wx, &wy);
    left_index = template_find_left_position(io, t->contig, t->num_contigs,
					  t->contig_offset, wx);

    for (i = 0; i < num_contigs; i++) {

	for (j = 0; j < t->num_contigs; j++) {
	    if (t->contig[j] == contig_list[i])
		break;
	}
	if (t->num_contigs > 1) {
	    ReOrderContigs(t->contig, t->cursor, t->cursor_visible, j,
			   left_index);
	}
/*
	sprintf(cmd, "UpdateTemplateContigList %d %s %d %d %d %d \n", *handle_io(io),
		t->frame, t->contig[j], t->contig[left_index], j, left_index);

	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    printf("update_template_contig_order: %s\n", interp->result);
*/
    }

    /*
     * update template display - MUST do before other notifications because
     * the contig offsets are generated during the template display
     */
    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_REDRAW;
    recalc = 1;
    gen.data = (void *)&recalc;
    result_notify(io, t->id, (reg_data *)&gen, 0);

    for (i = 0; i < t->num_wins; i++) {
	if (t->win_list[i]->id != t->id) {
	    result_notify(io, t->win_list[i]->id, (reg_data *)&gen, 0);
	}
    }
}

void
refresh_contig_order(Tcl_Interp *interp,
		     GapIO *io,
		     int template_id)
{
    obj_template_disp *t;
    GCardinal *order = ArrayBase(GCardinal, io->contig_order);
    int i, j;
    reg_order ro;
    reg_buffer_start rs;
    reg_buffer_end re;
    int c_from, c_to;

    t = result_data(io, template_id, 0);

    for (i = 1; i < t->num_contigs; i++) {
	c_from = c_to = -1;
	for (j = 0; j < NumContigs(io); j++) {
	    if (order[j] == t->contig[i]) {
		c_from = j;
	    }
	    if (order[j] == t->contig[i-1]) {
		c_to = j;
	    }
	}
	if (c_from != -1 && c_to != -1)
	    ReOrder(io, order, c_from, c_to+1);
    }

    /* Notify of the start of the flurry of updates */
    rs.job = REG_BUFFER_START;
    for (i = 0; i < t->num_contigs; i++) {
	contig_notify(io, t->contig[i], (reg_data *)&rs);
    }

    ro.job = REG_ORDER;
    ro.pos = t->contig[0];

    for (i = 0; i < t->num_contigs; i++)
	contig_notify(io, t->contig[i], (reg_data *)&ro);

    /* Notify the end of our updates */
    re.job = REG_BUFFER_END;
    for (i = 0; i < t->num_contigs; i++) {
	contig_notify(io, t->contig[i], (reg_data *)&re);
    }

    /* save contig order */
    ArrayDelay(io, io->db.contig_order, io->db.Ncontigs, io->contig_order);
    flush2t(io);
}
