/*
 * File: readpair.c:
 * Version: 2.0 (1.0 == FORTRAN version)
 *
 * Author: Staden Package group
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: "Find read pair" code.
 *
 * Created: 17 March 1994
 * Updated:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "gap-dbstruct.h"
#include "IO.h"
#include "fortran.h"
#include "xalloc.h"
#include "io-reg.h"
#include "cs-object.h"
#include "contig_selector.h"
#include "newgap_cmds.h"
#include "contigEditor.h"
#include "gap_globals.h"
#include "misc.h"
#include "template.h"
#include "text_output.h"
#include "tclXkeylist.h"
#include "complement.h"
#include "tcl_utils.h"
#include "readpair.h"

/*
 * Match callback.
 * 'obj' is a match contained within the 'repeat' list.
 */
void *readpair_obj_func(int job, void *jdata, obj_read_pair *obj,
			mobj_template *template) {
    static char buf[80];
    char cmd[1024];
    f_int *handle;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(template->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(template->io, cs_id, 0);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	if (io_rdonly(template->io) && ((obj->c1 > 0 && obj->c2 < 0) ||
					(obj->c1 < 0 && obj->c2 > 0))) {
	    return "Information\0Hide\0IGNORE\0"
		"IGNORE\0SEPARATOR\0Remove\0";
	} else {
	    return "Information\0Hide\0Invoke template display *\0"
		"Invoke join editor\0SEPARATOR\0Remove\0";
	}
	break;

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */
	    start_message();

	    vmessage("Read pair:\n");
	    vmessage("    From contig %s(#%d) at %d reading %s(#%d)\n",
		     get_contig_name(template->io, ABS(obj->c1)),
		     io_clnbr(template->io, ABS(obj->c1)), obj->pos1,
		     get_read_name(template->io, obj->read1), obj->read1);
	    vmessage("    With contig %s(#%d) at %d reading %s(#%d)\n",
		     get_contig_name(template->io, ABS(obj->c2)),
		     io_clnbr(template->io, ABS(obj->c2)), obj->pos2,
		     get_read_name(template->io, obj->read2), obj->read2);
	    {
		GReadings r1, r2;

		gel_read(template->io, obj->read1, r1);
		gel_read(template->io, obj->read2, r2);

		vmessage("    Direction of first read is %swards\n",
			 STRAND(r1) == GAP_STRAND_FORWARD ? "for" : "back");
		vmessage("    Direction of second read is %swards\n",
			 STRAND(r2) == GAP_STRAND_FORWARD ? "for" : "back");
	    }
	    vmessage("    Length %d\n\n", obj->length);
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, (obj_match *)obj,
		     (mobj_repeat *)template, csplot_hash);
	    break;

	case -2: /* default */
	case 2: /* Invoke template display */
	    {
                int cnum[2];
		char c_list[1024];
		char r_list[1024];
		char c1_name[100];
		char c2_name[100];
		GReadings r1, r2;

		obj->flags |= OBJ_FLAG_VISITED;
		template->current = obj - (obj_read_pair *)template->match;
		Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(template),
			    NULL);

                cnum[0] = ABS(obj->c1);
                cnum[1] = ABS(obj->c2);

                /* Complement a contig if needed */
		/* if it fails to complement the first contig for any reason,
		 * eg the contig editor is up, then it tries to complement
		 * the second
		 */
                if ((obj->c1 > 0) != (obj->c2 > 0)) {
		    if (io_rdonly(template->io)) {
			bell();
			break;
		    }
		    if (io_clength(template->io, ABS(obj->c1)) <
			io_clength(template->io, ABS(obj->c2))) {
			if (-1 == complement_contig(template->io,
						    ABS(obj->c1)))
			    if (-1 == complement_contig(template->io,
							ABS(obj->c2)))
				return NULL;
		    } else {
			if (-1 == complement_contig(template->io,
						    ABS(obj->c2)))
			    if (-1 == complement_contig(template->io,
							ABS(obj->c1)))
				return NULL;
		    }
                }

		/* FIXME ??? */
		if ((handle = handle_io(template->io)) == NULL) {
		    verror(ERR_FATAL, "readpairs", "invalid io handle");
		}

		gel_read(template->io, obj->read1, r1);
		gel_read(template->io, obj->read2, r2);

		/*
		 * determine the sense of each reading to find the correct
		 * orientation to plot the contigs in the template display.
		 * The reading to the left should have sense = 0; the reading
		 * to the right should have sense = 1
		 */
		if (r1.sense) {
		    strcpy(c2_name, get_contig_name(template->io,
						    ABS(obj->c1)));
		    strcpy(c1_name, get_contig_name(template->io,
						    ABS(obj->c2)));
		} else {
		    strcpy(c1_name, get_contig_name(template->io,
						    ABS(obj->c1)));
		    strcpy(c2_name, get_contig_name(template->io,
						    ABS(obj->c2)));
		}
		sprintf(c_list, "%s %s", c1_name, c2_name);
		sprintf(r_list, "%d %d", obj->read1, obj->read2);

		sprintf(cmd, "CreateTemplateDisplay %d {%s}",
			*handle, c_list);
		if (TCL_ERROR == Tcl_Eval(GetInterp(), cmd))
		    printf("%s\n", GetInterp()->result);
	    }
	    break;

	case 3: /* Join editor */
	    {
	        int cnum[2], llino[2], pos[2];
		GReadings r1, r2;
		int pt1, pt2;

		cnum[0] = ABS(obj->c1);
		cnum[1] = ABS(obj->c2);

		/* Complement a contig if needed */
		if ((obj->c1 > 0) != (obj->c2 > 0)) {
		    if (cnum[0] == cnum[1]) {
			verror(ERR_WARN, "join_editor",
			       "cannot display the same contig in two "
			       "different orientations");
			break;
		    }

		    if (-1 == complement_contig(template->io, cnum[0]))
			if (-1 == complement_contig(template->io, cnum[1]))
			    return NULL;
		}

		/* Get strand information */
		gel_read(template->io, obj->read1, r1);
		gel_read(template->io, obj->read2, r2);

		pt1 = STRAND(r1);
		pt2 = STRAND(r2);

		/* Set position based on strand */
		if (pt1 != pt2) {
		    if (r1.sense == 0) {
		        pos[0] = io_clength(template->io, ABS(obj->c1)) + 1;
			pos[1] = 1;
		    } else {
			pos[0] = 1;
		        pos[1] = io_clength(template->io, ABS(obj->c2)) + 1;
		      }
		} else {
		    if (pt1 == r1.sense) {
		        pos[0] = 1;
			pos[1] = 1;
		    } else {
		      pos[0] = io_clength(template->io, ABS(obj->c1)) + 1;
		      pos[1] = io_clength(template->io, ABS(obj->c2)) + 1;
		    }
		}

		llino[0] = io_clnbr(template->io, cnum[0]);
		llino[1] = io_clnbr(template->io, cnum[1]);

		join_contig(GetInterp(), template->io, cnum, llino, pos,
			    consensus_cutoff, quality_cutoff);
		break;
	    }

	case 4: /* Remove */
	    obj_remove(GetInterp(), cs->window, (obj_match *)obj,
		       (mobj_repeat *)template, csplot_hash);
	    break;

	}
	break;

    case OBJ_GET_BRIEF:
	sprintf(buf, "Read pair: %c#%d@%d with %c#%d@%d, len %d",
		obj->c1 > 0 ? '+' : '-', obj->read1, obj->pos1,
		obj->c2 > 0 ? '+' : '-', obj->read2, obj->pos2,
		obj->length);
	return buf;
    }

    return NULL;
}

void readpair_callback(GapIO *io, int contig, void *fdata, reg_data *jdata) {
    mobj_template *r = (mobj_template *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id, 0);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Find read pairs");
	break;


    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, (mobj_repeat *)r,
			csplot_hash, cs->window);
	break;


    case REG_COMPLEMENT:

	csmatch_complement(io, contig, r, csplot_hash, cs->window);
	break;


    case REG_GET_OPS:

	if (r->all_hidden)
	    jdata->get_ops.ops = "PLACEHOLDER\0PLACEHOLDER\0Information\0"
		"PLACEHOLDER\0Hide all\0Reveal all\0SEPARATOR\0Remove\0";
	else
	    jdata->get_ops.ops = "Use for 'Next'\0Reset 'Next'\0Information\0"
		"Configure\0Hide all\0Reveal all\0SEPARATOR\0Remove\0";
	break;


    case REG_INVOKE_OP:

	switch (jdata->invoke_op.op) {
	case 0: /* Next */
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(r), NULL);
	    break;
	case 1: /* Reset Next */
	    csmatch_reset_next((mobj_repeat *)r);
	    break;
	case 2: /* Information */
	    csmatch_info((mobj_repeat *)r, "Read pair");
	    break;
	case 3: /* Configure */
	    csmatch_configure(io, cs->window, (mobj_repeat *)r);
	    break;
	case 4: /* Hide all */
	    csmatch_hide(GetInterp(), cs->window, (mobj_repeat *)r,
			 csplot_hash);
	    break;
	case 5: /* Reveal all */
	    csmatch_reveal(GetInterp(), cs->window, (mobj_repeat *)r,
			   csplot_hash);
	    break;
	case 6: /* Remove */
	    csmatch_remove(io, cs->window,
			   (mobj_repeat *)r,
			   csplot_hash);
	    break;
	}
	break;


    case REG_PARAMS:

	jdata->params.string = r->params;
	break;


    case REG_NUMBER_CHANGE:

	csmatch_renumber(io, contig, jdata->number.number,
			 (mobj_repeat *)r, csplot_hash, cs->window);
	break;


    case REG_ORDER:

	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs->window,
		       (mobj_repeat *)r,
		       csplot_hash);
	break;


    case REG_DELETE:

	csmatch_contig_delete(io, (mobj_repeat *)r, contig,
			      cs->window, csplot_hash);
	break;

    case REG_LENGTH:
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;
    }
}

/*
 * Plot the templates which span more than one contig on contig selector plot.
 * The direction of the line is governed by the "direction" of the template
 * A line going from top right to bottom left (/) indicates the second contig
 * needs completing
 * A line going from top left to bottom right (\) indicates the second contig
 * does not need completing
 */
int PlotTempMatches(GapIO *io, template_c **tarr) {
    int i, j, k, id;
    int *contig, *pos, *dir, *length, *readnum;
    mobj_template *template;
    obj_read_pair *matches;
    int array_size = Ncontigs(io), array_inc = 10;
    int count;
    int n_matches = 0;
    int max_matches = 64; /* dynamically grows */
    gel_cont_t *gc;
    item_t *item, *it;
    char *val;

    if (NULL == (contig = (int *)xmalloc(array_size * sizeof(int))))
	return -1;
    if (NULL == (pos = (int *)xmalloc(array_size * sizeof(int))))
	return -1;
    if (NULL == (dir = (int *)xmalloc(array_size * sizeof(int))))
	return -1;
    if (NULL == (length = (int *)xmalloc(array_size * sizeof(int))))
	return -1;
    if (NULL == (readnum = (int *)xmalloc(array_size * sizeof(int))))
	return -1;
    if (NULL == (template = (mobj_template *)xmalloc(sizeof(mobj_template))))
	return -1;
    if (NULL == (matches = (obj_read_pair *)xmalloc(max_matches *
						    sizeof(obj_read_pair))))
	return -1;

    for (i = 1; i <= Ntemplates(io); i++) {
	if (tarr[i]) {
	    count = 0;

	    for (item = head(tarr[i]->gel_cont); item; item = item->next) {
		int contig_tmp, strand;
		GReadings r;

		gc = (gel_cont_t *)(item->data);

		if (gc->contig < 0)
		    continue;

		gel_read(io, gc->read, r);
		contig[count] = gc->contig;
		readnum[count] = gc->read;
		pos[count]    = r.position;
		strand = (PRIMER_TYPE_GUESS(r) == GAP_PRIMER_FORWARD ||
			  PRIMER_TYPE_GUESS(r) == GAP_PRIMER_CUSTFOR) ? 0 : 1;
		dir[count]    = r.sense == r.strand ? 1 : -1;
		length[count] = r.sequence_length;

		contig_tmp = gc->contig;
		gc->contig = -gc->contig;
		for (it = item->next; it; it = it->next) {
		    gc = (gel_cont_t *)(item->data);
		    if (gc->contig == contig_tmp)
			gc->contig = -gc->contig;
		}

		/*
		 * check that count is less than Ncontigs
		 * if it is larger, then need to allocate more space to
		 * arrays
		 */
		count++;
	        if (count >= array_size) {
	            /* set array_size to be count plus arbitary increment */
	            array_size = count + array_inc;
	            if (NULL == (contig = (int *)realloc(contig, array_size *
							 sizeof(int))))
	                return -1;
	            if (NULL == (pos = (int *)realloc(pos, array_size *
						      sizeof(int))))
	                return -1;
	            if (NULL == (dir = (int *)realloc(dir, array_size *
						      sizeof(int))))
	                return -1;
	            if (NULL == (length = (int *)realloc(length, array_size *
							 sizeof(int))))
	                return -1;
	            if (NULL == (readnum = (int *)realloc(readnum,array_size *
							  sizeof(int))))
	                return -1;
	        }
	    }

	    /* compare each contig with every other */
	    for (j = 0; j < count-1; j++) {
		for (k = j + 1; k < count; k++) {
		    /*
		     * check that individual pairs of contigs are in different
		     * contigs
		     */
		    if (contig[j] == contig[k]) {
			continue;
		    }
		    matches[n_matches].func =
			(void *(*)(int, void *, struct obj_match_t *,
				   struct mobj_repeat_t *))readpair_obj_func;
		    matches[n_matches].data = template;
		    matches[n_matches].c1 = dir[j] * contig[j];
		    matches[n_matches].c2 = dir[k] * contig[k];
		    matches[n_matches].pos1 = pos[j];
		    matches[n_matches].pos2 = pos[k];
		    matches[n_matches].length = (length[j] + length[k]) /2;
		    matches[n_matches].read1 = readnum[j];
		    matches[n_matches].read2 = readnum[k];
		    matches[n_matches].flags = 0;
		    if (++n_matches >= max_matches) {
			obj_read_pair *tmp;

			max_matches *= 2;
			tmp = (obj_read_pair *)
			    xrealloc(matches,
				     max_matches * sizeof(obj_read_pair));
			if (NULL == tmp) {
			    xfree(contig);
			    xfree(pos);
			    xfree(dir);
			    xfree(length);
			    xfree(readnum);
			    xfree(template);
			    xfree(matches);
			    return -1;
			} else {
			    matches = tmp;
			}
		    }
		}
	    }
	}
    }

    /*
     * Free memory and return if we've got no matches.
     */
    if (0 == n_matches) {
	xfree(contig);
	xfree(pos);
	xfree(dir);
	xfree(length);
	xfree(readnum);
	xfree(template);
	xfree(matches);

	return 0;
    }

    template->num_match = n_matches;
    template->match = (obj_match *)matches;
    template->io = io;
    strcpy(template->tagname, CPtr2Tcl(template));

    val = get_default_string(GetInterp(), gap_defs, "READPAIR.COLOUR");
    strcpy(template->colour, val);

    template->linewidth = get_default_int(GetInterp(), gap_defs,
					  "READPAIR.LINEWIDTH");

    template->params = (char *)xmalloc(10);
    if (template->params)
	sprintf(template->params, "none");
    template->all_hidden = 0;
    template->current = -1;
    template->reg_func = readpair_callback;
    template->match_type = REG_TYPE_READPAIR;

    PlotRepeats(io, template);
    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(template), NULL);

    xfree(contig);
    xfree(pos);
    xfree(dir);
    xfree(length);
    xfree(readnum);

    /*
     * Register the repeat search with each of the contigs used.
     * Currently we assume that this is all.
     */
    id = register_id();
    for (i = 1; i <= NumContigs(io); i++) {
	contig_register(io, i, readpair_callback, (void *)template, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ORDER, REG_TYPE_READPAIR);
    }

    return(0);
}

static void do_report(GapIO *io, template_c **tarr, int *order) {
    int i, listed_problems = 0;
    char name[DB_NAMELEN + 1];
    GTemplates t;
    GReadings r;
    item_t *item;
    char *mode;
    int length;

    for (i = 0; order[i]; i++) {
	template_c *tc = tarr[order[i]];
	template_read(io, tc->num, t);
	TextRead(io, t.name, name, sizeof(name)-1);

	if (tc->consistency && !listed_problems){
	    vmessage("*** Possibly problematic templates listed below ***\n");
	    listed_problems = 1;
	}

	/*
	 * Special case: when we are spanning contigs, and we have a
	 * forward read in one and a reverse read in the other, then we
	 * computed the possible distance based upon the distance to the
	 * end of each contig.
	 */
	mode = NULL;
	if (tc->flags & TEMP_FLAG_SPANNING) {
	    int c1 = 0, distf = 0, distr = 0;
	    GReadings r;

	    for (item = head(tc->gel_cont); item; item = item->next) {
		gel_cont_t *gc = (gel_cont_t *)item->data;

		if (!c1)
		    c1 = gc->contig;
		else if (gc->contig == c1)
		    continue;

		gel_read(io, gc->read, r);

		switch(PRIMER_TYPE(r)) {
		case GAP_PRIMER_UNKNOWN:
		    continue;

		case GAP_PRIMER_FORWARD:
		case GAP_PRIMER_CUSTFOR:
		    distf = r.sense
			? r.position + r.sequence_length - 1
			: io_clength(io, gc->contig) - r.position + 1;
		    break;

		case GAP_PRIMER_REVERSE:
		case GAP_PRIMER_CUSTREV:
		    distr = r.sense
			? r.position + r.sequence_length - 1
			: io_clength(io, gc->contig) - r.position + 1;
		}
	    }

	    if (distf && distr) {
		mode = "computed";
		length = distf + distr;

		if (length < t.insert_length_min ||
		    length > t.insert_length_max) {
		    tc->consistency |= TEMP_CONSIST_DIST;
		}
	    }
	}

	if (mode == NULL) {
	    mode = (tc->flags & TEMP_FLAG_EXPECTED) ? "expected" : "observed";
	    length = ABS(tc->end - tc->start);
	}

	/* Template line */
	vmessage("Template %12s(%4d), length %d-%d(%s %d) score %f\n",
		 name, tc->num, t.insert_length_min, t.insert_length_max,
		 mode, length, tc->score);

	/* Reading lines */
	for (item = head(tc->gel_cont); item; item = item->next) {
	    gel_cont_t *gc = (gel_cont_t *)item->data;

	    gel_read(io, gc->read, r);
	    strcpy(name, io_rname(io, gc->read));

	    vmessage("%c%c%c%c Reading %*s(%+5d%c), pos %6d%+5d, contig %4d\n",
		     (tc->consistency & TEMP_CONSIST_UNKNOWN) ? '?' : ' ',
		     (tc->consistency & TEMP_CONSIST_DIST)    ? 'D' : ' ',
		     (tc->consistency & TEMP_CONSIST_PRIMER)  ? 'P' : ' ',
		     (tc->consistency & TEMP_CONSIST_STRAND)  ? 'S' : ' ',
		     DB_NAMELEN, name,
		     gc->read * (r.sense ? -1 : 1),
		     "?FRfr"[PRIMER_TYPE_GUESS(r)],
		     r.position, r.end - r.start - 1,
		     chain_left(io, gc->read));
	}

	vmessage("\n");
    }
}

/*
 * External callable functions - the C interface
 * ---------------------------------------------
 */
int find_read_pairs(GapIO *io, int num_contigs, int *contig_array) {
    int *order;
    template_c **tarr;

    /* Initialise templates to all */
    if (NULL == (tarr = init_template_checks(io, num_contigs, contig_array,0)))
	return -1;

    /*
     * Strip down to only those with multiple readings. Then sort on
     * problem and report.
     */
    remove_single_templates(io, tarr);
    check_all_templates(io, tarr);
    if (NULL == (order = sort_templates(io, tarr))) {
	uninit_template_checks(io, tarr);
	return -1;
    }
    do_report(io, tarr, order);

    /* Find only those templates spanning multiple contigs and plot them. */
    contig_spanning_templates(io, tarr);
    PlotTempMatches(io, tarr);

    /* Tidy up */
    uninit_template_checks(io, tarr);
    xfree(order);

    return 0;
}
