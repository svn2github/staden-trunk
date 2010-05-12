#include <limits.h>
#include <string.h>
#include <os.h>
#include <io_lib/Read.h>

#include "trace_print.h"
#include "tk.h"
#include "split.h"

int ps_configure_trace(DNATrace *t, int argc, char **argv) {
    cli_args	a[] = {
	{"-title",      ARG_STR,    1, "",    offsetof(PS_TRACE, title)},
	{"-scale_y",    ARG_DOUBLE, 1, "1.0", offsetof(PS_TRACE, scale_y)},
	{"-scale_x",    ARG_DOUBLE, 1, "1.0", offsetof(PS_TRACE, scale_x)},
	{"-first_base", ARG_INT,    1, "0",   offsetof(PS_TRACE, first_base)},
	{"-last_base",  ARG_INT,    1, "-1",  offsetof(PS_TRACE, last_base)},
	{NULL,          0,          0, NULL,  0}
    };

    if(NULL == t->read) {
	return TCL_ERROR;
    }

    if(-1 == parse_args(a, &(t->ps_trace), argc, argv)) {
	return TCL_ERROR;
    }

    t->ps_trace.title = strdup(t->ps_trace.title);

    t->ps_trace.basePos_lookup =
	trace_index_to_basePos(t->read->basePos, t->read->NBases, t->read->NPoints);

    return TCL_OK;
}

int ps_configure_trace_line(DNATrace *t, int argc, char **argv) {
    LINE_OPTIONS	*options;

    switch(**argv) {
    case 'a': case 'A':
	options = &(t->ps_trace.options_A);
	break;
    case 'c': case 'C':
	options = &(t->ps_trace.options_C);
	break;
    case 'g': case 'G':
	options = &(t->ps_trace.options_G);
	break;
    case 't': case 'T':
	options = &(t->ps_trace.options_T);
	break;
    default:
	options = &(t->ps_trace.options_N);
	break;
    }

    return ps_configure_line(t->interp, t->tkwin, options, argc, argv);
}

int free_ps_trace(PS_TRACE *ps_trace) {
    if(ps_trace->basePos_lookup != NULL) {
	xfree(ps_trace->basePos_lookup);
	ps_trace->basePos_lookup = NULL;
    }
    free_ps_line(&(ps_trace->options_A));
    free_ps_line(&(ps_trace->options_C));
    free_ps_line(&(ps_trace->options_G));
    free_ps_line(&(ps_trace->options_T));
    free_ps_line(&(ps_trace->options_N));

    return 0;
}

int trace_print(DNATrace *t, char *file_name) {
    /*
     * ps_configure() should be called at
     * some point before this function.
     */

    FILE	*ps_file;

    ps_file = ps_fopen(file_name, t->ps_options);
    if(NULL == ps_file) {
	return -1;
    }

    /* Set things common to all four traces. */
    t->ps_trace.trace_height = t->ps_options.panel_height - (2.1 * t->ps_options.font_size);
    t->ps_trace.seq_y_pos = t->ps_options.panel_height - (2 * t->ps_options.font_size);
    t->ps_trace.num_y_pos = t->ps_options.panel_height - t->ps_options.font_size;
    t->ps_trace.scale_x = t->scale_x;

    if (t->read->maxTraceVal) {
	t->ps_trace.scale_y = t->scale_y * ((double) t->ps_trace.trace_height)
	    / (double)t->read->maxTraceVal;
    } else {
	t->ps_trace.scale_y = 0;
    }

    if(-1 == ps_trace_draw_trace(t, ps_file))
	return -1;

    if(fclose(ps_file) == EOF)
	return -1;

    return 0;
}

int ps_trace_draw_trace(DNATrace *t, FILE *ps_file) {
    int		y_trans, n_translations, n_panels, page_n;
    int		i;
    PS_POINT	*p;
    PS_TEXT	*A, *C, *G, *T, *gap, *numbers;
    int		np, np_max, nA, nC, nG, nT, nGap, n_num;
    int		first_base, last_base;
    int		first_point, last_point;

    first_base = ((t->ps_trace.first_base <= 0) || (t->ps_trace.first_base >= t->read->NBases))
	? 0
	: t->ps_trace.first_base - 1;
    last_base = ((t->ps_trace.last_base < first_base) || (t->ps_trace.last_base >= t->read->NBases))
	? (t->read->NBases - 1)
	: t->ps_trace.last_base - 1;
    first_point = t->read->basePos[first_base];
    last_point = t->read->basePos[last_base];

    y_trans = -1 * (t->ps_options.panel_height + t->ps_options.panel_sep);
    n_translations = 0;
    n_panels = 0;
    page_n = 1;
    np_max = (int) (((double) t->ps_options.panel_width) / t->ps_trace.scale_x);
    for(i = first_point; i <= last_point; i += np_max) {
	/* Allow for the fact that the last section of trace  */
	/* may have less points than will fit into one width. */
	np = ((last_point - i + 1) < np_max) ? last_point - i + 1: np_max;

	/* Start a new page. */
	if(n_panels == 0) {
	    ps_newpage(ps_file, t->ps_trace.title, page_n, t->ps_options);
	    n_translations = 0;
	}

	/* Translate to origin of panel. */
	fprintf(ps_file, "%d %d t\n", 0, y_trans);
	n_translations++;

	/* Calculate and draw the trace segment. */
	if(NULL == (p = ps_trace_segment(t->read->traceA, i, np, t->ps_trace.scale_x,
					 t->ps_trace.scale_y, t->ps_trace.trace_height)))
	    return -1;
	ps_draw_lines(ps_file, t->ps_trace.options_A, p, np);
	xfree(p);

	if(NULL == (p = ps_trace_segment(t->read->traceC, i, np, t->ps_trace.scale_x,
					 t->ps_trace.scale_y, t->ps_trace.trace_height)))
	    return -1;
	ps_draw_lines(ps_file, t->ps_trace.options_C, p, np);
	xfree(p);

	if(NULL == (p = ps_trace_segment(t->read->traceG, i, np, t->ps_trace.scale_x,
					 t->ps_trace.scale_y, t->ps_trace.trace_height)))
	    return -1;
	ps_draw_lines(ps_file, t->ps_trace.options_G, p, np);
	xfree(p);

	if(NULL == (p = ps_trace_segment(t->read->traceT, i, np, t->ps_trace.scale_x,
					 t->ps_trace.scale_y, t->ps_trace.trace_height)))
	    return -1;
	ps_draw_lines(ps_file, t->ps_trace.options_T, p, np);
	xfree(p);

	/* Now find the bases identified in this segment, and draw them. */
       	if(-1 == ps_sequence_segment(t, i, np, &A, &C, &G, &T, &gap, &nA, &nC, &nG, &nT, &nGap))
	    return -1;
	ps_draw_text(ps_file, A, nA, t->ps_trace.options_A.rgb, 'c');
	ps_draw_text(ps_file, C, nC, t->ps_trace.options_C.rgb, 'c');
	ps_draw_text(ps_file, G, nG, t->ps_trace.options_G.rgb, 'c');
	ps_draw_text(ps_file, T, nT, t->ps_trace.options_T.rgb, 'c');
	ps_draw_text(ps_file, gap, nGap, t->ps_trace.options_N.rgb, 'c');
	xfree_ps_text(A, nA);
	xfree_ps_text(C, nC);
	xfree_ps_text(G, nG);
	xfree_ps_text(T, nT);
	xfree_ps_text(gap, nGap);

	/* Now number every 10th base */
	if(-1 == ps_numbers_segment(t, i, np, &numbers, &n_num))
	    return -1;
	ps_draw_text(ps_file, numbers, n_num, t->ps_trace.options_N.rgb, 'e');
	xfree_ps_text(numbers, n_num);

	/* Check the number of panels against the maximum per page. */
	n_panels++;
	if(n_panels >= t->ps_options.n_panel) {
	    /* Translate back to start. */
	    fprintf(ps_file, "%d %d t\n", 0, -1 * (n_translations * y_trans));
	    n_translations = 0;

	    n_panels = 0;
	    ps_finishpage(ps_file);
	    page_n++;
	}
    }
    if((n_panels > 0) && (n_panels < t->ps_options.n_panel)) {
	    ps_finishpage(ps_file);
    }

    return 0;
}

PS_POINT *ps_trace_segment(TRACE *trace, int first, int NPoints, double x_scale, double y_scale, int y_max) {
    PS_POINT	*p;
    int		i;

    p = (PS_POINT *) xmalloc(NPoints * sizeof(PS_POINT));
    if(NULL == p) {
	return NULL;
    }

    for(i = 0; i < NPoints; i++) {
	p[i].x = (int) (i * x_scale);
	p[i].y = (int) (trace[i + first] * y_scale);
	if(p[i].y > y_max) {
	    p[i].y = y_max;
	}
    }

    return p;
}

int ps_sequence_segment(DNATrace *t, int first, int NPoints,
			  PS_TEXT **A, PS_TEXT **C, PS_TEXT **G, PS_TEXT **T, PS_TEXT **gap,
			  int *nA, int *nC, int *nG, int *nT, int *nGap) {
    /* Note that 'first' and 'NPoints' refer to the index of the */
    /* first item in TRACE and the number of TRACE items to look */
    /* at, not to values in the base and basePos arrays.         */
    int		i, j;

    /* Get to the first base that appears on the segment of TRACE data. */
    j = first;
    while((t->ps_trace.basePos_lookup[j] == -1) && (j < first + NPoints)) {j++;}
    i = t->ps_trace.basePos_lookup[j];

    /* Enter all those bases contained in this segment of the trace */
    /* into an array, along with their positions on the segment.    */
    *nA = *nC = *nG = *nT = *nGap = 0;
    if(NULL == (*A = (PS_TEXT *) xmalloc(NPoints * sizeof(PS_TEXT))))
	return -1;
    if(NULL == (*C = (PS_TEXT *) xmalloc(NPoints * sizeof(PS_TEXT))))
	return -1;
    if(NULL == (*G = (PS_TEXT *) xmalloc(NPoints * sizeof(PS_TEXT))))
	return -1;
    if(NULL == (*T = (PS_TEXT *) xmalloc(NPoints * sizeof(PS_TEXT))))
	return -1;
    if(NULL == (*gap = (PS_TEXT *) xmalloc(NPoints * sizeof(PS_TEXT))))
	return -1;

    while((t->read->basePos[i] < (first + NPoints)) && (i < t->read->NBases)) {
	switch(t->read->base[i]){
	case 'A':
	case 'a':
	    char_to_ps_text(&((*A)[*nA]), t->read->base[i],
			    (int) ((t->read->basePos[i] - first) * t->ps_trace.scale_x),
			    t->ps_trace.seq_y_pos);
	    (*nA)++;
	    break;
	case 'C':
	case 'c':
	    char_to_ps_text(&((*C)[*nC]), t->read->base[i],
			    (int) ((t->read->basePos[i] - first) * t->ps_trace.scale_x),
			    t->ps_trace.seq_y_pos);
	    (*nC)++;
	    break;
	case 'G':
	case 'g':
	    char_to_ps_text(&((*G)[*nG]), t->read->base[i],
			    (int) ((t->read->basePos[i] - first) * t->ps_trace.scale_x),
			    t->ps_trace.seq_y_pos);
	    (*nG)++;
	    break;
	case 'T':
	case 't':
	    char_to_ps_text(&((*T)[*nT]), t->read->base[i],
			    (int) ((t->read->basePos[i] - first) * t->ps_trace.scale_x),
			    t->ps_trace.seq_y_pos);
	    (*nT)++;
	    break;
	default:
	    char_to_ps_text(&((*gap)[*nGap]), t->read->base[i],
			    (int) ((t->read->basePos[i] - first) * t->ps_trace.scale_x),
			    t->ps_trace.seq_y_pos);
	    (*nGap)++;
	    break;
	}
	i++;
    }

    if(NULL == (*A = (PS_TEXT *) xrealloc(*A, *nA * sizeof(PS_TEXT) + 1)))
	return -1;
    if(NULL == (*C = (PS_TEXT *) xrealloc(*C, *nC * sizeof(PS_TEXT) + 1)))
	return -1;
    if(NULL == (*G = (PS_TEXT *) xrealloc(*G, *nG * sizeof(PS_TEXT) + 1)))
	return -1;
    if(NULL == (*T = (PS_TEXT *) xrealloc(*T, *nT * sizeof(PS_TEXT) + 1)))
	return -1;
    if(NULL == (*gap = (PS_TEXT *) xrealloc(*gap, *nGap * sizeof(PS_TEXT) + 1)))
	return -1;

    return 0;
}

int ps_numbers_segment(DNATrace *t, int first, int NPoints, PS_TEXT **numbers, int *n_num) {
    int	i, j, first_base, last_base;

    /* Find the index of the first base in basePos that appears on the segment of TRACE data. */
    i = first;
    while((t->ps_trace.basePos_lookup[i] == -1) && (i < first + NPoints)) {i++;}
    first_base = t->ps_trace.basePos_lookup[i];

    /* Find the index of the last base in basePos that appears on the segment of TRACE data. */
    i = first + NPoints - 1;
    while((t->ps_trace.basePos_lookup[i] == -1) && (i >= first)) {i--;}
    last_base = t->ps_trace.basePos_lookup[i];

    if(NULL == (*numbers = (PS_TEXT *) xmalloc(NPoints * sizeof(PS_TEXT))))
	return -1;

    *n_num = 0;
    i = 0;
    while(i <= last_base - first_base) {
	/* If sequence is complemented, then count down. */
	j = t->comp ? last_base - i: first_base + i;
	/* Add one because sequence numbering should start at one, not zero.   */
	if((j + 1) % 10 == 0) {
	    int_to_ps_text(&((*numbers)[*n_num]), (j + 1),
			   (int) ((t->read->basePos[j] - first) * t->ps_trace.scale_x),
			   t->ps_trace.num_y_pos);
	    (*n_num)++;
	}
	i++;
    }
    if(NULL == (*numbers = (PS_TEXT *) xrealloc(*numbers, *n_num * sizeof(PS_TEXT) + 1)))
	return -1;

    return 0;
}

