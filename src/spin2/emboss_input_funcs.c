#include <tcl.h>
#include <unistd.h>
#include <string.h>
#include <float.h>

#include "xalloc.h"
#include "text_output.h"
#include "seq_results.h"
#include "seq_plot_funcs.h"
#include "emboss_input_funcs.h"
#include "tcl_utils.h"
#include "spin_globals.h"
#include "sequence_pair_display.h"
#include "tclCanvGraph.h"
#include "seq_element.h"

/* emboss colours taken from ajax/ajgraph.c */
static char *colournum[] = { "BLACK", "RED", "YELLOW", "GREEN", "AQUAMARINE",
			     "PINK", "WHEAT", "GREY", "BROWN", "BLUE", 
			     "BLUEVIOLET", "CYAN", "TURQUOISE", "MAGENTA", 
			     "SALMON", "WHITE"};

void emboss_graph_callback(int seq_num, void *obj, seq_reg_data *jdata);

void emboss_graph_shutdown(Tcl_Interp *interp,
			   seq_result *s_result,
			   element *e)
{
    in_emboss *input = s_result->input;
    seq_reg_key_name info;
    static char buf[80];

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    seq_deregister(GetSeqNum(s_result->seq_id[HORIZONTAL]), 
		   emboss_graph_callback, (seq_result *)s_result);

    if (s_result->seq_id[VERTICAL] != -1) {
	seq_deregister(GetSeqNum(s_result->seq_id[VERTICAL]), 
		       emboss_graph_callback, (seq_result *)s_result);
    }

    if (e->num_results > 0) {
	e->replot_func(e);
    }
    
    if (s_result->graph == SEQ_E_DOT) {
	DestroySequencePairDisplay(interp, s_result->id);
    }
    free(input->params);
    xfree(s_result->input);
    xfree(s_result);
}

void emboss_graph_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_emboss *input = s_result->input;
    int result_id = s_result->id;
    char cmd[1024];
    element *e = s_result->e;
    plot_data *result;
    Tcl_Interp *interp;
    text_emboss *text_data = (text_emboss *)s_result->text_data;

    if (e) {
	result = find_plot_data(e, result_id);
	interp = e->c->interp;
    }

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Emboss graph plot");
	break;    

    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "%s #%d", result->name, result_id);
	break;

    case SEQ_GET_BRIEF:
	{
	    jdata->name.line[0] = '\0';
	    if (text_data->title)
		strcat(jdata->name.line, text_data->title);

	    if (text_data->maintitle)
		strcat(jdata->name.line, text_data->maintitle);
	}
	break;
	
    case SEQ_GET_OPS:
	if (result->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0"
		"PLACEHOLDER\0PLACEHOLDER\0PLACEHOLDER\0Reveal\0SEPARATOR\0Remove\0";
	} else {
	    if (s_result->graph == SEQ_E_DOT) {
		jdata->get_ops.ops = "Information\0List results\0Configure\0"
		    "Display sequences\0Hide\0PLACEHOLDER\0SEPARATOR\0Remove\0";
	    } else {
		jdata->get_ops.ops = "Information\0List results\0Configure\0"
		    "PLACEHOLDER\0Hide\0PLACEHOLDER\0SEPARATOR\0Remove\0";
	    }
	}
	break;
    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* information */
	    vfuncheader("input parameters");
	    vmessage("%s\n", input->params);
	    break;
	case 1: /* results */
	    Tcl_Eval(interp, "SetBusy");
	    vfuncheader("results");
	    s_result->txt_func(s_result);
	    Tcl_Eval(interp, "ClearBusy");
	    break;
	case 2: /* configure */
	    sprintf(cmd, "result_config %d %d %s %d %s", 
		    result_id, result->line_width, result->colour, e->id, 
		    e->win);
	    if (TCL_OK != Tcl_Eval(interp, cmd)){
		puts(interp->result);
	    }
	    break;
	case 3: /* display sequences */
	    SequencePairDisplay(interp, e->win, result_id, 
				s_result->seq_id[HORIZONTAL], 
				s_result->seq_id[VERTICAL]);
	    break;
	case 4: /* hide all */
	    result->hidden = 1;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 5: /* reveal all */
	    result->hidden = 0;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 6: /* remove */ 
	    {
		emboss_graph_shutdown(interp, s_result, e);
		remove_result_from_element(e, result_id);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	s_result->pr_func(s_result, (seq_reg_plot *)jdata);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case RESULT:
	    jdata->info.result = (void *)s_result;
	    break;
	case WIN_NAME:
	    {
		char *r_win = e->win;
		jdata->info.result = (void *)r_win;

		break;
	    }
	case WIN_SIZE:
	    {
		static d_point pt;
		if (s_result->type == SEQ_TYPE_GRAPH_PLOT) {
		    pt.x = get_default_int(interp, tk_utils_defs, 
					    w("ELEMENT.PLOT_WIDTH"));
		    pt.y = get_default_double(interp, tk_utils_defs,
					       w("ELEMENT.PLOT_HEIGHT"));
		} else if (s_result->type == SEQ_TYPE_DOT_PLOT) {
		    pt.x = get_default_int(interp, tk_utils_defs, 
					    w("DOT.PLOT_WIDTH"));
		    pt.y = get_default_double(interp, tk_utils_defs,
					       w("DOT.PLOT_HEIGHT"));

		}
		jdata->info.result = (void *)&pt;
		break;
	    } /* WIN_SIZE */
	} /* switch SEQ_RESULT_INFO */
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    emboss_graph_shutdown(interp, s_result, e);
	    break;
	} /* SEQ_DELETE, SEQ_QUIT */
    } /* switch */
}

void emboss_graph_text_func(void *obj)
{

}

int store_emboss_graph(int seq_num,
		       int start,
		       int end,
		       e_graph *data,
		       in_emboss *input,
		       text_emboss *text_data,
		       Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int id;

#ifdef DEBUG
    for (i = 0; i < data->n_pts; i++) {
	printf("i %d pos %d score %f\n", i, data->p_array[i].pos, 
	       data->p_array[i].score);
    }
#endif

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    id = get_reg_id();

    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    s_result->seq_id[VERTICAL] = -1;

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->type = SEQ_TYPE_GRAPH_PLOT;
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj;
    s_result->text_data = text_data;
    s_result->e = NULL;

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = emboss_graph_callback;
    s_result->txt_func = emboss_graph_text_func;

    seq_register(seq_num, emboss_graph_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int read_emboss_data_file(char *filename,
			  e_graph **out_data,
			  Tcl_Obj **graph_obj,
			  text_emboss **text_data)
{
    FILE *fp;
    int i, j;
    double pos;
    char line[1024];
    char dummy[100];
    char plot_type[50];
    char title[1025];
    int num_graphs;
    int graph_num;
    double xmina, xmaxa, ymina, ymaxa;
    double xmin, xmax, ymin, ymax;
    double scale_xmin, scale_xmax, scale_ymin, scale_ymax;
    char maintitle[1024];
    char subtitle[1024];
    char xtitle[1024];
    char ytitle[1024];
    e_graph *data;
    e_obj *d_obj = NULL;
    e_obj *g_obj = NULL;
    p_score *p_array = NULL;
    char obj_type[100];
    int colour;
    Graph *graph;

    if (NULL == (fp = fopen(filename, "r"))) {
	verror(ERR_WARN, "emboss", "unable to open file %s \n", filename);
	return -1;
    }

    if (NULL == (data = (e_graph *)xmalloc(sizeof(e_graph))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;

    if (NULL == (*text_data = (text_emboss *)xmalloc(sizeof(text_emboss))))
	return -1;

    Tcl_InitGraph(&graph);

    data->n_pts = 0;
    data->n_data_obj = 0;
    data->n_graph_obj = 0;
    data->title = NULL;
    data->maintitle = NULL;
    data->subtitle = NULL;
    data->xtitle = NULL;
    data->ytitle = NULL;

    (*text_data)->title = NULL;
    (*text_data)->maintitle = NULL;
    (*text_data)->subtitle = NULL;
    (*text_data)->xtitle = NULL;
    (*text_data)->ytitle = NULL;

    /* read in dat file */
    fgets(line, 1024, fp);
    strcpy(plot_type, &line[2]);
    
    while (fgets(line, 1024, fp)) {
	if (strncmp(line, "##Title", 7) ==0) {
	    strcpy(title, &line[8]);
	} else if (strncmp(line, "##Graphs", 8) == 0) {
	    sscanf(line, "##Graphs %d\n", &num_graphs);
	} else if (strncmp(line, "##Number", 8) == 0) {
	    sscanf(line, "##Number %d\n", &graph_num);
	} else if (strncmp(line, "##Points", 8) == 0) {
	    sscanf(line, "##Points %d\n", &data->n_pts);
#ifdef DEBUG
	    printf("NUM_PTS %d\n", data->n_pts);
#endif

	} else if (strncmp(line, "##XminA", 7) == 0) {
	    sscanf(line, "##XminA %lf XmaxA %lf YminA %lf YmaxA %lf\n",
		   &xmina, &xmaxa, &ymina, &ymaxa);
	} else if (strncmp(line, "##Xmin", 6) == 0) {
	    sscanf(line, "##Xmin %lf Xmax %lf Ymin %lf Ymax %lf\n",
		   &xmin, &xmax, &ymin, &ymax);
	} else if (strncmp(line, "##ScaleXmin", 11) == 0) {
	    sscanf(line, "##ScaleXmin %lf ScaleXmax %lf ScaleYmin %lf ScaleYmax %lf\n", 
		   &scale_xmin, &scale_xmax, &scale_ymin, &scale_ymax);
	} else if (strncmp(line, "##Maintitle", 11) == 0) {
	    strcpy(maintitle, &line[11]);
	} else if (strncmp(line, "##Subtitle", 10) == 0) {
	    strcpy(subtitle, &line[10]);
	} else if (strncmp(line, "##Xtitle", 8) == 0) {
	    strcpy(xtitle, &line[8]);
	} else if (strncmp(line, "##Ytitle", 8) == 0) {
	    strcpy(ytitle, &line[8]);
	} else if (strncmp(line, "##DataObjects", 13) == 0) {
	  fscanf(fp, "##Number %d\n", &data->n_data_obj);
	  if (NULL == (d_obj = (e_obj *)xmalloc((data->n_data_obj+1) * 
						sizeof(e_obj))))
	      return -1;
	  
	  for (i = 0; i < data->n_data_obj; i++) {
	      /* have to deal with Filled Rectangle */
	      fgets(line, 1024, fp);
	      if (strncmp(line, "Filled Rectangle", 16) == 0) {
		  sscanf(line, "%s %s x1 %lf y1 %lf x2 %lf y2%lf colour %d\n", 
			 dummy, dummy, &d_obj[i].pos.x0, &d_obj[i].pos.y0, 
			 &d_obj[i].pos.x1, &d_obj[i].pos.y1, &colour);
		  d_obj[i].type = G_RECTANGLEFILL;
	      } else {
		  sscanf(line, "%s x1 %lf y1 %lf x2 %lf y2%lf colour %d\n", 
			 obj_type, &d_obj[i].pos.x0, &d_obj[i].pos.y0, 
			 &d_obj[i].pos.x1, &d_obj[i].pos.y1, &colour);
	      }
	      strcpy(d_obj[i].colour, colournum[colour]);

	      if (strcmp(obj_type, "Line") == 0) {
		  d_obj[i].type = G_LINES;
	      } else if (strcmp(obj_type, "Rectangle") == 0) {
		  d_obj[i].type = G_RECTANGLE;
	      } else if (strcmp(obj_type, "Text") == 0) {
		  d_obj[i].type = G_TEXT;
	      }
#ifdef DEBUG
		printf("i %d type %d %f %f %f %f %s\n", i, d_obj[i].type,
		       d_obj[i].pos.x0, d_obj[i].pos.y0, d_obj[i].pos.x1, 
		       d_obj[i].pos.y1, d_obj[i].colour);
#endif
	  }
	} else if (strncmp(line, "##GraphObjects", 14) == 0) {
	    fscanf(fp, "##Number %d\n", &data->n_graph_obj);
	    if (NULL == (g_obj = (e_obj *)xmalloc((data->n_graph_obj+1) * 
						  sizeof(e_obj))))
		return -1;
	    
	    for (i = 0; i < data->n_graph_obj; i++) {
		/* have to deal with Filled Rectangle */
		fgets(line, 1024, fp);

		if (strncmp(line, "Filled Rectangle", 16) == 0) {
		    sscanf(line, "%s %s x1 %lf y1 %lf x2 %lf y2 %lf colour %d\n", 
			   dummy, dummy, &g_obj[i].pos.x0, &g_obj[i].pos.y0, 
			   &g_obj[i].pos.x1, &g_obj[i].pos.y1, &colour);
		    g_obj[i].type = G_RECTANGLEFILL;
		} else {
		    sscanf(line, "%s x1 %lf y1 %lf x2 %lf y2 %lf colour %d\n", 
		       obj_type, &g_obj[i].pos.x0, &g_obj[i].pos.y0, 
		       &g_obj[i].pos.x1, &g_obj[i].pos.y1, &colour);
		}
		strcpy(g_obj[i].colour, colournum[colour]);

		if (strncmp(obj_type, "Line", 4) == 0) {
		    g_obj[i].type = G_LINES;
		} else if (strncmp(obj_type, "Rectangle", 9) == 0) {
		    g_obj[i].type = G_RECTANGLE;
		} else if (strncmp(obj_type, "Text", 4) == 0) {
		    g_obj[i].type = G_TEXT;
		}
		
#ifdef DEBUG
		printf("i %d type %d %f %f %f %f %s\n", i, g_obj[i].type,
		       g_obj[i].pos.x0, g_obj[i].pos.y0, g_obj[i].pos.x1, 
		       g_obj[i].pos.y1, g_obj[i].colour);
#endif
	    }
	} else if (strncmp(line, "##", 2) != 0) {
	    if (NULL == (p_array = (p_score *)xmalloc((data->n_pts+1) * 
						      sizeof(p_score))))
		return -1;
	    
	    
	    /* do first point */
	    if (0 == sscanf(line, "%lf %lf\n", &pos, &p_array[0].score)) {
		return -1;
	    }
	    p_array[0].pos = (int)pos;
	    for (i = 1; i < data->n_pts; i++) {
		fscanf(fp, "%lf %lf\n", &pos, &p_array[i].score);
		p_array[i].pos = (int)pos;
#ifdef DEBUG
		printf("i %d pos %d score %f\n", i, p_array[i].pos, 
		       p_array[i].score);
#endif
	    }
	}
    }

    data->p_array = p_array;
    data->d_obj = d_obj;
    data->g_obj = g_obj;

    data->dim.x0 = scale_xmin;
    data->dim.x1 = scale_xmax;
    data->dim.y0 = scale_ymin;
    data->dim.y1 = scale_ymax;

    /* convert data structure into Graph structure */
    if (data->n_pts > 0) {
	if (NULL == (graph->p_arrays = (parray *)xmalloc(sizeof(parray))))
	    return -1;

	graph->n_parrays = 1;
	graph->p_arrays[0].n_pts = data->n_pts;
	if (NULL == (graph->p_arrays[0].p_array = 
		     (g_pt *)xmalloc(sizeof(g_pt) * 
				     graph->p_arrays[0].n_pts)))
	    return -1;
	    
	graph->p_arrays[0].type = G_LINE;
	for (i = 0; i < data->n_pts; i++) {
	    graph->p_arrays[0].p_array[i].x = p_array[i].pos;
	    graph->p_arrays[0].p_array[i].y = p_array[i].score;
	}
    }

    if (data->n_graph_obj > 0 || data->n_data_obj > 0) {
	if (NULL==(graph->d_arrays=(darray *)xmalloc((data->n_data_obj+
						      data->n_graph_obj + 1) * 
						     sizeof(darray))))
	    return -1;
	
	graph->n_darrays = data->n_graph_obj + data->n_data_obj;

	for (i = 0; i < data->n_data_obj; i++) {
	    if (NULL == (graph->d_arrays[i].d_array = 
			 (gd_line *)xmalloc(sizeof(gd_line))))
		return -1;
	    graph->d_arrays[i].n_dlines = 1;
	    graph->d_arrays[i].d_array[0].x0 = d_obj[i].pos.x0;
	    graph->d_arrays[i].d_array[0].x1 = d_obj[i].pos.x1;
	    graph->d_arrays[i].d_array[0].y0 = d_obj[i].pos.y0;
	    graph->d_arrays[i].d_array[0].y1 = d_obj[i].pos.y1;
	    graph->d_arrays[i].type = d_obj[i].type;
	    
	} 
	for (j = 0; j < data->n_graph_obj; j++, i++) {
	    if (NULL == (graph->d_arrays[i].d_array = 
			 (gd_line *)xmalloc(sizeof(gd_line))))
		return -1;
	    graph->d_arrays[i].n_dlines = 1;
	    graph->d_arrays[i].d_array[0].x0 = g_obj[j].pos.x0;
	    graph->d_arrays[i].d_array[0].x1 = g_obj[j].pos.x1;
	    graph->d_arrays[i].d_array[0].y0 = g_obj[j].pos.y0;
	    graph->d_arrays[i].d_array[0].y1 = g_obj[j].pos.y1;
	    graph->d_arrays[i].type = g_obj[j].type;
	}
    }

    graph->dim.x0 = scale_xmin;
    graph->dim.x1 = scale_xmax;
    graph->dim.y0 = scale_ymin;
    graph->dim.y1 = scale_ymax;
    *graph_obj = Tcl_NewGraphObj(graph);

    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    if (data->n_pts > 0) {
	xfree(graph->p_arrays[0].p_array);
	xfree(graph->p_arrays);
    }
    if (data->n_graph_obj > 0 || data->n_data_obj > 0) {
	for (i = 0; i < data->n_data_obj; i++) {
	    xfree(graph->d_arrays[i].d_array);
	}
	xfree(graph->d_arrays);
    }
    xfree(graph);

    if (strcmp(title, " \n") != 0) { 
	data->title = strdup(title);
	(*text_data)->title = strdup(title);
    }
    if (strcmp(maintitle, " \n") != 0) {
	data->maintitle = strdup(maintitle);
	(*text_data)->maintitle = strdup(maintitle);
    }
    if (strcmp(subtitle, " \n") != 0) {
	data->subtitle = strdup(subtitle);
	(*text_data)->subtitle = strdup(subtitle);
    }
    if (strcmp(xtitle, " \n") != 0) {
	data->xtitle = strdup(xtitle);
	(*text_data)->xtitle = strdup(xtitle);
    }
    if (strcmp(ytitle, " \n") != 0) {
	data->ytitle = strdup(ytitle);
	(*text_data)->ytitle = strdup(ytitle);
    }

#ifdef DEBUG
    printf("plot_type %s\n", plot_type);
    printf("title %s\n", title);
    printf("graphs %d %d %d\n", num_graphs, graph_num, data->n_pts);
    printf("xmina %f %f %f %f\n", xmina, xmaxa, ymina, ymaxa);
    printf("xmin %f %f %f %f\n", xmin, xmax, ymin, ymax);
    printf("scale %f %f %f %f\n", scale_xmin, scale_xmax, scale_ymin, scale_ymax);
    if (num_graphs > 1) { 
	printf("%s\n %s\n %s\n %s\n", maintitle, subtitle, xtitle, ytitle);
    } else {
	printf("%s\n %s\n %s\n", maintitle, xtitle, ytitle);
    }
#endif

    fclose (fp);
    *out_data = data;
    return 0;
}

int init_emboss_graph_create(Tcl_Interp *interp, 
			     int seq_id, 
			     int start, 
			     int end,
			     char *filename,
			     Tcl_Obj **graph_obj,
			     int *id)
{
    int seq_num, seq_len;
    in_emboss *input;
    Tcl_DString input_params;
    text_emboss *text_data;

    e_graph *data = NULL;

    seq_num = GetSeqNum(seq_id);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }
    seq_len = end - start + 1;

    if (NULL == (input = (in_emboss *)xmalloc (sizeof(in_emboss))))
	return -1;

    read_emboss_data_file(filename, &data, graph_obj, &text_data);
    if (!data) {
	verror(ERR_FATAL,"emboss", "error in reading results\n");
	return -1;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n",
		       GetSeqName(seq_num), start, end);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_emboss_graph(seq_num, start, end, data, input, 
					text_data, graph_obj))) {
	verror(ERR_FATAL,"emboss", "error in saving results\n");
	return -1;
    }
    xfree(data);
    return 0;
}

int init_emboss_stick_create(Tcl_Interp *interp, 
			     int seq_id, 
			     int start, 
			     int end,
			     char *filename,
			     int *id)
{
    FILE *fp;
    int i;
    int *pos = NULL;
    int *res = NULL;
    int seq_num, seq_len;
    in_emboss *input;
    Tcl_DString input_params;
    char dummy1[100];
    char dummy2;

    seq_num = GetSeqNum(seq_id);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }
    seq_len = end - start + 1;

    if (NULL == (pos = (int *)xmalloc((seq_len+1) * sizeof(int))))
	return -1;
    if (NULL == (res = (int *)xmalloc((seq_len+1) * sizeof(res))))
	return -1;
    if (NULL == (input = (in_emboss *)xmalloc (sizeof(in_emboss))))
	return -1;

    if (NULL == (fp = fopen(filename, "r"))) {
	printf("unable to open file\n");
	return -1;
    }

    if (fgetc(fp) == 'P') {
	printf("first char\n");
	fgets(dummy1, 100, fp);
	fgets(dummy1, 100, fp);
	fgets(dummy1, 100, fp);
    } else {
	rewind(fp);
    }

    i = 0;
    while (fscanf(fp, "%d %c %d\n", &pos[i], &dummy2, &res[i]) != EOF) { 
#ifdef DEBUG
	printf("pos %d result %d\n", pos[i], res[i]);
#endif
	i++;
    }
    fclose (fp);

    printf("num points %d\n", i);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n",
		       GetSeqName(seq_num), start, end);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

#ifdef TODO
    if (-1 == (*id = store_emboss_graph(seq_num, start, end, pos, score, 
					i, input))) {
	verror(ERR_FATAL,"emboss", "error in saving results\n");
	return -1;
    }
#endif
    if (pos)
	xfree(pos);
    if (res)
	xfree(res);

    return 0;
}

int init_emboss_graph_plot(Tcl_Interp *interp,
			   int seq_id_h,
			   int seq_id_v,
			   int result_id,
			   char *e_win, 
			   int element_id,
			   char *c_win, 
			   int container_id,
			   Tcl_Obj *results,
			   int line_width,
			   char *colour,
			   char *element_type,
			   char *name)
{
    seq_result *s_result;
    Graph *graph;
    plot_data *result;
    configs *configure;

    s_result = seq_id_to_result(result_id);

    if (NULL == (result = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    if (NULL == (result->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;
    
    if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;

    configure->position = -1.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    result->configure[0] = configure;
    result->n_configure = 1;
    result->sf_m = 1.0;
    result->sf_c = 0.0;
    result->result_id = result_id;
    result->scale = SCALE_X | SCALE_Y;
    result->hidden = 0;
    result->line_width = line_width;
    result->colour = strdup(colour);
    result->len_ruler = 1;
    result->amp_ruler = 0;
    result->name = strdup(name);
    sprintf(result->tags, "id%d", result_id);

    graph = Tcl_GetGraphFromObj(results);

    init_seq_element(interp, s_result, seq_id_h, seq_id_v, result, 
		     container_id, element_id, e_win, c_win, 
		     HORIZONTAL, HORIZONTAL|VERTICAL, graph, 
		     element_type);
    return 0;
}

int store_emboss_dot(int seq_num_h, 
		     int start_h, 
		     int end_h,
		     int seq_num_v, 
		     int start_v, 
		     int end_v,
		     e_graph *data,
		     in_emboss *input,
		     text_emboss *text_data,
		     Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int id;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    id = get_reg_id();

    /* data assignment */
    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num_h);
    s_result->seq_id[VERTICAL] = GetSeqId(seq_num_v);

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->type = SEQ_TYPE_DOT_PLOT;
    s_result->gr_type = SEQ_TYPE_DOT_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj;
    s_result->text_data = text_data;
    s_result->e = NULL;

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = emboss_graph_callback;
    s_result->txt_func = emboss_graph_text_func;

    seq_register(seq_num_h, emboss_graph_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    seq_register(seq_num_v, emboss_graph_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
     return id;
}

int init_emboss_dot_create(Tcl_Interp *interp, 
			   int seq_id_h, 
			   int start_h, 
			   int end_h,
			   int seq_id_v, 
			   int start_v, 
			   int end_v,
			   char *filename,
			   Tcl_Obj **graph_obj,
			   int *id)
{
    int seq_num_h, seq_len_h, seq_num_v, seq_len_v;
    in_emboss *input;
    Tcl_DString input_params;
    int *seq1_match = NULL;
    int *seq2_match = NULL; 
    int *len_match = NULL;
    e_graph *data = NULL;
    text_emboss *text_data;

    seq_num_h = GetSeqNum(seq_id_h);
    seq_num_v = GetSeqNum(seq_id_v);

    seq_len_h = GetSeqLength(seq_num_h);
    seq_len_v = GetSeqLength(seq_num_v);

    /* if the end has not been defined, set it to be the sequence length */
   if (end_h == -1)
	end_h = seq_len_h;

    if (end_v == -1)
	end_v = seq_len_v;

    seq_len_h = end_h - start_h + 1;
    seq_len_v = end_v - start_v + 1;

    read_emboss_data_file(filename, &data, graph_obj, &text_data);
    if (!data) {
	verror(ERR_FATAL,"emboss", "error in reading results\n");
	return -1;
    }

    if (NULL == (seq1_match = (int *)xmalloc((data->n_data_obj+1) * sizeof(int))))
	return -1;
    if (NULL == (seq2_match = (int *)xmalloc((data->n_data_obj+1) * sizeof(int))))
	return -1;
    if (NULL == (len_match = (int *)xmalloc((data->n_data_obj+1) * sizeof(int))))
	return -1;

    if (NULL == (input = (in_emboss *)xmalloc (sizeof(in_emboss))))
	return -1;

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\nsequence %s: from %d to %d\n",
		       GetSeqName(seq_num_h), start_h, end_h,
		       GetSeqName(seq_num_v), start_v, end_v);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_emboss_dot(seq_num_h, start_h, end_h,
				      seq_num_v, start_v, end_v, 
				      data, input, text_data, graph_obj))) {
	verror(ERR_FATAL,"emboss", "error in saving results\n");
	return -1;
    }

    xfree(seq1_match);
    xfree(seq2_match);
    xfree(len_match);
    xfree(data);
    return 0;
}

int init_emboss_dot_plot(Tcl_Interp *interp,
			 int seq_id_h,
			 int seq_id_v,
			 int result_id,
			 char *e_win, 
			 int element_id,
			 char *c_win, 
			 int container_id,
			 Tcl_Obj *results,
			 int line_width,
			 char *colour,
			 char *element_type,
			 char *name)
{
    seq_result *s_result;
    Graph *graph;
    plot_data *result;
    configs *configure;

    s_result = seq_id_to_result(result_id);

    if (NULL == (result = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    if (NULL == (result->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;
    
    if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;

    configure->position = -1.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    result->configure[0] = configure;
    result->n_configure = 1;
    result->sf_m = 1.0;
    result->sf_c = 0.0;
    result->result_id = result_id;
    result->scale = SCALE_X | SCALE_Y;
    result->hidden = 0;
    result->line_width = line_width;
    result->colour = strdup(colour);
    result->len_ruler = 1;
    result->amp_ruler = 0;
    result->name = strdup(name);
    sprintf(result->tags, "id%d", result_id);

    graph = Tcl_GetGraphFromObj(results);

    init_seq_element(interp, s_result, seq_id_h, seq_id_v, result, 
		     container_id, element_id, e_win, c_win, 
		     HORIZONTAL|VERTICAL, HORIZONTAL|VERTICAL, graph, 
		     element_type);
    return 0;
}
