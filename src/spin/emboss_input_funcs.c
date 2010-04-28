#include <staden_config.h>

#include <tcl.h>
#include <unistd.h>
#include <string.h>
#include <float.h>

#include "xalloc.h"
#include "text_output.h"
#include "seq_results.h"
#include "seq_plot_funcs.h"
#include "seq_raster.h"
#include "emboss_input_funcs.h"
#include "tcl_utils.h"
#include "spin_globals.h"
#include "sequence_pair_display.h"

/* emboss colours taken from ajax/ajgraph.c */
static char *colournum[] = { "BLACK", "RED", "YELLOW", "GREEN", "AQUAMARINE",
			     "PINK", "WHEAT", "GREY", "BROWN", "BLUE", 
			     "BLUEVIOLET", "CYAN", "TURQUOISE", "MAGENTA", 
			     "SALMON", "WHITE"};

void emboss_graph_callback(int seq_num, void *obj, seq_reg_data *jdata);

void emboss_graph_shutdown(Tcl_Interp *interp,
			   seq_result *result,
			   char *raster_win,
			   int seq_num)
{
    in_emboss *input = result->input;
    out_raster *output = result->output;
    char *tmp;
    seq_reg_key_name info;
    static char buf[80];
    int raster_id;
    double wx0, wy0, wx1, wy1;
    Tk_Raster *raster;
    Tcl_CmdInfo info1;
    RasterResult *raster_result;
    e_graph *g_data;

    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    seq_deregister(GetSeqNum(result->seq_id[HORIZONTAL]), 
		   emboss_graph_callback, (seq_result *)result);

    if (result->seq_id[VERTICAL] != -1) {
	seq_deregister(GetSeqNum(result->seq_id[VERTICAL]), 
		       emboss_graph_callback, (seq_result *)result);
    }

    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
	ReplotAllCurrentZoom(interp, raster_win);
	tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
    
	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", 
				  tmp, NULL)){
	    verror(ERR_WARN, "emboss", "graph shutdown %s \n", 
		   Tcl_GetStringResult(interp));
	}

	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "emboss", "graph remove %s \n", Tcl_GetStringResult(interp));
	
	Tcl_GetCommandInfo(interp, raster_win, &info1);
	raster = (Tk_Raster*)info1.clientData;

	/* find original y before reset size */
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

	/* update the scale of the remaining results in the raster window */
	SeqReSetRasterWindowSize(interp, raster_win, result->graph);
	ReSetRasterWindowWorld(interp, raster_win, wy1, result->graph);
	ReplotAllRasterWindow(interp, raster_win);
    }
    if (result->graph == SEQ_E_DOT) {
	DestroySequencePairDisplay(interp, result->id);
    }

    g_data = result->data;

    if (g_data->p_array) 
	xfree(g_data->p_array);
    if (g_data->d_obj)
	xfree(g_data->d_obj);
    if (g_data->g_obj)
	xfree(g_data->g_obj);

    if (g_data->title)
	free(g_data->title);
    if (g_data->maintitle)
	free(g_data->maintitle);
    if (g_data->subtitle)
	free(g_data->subtitle);
    if (g_data->xtitle)
	free(g_data->xtitle);
    if (g_data->ytitle)
	free(g_data->ytitle);


    if (output->configure) {
	xfree(output->configure[0]);
	if (output->n_configure == 2) 
	    xfree(output->configure[1]);
	xfree(output->configure);
    }
    free(output->name);

    xfree(result->data);

    free(input->params);
    xfree(result->input);
    xfree(result->output);
    xfree(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);

}

void emboss_graph_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_emboss *input = result->input;
    out_raster *output = result->output;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Emboss graph plot");
	break;    

    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "%s #%d", output->name, result->id);
	break;

    case SEQ_GET_BRIEF:
	{
	    e_graph *g_data;
	    g_data = result->data;
	
	    /*
	    if (g_data->maintitle) {
		sprintf(jdata->name.line, "%s", g_data->maintitle);
	    } else if (g_data->title){
		sprintf(jdata->name.line, "%s", g_data->title);
	    }
	    */
	    jdata->name.line[0] = '\0';
	    if (g_data->title)
		strcat(jdata->name.line, g_data->title);

	    if (g_data->maintitle)
		strcat(jdata->name.line, g_data->maintitle);

#ifdef REMOVE
	    if (result->seq_id[VERTICAL] == -1) {
		sprintf(jdata->name.line, "%s: seq=%s", 
			output->name, GetSeqName(GetSeqNum(result->seq_id[HORIZONTAL])));
	    } else {
		sprintf(jdata->name.line, "%s: seq_h=%s seq_v=%s", 
		    output->name, 
			GetSeqName(GetSeqNum(result->seq_id[HORIZONTAL])), 
			GetSeqName(GetSeqNum(result->seq_id[VERTICAL])));
	    }
#endif
	}
	break;
	
    case SEQ_GET_OPS:
	if (output->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0"
		"PLACEHOLDER\0PLACEHOLDER\0PLACEHOLDER\0Reveal\0SEPARATOR\0Remove\0";
	} else {
	    if (result->graph == SEQ_E_DOT) {
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
	    Tcl_Eval(output->interp, "SetBusy");
	    vfuncheader("results");
	    result->txt_func(result);
	    Tcl_Eval(output->interp, "ClearBusy");
	    break;
	case 2: /* configure */
	    sprintf(cmd, "RasterConfig %d", id);
	    if (TCL_OK != Tcl_Eval(output->interp, cmd)){
		puts(Tcl_GetStringResult(output->interp));
	    }
	    break;
	case 3: /* display sequences */
	    SequencePairDisplay(output->interp, output->raster_win, id, 
				result->seq_id[HORIZONTAL], 
				result->seq_id[VERTICAL]);
	    break;
	case 4: /* hide all */
	    output->hidden = 1;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 5: /* reveal all */
	    output->hidden = 0;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 6: /* remove */ 
	    {
		Tcl_Interp *interp = output->interp;
		emboss_graph_shutdown(interp, result, output->raster_win, 
				      seq_num);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	result->pr_func(result, (seq_reg_plot *)jdata);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case OUTPUT:
	    jdata->info.result = (void *)output;
	    break;
	case INPUT: 
	    jdata->info.result = (void *)input;
	    break;
	case DIMENSIONS: 
	    {
		e_graph *g_data;

		g_data = result->data;
		jdata->info.result = (void *)&g_data->dim;

		break;
	    }
	case INDEX: 
	    jdata->info.result = (void *)id;
	    break;
	case RESULT:
	    jdata->info.result = (void *)result;
	    break;
	case WIN_NAME:
	    {
		char *r_win = output->raster_win;
		jdata->info.result = (void *)r_win;
		break;
	    }
	case WIN_SIZE:
	    {
		d_point *pt;
		Tcl_Interp *interp = output->interp;
		if (NULL == (pt = (d_point *)xmalloc(sizeof(d_point))))
		    return;

		if (result->graph == SEQ_GRAPH) {
		    pt->x = get_default_int(interp, spin_defs, 
					    w("EMBOSS.RASTER.GRAPH.PLOT_WIDTH"));
		    pt->y = get_default_double(interp, spin_defs,
					       w("EMBOSS.RASTER.GRAPH.PLOT_HEIGHT"));
		} else if (result->graph == SEQ_E_DOT) {
		    pt->x = get_default_int(interp, spin_defs, 
					    w("EMBOSS.RASTER.DOT.PLOT_WIDTH"));
		    pt->y = get_default_double(interp, spin_defs,
					       w("EMBOSS.RASTER.DOT.PLOT_HEIGHT"));

		}
		jdata->info.result = (void *)pt;
		break;
	    } /* WIN_SIZE */
	} /* switch SEQ_RESULT_INFO */
	break;
    case SEQ_HIDE: 
	output->hidden = 1;
	break;
    case SEQ_REVEAL: 
	output->hidden = 0;
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    Tcl_Interp *interp = output->interp;
	    emboss_graph_shutdown(interp, result, output->raster_win,
				  seq_num);
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
		       in_emboss *input)
{
    seq_result *result;
    int id;

#ifdef DEBUG
    for (i = 0; i < data->n_pts; i++) {
	printf("i %d pos %d score %f\n", i, data->p_array[i].pos, 
	       data->p_array[i].score);
    }
#endif

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    result->data = (e_graph *)data;

    id = get_reg_id();

    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    result->id = id; 
    result->input = (void *)input; 
    result->output = NULL;
    result->type = SEQ_TYPE_GRAPH_PLOT;
    result->frame = 0;
    result->graph = SEQ_GRAPH;

    result->pr_func = emboss_graph_plot_func;
    result->op_func = emboss_graph_callback;
    result->txt_func = emboss_graph_text_func;

    seq_register(seq_num, emboss_graph_callback, (void *)result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int read_emboss_data_file(char *filename,
			  e_graph **out_data)
{
    FILE *fp;
    int i;
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
    
    if (NULL == (fp = fopen(filename, "r"))) {
	verror(ERR_WARN, "emboss", "unable to open file %s \n", filename);
	return -1;
    }

    if (NULL == (data = (e_graph *)xmalloc(sizeof(e_graph))))
	return -1;
    data->n_pts = 0;
    data->n_data_obj = 0;
    data->n_graph_obj = 0;
    data->title = NULL;
    data->maintitle = NULL;
    data->subtitle = NULL;
    data->xtitle = NULL;
    data->ytitle = NULL;

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
		  d_obj[i].type = E_RECTANGLEFILL;
	      } else {
		  sscanf(line, "%s x1 %lf y1 %lf x2 %lf y2%lf colour %d\n", 
			 obj_type, &d_obj[i].pos.x0, &d_obj[i].pos.y0, 
			 &d_obj[i].pos.x1, &d_obj[i].pos.y1, &colour);
	      }
	      strcpy(d_obj[i].colour, colournum[colour]);

	      if (strcmp(obj_type, "Line") == 0) 
		  d_obj[i].type = E_LINE;
	      else if (strcmp(obj_type, "Rectangle") == 0) 
		  d_obj[i].type = E_RECTANGLE;
	      else if (strcmp(obj_type, "Text") == 0) 
		  d_obj[i].type = E_TEXT;
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
		    g_obj[i].type = E_RECTANGLEFILL;
		} else {
		    sscanf(line, "%s x1 %lf y1 %lf x2 %lf y2 %lf colour %d\n", 
		       obj_type, &g_obj[i].pos.x0, &g_obj[i].pos.y0, 
		       &g_obj[i].pos.x1, &g_obj[i].pos.y1, &colour);
		}
		strcpy(g_obj[i].colour, colournum[colour]);

		if (strncmp(obj_type, "Line", 4) == 0) {
		    g_obj[i].type = E_LINE;
		} else if (strncmp(obj_type, "Rectangle", 9) == 0) {
		    g_obj[i].type = E_RECTANGLE;
		} else if (strncmp(obj_type, "Text", 4) == 0) {
		    g_obj[i].type = E_TEXT;
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

    if (strcmp(title, " \n") != 0) 
	data->title = strdup(title);
    
    if (strcmp(maintitle, " \n") != 0) 
	data->maintitle = strdup(maintitle);

    if (strcmp(subtitle, " \n") != 0) 
	data->subtitle = strdup(subtitle);

    if (strcmp(xtitle, " \n") != 0) 
	data->xtitle = strdup(xtitle);

    if (strcmp(ytitle, " \n") != 0) 
	data->ytitle = strdup(ytitle);

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

    *out_data = data;
    fclose (fp);
    return 0;
}

int init_emboss_graph_create(Tcl_Interp *interp, 
			     int seq_id, 
			     int start, 
			     int end,
			     char *filename,
			     int *id)
{
    int seq_num, seq_len;
    in_emboss *input;
    Tcl_DString input_params;

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

    read_emboss_data_file(filename, &data);
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

    if (-1 == (*id = store_emboss_graph(seq_num, start, end, data, input))) {
	verror(ERR_FATAL,"emboss", "error in saving results\n");
	return -1;
    }
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
			   int seq_id, 
			   int result_id,
			   char *name,
			   char *raster_win,
			   int raster_id,
			   char *colour, 
			   int line_width)
{
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    out_raster *output; 
    char *opts[5];
    int seq_num;
    seq_cursor_notify cn;
    RasterResult *raster_result;
    cursor_t *cursor;
    int superimpose = 1;
    seq_result *result;
    config *configure;
    in_emboss *input;
    e_graph *data;

    if (NULL == (output = (out_raster *)xmalloc(sizeof(out_raster))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    result = result_data(result_id, seq_num);
    input = result->input;
    result->output = output;
    data = result->data;

    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	return -1;

    raster = (Tk_Raster*)info.clientData;

    RasterInitPlotFunc(raster, SeqRasterPlotFunc);
    
    /* need to check if superimposing result on another plot */
    raster_result = raster_id_to_result(raster_id);
    if (raster_result->num_results == 0) {
	superimpose = 0;
    }

    if (NULL == (opts[1] = (char *)xmalloc((strlen(colour)+1) * 
					   sizeof(char))))
	return -1;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;

    output->name = strdup(name);
    strcpy(output->raster_win, raster_win);
    output->raster_id = raster_id;
    output->interp = interp;
    output->hidden = 0;

    opts[0] = "-fg";
    strcpy(opts[1], colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = NULL;

    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);

    if (NULL == (configure = (config *)xmalloc(sizeof(config))))
	return -1;

    if (NULL == (output->configure = (config **)xmalloc(sizeof(config*))))
	return -1;

    configure->position = 0.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    output->configure[0] = configure;
    output->n_configure = 1;
    output->scroll = 'b';
    output->sf_m = 1.0;
    output->sf_c = 0.0;

    if (superimpose) {
	SeqSuperimposeResult(interp, output->raster_win, result_id, 
			     data->dim.x0, data->dim.y0, 
			     data->dim.x1, data->dim.y1);
    } else {
      /* 
       * update world scroll region for this raster, if necessary to be the 
       * largest range 
       */
	RasterSetWorldScroll(raster, data->dim.x0, data->dim.y0, 
			     data->dim.x1, data->dim.y1);
	SeqAddRasterToWindow(interp, raster_win, result->graph);
    }

    raster_result = raster_id_to_result(raster_id);
    cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);

    /* move cursor to start position if no cursor is yet present */
    if (raster_result->cursor_array[cursor->id].prev_pos == -1 
	&& data->dim.x0 > raster_result->cursor_array[cursor->id].prev_pos) {
	cursor->abspos = data->dim.x0;
    }

    AddResultToRaster(raster_result);
    ReplotAllCurrentZoom(interp, output->raster_win);
    xfree(opts[1]);
    xfree(opts[3]);

    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(seq_num, (seq_reg_data *)&cn); 

    return 0;

}

int store_emboss_dot(int seq_num_h, 
		     int start_h, 
		     int end_h,
		     int seq_num_v, 
		     int start_v, 
		     int end_v,
		     e_graph *data,
		     in_emboss *input)
{
    seq_result *result;
    int id;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    result->data = data;

    id = get_reg_id();

    /* data assignment */
    result->seq_id[HORIZONTAL] = GetSeqId(seq_num_h);
    result->seq_id[VERTICAL] = GetSeqId(seq_num_v);

    result->id = id; 
    result->input = (void *)input; 
    result->output = NULL;
    result->type = SEQ_TYPE_DOT_PLOT;
    result->frame = 0;
    result->graph = SEQ_E_DOT;

    result->pr_func = emboss_graph_plot_func;
    result->op_func = emboss_graph_callback;
    result->txt_func = emboss_graph_text_func;

    seq_register(seq_num_h, emboss_graph_callback, (void *)result, 
		 SEQ_PLOT_PERM, id);
    seq_register(seq_num_v, emboss_graph_callback, (void *)result, 
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
			   int *id)
{
    int seq_num_h, seq_len_h, seq_num_v, seq_len_v;
    in_emboss *input;
    Tcl_DString input_params;
    int *seq1_match = NULL;
    int *seq2_match = NULL; 
    int *len_match = NULL;
    e_graph *data = NULL;

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

    read_emboss_data_file(filename, &data);
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
				      data, input))) {
	verror(ERR_FATAL,"emboss", "error in saving results\n");
	return -1;
    }

    xfree(seq1_match);
    xfree(seq2_match);
    xfree(len_match);
    return 0;
}

int init_emboss_dot_plot(Tcl_Interp *interp, 
			 int seq_id_h,
			 int seq_id_v,
			 int result_id,
			 char *name,
			 char *raster_win,
			 int raster_id,
			 char *colour, 
			 int line_width)
{
     char *opts[7];
     e_graph *data;
     seq_result *result;

     if (NULL == (opts[1] = (char *)xmalloc((strlen(colour)+1) * 
					    sizeof(char))))
	 return -1;
     if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	 return -1;
     if (NULL == (opts[5] = (char *)xmalloc(15 * sizeof(char))))
	 return -1;
     
     opts[0] = "-fg";
     strcpy(opts[1], colour);
     opts[2] = "-linewidth";
     sprintf(opts[3], "%d", line_width);
     opts[4] = "-capstyle";
     strcpy(opts[5], "round");
     opts[6] = NULL;
     
     result = result_data(result_id, GetSeqNum(seq_id_h));
     data = result->data;    
     
     init_dot_plot(interp, seq_id_h, seq_id_v, result_id, name, raster_win,
		   raster_id, opts, 6, LINE, data->dim);
     
     xfree(opts[1]);
     xfree(opts[3]);
     xfree(opts[5]);
     
     return 0;
}
