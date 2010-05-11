#include <assert.h>
#include <string.h>
#include <os.h>
#include <io_lib/Read.h>
#include <io_lib/traceType.h>
#include <io_lib/expFileIO.h>

#include "misc.h"
#include "tkTrace.h"
#include "xalloc.h"
#include "postscript.h"
#include "trace_print.h"

#define PYRO_SCALE 1000

static int trace_load_private(DNATrace *t);
void trace_pyroalign(Read *r);

void trace_unload(DNATrace *t) {
    if (t->read)
	read_deallocate(t->read);

    if (t->tracePos)
	xfree(t->tracePos);

    if (t->tracePosE)
	xfree(t->tracePosE);

    if (t->edBases)
	xfree(t->edBases);

    if (t->edPos)
	xfree(t->edPos);

    if (t->edConf)
	xfree(t->edConf);

    free_ps_options(&(t->ps_options));
    free_ps_trace(&(t->ps_trace));

    t->read = NULL;
    t->tracePos = NULL;
    t->tracePosE = NULL;
    t->edBases = NULL;
    t->edPos = NULL;
    t->edConf = NULL;
}

/*
 * Initialises the tracePos and tracePosE arrays.
 */
void trace_init_pos(DNATrace *t) {
    register int i,j,k;

    if (t->read->NPoints) {
	/* A real trace */

	/* For tracePos - the mapping to original bases */
	for (j = i = 0; i < t->read->NBases; i++) {
	    k = t->read->basePos[i];
	    if (k >= t->read->NPoints)
		k = t->read->NPoints-1;
	    for (; j <= k; j++)
		t->tracePos[j] = i;
	}

	for (; j < t->read->NPoints; j++)
	    t->tracePos[j] = i-1;

	/* For tracePosE - the mapping to edited bases */
	if (t->comp) {
	    for (j = i = 0; i < t->Ned; i++) {
		while (!(k = t->edPos[i]))
		    if (++i == t->Ned) {
			k = t->read->NPoints;
			break;
		    }
		k = t->Ned-k;
		if (k >= t->read->NBases) {
		    printf("Reading past end of array. Ned=%d bases=%d\n",
			   t->Ned, t->read->NBases);
		    k = t->read->NBases - 1;
		}
		if (k < 0) {
		    printf("Reading past start of array\n");
		    k = 0;
		}
		k = t->read->basePos[k];
		if (k >= t->read->NPoints)
		    k = t->read->NPoints-1;
		for (; j <= k; j++)
		    t->tracePosE[j] = i;
	    }

	    for (; j < t->read->NPoints; j++)
		t->tracePosE[j] = i-1;
	} else {
	    for (j = i = 0; i < t->Ned; i++) {
		while (!(k = t->edPos[i]))
		    if (++i == t->Ned) {
			k = t->read->NPoints;
			break;
		    }
		k--;
		if (k >= t->read->NBases) {
		    printf("Reading past end of array. Ned=%d bases=%d\n",
			   t->Ned, t->read->NBases);
		    k = t->read->NBases - 1;
		}
		if (k < 0) {
		    printf("Reading past start of array\n");
		    k = 0;
		}
		k = t->read->basePos[k];
		if (k >= t->read->NPoints)
		    k = t->read->NPoints-1;
		for (; j <= k; j++)
		    t->tracePosE[j] = i;
	    }

	    for (; j < t->read->NPoints; j++)
		t->tracePosE[j] = i-1;
	}

    } else {
	/* A text file - invent some values */
	for (i = 0; i < t->read->NBases; i++) {
	    t->read->basePos[i] = (i+1) * 8;
	}

	for (i = 0; i < t->read->NBases * 8; i++) {
	    t->tracePos[i] = t->tracePosE[i] = i / 8;
	}

	for (i = t->read->NBases * 8; i < (t->read->NBases+1) * 8; i++) {
	    t->tracePos[i] = t->read->NBases - 1;
	}

	t->read->NPoints = (t->read->NBases+1) * 8;
	t->read->maxTraceVal = 0;

	t->read->traceA = (TRACE *)xrealloc(t->read->traceA,
					    t->read->NPoints * sizeof(TRACE));
	t->read->traceC = (TRACE *)xrealloc(t->read->traceC,
					    t->read->NPoints * sizeof(TRACE));
	t->read->traceG = (TRACE *)xrealloc(t->read->traceG,
					    t->read->NPoints * sizeof(TRACE));
	t->read->traceT = (TRACE *)xrealloc(t->read->traceT,
					    t->read->NPoints * sizeof(TRACE));

	memset(t->read->traceA, 0, t->read->NPoints * sizeof(TRACE));
	memset(t->read->traceC, 0, t->read->NPoints * sizeof(TRACE));
	memset(t->read->traceG, 0, t->read->NPoints * sizeof(TRACE));
	memset(t->read->traceT, 0, t->read->NPoints * sizeof(TRACE));
    }
}

/*
 * Some solexa data has extreme outliers and so an entire tile may get encoded
 * with a DC offset which is inappropriate for each trace as a whole.
 * Here we re-encode using the minimum range needed.
 */
void trace_recalc_baseline(DNATrace *t) {
    Read *r = t->read;
    int i;
    int min_val = INT_MAX, max_val = 0;

    /* Identify the real dynamic range and offset */
    for (i = 0; i < r->NPoints; i++) {
	if (min_val > r->traceA[i]) min_val = r->traceA[i];
	if (min_val > r->traceC[i]) min_val = r->traceC[i];
	if (min_val > r->traceG[i]) min_val = r->traceG[i];
	if (min_val > r->traceT[i]) min_val = r->traceT[i];
	if (max_val < r->traceA[i]) max_val = r->traceA[i];
	if (max_val < r->traceC[i]) max_val = r->traceC[i];
	if (max_val < r->traceG[i]) max_val = r->traceG[i];
	if (max_val < r->traceT[i]) max_val = r->traceT[i];
    }

    /* Shift it down */
    for (i = 0; i < r->NPoints; i++) {
	r->traceA[i] -= min_val;
	r->traceC[i] -= min_val;
	r->traceG[i] -= min_val;
	r->traceT[i] -= min_val;
    }
    r->baseline    -= min_val;
    r->maxTraceVal -= min_val;
}

/*
 * Loads a Trace structure from a file
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int trace_load(DNATrace *t, char *file, char *format) {
    int form = trace_type_str2int(format);
    int tmp;

    /* Load the Read structure */
    if (t->read)
	trace_unload(t);

    tmp = read_experiment_redirect(2);
    if (NULLRead == (t->read = read_reading(file, form))) {
	read_experiment_redirect(tmp);
	return -1;
    }
    read_experiment_redirect(tmp);

    /* Auto-detect pyrosequencing by presence of flows */
    if (t->read->flow_order && t->read->flow && t->read->nflows) {
	t->style = STYLE_PYRO;
    }

    if (t->style == STYLE_PYRO) {
	trace_pyroalign(t->read);
	t->yticks = PYRO_SCALE;
    }


    /* Auto-detect solexa style sequences by 1 base per sample */
    if (t->read->NBases == t->read->NPoints && t->read->NBases != 0) {
	t->style = STYLE_STICK;
	trace_recalc_baseline(t);
    }

    return trace_load_private(t);
}

/*
 * Creates a Trace structure from a memory copy of a Read structure
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int trace_memory_load(DNATrace *t, Read *r) {
    t->read = r;

    return trace_load_private(t);
}

/*
 * Returns:
 *   0 for success
 *  -1 for failure
 */
static int trace_load_private(DNATrace *t) {
    register int i;
    Exp_info *e = NULL;

    /*
     * Allocate and initialise the edited base values.
     * These come directly from the Read structure. When dealing with an
     * experiment file we've requested that we use the original SCF data,
     * rather than those held in the experiment file. This information
     * can still be gleaned from the Read.orig_trace field.
     */
    t->MaxNed = t->read->NBases + TRACE_EDITS;
    t->Ned = t->read->NBases;
    if (NULL == (t->edBases = (char *)xmalloc(t->MaxNed+1))) {
	trace_unload(t);
	return -1;
    }
    if (NULL == (t->edPos = (int_2 *)xmalloc(t->MaxNed * sizeof(int_2)))) {
	trace_unload(t);
	return -1;
    }
    if (NULL == (t->edConf = (int1 *)xmalloc(t->MaxNed))) {
	trace_unload(t);
	return -1;
    }
    memcpy(t->edBases, t->read->base, t->read->NBases);
    t->edBases[t->read->NBases] = 0;

    if (t->read->prob_A && t->read->prob_C &&
	t->read->prob_G && t->read->prob_T) {
	for (i = 0; i < t->read->NBases; i++) {
	    switch (t->read->base[i]) {
	    case 'a':
	    case 'A':
		t->edConf[i] = t->read->prob_A[i];
		break;
	    case 'c':
	    case 'C':
		t->edConf[i] = t->read->prob_C[i];
		break;
	    case 'g':
	    case 'G':
		t->edConf[i] = t->read->prob_G[i];
		break;
	    case 't':
	    case 'T':
		t->edConf[i] = t->read->prob_T[i];
		break;
	    default:
		t->edConf[i] = (t->read->prob_A[i] +
				t->read->prob_C[i] +
				t->read->prob_G[i] +
				t->read->prob_T[i]) / 4;
		break;
	    }
	    t->edPos[i] = i+1;
	}
    } else {
	for (i = 0; i < t->read->NBases; i++) {
	    t->edConf[i] = 0;
	    t->edPos[i] = i+1;
	}
    }

    t->leftVector = -1;
    t->rightVector = -1;
    /*
     * If we've loaded an experiment file, then we'll use this new sequence as
     * the edited sequence. This requires updating both the edPos, edBases
     * and edConf arrays as well as their sizes.
     */
    if (t->read->orig_trace_format == TT_EXP && t->read->orig_trace) {
	char *str;

	e = (Exp_info *)t->read->orig_trace;
	if (exp_Nentries(e, EFLT_SQ) && (str = exp_get_entry(e, EFLT_SQ))) {
	    t->Ned = strlen(str);
	    t->MaxNed = t->Ned + TRACE_EDITS;

	    t->edBases = (char *)xrealloc(t->edBases, t->MaxNed+1);
	    strcpy(t->edBases, str);
	}

	if (exp_Nentries(e, EFLT_ON) && (str = exp_get_entry(e, EFLT_ON))) {
	    t->edPos = (int_2 *)xrealloc(t->edPos, t->MaxNed * 2);
	    str2opos(t->edPos, t->MaxNed, str);
	}

	if (exp_Nentries(e, EFLT_AV) && (str = exp_get_entry(e, EFLT_AV))) {
	    t->edConf = (int1 *)xrealloc(t->edConf, t->MaxNed);
	    str2conf(t->edConf, t->MaxNed, str);
	}

	if (exp_Nentries(e, EFLT_SL)) {
	    exp_get_int(e, EFLT_SL, &t->leftVector);
	}

	if (exp_Nentries(e, EFLT_SR)) {
	    exp_get_int(e, EFLT_SR, &t->rightVector);
	}
    }

    t->comp = 0;
    /* Allocate and initialise the sample-base coordinate translations */
    if (t->read->NPoints) {
	if (NULL == (t->tracePos = (uint_2 *)xmalloc(t->read->NPoints *
						    sizeof(uint_2)))) {
	    trace_unload(t);
	    return -1;
	}

	if (NULL == (t->tracePosE = (uint_2 *)xmalloc(t->read->NPoints *
						     sizeof(uint_2)))) {
	    trace_unload(t);
	    return -1;
	}

	trace_init_pos(t);

    } else {
	/* invent some tracePos, tracePosE and basePos numbers */
	if (NULL == (t->tracePos = (uint_2 *)xmalloc(8 * (t->read->NBases+1) *
						    sizeof(uint_2))) ||
	    NULL == (t->tracePosE = (uint_2 *)xmalloc(8 * (t->read->NBases+1) *
						      sizeof(uint_2))) ||
            NULL == (t->read->basePos =
                     (uint_2 *)xmalloc((t->read->NBases+1) * sizeof(uint_2)))){
	    trace_unload(t);
	    return -1;
	}

	trace_init_pos(t);
    }

    return 0;
}



/*-----------------------------------------*/
/* Private helper methods for trace_save() */
/*-----------------------------------------*/

static int trace_save_as_experiment_file( DNATrace* t, char* file )
{
    int       i;
    Exp_info* e;
    mFILE*     fp;
    char*     _base;
    int       _NBases;
    char*     buf = 0;
    int       ret = -1;



    /* Switch edited and original sequences */
    _base           = t->read->base;
    _NBases         = t->read->NBases;
    t->read->base   = t->edBases;
    t->read->NBases = t->Ned;



    /* Convert read to experiment file */
    e = read2exp( t->read, file );



    /* Restore original read structure */
    t->read->base   = _base;
    t->read->NBases = _NBases;
    if( !e )
	return -1;



    /* Create the output file */
    fp = mfopen( file, "w" );
    if( !fp )
	goto cleanup_and_exit;
    e->fp = fp;



    /* Allocate some working storage */
    _NBases = t->Ned;
    i = 5 * _NBases + 2;
    if( i < 256 )
	i = 256;
    buf = (char*) xmalloc( i );
    if( !buf )
	goto cleanup_and_exit;



    /* Add original position information */
    for( i=0; i<_NBases; i++ )
	if( i+1 != t->edPos[i] )
	    break;
    if( i != _NBases )
    {
	Array a;

	/* Lop off old AV and ON lines */
	a = e->entries[EFLT_ON];
	for (i = 0; i < e->Nentries[EFLT_ON]; i++) {
	    if (arr(char *,a,i) != NULL)
		xfree(arr(char *,a,i));
	    e->Nentries[EFLT_ON] = 0;
	}
	a = e->entries[EFLT_AV];
	for (i = 0; i < e->Nentries[EFLT_AV]; i++) {
	    if (arr(char *,a,i) != NULL)
		xfree(arr(char *,a,i));
	    e->Nentries[EFLT_AV] = 0;
	}

	exp_print_mfile( fp, e );
	opos2str( (int_2*) t->edPos, _NBases, buf );
	exp_put_str( e, EFLT_ON, buf, strlen(buf) );
	conf2str( t->edConf, _NBases, buf );
	exp_put_str( e, EFLT_AV, buf, strlen(buf) );
    }
    else
    {
	exp_print_mfile( fp, e );
    }



    /* Add cutoff info. */
    if( t->leftVector != -1 ) 
    {
	sprintf( buf, "%d", t->leftVector );
	exp_put_str( e, EFLT_SL, buf, strlen(buf) );
    }
    if( t->rightVector != -1 ) 
    {
	sprintf( buf, "%d", t->rightVector );
	exp_put_str( e, EFLT_SR, buf, strlen(buf) );
    }
    ret = 0;



    /* Exit */
    cleanup_and_exit:
    if(buf)   xfree( buf );
    if(e->fp) mfclose( e->fp );
    e->fp = NULL;
    exp_destroy_info( e );
    return ret;
}



static int trace_save_as_plain_text_file( DNATrace* t, char* file )
{
    int ret;


    /* Save and swap values */
    char* _base     = t->read->base;
    int   _NBases   = t->read->NBases;
    t->read->base   = t->edBases;
    t->read->NBases = t->Ned;



    /* Set clip points */
    t->read->base += t->read->leftCutoff;
    if( t->read->rightCutoff )
	t->read->NBases = t->read->rightCutoff - t->read->leftCutoff - 1;
    else
	t->read->NBases = t->read->NBases - t->read->leftCutoff;



    /* Save sequence */
    ret = write_reading( file, t->read, TT_PLN );



   /* Restore values */
   t->read->base   = _base;
   t->read->NBases = _NBases;
   return ret;
}



static int trace_save_as_trace_file( DNATrace* t, char* file, int format )
{
    uint_2* _basePos;
    int     _NBases;
    char*   _base;
    char*   _prob_A;
    char*   _prob_C;
    char*   _prob_G;
    char*   _prob_T;
    char*   prob_A  = 0;
    char*   prob_C  = 0;
    char*   prob_G  = 0;
    char*   prob_T  = 0;
    uint_2* basePos = 0;
    int     ret     = -1;



    /* Memory allocation */
    _NBases = t->Ned;
    prob_A  = (char*)   xmalloc( _NBases );
    prob_C  = (char*)   xmalloc( _NBases );
    prob_G  = (char*)   xmalloc( _NBases );
    prob_T  = (char*)   xmalloc( _NBases );
    basePos = (uint_2*) xmalloc( _NBases * sizeof(uint_2) );
    if( !basePos || !prob_A || !prob_C || !prob_G || !prob_T )
	goto cleanup_and_exit;
    

    
    /* Fill the new base position & confidence arrays with edited data */
    read_update_base_positions( t->read, t->comp, t->Ned, t->edBases, t->edPos, basePos );
    read_update_confidence_values( t->Ned, t->edBases, t->edConf, prob_A, prob_C, prob_G, prob_T );



    /* Save the edited trace without altering the original using the ole swap trick */
    _prob_A  = t->read->prob_A;
    _prob_C  = t->read->prob_C;
    _prob_G  = t->read->prob_G;
    _prob_T  = t->read->prob_T;
    _base    = t->read->base;
    _NBases  = t->read->NBases;
    _basePos = t->read->basePos;
    t->read->prob_A  = prob_A;
    t->read->prob_C  = prob_C;
    t->read->prob_G  = prob_G;
    t->read->prob_T  = prob_T;
    t->read->base    = t->edBases;
    t->read->NBases  = t->Ned;
    t->read->basePos = basePos;



    /* Save the trace */
    ret = write_reading( file, t->read, format );



    /* Restore original fields */
    t->read->prob_A  = _prob_A;
    t->read->prob_C  = _prob_C;
    t->read->prob_G  = _prob_G;
    t->read->prob_T  = _prob_T;
    t->read->base    = _base;
    t->read->NBases  = _NBases;
    t->read->basePos = _basePos;



    /* Exit */
    cleanup_and_exit:
    if(basePos) xfree(basePos);
    if(prob_A)  xfree(prob_A);
    if(prob_C)  xfree(prob_C);
    if(prob_G)  xfree(prob_G);
    if(prob_T)  xfree(prob_T);
    return ret;
}



/*
 * Saves a given trace file in the specified output format. Just a stub
 * that calls helpers to do the real work.
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */

int trace_save( DNATrace* t, char* file, char *format )
{
    int ret;
    int fmt;


    /* Check input */
    assert(t);
    assert(file);
    assert(*file);
    assert(format);
    assert(*format);
    assert(t->read);
    if( !t || !file || !format || !t->read || !(*file) || !(*format) )
	return -1;


    /* Deduce format */
    fmt = trace_type_str2int( format );


    /* Do the save */
    switch( fmt )
    {
	case TT_EXP: ret = trace_save_as_experiment_file( t, file );  break;
	case TT_PLN: ret = trace_save_as_plain_text_file( t, file );  break;
        default:     ret = trace_save_as_trace_file( t, file, fmt );  break;
    }
    return ret;
}



/*
 * Starting as position 'pos' in edited sequence, we search backwards for the
 * first base corresponding to an original base.
 */
int trace_find_prev_orig(DNATrace *t, int pos) {
    for (; pos >=0; pos--)
	if (t->edPos[pos] != 0)
	    return pos;

    return 0;
}

/*
 * Insert a base leftwards of 'pos'. Position 0 is left of first base.
 */
void trace_insert(DNATrace *t, int pos, char base) {
    int len = t->Ned - pos + 1;
    int pos2;

    if (pos + len > t->MaxNed)
	len = t->MaxNed - (pos+1);

    /* original positions */
    memmove(&t->edPos[pos+1], &t->edPos[pos], len * sizeof(*t->edPos));
    t->edPos[pos] = 0;

    /* Confidence */
    memmove(&t->edConf[pos+1], &t->edConf[pos], len * sizeof(*t->edConf));
    t->edConf[pos] = 100;

    /* sequence */
    memmove(&t->edBases[pos+1], &t->edBases[pos], len);
    t->edBases[pos] = base;

    /* sample points */
    pos2 = t->read->basePos[t->edPos[trace_find_prev_orig(t, pos - 1)]] + 1;
    while (t->tracePosE[pos2] < pos)
	pos2++;

    /* NPoints is 0 for plain, so this works */
    for (; pos2 < t->read->NPoints; pos2++)
	t->tracePosE[pos2]++;

    if (t->read->leftCutoff && t->read->leftCutoff >= pos) {
	t->read->leftCutoff++;
    }

    if (t->leftVector && t->leftVector >= pos) {
	t->leftVector++;
    }

    if (t->read->rightCutoff && t->read->rightCutoff >= pos) {
	t->read->rightCutoff++;
    }

    if (t->rightVector && t->rightVector >= pos) {
	t->rightVector++;
    }

    t->Ned++;
    t->cursor_pos++;
}

/*
 * Delete a base leftwards of pos
 */
void trace_delete(DNATrace *t, int pos) {
    int len = t->Ned - pos, pos2;

    if (pos < 1)
	return;

    /* sample points */
    pos2 = t->read->basePos[trace_find_prev_orig(t, pos - 1)] + 1;

    /* original positions */
    memmove(&t->edPos[pos-1], &t->edPos[pos], len * sizeof(*t->edPos));

    /* confidence */
    memmove(&t->edConf[pos-1], &t->edConf[pos], len * sizeof(*t->edConf));

    /* sequence */
    memmove(&t->edBases[pos-1], &t->edBases[pos], len);

    while (t->tracePosE[pos2] < pos)
	pos2++;

    /* NPoints is 0 for plain, so this works */
    for (; pos2 < t->read->NPoints; pos2++)
	t->tracePosE[pos2]--;

    if (t->read->leftCutoff >= pos) {
	t->read->leftCutoff--;
    }

    if (t->leftVector >= pos) {
	t->leftVector--;
    }

    if (t->read->rightCutoff >= pos) {
	t->read->rightCutoff--;
    }

    if (t->rightVector >= pos) {
	t->rightVector--;
    }

    t->Ned--;
    t->cursor_pos--;
}

/*
 * Change a base at pos
 */
void trace_change(DNATrace *t, int pos, char base) {
    t->edPos[pos] = 0;
    t->edConf[pos] = 100;
    t->edBases[pos] = base;
}


/*
 * Takes a pyrosequencing style trace with a series of spikes and multiple
 * bases positioned at the same spike, and realigns it so that basecalls are
 * the primary spacing with one "x unit" per base-call or spike.
 *
 * The traceA, traceC, traceG and traceT channels are created from an
 * alignment of the sequencing to the flow/flow_order data.
 *
 * This means that the standard trace viewing code works as expected.
 */
void trace_pyroalign(Read *r) {
    int i, k, len = 0, last;
    TRACE *out[4];
    int ip, op;
    int lookup[256];

    /* Compute size */
    for (last = -1, i = 0; i < r->NBases; i++) {
	if (r->basePos[i] == last)
	    len++;
	else
	    len += r->basePos[i] - last;
	last = r->basePos[i];
    }
    len += r->nflows - last;
    
    out[0] = (TRACE *)xcalloc(len+1, sizeof(TRACE));
    out[1] = (TRACE *)xcalloc(len+1, sizeof(TRACE));
    out[2] = (TRACE *)xcalloc(len+1, sizeof(TRACE));
    out[3] = (TRACE *)xcalloc(len+1, sizeof(TRACE));
    memset(lookup, 0, 256*sizeof(lookup[0]));
    lookup['A'] = lookup['a'] = 0;
    lookup['C'] = lookup['c'] = 1;
    lookup['G'] = lookup['g'] = 2;
    lookup['T'] = lookup['t'] = 3;

    r->maxTraceVal = 1;
    for (k = ip = 0, op = 1; ip < r->nflows || k < r->NBases; ip++) {
	int v;
	v = out[lookup[(r->flow_order[ip])]][op] =
	    MAX(r->flow[ip] * PYRO_SCALE, 1);
	if (r->maxTraceVal < v)
	    r->maxTraceVal = v;
	if (k < r->NBases) {
	    if (r->basePos[k] == ip+1) {
		r->basePos[k] = op;
		for(k++; k < r->NBases && r->basePos[k] == ip+1; k++) {
		    r->basePos[k] = ++op;
		}
	    }
	}
	op++;
    }

    if (r->traceA) xfree(r->traceA); r->traceA = out[0];
    if (r->traceC) xfree(r->traceC); r->traceC = out[1];
    if (r->traceG) xfree(r->traceG); r->traceG = out[2];
    if (r->traceT) xfree(r->traceT); r->traceT = out[3];
    r->NPoints = op;

    r->maxTraceVal = PYRO_SCALE * (int)((r->maxTraceVal+PYRO_SCALE-1) /
					PYRO_SCALE);

#if 0
    /*
     * This works as a way to inline edit the flow and flow_order arrays,
     * but it breaks lots of code elsewhere that looks in NPoints, traceA,
     * etc.
     */
    /* Regenerate data */
    out_flow = (char *)xcalloc(len, sizeof(float));
    in_flow  = r->flow_order;
    out_proc = (float *)xcalloc(len, sizeof(float));
    in_proc  = r->flow;

    for (k = op = ip = 0; ip < r->nflows && k < r->NBases; ip++) {
	out_proc[op] = in_proc[ip];
	out_flow[op] = in_flow[ip];
	if (r->basePos[k] == ip+1) {
	    printf("1:r->basePos[%d] = %d -> %d\n",
		   k, ip+1, op+1);
	    r->basePos[k] = op+1;
	    for(k++; k < r->NBases && r->basePos[k] == ip+1; k++) {
		printf("2:r->basePos[%d] = %d -> %d\n",
		       k, ip+1, 1+op+1);
		r->basePos[k] = ++op+1;
		out_proc[op] = 0.5;
		out_flow[op] = '-';
	    }
	}
	op++;
    }

    xfree(r->flow); r->flow = out_proc;
    xfree(r->flow_order); r->flow_order = out_flow;
    r->nflows = op;
#endif
}
