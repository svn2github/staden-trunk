#include <tcl.h>
#include <float.h>
#include <string.h>
#include <math.h>

#include "misc.h"
#include "tcl_utils.h"
#include "sequence_formats.h"
#include "seq_results.h"
#include "sip_results.h"
#include "sip_globals.h"
#include "xalloc.h"
#include "align.h"
#include "probs.h"
#include "dna_utils.h"
#include "readpam.h"
#include "tclCanvGraph.h"
#include "sip_similar_spans.h"

static void free_mat_name(mat_name *mat);
static mat_name *alloc_mat_name(void);
static void copy_mat_name(mat_name *to, mat_name *from);

static mat_name *identity_prot_mat = NULL;
static mat_name *prot_mat = NULL;
static mat_name *dna_mat = NULL;
static mat_name *identity_dna_mat = NULL;
static int replot_temp;
static int remove_dup = -1;
static int max_matches;
static int def_matches;


void SipFreeResults(void)
{
    if (identity_prot_mat) {
	free_mat_name(identity_prot_mat);
    }
    if (identity_dna_mat) {
	free_mat_name(identity_dna_mat);
    }
    if (prot_mat) {
	free_mat_name(prot_mat);
    }
    if (dna_mat) {
	free_mat_name(dna_mat);
    }
}

/* setting and getting global variables */
void set_remove_dup(int yn) 
{
    if (yn) {
	remove_dup = 1;
    } else {
	remove_dup = 0;
    }
}

int get_remove_dup(void)
{
    return remove_dup;
}

void set_replot_temp(int response)
{
    replot_temp = response;
}

int get_replot_temp(void)
{
    return replot_temp;
}
void set_max_matches(int num)
{
    max_matches = num;
}

int get_max_matches(void)
{
    return max_matches;
}

void set_def_matches(int num)
{
    def_matches = num;
}

int get_def_matches(void)
{
    return def_matches;
}

/*
 * set the identity score matrix file for proteins or dna
 * called once for PROTEIN and DNA on startup
 * used by find identity routines
 */
int set_matrix_identity(int type)
{
    if (type == PROTEIN) {
	set_char_set(PROTEIN);

	/* initialise prot_mat */
	if (identity_prot_mat == NULL) {
	    if (NULL == (identity_prot_mat = alloc_mat_name()))
		return -1;
	}
	identity_prot_matrix(&identity_prot_mat->matrix);
    } else {
	set_char_set(DNA);

	/* initialise dna_mat */
	if (identity_dna_mat == NULL) {
	    if (NULL == (identity_dna_mat = alloc_mat_name()))
		return -1;
	}
	identity_dna_matrix(&identity_dna_mat->matrix);
    }
    return 0;
}

/*
 * Frees a mat_name structure as allocated by alloc_mat_name.
 * This copes with a variety of NULL inputs and partially allocated matrices.
 */
static void free_mat_name(mat_name *mat) {
    int i;

    if (!mat)
	return;

    if (mat->name)
	xfree(mat->name);

    if (mat->matrix) {
	for (i = 0; i < MAX_SCORE_MATRIX; i++) {
	    if (mat->matrix[i])
		free(mat->matrix[i]);
	}
	free(mat->matrix);
    }

    free(mat);
}

/*
 * Allocates an initialises (to blank) a mat_name structure.
 * This structure should be freed with a call to free_mat_name().
 * The matrix is returned initialised to contain all zeros and a null name.
 *
 * Returns pointer on success
 *         NULL on failure.
 */
static mat_name *alloc_mat_name(void) {
    int i;
    mat_name *mat = NULL;

    /* Allocate memory */
    if (NULL == (mat = (mat_name *)xcalloc(1, sizeof(*mat)))) {
	free_mat_name(mat);
	return NULL;
    }
    if (NULL == (mat->name = (char *)xcalloc(FILENAME_MAX, 1))) {
	free_mat_name(mat);
	return NULL;
    }
    if (NULL == (mat->matrix = (int **)xcalloc(MAX_SCORE_MATRIX,
					       sizeof(int *)))) {
	free_mat_name(mat);
	return NULL;
    }
    for (i = 0; i < MAX_SCORE_MATRIX; i++) {
	if (NULL == (mat->matrix[i] = (int *)xcalloc(MAX_SCORE_MATRIX,
						     sizeof(int)))) {
	    free_mat_name(mat);
	    return NULL;
	}
    }

    /* Already initialised to zero by the above calloc, so none needed */

    return mat;
}

/*
 * Copies one mat_name structure to another.
 */
static void copy_mat_name(mat_name *to, mat_name *from) {
    int i;

    if (!to || !from)
	return;

    if (to->name && from->name) {
	strncpy(to->name, from->name, FILENAME_MAX);
    }

    if (to->matrix && from->matrix) {
	for (i = 0; i < MAX_SCORE_MATRIX; i++) {
	    if (to->matrix[i] && from->matrix[i])
		memcpy(to->matrix[i], from->matrix[i],
		       MAX_SCORE_MATRIX * sizeof(int));
	}
    }
}

/*
 * set the score matrix file for proteins or dna
 * called once for PROTEIN and DNA on startup
 * called whenever the score matrix is changed (protein only at the present)
 */
int
set_matrix_file(char *file,
		int type)
{
    int i, j;
    mat_name *temp_mat = NULL;  

    if (type == PROTEIN) {
	set_char_set(PROTEIN);

	/* initialise prot_mat */
	if (prot_mat == NULL) {
	    if (NULL == (prot_mat = alloc_mat_name()))
		return -1;
	} else {
	    /* Take temporary backup */
	    if (NULL == (temp_mat = alloc_mat_name()))
		return -1;
	    copy_mat_name(temp_mat, prot_mat);
	}

	/* initialise matrix to 0 */
	for (i = 0; i < MAX_SCORE_MATRIX; i++) {
	    for (j = 0; j < MAX_SCORE_MATRIX; j++) {
		prot_mat->matrix[i][j] = 0;
	    }
	}

	/* if file is NULL, set up the protein identity matrix */
	if (file == NULL) {
	    identity_prot_matrix(&prot_mat->matrix);
	    if (prot_mat->name)
		free(prot_mat->name);
	    prot_mat->name = NULL;
	    free_mat_name(temp_mat);
	} else {
	    strcpy(prot_mat->name, file);
	    /*
	     * need to check if file is a valid matrix and reinstate previous
	     * matrix if failed to load new one
	     */
	    if (-1 == (create_pam_matrix(file, &prot_mat->matrix))) {
		copy_mat_name(prot_mat, temp_mat);
		free_mat_name(temp_mat);
	        return -1;
	    }
	    free_mat_name(temp_mat);
	}
    } else {
	set_char_set(DNA);

	/* initialise dna_mat */
	if (dna_mat == NULL) {
	    if (NULL == (dna_mat = alloc_mat_name()))
		return -1;
	}
	if (dna_mat->name)
	    free(dna_mat->name);
	dna_mat->name = NULL;

	/* if file is NULL, set up the dna identity matrix */
	if (file == NULL) {
	    identity_dna_matrix(&dna_mat->matrix);
	} else {
	    verror(ERR_WARN, "set score matrix", "Ability to change the DNA"
		   "score matrix is not supported at present");
	    /* 
	     * if (NULL == (dna_mat->name = (char *)xmalloc(FILENAME_MAX * 
	     * sizeof(char))))
	     * return -1;
	     * strcpy(dna_mat->name, file);
	     * create_pam_matrix(file, &dna_mat->matrix);
	     */
	    return -1;
	}
    }
    return 0;
}

/*
 * return current score matrix file for protein or dna
 */
int **
get_matrix_file(int type)
{
    if (type == PROTEIN) {
	return prot_mat->matrix;
    } else {
	return dna_mat->matrix;
    }
}

/*
 * return current identity score matrix file for protein or dna
 */
int **
get_matrix_identity(int type)
{
    if (type == PROTEIN) {
	return identity_prot_mat->matrix;
    } else {
	return identity_dna_mat->matrix;
    }
}

/*
 * return current score matrix file for protein or dna
 */
char *
get_matrix_name(int type)
{

    if (type == PROTEIN) {
	return prot_mat->name;
    } else {
	return dna_mat->name;
    }
}

void SipFreeResult(seq_result *result) {
    if (result) {
#ifdef FIXME
	d_plot *data = result->data;
	out_raster *output = result->output;

	xfree(data->p_array);
	xfree(data);
	free(output->name);

	xfree(output->configure[0]);
	xfree(output->configure);
	xfree(result->output);
#endif
	xfree(result->input);
	xfree(result);
    }
}

