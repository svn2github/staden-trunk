/*
 * File: gap-dbstruct.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: gap database file struct
 *
 * Created: 28 October 1992
 * Updated:
 *
 */

#include <stdio.h> /* IMPORT: NULL */
#include <string.h>

#include "g-error.h"
#include "gap-dbstruct.h"


char *file_list[] = {
    "",
};

size_t block_sizes[] = {
    16,				/* was 4 */
};	

GCardinal max_recs[] = {
    100000,
};




char *gap_construct_file(char *database,char *file, char *version, char *fillbuf)
/*
 * Given the database name, the file in the database and the version
 * construct an amalgam of all three components.
 *
 * currently this is:
 *	{database}.{file}{version}
 */
{
    static char fname[1024];
    char *buffer;

    if (fillbuf == NULL) {
	if (strlen(database) + strlen(file) + strlen(version) + 1 >= 1024) {
	    (void)gerr_set(GERR_NAME_TOO_LONG);
	    return NULL;
	}
	buffer = fname;
    } else
	buffer = fillbuf;

    sprintf(buffer,"%s.%s%s",database,file,version);

    return buffer;
}




/*************************************************************
 * Miscellaneous debugging routines
 *************************************************************/

void dumpGDatabase(GDatabase *r)
{
    printf("**\n");
    printf("** GDatabase\n");
    printf("**          version = %d\n", r->version);
    printf("**  maximum_db_size = %d\n", r->maximum_db_size);
    printf("**   actual_db_size = %d\n", r->actual_db_size);
    printf("**      max_gel_len = %d\n", r->max_gel_len);
    printf("**       data_class = %d\n", r->data_class);
    printf("**      num_contigs = %d\n", r->num_contigs);
    printf("**     num_readings = %d\n", r->num_readings);
    printf("**        Nfreerecs = %d\n", r->Nfreerecs);
    printf("**         freerecs = %d\n", r->freerecs);
    printf("**         Ncontigs = %d\n", r->Ncontigs);
    printf("**          contigs = %d\n", r->contigs);
    printf("**        Nreadings = %d\n", r->Nreadings);
    printf("**         readings = %d\n", r->readings);
    printf("**     Nannotations = %d\n", r->Nannotations);
    printf("**      annotations = %d\n", r->annotations);
    printf("** free_annotations = %d\n", r->free_annotations);
    printf("**       Ntemplates = %d\n", r->Ntemplates);
    printf("**        templates = %d\n", r->templates);
    printf("**          Nclones = %d\n", r->Nclones);
    printf("**           clones = %d\n", r->clones);
    printf("**         Nvectors = %d\n", r->Nvectors);
    printf("**          vectors = %d\n", r->vectors);
    printf("**\n");
}

void dumpGReadings(GReadings *r)
{
    printf("**\n");
    printf("** GReadings\n");
    printf("**            name = %8d\n", r->name);
    printf("**      trace_name = %8d\n", r->trace_name);
    printf("**      trace_type = %8d\n", r->trace_type);
    printf("**            left = %8d\n", r->left);
    printf("**           right = %8d\n", r->right);
    printf("**        position = %8d\n", r->position);
    printf("**          length = %8d\n", r->length);
    printf("**           sense = %8d\n", r->sense);
    printf("**        sequence = %8d\n", r->sequence);
    printf("**      confidence = %8d\n", r->confidence);
    printf("**  orig_positions = %8d\n", r->orig_positions);
    printf("**     annotations = %8d\n", r->annotations);
    printf("** sequence_length = %8d\n", r->sequence_length);
    printf("**           start = %8d\n", r->start);
    printf("**             end = %8d\n", r->end);
    printf("**        template = %8d\n", r->template);
    printf("**          strand = %8d\n", r->strand);
    printf("**          primer = %8d\n", r->primer);
    printf("**\n");
}


void dumpGContigs(GContigs *r)
{
    printf("**\n");
    printf("** GContigs\n");
    printf("**             left = %8d\n", r->left);
    printf("**            right = %8d\n", r->right);
    printf("**           length = %8d\n", r->length);
    printf("**	    annotations = %8d\n", r->annotations);
    printf("**\n");
}


