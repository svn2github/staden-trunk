#include <string.h>
#include "list.h"
#include "process.h"
#include "dapDB.h"
#include "bapDB.h"
#include "gapDB.h"
#include "flat_sd.h"
#include "misc.h"

void open_for_read(List *from)
{
    char *a;

    if ( (a = assoc(from,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	flat_sd_open_for_read(from);
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	xdap_middle_open_for_read(from);
    else if (strcmp(a,db_type_late_xdap) == 0)
	xdap_late_open_for_read(from);
    else if (strcmp(a,db_type_gap) == 0)
	gap_open_for_read(from);
    else
	crash("Source database type not supported\n");
}

void open_for_write(List *to)
{
    char *a;

    if ( (a = assoc(to,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	flat_sd_open_for_write(to);
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_late_xdap) == 0)
	xdap_late_open_for_write(to);
    else if (strcmp(a,db_type_gap) == 0)
	gap_open_for_write(to);
    else
	crash("Destination database type not supported\n");
}

void close_files(List *to)
{
    char *a;

    if ( (a = assoc(to,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	flat_sd_close(to);
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	xdap_middle_close(to);
    else if (strcmp(a,db_type_late_xdap) == 0)
	xdap_late_close(to);
    else if (strcmp(a,db_type_gap) == 0)
	gap_close(to);
    else
	crash("Destination database type not supported\n");
}





List *read_header(List *from)
{
    char *a;

    if ( (a = assoc(from,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	return (List *)flat_sd_read_header();
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	return (List *)xdap_middle_read_header();
    else if (strcmp(a,db_type_late_xdap) == 0)
	return (List *)xdap_late_read_header();
    else if (strcmp(a,db_type_gap) == 0)
	return (List *)gap_read_header();
    else
	crash("Source database type not supported\n");
    return (List *)0; /* stops warnings */
}


void write_header(List *to, List *l)
{
    char *a;

    if ( (a = assoc(to,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	flat_sd_write_header(l);
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_late_xdap) == 0)
	xdap_late_write_header(l);
    else if (strcmp(a,db_type_gap) == 0)
	gap_write_header(l);
    else
	crash("Destination database type not supported\n");
}


List *read_gel_data(List *from)
{
    char *a;

    if ( (a = assoc(from,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	return (List *)flat_sd_read_gel_data();
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	return (List *)xdap_middle_read_gel_data();
    else if (strcmp(a,db_type_late_xdap) == 0)
	return (List *)xdap_late_read_gel_data();
    else if (strcmp(a,db_type_gap) == 0)
	return (List *)gap_read_gel_data();
    else
	crash("Source database type not supported\n");
    return (List *)0; /* stops warnings */
}

void write_gel_data(List *to, List *l)
{
    char *a;

    if ( (a = assoc(to,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	flat_sd_write_gel_data(l);
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_late_xdap) == 0)
	xdap_late_write_gel_data(l);
    else if (strcmp(a,db_type_gap) == 0)
	gap_write_gel_data(l);
    else
	crash("Source database type not supported\n");
}

List *read_contig_data(List *from)
{
    char *a;

    if ( (a = assoc(from,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	return (List *)flat_sd_read_contig_data();
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	return (List *)xdap_middle_read_contig_data();
    else if (strcmp(a,db_type_late_xdap) == 0)
	return (List *)xdap_late_read_contig_data();
    else if (strcmp(a,db_type_gap) == 0)
	return (List *)gap_read_contig_data();
    else
	crash("Source database type not supported\n");
    return (List *)0; /* stops warnings */
}

void write_contig_data(List *to, List *l)
{
    char *a;

    if ( (a = assoc(to,db_type)) == NULL)
	crash("Type not specified for source database\n");
    else if (strcmp(a,db_type_RS_flat_file) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_SD_flat_file) == 0)
	flat_sd_write_contig_data(l);
    else if (strcmp(a,db_type_sap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_early_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_middle_xdap) == 0)
	crash("Source database type not supported\n");
    else if (strcmp(a,db_type_late_xdap) == 0)
	xdap_late_write_contig_data(l);
    else if (strcmp(a,db_type_gap) == 0)
	gap_write_contig_data(l);
    else
	crash("Source database type not supported\n");
}

void process(List *from, List *to)
{
    List *l;


    /*
    ** Initialise read
    */
    open_for_read(from);
    open_for_write(to);

    l = read_header(from);
    write_header(to,l);
    destroy_list(l);

    /*
    ** Process Gels
    */
    for (l = read_gel_data(from);
	 !isNil(l);
	 l = read_gel_data(from)) {
	write_gel_data(to,l);
	destroy_list(l);
    }

    /*
    ** Process Contigs
    */
    for (l = read_contig_data(from);
	 !isNil(l);
	 l = read_contig_data(from)) {
	write_contig_data(to,l);
	destroy_list(l);
    }


    /*
    ** Tidy up read
    */
    close_files(from);
    close_files(to);
	
}






