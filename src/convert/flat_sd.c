#include <staden_config.h>

#include <string.h>
#include "list.h"
#include "process.h"
#include "misc.h"

typedef char IOString[200];

static IOString flat_file;
static FILE *flat_fp;
static List *next_list;

static void set_file_names(char *name)
{
   strcpy(flat_file,name); strcat(flat_file,".flat");
}

void flat_sd_open_for_write(List *l)
{
    char *name;

    name = assoc(l,db_name);

    /*
    ** Create file name
    */
    set_file_names(name);


    /*
    ** Check for existance
    */
    if ( file_exists(flat_file) ) crash("Flat file %s already exists\n",flat_file);

    /*
    ** Open files
    */
    if ( (flat_fp = (fopen(flat_file,"w"))) == NULL) crash("Cannot open file %s\n",flat_file);

}

void flat_sd_write_header(List *l)
{
    fprint_list_f(flat_fp,l);
}

void flat_sd_write_gel_data(List *l)
{
    fprint_list_f(flat_fp,l);
}

void flat_sd_write_contig_data(List *l)
{
    fprint_list_f(flat_fp,l);
}


void flat_sd_close(List *l)
{
    fclose(flat_fp);
}

void flat_sd_open_for_read(List *l)
{
    char *name;

    name = assoc(l,db_name);

    /*
    ** Create file name
    */
    set_file_names(name);

    /*
    ** Open files
    */
    if ( (flat_fp = (fopen(flat_file,"r"))) == NULL) crash("Cannot open file %s\n",flat_file);

    next_list = nil;

}

List *flat_sd_read_header(void)
{
    List *l;

    if (isNil(next_list))
	next_list=read_list(flat_fp);

    if (isNil(next_list))
	l = nil;
    else {
	if (strcmp(db_from,atomVal(car(next_list))))
	    l = nil;
	else {
	    l = next_list;
	    next_list = nil;
	}
    }

    return l;
}

List *flat_sd_read_gel_data(void)
{
    List *l;

    if (isNil(next_list))
	next_list=read_list(flat_fp);

    if (isNil(next_list))
	l = nil;
    else {
	if (strcmp(gel_rec,atomVal(car(next_list))))
	    l = nil;
	else {
	    l = next_list;
	    next_list = nil;
	}
    }

    return l;
}

List *flat_sd_read_contig_data(void)
{
    List *l;

    if (isNil(next_list))
	next_list=read_list(flat_fp);

    if (isNil(next_list))
	l = nil;
    else {
	if (strcmp(contig_rec,atomVal(car(next_list))))
	    l = nil;
	else {
	    l = next_list;
	    next_list = nil;
	}
    }

    return l;
}
