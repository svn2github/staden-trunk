#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "process.h"
#include "bapIO.h"
#include "list.h"
#include "misc.h"

/*
** Bap IO
*/
static BapIO io;


static int cur_gel_index;
static int cur_contig_index;



static List *file_details(void)
{
    return
	build_list(
		   atom_str(db_files),
		   build_list(
			      atom_str(db_files_arch),
			      atom_str(io.ar_file),
			      nil),
		   build_list(
			      atom_str(db_files_rel),
			      atom_str(io.rl_file),
			      nil),
		   build_list(
			      atom_str(db_files_seq),
			      atom_str(io.sq_file),
			      nil),
		   build_list(
			      atom_str(db_files_tag),
			      atom_str(io.tg_file),
			      nil),

		   build_list(
			      atom_str(db_files_com),
			      atom_str(io.cc_file),
			      nil),
		   nil
	       );
}




static List *db_details(void)
{
    
    return
	build_list(
		   build_list(
			      atom_str(db_max_db_size),
			      atom_int(io.max_db_size),
			      nil),
		   build_list(
			      atom_str(db_max_gels),
			      atom_int(io.max_gels),
			      nil),
		   build_list(
			      atom_str(db_max_gel_length),
			      atom_int(io.max_gel_length),
			      nil),
		   build_list(
			      atom_str(db_data_class),
			      atom_int(io.data_class),
			      nil),
		   build_list(
			      atom_str(db_num_gels),
			      atom_int(io.num_gels),
			      nil),
		   build_list(
			      atom_str(db_num_contigs),
			      atom_int(io.num_contigs),
			      nil),
		   nil);

}




void xdap_late_open_for_read(List *l)
/*
**
*/
{
    char *name;
    char *version;

    name = assoc(l,db_name);
    if (! name)	crash("No database name specified\n");

    version = assoc(l,db_version);
    if (! version) crash("No version specified\n");


    bap_open_for_read(&io,name,version);

    cur_gel_index = 1;
    cur_contig_index = 1;
}


List *xdap_late_read_header(void)
/*
**
*/
{
    List *files;
    List *details;
    
    files = file_details();
    details = db_details();

    return
	join_list(
		  build_list(
			     atom_str(db_from),
			     nil),
		  details,
		  build_list(
			     files,
			     nil),
		  nil);


}



void xdap_late_close(List *l)
/*
** Close all relevant files
*/
{

    bap_close_files(&io);

}





List *xdap_late_read_gel_data(void)
{
    List *gel_details;

    if (cur_gel_index > io.num_gels)
	gel_details = nil;
    else {
	bap_ar_file_rec ar_line;
	bap_rl_file_rec rl_line;
	bap_sq_file_rec sq_line;
	bap_tg_file_rec tg_line;

	int index;
	char name[17];
	char *seq;
	int length;
	int comp;
	int pos;
	int l_nbr;
	int r_nbr;

	sq_line = (char *) malloc(io.max_gel_length+1);

	bap_read_ar(&io,cur_gel_index,&ar_line);
	bap_read_rl(&io,cur_gel_index,&rl_line);
	bap_read_tg(&io,cur_gel_index,&tg_line);
	bap_read_sq(&io,cur_gel_index,sq_line);


	index = cur_gel_index;
	f2cstr(ar_line.lines.name,BAP_FILE_NAME_LENGTH,name,
	       (size_t)BAP_FILE_NAME_LENGTH);
	length = abs(rl_line.lines.length);
	comp = (rl_line.lines.length < 0);
	seq = sq_line; seq[length] = '\0';
	pos = rl_line.lines.rel_pos;
	l_nbr = rl_line.lines.left_nbr;
	r_nbr = rl_line.lines.right_nbr;

	gel_details =
	    build_list(
		       atom_str(gel_rec),
		       build_list(
				  atom_str(gel_index),
				  atom_int(index),
				  nil),
		       build_list(
				  atom_str(gel_name),
				  atom_str(name),
				  nil),
		       build_list(
				  atom_str(gel_length),
				  atom_int(length),
				  nil),
		       build_list(
				  atom_str(gel_comp),
				  atom_int(comp),
				  nil),
		       build_list(
				  atom_str(gel_seq),
				  atom_str(seq),
				  nil),
		       build_list(
				  atom_str(gel_pos),
				  atom_int(pos),
				  nil),
		       build_list(
				  atom_str(gel_l_nbr),
				  atom_int(l_nbr),
				  nil),
		       build_list(
				  atom_str(gel_r_nbr),
				  atom_int(r_nbr),
				  nil),
		       nil);



	/* Get raw data details */
	if (tg_line.lines.comment) {
	    List *raw_data_details;

	    char *rd;
	    int rd_length;
	    int rd_cut;
	    int rd_ulen;
	    char rd_type[5];
	    char rd_file[19];
	    
	    rd = bap_read_comment(&io, tg_line.lines.comment);
	    sscanf(rd,"%6d%6d%6d%*s",&rd_length, &rd_cut, &rd_ulen);
	    f2cstr(&rd[18],4,rd_type,(size_t)4);
	    f2cstr(&rd[22],18,rd_file,(size_t)18);


	    raw_data_details =
		build_list(
			   build_list(
				      atom_str(gel_rd_length),
				      atom_int(rd_length),
				      nil),
			   build_list(
				      atom_str(gel_rd_cut),
				      atom_int(rd_cut),
				      nil),
			   build_list(
				      atom_str(gel_rd_ulen),
				      atom_int(rd_ulen),
				      nil),
			   build_list(
				      atom_str(gel_rd_type),
				      atom_str(rd_type),
				      nil),
			   build_list(
				      atom_str(gel_rd_file),
				      atom_str(rd_file),
				      nil),
			   nil);

	    join_list (gel_details, raw_data_details, nil);

	}



	/*
	** Process tags, maintaining separate lists for
	** (a) special tags
	** (b) annotation
	** (c) edits
	*/
	{
	    List *specials;
	    List *notes;
	    List *edits;

	    int_4 next;
	    specials = nil;
	    notes =
		build_list(
			   atom_str(gel_annotation),
			   nil);
	    edits =
		build_list(
			   atom_str(gel_edits),
			   nil);

	    while (tg_line.lines.next) {
		next = tg_line.lines.next;
		bap_read_tg(&io,next,&tg_line);

		if (strncmp(tg_line.lines.type.c,"*LC*",4)==0) {
		    if (tg_line.lines.comment) {
			List *lc;
			lc = build_list(
					atom_str(gel_l_cut_seq),
					atom_str(bap_read_comment(&io,tg_line.lines.comment)),
					nil);
			if (isNil(specials))
			    specials = build_list(lc,nil);
			else
			    specials = join_list(specials,build_list(lc,nil),nil);
		    }
		} else if (strncmp(tg_line.lines.type.c,"*RC*",4)==0) {
		    if (tg_line.lines.comment) {
			List *rc;
			rc = build_list(
					atom_str(gel_r_cut_seq),
					atom_str(bap_read_comment(&io,tg_line.lines.comment)),
					nil);
			if (isNil(specials))
			    specials = build_list(rc,nil);
			else
			    specials = join_list(specials,build_list(rc,nil),nil);
		    }
		} else if (strncmp(tg_line.lines.type.c,"*",1)==0) {
		    List *ed;
		    char base[2];
		    base[0] = tg_line.lines.type.c[3];
		    base[1] = '\0';
		    ed = build_list(
				    build_list(
					       atom_str(gel_ed_op),
					       atom_str( (strncmp(tg_line.lines.type.c,"*IN",3))==0
							? gel_ed_insert : gel_ed_delete),
					       nil),
				    build_list(
					       atom_str(gel_ed_base),
					       atom_str(base),
					       nil),
				    build_list(
					       atom_str(gel_ed_base_pos),
					       atom_int( tg_line.lines.position ),
					       nil),
				    nil);
		    edits = join_list(edits, build_list(ed,nil),nil);
		} else {
		    List *an;
		    char type[5];
		    char *com;
		    strncpy(type,tg_line.lines.type.c,4);
		    type[4]='\0';
		    com = bap_read_comment(&io,tg_line.lines.comment);
		    an = build_list(
				    build_list(
					       atom_str(gel_an_pos),
					       atom_int(tg_line.lines.position),
					       nil),
				    build_list(
					       atom_str(gel_an_len),
					       atom_int(tg_line.lines.length),
					       nil),
				    build_list(
					       atom_str(gel_an_type),
					       atom_str(type),
					       nil),
				    (com == NULL) ? nil :
				    build_list(
					       atom_str(gel_an_comment),
					       atom_str(com),
					       nil),
				    nil);
		    notes = join_list(notes,build_list(an,nil),nil);

		}
	    }
	    if (isNil(specials))
		gel_details = join_list(gel_details,build_list(edits,nil),build_list(notes,nil),nil);
	    else
		gel_details = join_list(gel_details,specials,build_list(edits,nil),build_list(notes,nil),nil);

	}



	cur_gel_index++;
	free(sq_line);

    }

    return gel_details;

}





List *xdap_late_read_contig_data(void)
{
    List *contig_details;

    if (cur_contig_index > io.num_contigs)
	contig_details = nil;
    else {
	bap_rl_file_rec rl_line;
	int length;
	int index;
	int left_end;
	int right_end;

	index = io.max_gels-cur_contig_index;
	bap_read_rl(&io,index,&rl_line);
	length = rl_line.clines.length;
	left_end = rl_line.clines.left_end;
	right_end = rl_line.clines.right_end;

	contig_details =
	    build_list(
		       atom_str(contig_rec),
		       build_list(
				  atom_str(contig_index),
				  atom_int(index),
				  nil),
		       build_list(
				  atom_str(contig_length),
				  atom_int(length),
				  nil),
		       build_list(
				  atom_str(contig_left_end),
				  atom_int(left_end),
				  nil),
		       build_list(
				  atom_str(contig_right_end),
				  atom_int(right_end),
				  nil),
		       nil);

	cur_contig_index++;
    }

    return contig_details;

}





void xdap_late_open_for_write(List *l)
/*
**
*/
{
    char *name;
    char *version;


    name = assoc(l,db_name);
    if (! name)	crash("No database name specified\n");

    version = assoc(l,db_version);
    if (! version) crash("No version specified\n");


    bap_open_for_write(&io,name,version);

    cur_gel_index = 1;
    cur_contig_index = 1;

}



void xdap_late_write_header(List *l)
{
    char *a;
    bap_rl_file_rec rl_dbheader;
    bap_rl_file_rec rl_header;
    bap_ar_file_rec ar_line;
    bap_sq_file_rec sq_line;
    bap_tg_file_rec tg_line;
    bap_cc_file_rec cc_line;

    if ( (a = assoc(l,db_data_class)) == NULL)
	crash("Sequence type (DNA or protein) not specified\n");
    else
	io.data_class = atoi(a);

    if ( (a = assoc(l,db_num_gels)) == NULL)
	crash("Number of gels not specified\n");
    else
	io.num_gels = atoi(a);

    if ( (a = assoc(l,db_num_contigs)) == NULL)
	crash("Number of contigs not specified\n");
    else
	io.num_contigs = atoi(a);

    if ( (a = assoc(l,db_max_gel_length)) == NULL)
	crash("Maximum length of a gel reading not specified\n");
    else
	io.max_gel_length = atoi(a);

    if ( (a = assoc(l,db_max_gels)) == NULL)
	crash("Maximum number of gels not specified\n");
    else
	io.max_gels = atoi(a);

    if ( (a = assoc(l,db_max_db_size)) == NULL ) {
	io.max_db_size = 1000;
    } else {
	io.max_db_size = atoi(a);
    }

    rl_dbheader.dbheader.maxdb = io.max_db_size;
    rl_dbheader.dbheader.idbsiz = io.max_gels;
    rl_dbheader.dbheader.maxgel = io.max_gel_length;
    rl_dbheader.dbheader.idm = io.data_class;
    bap_write_rl(&io,bap_rl_dbheader_rec(&io),&rl_dbheader);


    rl_header.header.num_gels = io.num_gels;
    rl_header.header.num_contigs = io.num_contigs;
    bap_write_rl(&io,bap_rl_header_rec(&io),&rl_header);

    bap_write_ar(&io,1,&ar_line);

    sq_line = (bap_sq_file_rec) malloc(io.max_gel_length);
    bap_write_sq(&io,1,sq_line);
    free(sq_line);


    tg_line.header.free_list = 0;
    tg_line.header.count = bap_tg_header_rec(&io);
    bap_write_tg(&io,bap_tg_header_rec(&io),&tg_line);

    cc_line.header.free_list = 0;
    cc_line.header.count = bap_cc_header_rec(&io);
    bap_write_cc(&io,bap_cc_header_rec(&io),&cc_line);
}



void xdap_late_write_gel_data(List *l)
{

    char *a;
    int i;
    bap_rl_file_rec rl_line;
    bap_ar_file_rec ar_line;
    bap_sq_file_rec sq_line;
    bap_tg_file_rec tg_line;
    List *edits;
    List *notes;

    sq_line = (bap_sq_file_rec) malloc(io.max_gel_length);

    /*
    ** Relationship line
    */
    if ( (a = assoc(l,gel_l_nbr)) == NULL)
	crash("No left neighbour for gel %d\n",cur_gel_index);
    else
	rl_line.lines.left_nbr = atoi(a);

    if ( (a = assoc(l,gel_r_nbr)) == NULL)
	crash("No right neighbour for gel %d\n",cur_gel_index);
    else
	rl_line.lines.right_nbr = atoi(a);

    if ( (a = assoc(l,gel_length)) == NULL)
	crash("Length of gel reading not specified for gel %d\n", cur_gel_index);
    else
	rl_line.lines.length = atoi(a);

    if ( (a = assoc(l,gel_comp)) == NULL)
	crash("Not known if gel %d complemented\n", cur_gel_index);
    else {
	i = atoi(a);
	if (i) rl_line.lines.length = -rl_line.lines.length;
    }

    if ( (a = assoc(l,gel_pos)) == NULL)
	crash("No position in contig specified for gel %d\n", cur_gel_index);
    else
	rl_line.lines.rel_pos = atoi(a);

    bap_write_rl(&io,cur_gel_index,&rl_line);



    /*
    ** Archive line
    */
    if ( (a = assoc(l,gel_name)) == NULL)
	crash("No gel name specified for gel %d\n", cur_gel_index);
    else
	c2fstr(a,BAP_FILE_NAME_LENGTH,ar_line.lines.name,
	       (size_t)BAP_FILE_NAME_LENGTH);
	    
    bap_write_ar(&io,cur_gel_index,&ar_line);



    /*
    ** Sequence
    */
    if ( (a = assoc(l,gel_seq)) == NULL)
	crash("No sequence for gel %d\n", cur_gel_index);
    else
	c2fstr(a,strlen(a), sq_line, (size_t)io.max_gel_length);

    bap_write_sq(&io,cur_gel_index,sq_line);


    /*
    ** Initialise tag fields
    */
    tg_line.lines.next = 0;
    tg_line.lines.length = 0;
    tg_line.lines.comment = 0;
    tg_line.lines.position = 0;
    bap_write_tg(&io, cur_gel_index, &tg_line);

    /*
    ** Raw data
    */
    if ( (a = assoc(l,gel_rd_length)) != NULL) {
	int length;
	int cut;
	int ulen;
	char type[5];
	char file[18];
	char s[41];

	length = atoi(a);

	if ( (a = assoc(l,gel_rd_cut)) == NULL)
	    crash ("No raw data left cutoff specified for gel %d\n",cur_gel_index);
	else
	    cut = atoi(a);

	if ( (a = assoc(l,gel_rd_ulen)) == NULL)
	    crash ("No raw data length specified for gel %d\n",cur_gel_index);
	else
	    ulen = atoi(a);

	if ( (a = assoc(l,gel_rd_type)) == NULL)
	    crash ("No raw data file type specified for gel %d\n",cur_gel_index);
	else
	    c2fstr(a,strlen(a),type,sizeof(type));

	if ( (a = assoc(l,gel_rd_file)) == NULL)
	    crash ("No raw data file specified for gel %d\n",cur_gel_index);
	else
	    c2fstr(a,strlen(a),file,sizeof(file));

	/*
	 * Digital Unix V4.0 has a bug with printf. %4.4s needs null
	 * term. strings, even when they're 4 or more characters long. This
	 * is contrary to the standard!.
	 */
	type[4] = 0;
	sprintf(s,"%6d%6d%6d%4.4s%18.18s",length,cut,ulen,type,file);

	tg_line.lines.comment = bap_write_comment(&io,s);
	bap_write_tg(&io, cur_gel_index, &tg_line);
    }






    /*
    ** Edits
    */
    edits = index_list_by_str(l,gel_edits);
    if ( ! isNil(edits)) {

	bap_tg_file_rec tg_last_ed;
	long last_ed_index;
	bap_tg_file_rec tg_ed;
	long ed_index;

	tg_last_ed = tg_ed = tg_line;
	last_ed_index = ed_index = cur_gel_index;

	for( edits = cdr(edits); ! isNil(edits) ; edits = cdr(edits)) {
	    List *e;
	    int pos;
	    char type[4];

	    e = car(edits);

	    if ( (a = assoc(e,gel_ed_op)) == NULL)
		crash ("No edit operation specified, gel %d\n",cur_gel_index);
	    else {
		if (strcmp(a,gel_ed_delete)==0)
		    strcpy(type,"*DE");
		else
		    strcpy(type,"*IN");
	    }

	    if ( (a = assoc(e,gel_ed_base)) == NULL)
		crash ("No base specified for edit, gel %d\n",cur_gel_index);
	    else
		type[3] = a[0];

	    if ( (a = assoc(e,gel_ed_base_pos)) == NULL)
		crash ("No base position specified for edit, gel %d\n",cur_gel_index);
	    else
		pos = atoi(a);

	    ed_index = bap_get_free_tag(&io);
	    tg_last_ed.lines.next = ed_index;

	    bap_write_tg(&io,last_ed_index,&tg_last_ed);
	
	    tg_ed.lines.position = pos;
	    tg_ed.lines.comment = 0;
	    tg_ed.lines.next = 0;
	    strncpy(tg_ed.lines.type.c,type,4);

	    tg_last_ed = tg_ed;
	    last_ed_index = ed_index;

	}

	bap_write_tg(&io,ed_index,&tg_ed);

    }

    /*
    ** Right cut offs
    */
    if ( (a = assoc(l, gel_r_cut_seq)) != NULL) {
	bap_tg_file_rec rcut;
	rcut.lines.position = 0;
	rcut.lines.length = 0;
	rcut.lines.next = 0;
	strncpy(rcut.lines.type.c,"*RC*",4);
	rcut.lines.comment = bap_write_comment(&io,a);
	bap_insert_tag(&io,cur_gel_index,rcut);
    }


    /*
    ** Left cut offs
    */
    if ( (a = assoc(l, gel_l_cut_seq)) != NULL) {
	bap_tg_file_rec lcut;
	lcut.lines.position = 0;
	lcut.lines.length = 0;
	lcut.lines.next = 0;
	strncpy(lcut.lines.type.c,"*LC*",4);
	lcut.lines.comment = bap_write_comment(&io,a);
	bap_insert_tag(&io,cur_gel_index,lcut);
    }


    /*
    ** Annotation
    */
    notes = index_list_by_str(l,gel_annotation);
    if (!isNil(notes)){

	for( notes = cdr(notes); ! isNil(notes) ; notes = cdr(notes)) {
	    List *n;
	    bap_tg_file_rec tg_rec;
	    
	    n = car(notes);

	    if ( (a = assoc(n,gel_an_pos)) == NULL)
		crash("No position for annotation, gel %d\n",cur_gel_index);
	    else
		tg_rec.lines.position = atoi(a);

	    if ( (a = assoc(n,gel_an_len)) == NULL)
		crash("No length for annotation, gel %d\n",cur_gel_index);
	    else
		tg_rec.lines.length = atoi(a);

	    if ( (a = assoc(n,gel_an_type)) == NULL)
		crash("No type for annotation, gel %d\n",cur_gel_index);
	    else
		c2fstr(a,strlen(a),tg_rec.lines.type.c,(size_t)4);

	    if ( (a = assoc(n,gel_an_comment)) == NULL)
		tg_rec.lines.comment = 0;
	    else
		tg_rec.lines.comment = bap_write_comment(&io,a);

	    bap_insert_tag(&io,cur_gel_index,tg_rec);
	}

    }

    cur_gel_index++;

    free(sq_line);

}





void xdap_late_write_contig_data(List *l)
{

    char *a;
    bap_rl_file_rec rl_line;
    long index;

    index = io.max_gels-cur_contig_index;
    
    if ( (a = assoc(l,contig_left_end)) == NULL)
	crash("No left end for contig %d\n", cur_contig_index);
    else
	rl_line.clines.left_end = atoi(a);

    if ( (a = assoc(l,contig_right_end)) == NULL)
	crash("No right end for contig %d\n",cur_contig_index);
    else
	rl_line.clines.right_end = atoi(a);

    if ( (a = assoc(l,contig_length)) == NULL)
	crash("No length for contig %d\n",cur_contig_index);
    else
	rl_line.clines.length = atoi(a);


    bap_write_rl(&io,index,&rl_line);

    cur_contig_index++;

}



