#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "process.h"
#include "list.h"
#include "dapIO.h"
#include "dapDB.h"
#include "misc.h"

/*
** For dap io
*/
static DapIO io;








int maxgel;


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







void xdap_middle_open_for_read(List *l)
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

    dap_open_for_read(&io, name, version);

    cur_gel_index = 1;
    cur_contig_index = 1;
}


List *xdap_middle_read_header(void)
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








List *xdap_middle_read_gel_data(void)
{
    List *gel_details;

    if (cur_gel_index > io.num_gels)
	gel_details = nil;
    else {
	dap_ar_file_rec ar_line;
	dap_rl_file_rec rl_line;
	dap_sq_file_rec sq_line;
	dap_tg_file_rec tg_line;

	int index;
	char name[13];
	char *seq;
	int length;
	int comp;
	int pos;
	int l_nbr;
	int r_nbr;

	sq_line = (char *) malloc(io.max_gel_length+1);

	dap_read_ar(&io,cur_gel_index,&ar_line);
	dap_read_rl(&io,cur_gel_index,&rl_line);
	dap_read_tg(&io,cur_gel_index,&tg_line);
	dap_read_sq(&io,cur_gel_index,sq_line);


	index = cur_gel_index;
	f2cstr(ar_line.lines.name,10,name,12);
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
	    
	    rd = dap_read_comment(&io, tg_line.lines.comment);
	    sscanf(rd,"%6d%6d%6d%*s",&rd_length, &rd_cut, &rd_ulen);
	    f2cstr(&rd[18],4,rd_type,4);
	    f2cstr(&rd[22],18,rd_file,18);


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
		dap_read_tg(&io,next,&tg_line);

		if (strncmp(tg_line.lines.type.c,"*LC*",4)==0) {
		    if (tg_line.lines.comment) {
			List *lc;
			lc = build_list(
					atom_str(gel_l_cut_seq),
					atom_str(dap_read_comment(&io,tg_line.lines.comment)),
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
					atom_str(dap_read_comment(&io,tg_line.lines.comment)),
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
		    com = dap_read_comment(&io,tg_line.lines.comment);
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





List *xdap_middle_read_contig_data(void)
{
    List *contig_details;

    if (cur_contig_index > io.num_contigs)
	contig_details = nil;
    else {
	dap_rl_file_rec rl_line;
	int length;
	int index;
	int left_end;
	int right_end;

	index = io.max_gels-cur_contig_index;
	dap_read_rl(&io,index,&rl_line);
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


void xdap_middle_close(List *l)
/*
** Close all relevant files
*/
{
    dap_close_files(&io);
}
