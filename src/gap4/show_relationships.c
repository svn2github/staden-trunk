#include <stdio.h>
#include <tcl.h>
#include <tk.h>

#include "io_utils.h"
#include "cli_arg.h"
#include "newgap_structs.h"
#include "misc.h"
#include "text_output.h"

int show_relationships(GapIO *io,
		       contig_list_t *contigs,
		       int num_contigs, 
		       int ordered)
{
    GReadings r;
    char *name;
    int i, j;
    int clen, cleft, cright;
    
    char *contig_line = " CONTIG LINES \n"
	" CONTIG            NUMBER   LENGTH                ENDS \n"
	    "                                              LEFT    RIGHT\n";
    
    char *gel_line =" GEL LINES \n"
	" %-*s   NUMBER POSITION LENGTH      NEIGHBOURS\n"
	    " %-*s                              LEFT    RIGHT\n";
    if (num_contigs == io->db.num_contigs) {
	/* all contigs */
	if (ordered == 1) {
	    for (i = 0; i < num_contigs; i++){ 
		clen   = io_clength(io, contigs[i].contig);
		cleft  = io_clnbr(io, contigs[i].contig);
		cright = io_crnbr(io, contigs[i].contig);
		
		/* print contig information */
		vmessage(contig_line);
		vmessage("%25d %8d %15d %8d\n", contigs[i].contig,
			 clen, cleft, cright);
	    
		/* print gel information */
		vmessage(gel_line, DB_NAMELEN, "NAME", DB_NAMELEN, "");
		
		for (j = cleft; j; j = r.right) {
		    
		    gel_read(io, j, r);  
		    name = io_rname(io, j);
		    
		    if (((r.position + r.sequence_length - 1) >= 
			 contigs[i].start) && (r.position <= contigs[i].end)) {

			if (r.sense)
			    r.sequence_length *= -1;
		    
			vmessage_tagged("SEQID", "%-*s",
					DB_NAMELEN+1, name);
			vmessage(" %8d %8d %6d %8d %8d\n",
				 j, r.position, r.sequence_length, 
				 r.left, r.right);
		    }
		    
		} /* end for */
		
	    } /* end for */

	} else if (ordered == 0) {
	
	    /* print contig information */
	    vmessage(contig_line);
	
	    for (i = 0; i < num_contigs; i++){ 
		clen   = io_clength(io, contigs[i].contig);
		cleft  = io_clnbr(io, contigs[i].contig);
		cright = io_crnbr(io, contigs[i].contig);
		
		vmessage("%25d %8d %15d %8d\n", contigs[i].contig,
			 clen, cleft, cright);
		
	    } /* end for */
	    
	    /* print gel information */
	    vmessage(gel_line, DB_NAMELEN, "NAME", DB_NAMELEN, "");
	    
	    for (j = 1; j <= NumReadings(io); j++){
		
		gel_read(io, j, r);
		name = io_rname(io, j);
		
		if (r.sense)
		    r.sequence_length *= -1;
		
		vmessage_tagged("SEQID", "%-*s",
				DB_NAMELEN+1, name);
		vmessage(" %8d %8d %6d %8d %8d\n", 
			 j, r.position, r.sequence_length, 
			 r.left, r.right);
		
	    } /* end for */
	}

    } else {
	/*  single or subset of contigs */
	for (i = 0; i < num_contigs; i++) {
	    clen   = io_clength(io, contigs[i].contig);
	    cleft  = io_clnbr(io, contigs[i].contig);
	    cright = io_crnbr(io, contigs[i].contig);
	    
	    /* print contig information */
	    vmessage(contig_line);
	    vmessage("%25d %8d %15d %8d\n",contigs[i].contig,
		     clen, cleft, cright);
	    
	    /* print gel information */
	    vmessage(gel_line, DB_NAMELEN, "NAME", DB_NAMELEN, "");
	    
	    for (j=cleft; j; j=r.right) {
		gel_read(io, j, r);  
	    
		if (((r.position + r.sequence_length - 1) >= 
		     contigs[i].start) && (r.position <= contigs[i].end)) {
		    
		    name = io_rname(io, j);
		    
		    if (r.sense)
			r.sequence_length *= -1;
		    
		    vmessage_tagged("SEQID", "%-*s",
				    DB_NAMELEN+1, name);
		    vmessage(" %8d %8d %6d %8d %8d\n", 
			     j, r.position, r.sequence_length, 
			     r.left, r.right);
		} /* end if */
		
	    } /* end for */
	} /* end if */
    }
    
    return TCL_OK;
}
