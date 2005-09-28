#ifndef _process_h
#define _process_h

#include "list.h"

#define db_from "From"
#define db_to "To"


#define db_files "Files"
#define db_files_seq "Sequence"
#define db_files_arch "Archive"
#define db_files_rel "Relationships"
#define db_files_raw "Raw-Data"
#define db_files_tag "Tag"
#define db_files_com "Comment"
#define db_files_flat "Flat-File"


#define db_rec "Database"
#define db_name "Name"
#define db_version "Version"
#define db_max_gels "Max-Gels"
#define db_max_db_size "Max-Database-Size"
#define db_max_gel_length "Max-Gel-Length"
#define db_max_contigs "Max-Contigs"
#define db_data_class "Data-Class"
#define db_num_gels "Num-Gels"
#define db_num_contigs "Num-Contigs"

#define db_type "Type"
#define db_type_RS_flat_file "RS-Flat-File"
#define db_type_SD_flat_file "SD-Flat-File"
#define db_type_sap "Sap"
#define db_type_early_xdap "Early-Xdap"
#define db_type_middle_xdap "Middle-Xdap"
#define db_type_late_xdap "Late-Xdap"
#define db_type_gap "Gap"
#define db_vtype_RS_flat_file "Flat file - created with sapf"
#define db_vtype_SD_flat_file "Flat file - created with this program"
#define db_vtype_sap "Sap database - the original format with three database files"
#define db_vtype_early_xdap "Early xdap database - with the raw data (RD) file"
#define db_vtype_middle_xdap "xdap database"
#define db_vtype_late_xdap "xbap database"
#define db_vtype_gap "xgap database"

#define contig_rec "Contigs"
#define contig_index "Contig-Index"
#define contig_length "Length"
#define contig_left_end "Left-End"
#define contig_right_end "Right-End"

#define gel_rec "Gels"
#define gel_index "Gel-Index"
#define gel_name "Name"
#define gel_seq "Sequence"
#define gel_pos "Pos-In-Contig"
#define gel_length "Length"
#define gel_comp "Complemented"
#define gel_l_nbr "Left-Nbr"
#define gel_r_nbr "Right-Nbr"
#define gel_rd_length "RD-Length"
#define gel_rd_cut "RD-Cut-Off-Position"
#define gel_rd_ulen "RD-Usable-Length"
#define gel_rd_type "Trace-File-Type"
#define gel_rd_file "Trace-File-Name"
#define gel_l_cut_seq "Left-Cutoff"
#define gel_r_cut_seq "Right-Cutoff"
#define gel_annotation "Annotation"
#define gel_an_pos "Position"
#define gel_an_len "Length"
#define gel_an_type "Type"
#define gel_an_comment "Comment"
#define gel_edits "Edits"
#define gel_ed_pos "Position"
#define gel_ed_type "Type"
#define gel_ed_char "Character"
#define gel_ed_delete "Delete"
#define gel_ed_insert "Insert"
#define gel_ed_op "Op"
#define gel_ed_base "Base"
#define gel_ed_base_pos "Pos"

#define tag_rec "TagTypes"
#define tag_name "Name"
#define tag_type "Type"
#define tag_fg "Foreground-Colour"
#define tag_bg "Background-Colour"
#define tag_dt "Default-Comment"



extern void write_header(List *to, List *l);
extern void write_gel_data(List *to, List *l);
extern void write_contig_data(List *to, List *l);
extern void process(List *from, List *to);

#endif /* _process_h */
