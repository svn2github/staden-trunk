#if !defined(FEATURE_TABLE_H)
#define FEATURE_TABLE_H

#include "parse_feature.h"
#include "read_sequence.h"

#define num_qual 70

FEATURE_TABLE *init_feature_table ( void );
FEATURE_TABLE *realloc_feature_table ( FEATURE_TABLE *ft, int num_entry );
void free_feature_table (FEATURE_TABLE *ft);
ft_range *copy_ft_range (ft_range *r_sou);
ft_entry *copy_ft_entry ( ft_entry *entry);
FEATURE_TABLE *copy_feature_table ( FEATURE_TABLE *ft );
int extend_deleted_feature_table (FEATURE_TABLE *ft, ft_entry *entry, int entry_id);
int compare_ft_single_range (ft_range *r1, ft_range *r2);
void set_single_range_null (ft_range *r);
void copy_single_ft_range (ft_range *r_des, ft_range *r_sou);
void change_single_ft_range (ft_range *r, int string_len);
void change_single_ft_range_insert (ft_range *r, int string_len);
int check_ft_range (ft_range *r);
int add_to_deleted_feature_table (FEATURE_TABLE *ft, ft_entry *e, ft_entry *e_start, int num_del);
void mv_right_to_left (ft_range *r);
FEATURE_TABLE *delete_string_modify_feature_table (SEQUENCE *sequence, int string_len, int position);
void insert_string_extend_ft_range (SEQUENCE *sequence, char *string, int position );
void insert_string_break_ft_range (SEQUENCE *sequence, char *string, int position);
FEATURE_TABLE *insert_string_delete_ft_range (SEQUENCE *sequence, int position);

int add_feature_to_feature_table (SEQUENCE *sequence, ft_entry *entry);
int compare_ft_entry (ft_entry *e1, ft_entry *e2);
int insert_sub_ft_in_feature_table ( SEQUENCE *seq_des, FEATURE_TABLE *ft_sou);
void change_feature_table_location (FEATURE_TABLE *ft, int position);
FEATURE_TABLE *get_feature_table (SEQUENCE *sequence, int string_len, int position);
int add_ft_in_sequence (SEQUENCE *s, FEATURE_TABLE *ft);

#endif
