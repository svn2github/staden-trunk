#ifndef _dapIO_h
#define _dapIO_h

#include "newtypes.h"

/*
** Definition of dap database files
*/



/*
** Archive file (*.AR?)
*/
#define DAP_FILE_NAME_LENGTH 12
typedef union {

    struct _ar_header {
	int_4 idbsiz;
	int_4 maxgel;
	int_4 idm;
    } header;

    struct _ar_lines {
	char name[DAP_FILE_NAME_LENGTH];
    } lines;

} dap_ar_file_rec;

#define dap_ar_header_rec(I)  ( 1000 )
#define dap_ar_byte_index(I,R) ( ((R)-1) * sizeof(dap_ar_file_rec) )







/*
** Relationships file (*.RL?)
*/
typedef union {

    struct _rl_header { 
	int_4 num_gels;
	int_4 num_contigs;
	int_4 spare1;
	int_4 spare2;
    } header;

    struct _rl_lines {
	int_4 rel_pos;
	int_4 length;
	int_4 left_nbr;
	int_4 right_nbr;
    } lines;

    struct _rl_clines {
	int_4 length;
	int_4 spare3;
	int_4 left_end;
	int_4 right_end;
    } clines;

} dap_rl_file_rec;

#define dap_rl_header_rec(I) ( (I)->max_gels )
#define dap_rl_byte_index(I,R) ( ((R)-1) * sizeof(dap_rl_file_rec) )






/*
** Sequence file (*.SQ?)
*/
typedef char *dap_sq_file_rec;
#define dap_sq_byte_index(I,R) ( ((R)-1) * ((I)->max_gel_length) )



/*
** Tag files (*.TG?)
*/
typedef union {
    int i;
    char c[4];
} dap_tag_type;



typedef union {

    struct _tg_header {
	int_4 count;
	int_4 spare1;
	int_4 spare2;
	dap_tag_type spare3;
	int_4 free_list;
    } header;

    struct _tg_lines {
	int_4 position;
	int_4 length;
	int_4 comment;
	dap_tag_type type;
	int_4 next;
    } lines;

} dap_tg_file_rec;

#define dap_tg_byte_index(I,R) ( ((R)-1) *sizeof(dap_tg_file_rec) )
#define dap_tg_header_rec(I) ( (I)->max_gels )

/*
** Comment files (*.CC?)
*/
#define DAP_COMMENT_SIZE 40
typedef union {

    struct _cc_header {
	int_4 free_list;
	int_4 count;
	char spare[DAP_COMMENT_SIZE - sizeof(int_4)];
    } header;

    struct _cc_lines {
	int_4 next;
	char comment[DAP_COMMENT_SIZE];
    } lines;

} dap_cc_file_rec;

#define dap_cc_byte_index(I,R) ( ((R)-1) * sizeof(dap_cc_file_rec) )
#define dap_cc_header_rec(I) (1)


/*
** Useful variables
*/
typedef char IOString[200];

typedef struct {
    int max_gels;
    int num_gels;
    int num_contigs;
    int max_gel_length;
    int data_class;
    
    FILE *ar_fp;
    FILE *rl_fp;
    FILE *sq_fp;
    FILE *tg_fp;
    FILE *cc_fp;

    IOString ar_file;
    IOString rl_file;
    IOString sq_file;
    IOString tg_file;
    IOString cc_file;
} DapIO;





extern void dap_read_tg(DapIO *io, int rec, dap_tg_file_rec *t);
extern void dap_write_tg(DapIO *io, int rec, dap_tg_file_rec *t);
extern void dap_read_ar(DapIO *io, int rec, dap_ar_file_rec *t);
extern void dap_write_ar(DapIO *io, int rec, dap_ar_file_rec *t);
extern void dap_read_rl(DapIO *io, int rec, dap_rl_file_rec *t);
extern void dap_write_rl(DapIO *io, int rec, dap_rl_file_rec *t);
extern void dap_read_cc(DapIO *io, int rec, dap_cc_file_rec *t);
extern void dap_write_cc(DapIO *io, int rec, dap_cc_file_rec *t);
extern void dap_read_sq(DapIO *io, int rec, dap_sq_file_rec t);
extern void dap_write_sq(DapIO *io, int rec, dap_sq_file_rec t);
extern char *dap_read_comment(DapIO *io, int_4 cp);
extern void dap_open_for_read(DapIO *io, char *name, char *version);
extern void dap_open_for_write(DapIO *io, char *name, char *version);
extern void dap_close_files(DapIO *io);


#endif /* _dapIO_h */
