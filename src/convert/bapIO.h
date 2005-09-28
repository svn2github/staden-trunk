#ifndef _bapIO_h
#define _bapIO_h

#include "newtypes.h"

/*
** Definition of dap database files
*/



/*
** Archive file (*.AR?)
*/
#define BAP_FILE_NAME_LENGTH 16
typedef union {

    struct _ar_lines {
	char name[BAP_FILE_NAME_LENGTH];
    } lines;

} bap_ar_file_rec;

#define bap_ar_byte_index(I,R) ( ((R)-1) *sizeof(bap_ar_file_rec) )







/*
** Relationships file (*.RL?)
*/
typedef union {

    struct _rl_dbheader { 
	int_4 maxdb;
	int_4 idbsiz;
	int_4 maxgel;
	int_4 idm;
    } dbheader;

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

} bap_rl_file_rec;

#define bap_rl_header_rec(I) ( (I)->max_gels )
#define bap_rl_byte_index(I,R) ( (R) * sizeof(bap_rl_file_rec) )
#define bap_rl_dbheader_rec(I) (0)






/*
** Sequence file (*.SQ?)
*/
typedef char *bap_sq_file_rec;
#define bap_sq_byte_index(I,R) ( ((R)-1) * ((I)->max_gel_length) )



/*
** Tag files (*.TG?)
*/
typedef union {
    int i;
    char c[4];
} bap_tag_type;



typedef union {

    struct _tg_header {
	int_4 count;
	int_4 spare1;
	int_4 spare2;
	bap_tag_type spare3;
	int_4 free_list;
    } header;

    struct _tg_lines {
	int_4 position;
	int_4 length;
	int_4 comment;
	bap_tag_type type;
	int_4 next;
    } lines;

} bap_tg_file_rec;

#define bap_tg_byte_index(I,R) ( ((R)-1) *sizeof(bap_tg_file_rec) )
#define bap_tg_header_rec(I) ( (I)->max_gels )

/*
** Comment files (*.CC?)
*/
#define BAP_COMMENT_SIZE 40
typedef union {

    struct _cc_header {
	int_4 free_list;
	int_4 count;
	char spare[BAP_COMMENT_SIZE - sizeof(int_4)];
    } header;

    struct _cc_lines {
	int_4 next;
	char comment[BAP_COMMENT_SIZE];
    } lines;

} bap_cc_file_rec;

#define bap_cc_byte_index(I,R) ( ((R)-1) *sizeof(bap_cc_file_rec) )
#define bap_cc_header_rec(I) (1)



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
    int max_db_size;
    
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
} BapIO;





extern void bap_read_tg(BapIO *io, int rec, bap_tg_file_rec *t);
extern void bap_write_tg(BapIO *io, int rec, bap_tg_file_rec *t);
extern void bap_read_ar(BapIO *io, int rec, bap_ar_file_rec *t);
extern void bap_write_ar(BapIO *io, int rec, bap_ar_file_rec *t);
extern void bap_read_rl(BapIO *io, int rec, bap_rl_file_rec *t);
extern void bap_write_rl(BapIO *io, int rec, bap_rl_file_rec *t);
extern void bap_read_cc(BapIO *io, int rec, bap_cc_file_rec *t);
extern void bap_write_cc(BapIO *io, int rec, bap_cc_file_rec *t);
extern void bap_read_sq(BapIO *io, int rec, bap_sq_file_rec t);
extern void bap_write_sq(BapIO *io, int rec, bap_sq_file_rec t);
extern char *bap_read_comment(BapIO *io, int_4 cp);
extern int_4 bap_write_comment(BapIO *io, char *c);
extern int_4 bap_get_free_tag(BapIO *io);
extern void bap_insert_tag(BapIO *io, int_4 gel, bap_tg_file_rec t);
extern void bap_open_for_read(BapIO *io, char *name, char *version);
extern void bap_open_for_write(BapIO *io, char *name, char *version);
extern void bap_close_files(BapIO *io);


#endif /* _bapIO_h */
