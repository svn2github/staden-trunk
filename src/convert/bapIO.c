#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "bapIO.h"
#include "misc.h"
#include "os.h"


/*
** Tag File IO
*/

void swap_tg_file_rec(bap_tg_file_rec *src, bap_tg_file_rec *dst)
{
    int i=1;

    if (*(char*)&i) {
	swap_int4(src->lines.position,dst->lines.position);
	swap_int4(src->lines.length,dst->lines.length);
	swap_int4(src->lines.comment,dst->lines.comment);
	dst->lines.type = src->lines.type;
	swap_int4(src->lines.next,dst->lines.next);
    } else {
	dst->lines.position = src->lines.position;
	dst->lines.length = src->lines.length;
	dst->lines.comment = src->lines.comment;
	dst->lines.type = src->lines.type;
	dst->lines.next = src->lines.next;
    }
    
}




void bap_read_tg(BapIO *io, int_4 rec, bap_tg_file_rec *t)
{
    FILE *f = io->tg_fp;
    bap_tg_file_rec tbuf;

    if ( fseeko(f,(off_t)bap_tg_byte_index(io,rec),0) )
	crash("Seek failure on tag file, record %d\n",rec);

    if ( fread(&tbuf, sizeof(bap_tg_file_rec), 1, f) != 1)
	crash("Read failure on tag file\n");

    swap_tg_file_rec(&tbuf,t);
}




void bap_write_tg(BapIO *io, int rec, bap_tg_file_rec *t)
{
    FILE *f = io->tg_fp;
    bap_tg_file_rec tbuf;

    if ( fseeko(f,(off_t)bap_tg_byte_index(io,rec),0) )
	crash("Seek failure on tag file, record %d\n",rec);


    swap_tg_file_rec(t,&tbuf);

    if ( fwrite(&tbuf, sizeof(bap_tg_file_rec), 1, f) != 1)
	crash("Write failure on tag file\n");
}



/*
** Archive File IO
*/
void bap_read_ar(BapIO *io, int rec, bap_ar_file_rec *t)
{
    FILE *f = io->ar_fp;
    if ( fseeko(f,(off_t)bap_ar_byte_index(io,rec),0) )
	crash("Seek failure on archive file, record %d\n",rec);

    if ( fread(t, sizeof(bap_ar_file_rec), 1, f) != 1)
	crash("Read failure on archive file\n");
}

void bap_write_ar(BapIO *io, int rec, bap_ar_file_rec *t)
{
    FILE *f = io->ar_fp;
    if ( fseeko(f,(off_t)bap_ar_byte_index(io,rec),0) )
	crash("Seek failure on archive file, record %d\n",rec);

    if ( fwrite(t, sizeof(bap_ar_file_rec), 1, f) != 1)
	crash("Write failure on archive file\n");
}




/*
** Relationship file IO
*/

void swap_rl_file_rec(bap_rl_file_rec *src, bap_rl_file_rec *dst)
{
    int i=1;

    if (*(char*)&i) {
	swap_int4(src->lines.rel_pos,dst->lines.rel_pos);
	swap_int4(src->lines.length,dst->lines.length);
	swap_int4(src->lines.left_nbr,dst->lines.left_nbr);
	swap_int4(src->lines.right_nbr,dst->lines.right_nbr);
    } else {
	dst->lines.rel_pos = src->lines.rel_pos;
	dst->lines.length = src->lines.length;
	dst->lines.left_nbr = src->lines.left_nbr;
	dst->lines.right_nbr = src->lines.right_nbr;
    }
    
}



void bap_read_rl(BapIO *io, int rec, bap_rl_file_rec *t)
{
    FILE *f = io->rl_fp;
    bap_rl_file_rec tbuf;

    if ( fseeko(f,(off_t)bap_rl_byte_index(io,rec),0) )
	crash("Seek failure on relationships file, record %d\n",rec);

    if ( fread(&tbuf, sizeof(bap_rl_file_rec), 1, f) != 1)
	crash("Read failure on relationships file\n");

    swap_rl_file_rec(&tbuf,t);
}

void bap_write_rl(BapIO *io, int rec, bap_rl_file_rec *t)
{
    FILE *f = io->rl_fp;
    bap_rl_file_rec tbuf;
    if ( fseeko(f,(off_t)bap_rl_byte_index(io,rec),0) )
	crash("Seek failure on relationships file, record %d\n",rec);

    swap_rl_file_rec(t,&tbuf);

    if ( fwrite(&tbuf, sizeof(bap_rl_file_rec), 1, f) != 1)
	crash("Write failure on relationships file\n");
}




/*
** Comment file IO
*/


void swap_cc_file_rec_header(bap_cc_file_rec *src, bap_cc_file_rec *dst)
{
    int i=1;
    if (*(char*)&i) {
	swap_int4(src->header.free_list,dst->header.free_list);
	swap_int4(src->header.count,dst->header.count);
    } else {
	dst->header.free_list = src->header.free_list;
	dst->header.count = src->header.count;
    }
    
}




void swap_cc_file_rec_lines(bap_cc_file_rec *src, bap_cc_file_rec *dst)
{
    int i=1;
    if (*(char*)&i) {
	swap_int4(src->lines.next,dst->lines.next);
    } else {
	dst->lines.next = src->lines.next;
    }
    memcpy(dst->lines.comment,src->lines.comment,BAP_COMMENT_SIZE);
    
}

void bap_read_cc(BapIO *io, int rec, bap_cc_file_rec *t)
{
    FILE *f = io->cc_fp;
    bap_cc_file_rec tbuf;
    if ( fseeko(f,(off_t)bap_cc_byte_index(io,rec),0) )
	crash("Seek failure on comment file, record %d\n",rec);

    if ( fread(&tbuf, sizeof(bap_cc_file_rec), 1, f) != 1)
	crash("Read failure on comment file, record %d\n",rec);

    if (rec == bap_cc_header_rec(io))
	swap_cc_file_rec_header(&tbuf,t);
    else
	swap_cc_file_rec_lines(&tbuf,t);

}



void bap_write_cc(BapIO *io, int rec, bap_cc_file_rec *t)
{
    FILE *f = io->cc_fp;
    bap_cc_file_rec tbuf;
    if ( fseeko(f,(off_t)bap_cc_byte_index(io,rec),0) )
	crash("Seek failure on comment file, record %d\n",rec);

    if (rec == bap_cc_header_rec(io))
	swap_cc_file_rec_header(t,&tbuf);
    else
	swap_cc_file_rec_lines(t,&tbuf);

    if ( fwrite(&tbuf, sizeof(bap_cc_file_rec), 1, f) != 1)
	crash("Crash failure on comment file, record %d\n",rec);
}



/*
** Sequence file IO
*/
void bap_read_sq(BapIO *io, int rec, bap_sq_file_rec t)
{
    FILE *f = io->sq_fp;
    if ( fseeko(f,(off_t)bap_sq_byte_index(io,rec),0) )
	crash("Seek failure on sequence file, record %d\n",rec);

    if ( fread(t, io->max_gel_length, 1, f) != 1)
	crash("Read failure on sequence file\n");
}

void bap_write_sq(BapIO *io, int rec, bap_sq_file_rec t)
{
    FILE *f = io->sq_fp;
    if ( fseeko(f,(off_t)bap_sq_byte_index(io,rec),0) )
	crash("Seek failure on sequence file, record %d\n",rec);

    if ( fwrite(t, io->max_gel_length, 1, f) != 1)
	crash("Write failure on sequence file\n");
}



/*
** Comment IO - Strings
*/
char *bap_read_comment(BapIO *io, int_4 cp)
{
    bap_cc_file_rec c;
    int count;
    int_4 nc;
    char *com,*comptr;

    if (!cp) return NULL;
    /* determine how long string is */
    count = 1;
    nc=cp;
    bap_read_cc(io, nc, &c);
    while (c.lines.next != 0) {
	nc = c.lines.next;
	count++;
        bap_read_cc(io, nc, &c);
    }

    com = comptr = (char *)malloc(count * BAP_COMMENT_SIZE+1);
    nc=cp;
    bap_read_cc(io, nc, &c);
    strncpy(com,c.lines.comment,BAP_COMMENT_SIZE); com+=BAP_COMMENT_SIZE;
    while (c.lines.next != 0) {
	nc = c.lines.next;
	count++;
        bap_read_cc(io, nc, &c);
        strncpy(com,c.lines.comment,BAP_COMMENT_SIZE); com+=BAP_COMMENT_SIZE;
    }

    *com = '\0';

    return comptr;
    
}



static int_4 get_free_comment(BapIO *io)
{
    bap_cc_file_rec head;
    bap_cc_file_rec freerec;
    int_4 free_id;
    bap_read_cc(io,bap_cc_header_rec(io),&head);
    if (head.header.free_list != 0) {
	/*
	** if a free slot somewhere, use it
	*/
	free_id = head.header.free_list;
	bap_read_cc(io,free_id,&freerec);
	head.header.free_list = freerec.lines.next;
	bap_write_cc(io,bap_cc_header_rec(io),&head);
    } else {
	/*
	** extend comment list file
	*/
	free_id = ++head.header.count;
	bap_write_cc(io,bap_cc_header_rec(io),&head);
	bap_write_cc(io,free_id,&freerec);
    }

    return free_id;
}

int_4 bap_get_free_tag(BapIO *io)
{
    bap_tg_file_rec head;
    bap_tg_file_rec freerec;
    int_4 free_id;
    bap_read_tg(io,bap_tg_header_rec(io),&head);
    if (head.header.free_list != 0) {
	/*
	** if a free slot somewhere, use it
	*/
	free_id = head.header.free_list;
	bap_read_tg(io,free_id,&freerec);
	head.header.free_list = freerec.lines.next;
	bap_write_tg(io,bap_tg_header_rec(io),&head);
    } else {
	/*
	** extend comment list file
	*/
	free_id = ++head.header.count;
	bap_write_tg(io,bap_tg_header_rec(io),&head);
	bap_write_tg(io,free_id,&freerec);
    }

    return free_id;
}





void bap_insert_tag(BapIO *io, int_4 gel, bap_tg_file_rec t)
{
    int_4 next, last;
    int_4 free;
    bap_tg_file_rec tg,last_tg;

    last = gel;
    bap_read_tg(io,last,&last_tg);

    next = last_tg.lines.next;
    if (next) bap_read_tg(io,next,&tg);

    while (next && tg.lines.position <= t.lines.position) {
	last = next;
	last_tg = tg;
	next = tg.lines.next;
	if (next) bap_read_tg(io,next,&tg);
    }

    /* insert after last */
    free = bap_get_free_tag(io);
    t.lines.next = next;
    last_tg.lines.next = free;
    bap_write_tg(io,last,&last_tg);
    bap_write_tg(io,free,&t);

}


int_4 bap_write_comment(BapIO *io, char *c)
{
    bap_cc_file_rec com;
    int_4 cur,next,this_comment;
    int clen = strlen(c);
    int piece;

    /* write out first block of BAP_COMMENT_SIZE */
    this_comment=cur=get_free_comment(io);
    if (clen>BAP_COMMENT_SIZE)
	piece = BAP_COMMENT_SIZE;
    else
	piece = clen;

    {int i; for(i=0;i<BAP_COMMENT_SIZE;i++)com.lines.comment[i]=' ';}
    strncpy(com.lines.comment,c,piece);

    c+= piece;
    clen -= piece;
    while (clen > 0) {
	next = get_free_comment(io);
	com.lines.next = next;
	bap_write_cc(io,cur,&com);
	cur = next;
	if (clen<BAP_COMMENT_SIZE)
	    piece = clen;

	{int i; for(i=0;i<BAP_COMMENT_SIZE;i++)com.lines.comment[i]=' ';}
	strncpy(com.lines.comment,c,piece);

	c+= piece;
	clen -= piece;
    }
    com.lines.next = 0;
    if (piece!=BAP_COMMENT_SIZE)
	com.lines.comment[piece]='\0';
    bap_write_cc(io,cur,&com);

    return this_comment;
}










static void set_file_names(BapIO *io, char *name, char *version)
{
    strcpy(io->ar_file,name); strcat(io->ar_file,".AR"); strcat(io->ar_file,version);
    strcpy(io->rl_file,name); strcat(io->rl_file,".RL"); strcat(io->rl_file,version);
    strcpy(io->sq_file,name); strcat(io->sq_file,".SQ"); strcat(io->sq_file,version);
    strcpy(io->tg_file,name); strcat(io->tg_file,".TG"); strcat(io->tg_file,version);
    strcpy(io->cc_file,name); strcat(io->cc_file,".CC"); strcat(io->cc_file,version);
}

static void bap_open_files(BapIO *io, char *name, char *version, char *mode)
/*
**
*/
{
    /*
    ** Create file names
    */
    set_file_names(io,name,version);

    /*
    ** Open files
    */
    if ( ( io->ar_fp = fopen(io->ar_file,mode) ) == NULL )
	crash("Error opening archive file %s\n",io->ar_file);
    if ( ( io->rl_fp = fopen(io->rl_file,mode) ) == NULL )
	crash("Error opening relationships file %s\n",io->rl_file);
    if ( ( io->sq_fp = fopen(io->sq_file,mode) ) == NULL )
	crash("Error opening sequence file %s\n",io->sq_file);
    if ( ( io->tg_fp = fopen(io->tg_file,mode) ) == NULL )
	crash("Error opening tag file %s\n",io->tg_file);
    if ( ( io->cc_fp = fopen(io->cc_file,mode) ) == NULL )
	crash("Error opening tag-comment file %s\n",io->cc_file);

}


void bap_open_for_read(BapIO *io, char *name, char *version)
{
    bap_rl_file_rec rl_header;

    bap_open_files(io,name,version,"rb");
    
    bap_read_rl(io,bap_rl_dbheader_rec(io),&rl_header);
    io->max_gels = rl_header.dbheader.idbsiz;
    io->max_gel_length = rl_header.dbheader.maxgel;
    io->data_class = rl_header.dbheader.idm;

    bap_read_rl(io,bap_rl_header_rec(io),&rl_header);
    io->num_gels = rl_header.header.num_gels;
    io->num_contigs = rl_header.header.num_contigs;
}

void bap_open_for_write(BapIO *io, char *name, char *version)
{
    bap_open_files(io,name,version,"w+b");
}

void bap_open_for_update(BapIO *io, char *name, char *version)
{
    bap_open_files(io,name,version,"r+b");
}



void bap_close_files(BapIO *io)
/*
** Close all relevant files
*/
{

    fclose(io->ar_fp);
    fclose(io->rl_fp);
    fclose(io->sq_fp);
    fclose(io->tg_fp);
    fclose(io->cc_fp);

}

