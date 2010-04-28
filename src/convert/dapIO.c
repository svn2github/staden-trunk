#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "dapIO.h"
#include "misc.h"

/*
** Tag File IO
*/
void dap_read_tg(DapIO *io, int rec, dap_tg_file_rec *t)
{
    FILE *f = io->tg_fp;
    if ( fseeko(f,(off_t)dap_tg_byte_index(io,rec),0) )
	crash("Seek failure on tag file, record %d\n",rec);

    if ( fread(t, sizeof(dap_tg_file_rec), 1, f) != 1)
	crash("Read failure on tag file\n");
}

void dap_write_tg(DapIO *io, int rec, dap_tg_file_rec *t)
{
    FILE *f = io->tg_fp;
    if ( fseeko(f,(off_t)dap_tg_byte_index(io,rec),0) )
	crash("Seek failure on tag file, record %d\n",rec);

    if ( fwrite(t, sizeof(dap_tg_file_rec), 1, f) != 1)
	crash("Write failure on tag file\n");
}



/*
** Archive File IO
*/
void dap_read_ar(DapIO *io, int rec, dap_ar_file_rec *t)
{
    FILE *f = io->ar_fp;
    if ( fseeko(f,(off_t)dap_ar_byte_index(io,rec),0) )
	crash("Seek failure on archive file, record %d\n",rec);

    if ( fread(t, sizeof(dap_ar_file_rec), 1, f) != 1)
	crash("Read failure on archive file\n");
}

void dap_write_ar(DapIO *io, int rec, dap_ar_file_rec *t)
{
    FILE *f = io->ar_fp;
    if ( fseeko(f,(off_t)dap_ar_byte_index(io,rec),0) )
	crash("Seek failure on archive file, record %d\n",rec);

    if ( fwrite(t, sizeof(dap_ar_file_rec), 1, f) != 1)
	crash("Write failure on archive file\n");
}




/*
** Relationship file IO
*/
void dap_read_rl(DapIO *io, int rec, dap_rl_file_rec *t)
{
    FILE *f = io->rl_fp;
    if ( fseeko(f,(off_t)dap_rl_byte_index(io,rec),0) )
	crash("Seek failure on relationships file, record %d\n",rec);

    if ( fread(t, sizeof(dap_rl_file_rec), 1, f) != 1)
	crash("Read failure on relationships file\n");
}

void dap_write_rl(DapIO *io, int rec, dap_rl_file_rec *t)
{
    FILE *f = io->rl_fp;
    if ( fseeko(f,(off_t)dap_rl_byte_index(io,rec),0) )
	crash("Seek failure on relationships file, record %d\n",rec);

    if ( fwrite(t, sizeof(dap_rl_file_rec), 1, f) != 1)
	crash("Write failure on relationships file\n");
}





/*
** Comment file IO
*/
void dap_read_cc(DapIO *io, int rec, dap_cc_file_rec *t)
{
    FILE *f = io->cc_fp;
    if ( fseeko(f,(off_t)dap_cc_byte_index(io,rec),0) )
	crash("Seek failure on comment file, record %d\n",rec);

    if ( fread(t, sizeof(dap_cc_file_rec), 1, f) != 1)
	crash("Read failure on comment file, record %d\n",rec);
}

void dap_write_cc(DapIO *io, int rec, dap_cc_file_rec *t)
{
    FILE *f = io->cc_fp;
    if ( fseeko(f,(off_t)dap_cc_byte_index(io,rec),0) )
	crash("Seek failure on comment file, record %d\n",rec);

    if ( fwrite(t, sizeof(dap_cc_file_rec), 1, f) != 1)
	crash("Crash failure on comment file, record %d\n",rec);
}



/*
** Sequence file IO
*/
void dap_read_sq(DapIO *io, int rec, dap_sq_file_rec t)
{
    FILE *f = io->sq_fp;
    if ( fseeko(f,(off_t)dap_sq_byte_index(io,rec),0) )
	crash("Seek failure on sequence file, record %d\n",rec);

    if ( fread(t, io->max_gel_length, 1, f) != 1)
	crash("Read failure on sequence file\n");
}

void dap_write_sq(DapIO *io, int rec, dap_sq_file_rec t)
{
    FILE *f = io->sq_fp;
    if ( fseeko(f,(off_t)dap_sq_byte_index(io,rec),0) )
	crash("Seek failure on sequence file, record %d\n",rec);

    if ( fwrite(t, io->max_gel_length, 1, f) != 1)
	crash("Write failure on sequence file\n");
}



/*
** Comment IO - Strings
*/
char *dap_read_comment(DapIO *io, int_4 cp)
{
    dap_cc_file_rec c;
    int count;
    int_4 nc;
    char *com,*comptr;

    if (!cp) return NULL;
    /* determine how long string is */
    count = 1;
    nc=cp;
    dap_read_cc(io, nc, &c);
    while (c.lines.next != 0) {
	nc = c.lines.next;
	count++;
        dap_read_cc(io, nc, &c);
    }

    com = comptr = (char *)malloc(count * DAP_COMMENT_SIZE+1);
    nc=cp;
    dap_read_cc(io, nc, &c);
    strncpy(com,c.lines.comment,DAP_COMMENT_SIZE); com+=DAP_COMMENT_SIZE;
    while (c.lines.next != 0) {
	nc = c.lines.next;
	count++;
        dap_read_cc(io, nc, &c);
        strncpy(com,c.lines.comment,DAP_COMMENT_SIZE); com+=DAP_COMMENT_SIZE;
    }

    *com = '\0';

    return comptr;
    
}

#ifdef nodef
static int_4 get_free_comment(DapIO *io)
{
    dap_cc_file_rec head;
    dap_cc_file_rec freerec;
    int_4 free_id;
    dap_read_cc(io,dap_cc_header_rec(io),&head);
    if (head.header.free_list != 0) {
	/*
	** if a free slot somewhere, use it
	*/
	free_id = head.header.free_list;
	dap_read_cc(io,free_id,&freerec);
	head.header.free_list = freerec.lines.next;
	dap_write_cc(io,dap_cc_header_rec(io),&head);
    } else {
	/*
	** extend comment list file
	*/
	free_id = ++head.header.count;
	dap_write_cc(io,dap_cc_header_rec(io),&head);
	dap_write_cc(io,free_id,&freerec);
    }

    return free_id;
}

static int_4 get_free_tag(DapIO *io)
{
    dap_tg_file_rec head;
    dap_tg_file_rec freerec;
    int_4 free_id;
    dap_read_tg(io,dap_tg_header_rec(io),&head);
    if (head.header.free_list != 0) {
	/*
	** if a free slot somewhere, use it
	*/
	free_id = head.header.free_list;
	dap_read_tg(io,free_id,&freerec);
	head.header.free_list = freerec.lines.next;
	dap_write_tg(io,dap_tg_header_rec(io),&head);
    } else {
	/*
	** extend comment list file
	*/
	free_id = ++head.header.count;
	dap_write_tg(io,dap_tg_header_rec(io),&head);
	dap_write_tg(io,free_id,&freerec);
    }

    return free_id;
}





static void insert_tag(DapIO *io, int_4 gel, dap_tg_file_rec t)
{
    int_4 next, last;
    int_4 free;
    dap_tg_file_rec tg,last_tg;

    last = gel;
    dap_read_tg(io,last,&last_tg);

    next = last_tg.lines.next;
    if (next) dap_read_tg(io,next,&tg);

    while (next && tg.lines.position <= t.lines.position) {
	last = next;
	last_tg = tg;
	next = tg.lines.next;
	if (next) dap_read_tg(io,next,&tg);
    }

    /* insert after last */
    free = get_free_tag(io);
    t.lines.next = next;
    last_tg.lines.next = free;
    dap_write_tg(io,last,&last_tg);
    dap_write_tg(io,free,&t);

}


static int_4 write_comment(DapIO *io, char *c)
{
    dap_cc_file_rec com;
    int_4 cur,next,this_comment;
    int clen = strlen(c);
    int piece;

    /* write out first block of DAP_COMMENT_SIZE */
    this_comment=cur=get_free_comment(io);
    if (clen>DAP_COMMENT_SIZE)
	piece = DAP_COMMENT_SIZE;
    else
	piece = clen;

    {int i; for(i=0;i<DAP_COMMENT_SIZE;i++)com.lines.comment[i]=' ';}
    strncpy(com.lines.comment,c,piece);

    c+= piece;
    clen -= piece;
    while (clen > 0) {
	next = get_free_comment(io);
	com.lines.next = next;
	dap_write_cc(io,cur,&com);
	cur = next;
	if (clen<DAP_COMMENT_SIZE)
	    piece = clen;

	{int i; for(i=0;i<DAP_COMMENT_SIZE;i++)com.lines.comment[i]=' ';}
	strncpy(com.lines.comment,c,piece);

	c+= piece;
	clen -= piece;
    }
    com.lines.next = 0;
    if (piece!=DAP_COMMENT_SIZE)
	com.lines.comment[piece]='\0';
    dap_write_cc(io,cur,&com);

    return this_comment;
}
#endif /*nodef*/

static void set_file_names(DapIO *io, char *name, char *version)
{
    strcpy(io->ar_file,name); strcat(io->ar_file,".AR"); strcat(io->ar_file,version);
    strcpy(io->rl_file,name); strcat(io->rl_file,".RL"); strcat(io->rl_file,version);
    strcpy(io->sq_file,name); strcat(io->sq_file,".SQ"); strcat(io->sq_file,version);
    strcpy(io->tg_file,name); strcat(io->tg_file,".TG"); strcat(io->tg_file,version);
    strcpy(io->cc_file,name); strcat(io->cc_file,".CC"); strcat(io->cc_file,version);
}

static void dap_open_files(DapIO *io, char *name, char *version, char *mode)
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
	crash("No archive file %s\n",io->ar_file);
    if ( ( io->rl_fp = fopen(io->rl_file,mode) ) == NULL )
	crash("No relationships file %s\n",io->rl_file);
    if ( ( io->sq_fp = fopen(io->sq_file,mode) ) == NULL )
	crash("No sequence file %s\n",io->sq_file);
    if ( ( io->tg_fp = fopen(io->tg_file,mode) ) == NULL )
	crash("No tag file %s\n",io->tg_file);
    if ( ( io->cc_fp = fopen(io->cc_file,mode) ) == NULL )
	crash("No tag-comment file %s\n",io->cc_file);

}


void dap_open_for_read(DapIO *io, char *name, char *version)
{
    dap_ar_file_rec ar_header;
    dap_rl_file_rec rl_header;

    dap_open_files(io,name,version,"rb");
    
    dap_read_ar(io,dap_ar_header_rec(io),&ar_header);
    io->max_gels = ar_header.header.idbsiz;
    io->max_gel_length = ar_header.header.maxgel;
    io->data_class = ar_header.header.idm;

    dap_read_rl(io,dap_rl_header_rec(io),&rl_header);
    io->num_gels = rl_header.header.num_gels;
    io->num_contigs = rl_header.header.num_contigs;
}

void dap_open_for_write(DapIO *io, char *name, char *version)
{
    dap_open_files(io,name,version,"r+b");
}



void dap_close_files(DapIO *io)
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

