/*
 * File: gap-create.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: routine to create a new gap database
 *
 * Created: prior to 18 September 1992
 * Updated:
 *
 */

#include <staden_config.h>

#include <stdio.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h> /* import: creat */
#include <unistd.h> /* import: access */
#include <string.h> /* import: strcat */
#include <errno.h>

#if 0
/* jkb - why do we need DLL_IMPORT anyway? we're using errno.h so who cares. */

/* 6/1/99 johnt - need to explicitly import golbals with Visula C++ */
#ifdef _MSC_VER
# define DLL_IMPORT __declspec(dllimport)
#else
# define DLL_IMPORT
#endif

extern DLL_IMPORT int errno;
#endif

/*
 * Windows requires O_BINARY for open()ing binary files.
 * However Unix doesn't have this, so we'll define it to be 0 (ie no effect).
 */
#ifndef O_BINARY
#   define O_BINARY 0
#endif


#include "bitmap.h"

#include "g-os.h"
#include "g-io.h"
#include "g-filedefs.h"
#include "g-error.h"
#include "g-defs.h"
#include "g-misc.h" /* IMPORT: G_Number */


#include "gap-error.h"
#include "gap-if.h"
#include "gap-create.h"
#include "gap-dbstruct.h"
#include "gap-init.h"
#include "gap-io.h"
#include "gap-error.h"
#include "gap-defaults.h"

#include "FtoC.h"

int maxdb = 8000;
static int bitsize = G_32BIT;

static int create_new_files(char *fn, size_t s, GCardinal max_rec)
/*
 * create a new file header
 * returns:    0 - ok
 *          else - error
 */
{
    char auxfn[1024];
    int fd;
    AuxHeader auxheader;
    GCardinal i;
    int (*(*low_level_vector))(int fd, void *x, int num);
    int endian = 1;

    /*
     * Determine the default vectors for creating a database.
     * This is temporary as from here on the g-library will auto-sense when
     * it opens the database files.
     */
    if ( *(char *)&endian ) {
	low_level_vector = (bitsize == G_64BIT)
	    ? low_level_vectors_swapped64
	    : low_level_vectors_swapped32;
    } else {
	low_level_vector = (bitsize == G_64BIT)
	    ? low_level_vectors64
	    : low_level_vectors32;
    }

    /* check file name isn't too long */
    if ( strlen(fn) + strlen(G_AUX_SUFFIX) >= sizeof(auxfn) ) return gerr_set(GERR_NAME_TOO_LONG);
	
    strcpy(auxfn,fn);
    strcat(auxfn,G_AUX_SUFFIX);

    /* check files don't already exist */
    if (access(fn,F_OK)==0 || errno != ENOENT) return gerr_set(GERR_FILE_EXISTS);
    if (access(auxfn,F_OK)==0 || errno != ENOENT) return gerr_set(GERR_FILE_EXISTS);

    /* create files */
    /* LOW LEVEL IO HERE */
    if ( (fd = creat(fn,G_DEF_PERMS)) == -1 ) return gerr_set(GERR_CANT_CREATE);
    /* LOW LEVEL IO HERE */
    close(fd);

    /* LOW LEVEL IO HERE */
    if ( (fd = creat(auxfn,G_DEF_PERMS)) == -1 ) return gerr_set(GERR_CANT_CREATE);

    /* initialise header */
    auxheader.file_size = 0;
    auxheader.block_size = (int4) s;
    auxheader.num_records = 0;
    auxheader.max_records = max_rec;
    auxheader.last_time = G_YEAR_DOT;
    auxheader.flags = (GFlags) 0;
    auxheader.spare1 = (GFlags) 0;
    auxheader.free_time = G_YEAR_DOT;
    for (i=G_Number(auxheader.spare)-1;i>=0;i--) auxheader.spare[i]=0;
    auxheader.format = bitsize;

    /* write(fd,&auxheader,sizeof(auxheader)); */
    (void) (low_level_vector[GOP_WRITE_AUX_HEADER])(fd,&auxheader,1);

    /* LOW LEVEL IO HERE */
    close(fd);
    return 0;
}


static int gap_create_db(char *project,char *version)
/*
 *
 */
{
    int err;
    GCardinal i;

    for (i=0;i<GAP_FILES;i++) {
	char *file = gap_construct_file(project,file_list[i],version,NULL);
	if (file==NULL) return 1;
	err = create_new_files(file,block_sizes[i],max_recs[i]);
        if (err) return err;
    }

    return gerr_set(GERR_NONE);
}





static int gap_init_db(char *project,char *version,int read_only)
/*
 *
 */
{
    int err;
    GapServer *s;
    GapClient *client;

    GDatabase db;
    GVectors vec;
    GClones clone;
    GTemplates temp;
    GCardinal VEC = GR_Vectors_Default;
    GCardinal TEMP = GR_Templates_Default;
    GCardinal CLONE = GR_Clones_Default;

    BASE_TYPE BITS = ~(BASE_TYPE)0; /* records 0 to BIT_ELE-1 */
    int i;
    GView view[BIT_ELE];


    /*
     * open it
     */
    if ( (s = gap_open_server(project,version,read_only)) == NULL ) {
	GAP_ERROR("cannot open database");
	return(1);
    }

    /*
     * Connect client
     */
    client = g_connect_client(s,G_LOCK_EX);
    if (client == NULL) {
	GAP_ERROR("cannot connect client (gap-create-db)");
	return(1);
    }


    /*
     * Lock the records we need
     */
    for(i=0;i<BIT_ELE;i++) {
	view[i] = g_lock_N(client,GAP_DATABASE_FILE,i,G_LOCK_EX);
	if (view[i]<0) {
	    GAP_ERROR("Could not lock record");
	    return(1);
	}
	/*
	 * write something
	 */
	err = GAP_WRITE(client,view[i],NULL,0,GT_Unknown,sizeof(GCardinal));
	if (err) {
	    GAP_ERROR("could not write initial 'something'");
	    return 1;
	}
    }


    
    /*
     * initialise database header record
     */
    db.Nfreerecs        = 4; /* number of bytes set */
    db.freerecs         = GR_Freerecs;
    db.Ncontigs         = 0;
    db.contigs          = GR_Contigs;
    db.contig_order	= GR_Contig_Order;
    db.Nreadings        = 0;
    db.readings         = GR_Readings;
    db.Nannotations     = 0;
    db.annotations      = GR_Annotations;
    db.maximum_db_size  = maxdb;
    db.actual_db_size   = maxdb;
    db.max_gel_len      = GAP_READ_LEN;
    db.data_class       = GAP_DNA;
    db.version          = GAP_DB_VERSION;
    db.num_contigs      = 0;
    db.num_readings     = 0;
    db.free_annotations = 0;
    db.Ntemplates	= 0;
    db.templates	= GR_Templates;
    db.Nclones		= 0;
    db.clones		= GR_Clones;
    db.Nvectors		= 0;
    db.vectors		= GR_Vectors;
    db.Nnotes 		= 0;
    db.notes_a		= 0;
    db.notes		= 0;
    db.free_notes	= 0;

    /*
     * Default vector
     */
    vec.name		= GR_Unknown;
    vec.level		= GAP_LEVEL_UNKNOWN;

    /*
     * Default template
     */
    temp.name	        = GR_Unknown;
    temp.strands	= 0;
    temp.vector		= 0;
    temp.clone		= 0;
    temp.insert_length_min = DEFAULT_SI_from;
    temp.insert_length_max = DEFAULT_SI_to;

    /*
     * Default clone
     */
    clone.name		= GR_Unknown;;
    clone.vector	= 0;
    
    
    /*
     * write enough to keep gap happy
     */
    err = 0;
    err |= GAP_WRITE(client,view[GR_Database],&db,sizeof(db),GT_Database,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Freerecs],&BITS,sizeof(BITS),GT_Bitmap,sizeof(BITS));
    err |= GAP_WRITE(client,view[GR_Contigs],NULL,0,GT_Array,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Readings],NULL,0,GT_Array,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Annotations],NULL,0,GT_Array,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Templates],&TEMP,sizeof(TEMP),GT_Array,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Clones],&CLONE,sizeof(CLONE),GT_Array,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Vectors],&VEC,sizeof(VEC),GT_Array,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Templates_Default],&temp,sizeof(temp),GT_Templates,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Clones_Default],&clone,sizeof(clone),GT_Clones,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Vectors_Default],&vec,sizeof(vec),GT_Vectors,sizeof(GCardinal));
    err |= GAP_WRITE(client,view[GR_Unknown],UNKNOWN,sizeof(UNKNOWN),GT_Text,sizeof(char));
    err |= GAP_WRITE(client,view[GR_Contig_Order],NULL,0,GT_Array,sizeof(GCardinal));

    if (err) {
	GAP_ERROR("Could not write initial records");
	return 1;
    }

    /*
     * Unlock
     */
    for(i=0;i<BIT_ELE;i++) {
	err = g_unlock(client,view[i]);
	if (err) {
	    GAP_ERROR("Could not lock record");
	    return(1);
	}
    }

	
    if (err=g_disconnect_client(client)) {
	GAP_ERROR("problem disconnecting (gap-create-db)");
	return(1);
    }

    gap_shutdown_server(s);

    return 0;
}






int gap_new_db(char *project,char *version,int read_only)
/*
 *
 */
{
    int err;


    /*
     * create empty db
     */
    if ( err = gap_create_db(project,version) ) {
	GAP_ERROR("cannot create database");
	return(1);
    }


    /*
     * fill it with something
     */
    if ( err = gap_init_db(project,version,read_only) ) {
	GAP_ERROR("cannot initialise database");
	return(1);
    }



    return 0;
}


/*
 * Copies database <base> from version <from> to version <to>.
 * We assume that the files do not already exist. Safe when called from DBCOPY
 * as this already performs this check for us.
 * In the event of an error (eg file system full) files maybe left around.
 *
 * For more complex copies (such as altering the database size) we should
 * open our new DB and use writr0_ etc. This is best implemented with open2t_
 * and an extra function to switch 'active' GapIO.
 *
 * NB: The current method preserves time stamps, but does not aid with any
 * garbage collection.
 */

int cpdb(char *base, char *from, char *to) {
    int ifd, ofd;
    int i;
    char buf[BUFSIZ];

    /* Loop around copying each file pair for each GAP_FILE */
    for (i = 0; i < GAP_FILES; i++) {
	char fn_f[256], fn_t[256];
	int e, j, x, c;

	(void)gap_construct_file(base, file_list[i], from, fn_f);
	(void)gap_construct_file(base, file_list[i], to, fn_t);

	for (j = 0; j < 2; j++) {
	    /* do the opening */
	    if (-1 == (ifd = open(fn_f, O_RDONLY | O_BINARY))) {
		perror(fn_f);
		return -1;
	    }
	    if (-1 == (ofd = open(fn_t, O_RDWR | O_TRUNC | O_CREAT | O_BINARY,
				  0666))) {
		perror(fn_t);
		return -1;
	    }

	    /* now do the copying - LOW LEVEL IO HERE */
	    while ((e = read(ifd, buf, BUFSIZ)) > 0) {
		c = 0;
		while (e > 0) {
		    x = write(ofd, buf+c, e);
		    if (x == -1) {
			e = -2;
			break;
		    }
		    e -= x;
		    c += x;
		}
		if (e < 0)
		    break;
	    }
	    if (e < 0) {
		perror(e == -1 ? "read" : "write");
		return -1;
	    }

	    /* close the file */
	    close(ifd);
	    close(ofd);

	    /* repeat for AUX file */
	    strcat(fn_f, G_AUX_SUFFIX);
	    strcat(fn_t, G_AUX_SUFFIX);
	}
    } 

    return 0;
}

f_proc_ret cpdb_(char *base, char *from_, char *to_,
		 f_implicit base_l, f_implicit from_l, f_implicit to_l) {
    char from[256], to[256];

    /* NULL terminate base. We can assume there is room to place the NULL as
     * <to> is infact <base + base_l + 1> when called from the fortran DBCOPY
     * routine.
     */
    base[base_l] = '\0';
    /*
    strncpy(from, from_, from_l); from[from_l] = '\0';
    strncpy(to, to_, to_l); to[to_l] = '\0';
    */

    Fstr2Cstr(from_, from_l, from, 256);
    Fstr2Cstr(to_,   to_l,   to,   256);

    cpdb(base, from, to);

    f_proc_return();
}

void set_db_bitsize(int bs) {
    bitsize = bs;
}
