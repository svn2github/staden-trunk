/* alfsplit.c  
   Written by Richard Durbin 12/28/90.
   Takes big combined alf results file, and splits it into
   separate files for each clone.
   Only keep processed data, sequence data and experimental notes.
   Although the format of the small files is based on that of an
   ALF file officially split on the PC, they are unfortunately
   not reaadable by ALFManager software on the PC.
   */

/* first give full function prototypes for system functions
   these are incomplete in the Sun /usr/include/...
   
   Modified by Simon Dear 21 August 1991.
   Ignore value of s3 in readIndexEntry in check on sensible values
   This value is DirEntry.fType in Pharmacia documentation
   
   Modified by Simon Dear 25 October 1991.
   Machine independant I/O

   24 August 1992 [Simon Dear]
   MAJOR HACK - to allow for readings which have no clone name

   22 Sept 1994, James Bonfield
   Tidy up of output - only display index entries that contain data
   */

#include <staden_config.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdarg.h>  /* varargs needed for v*printf() prototypes */
#include <ctype.h>

#include "os.h"

typedef int BOOL;
#define TRUE 1
#define FALSE 0

/********** routines to read and write Index entries ***********/

static char junk[512] ;	/* for when we want to read/write junk */

/***** architecture independant reads ******/
static int_4 read_int_4(FILE *fp)

{
    unsigned char buf[sizeof(int_4)];
    
    if (fread(buf, sizeof(buf), 1, fp) != 1) return 0;
    return (int_4)
        (((uint_4)buf[0]) +
         ((uint_4)buf[1]<<8) +
         ((uint_4)buf[2]<<16) +
         ((uint_4)buf[3]<<24));
}

static int_2 read_int_2(FILE *fp)

{
    unsigned char buf[sizeof(int_2)];
    
    if (fread(buf, sizeof(buf), 1, fp) != 1) return 0;
    return (int_2)
        (((uint_2)buf[0]) +
         ((uint_2)buf[1]<<8));
}

static void write_int_4(FILE *fp, int_4 l)
{
    unsigned char buf[sizeof(int_4)];
    
    buf[0] = (unsigned char)(l&255);
    buf[1] = (unsigned char)(l>>8&255);
    buf[2] = (unsigned char)(l>>16&255);
    buf[3] = (unsigned char)(l>>24&255);
    
    fwrite(buf, sizeof(buf), 1, fp);
}

static void write_int_2(FILE *fp, int_2 l)
{
    unsigned char buf[sizeof(int_2)];
    
    buf[0] = (unsigned char)(l&255);
    buf[1] = (unsigned char)(l>>8&255);
    
    fwrite(buf, sizeof(buf), 1, fp);
}







typedef struct IndexEntryStruct
{ int_4 isTraces ;
  char  label[40] ;
  int_4 dataLen ;
  int_4 blockLen ;
  int_4 offset ;
  int   valid ;
} *IndexEntry ;

static BOOL readIndexEntry (FILE *fil, IndexEntry ent)
{
    short s1,s2,s3 ;
    
#define readInt() (read_int_4(fil))
#define readShort() (read_int_2(fil))
    
    clearerr (fil) ;
    
    s1 = readShort() ; /* flag: bit 0: 1==used, 0==free */
    s2 = readShort() ; /* subfile class */
    s3 = readShort() ; /* subfile type */

    /* Was this, but s3 can be 0 now
       if ((s1 != 1 || s3 != 1) && (s1 || s2 || s3))
       */
    if ((s1 != 1) && (s1 || s2))
	return FALSE ;
    ent->isTraces = (s2 == 4) ;
    
    fread (ent->label,40,1,fil) ;
    ent->dataLen = readInt() ;
    ent->blockLen = readInt() ;
    ent->offset = readInt() ;

    ent->valid = ent->dataLen ? 1 : 0;
    fread (junk,70,1,fil) ;
    
    return !ferror (fil) ;
}

static BOOL writeIndexEntry (FILE *fil, IndexEntry ent)
{
    
#define writeInt(xx) (write_int_4(fil,xx))
#define writeShort(xx) (write_int_2(fil,xx))
    
    clearerr (fil) ;
    
    writeShort(1) ;
    if (ent->isTraces)
	writeShort(4) ;
    else
	writeShort(2) ;
    writeShort(1) ;
    
    fwrite (ent->label,40,1,fil) ;
    writeInt(ent->dataLen) ;
    writeInt(ent->blockLen) ;
    writeInt(ent->offset) ;
    
    fwrite (junk,70,1,fil) ;
    
    return !ferror (fil) ;
}

/************************************************************/

void crash (char* format,...)
{
    va_list args ;
    
    va_start (args,format) ;
    vfprintf (stderr,format,args) ;
    va_end (args) ;
    
    exit (1) ;
}

/*****************/

static void readLine (FILE *fil, char* cp)
{
    while ((*cp = fgetc(fil)) && *cp != EOF && *cp != '\n')
	++cp ;
    *cp = 0 ;
}

/*****************/
#define MAXCLONES 10
int main (int argc, char* *argv)
{
    FILE *inEnt, *inData, *outEnt[MAXCLONES], *outData[MAXCLONES] ;
    /* open two pointers in each file - index and data */
    IndexEntry EN ;
    IndexEntry ent ;
    char expLine[4][20],name[MAXCLONES][256],note[MAXCLONES][80],fname[25];
    off_t seqOffset[MAXCLONES], dataOffset[MAXCLONES];
    int_4 seqDataLen[MAXCLONES], dataDataLen[MAXCLONES];
    int_4 seqBlockLen[MAXCLONES], dataBlockLen[MAXCLONES];
    char buf[512] ;
    int i,j,len ;
    size_t lastDot,lastSlash;
    
    if (argc != 2)
	crash ("Usage: alfsplit rawfilename\n") ;
    
    inData = fopen (argv[1],"r") ;
    if (!(inEnt = fopen (argv[1],"r")))
	crash ("Could not open file '%s'\n",argv[1]) ;
    
    /* first find the experimental notes entry and extract file names */
    
    ent = (IndexEntry) malloc (sizeof (struct IndexEntryStruct)) ;
    EN = (IndexEntry) malloc (sizeof (struct IndexEntryStruct)) ;
    if (fseeko (inEnt,(off_t)512,0))
	crash ("Could not seek to index in raw file\n") ;
    while (TRUE)
	{ if (!readIndexEntry (inEnt,EN))
	      crash ("Can't find Experimental Notes index entry\n") ;
	  if (!strcmp (EN->label,"ALF Experimental notes"))
	      break ;
      }
    
    if (fseeko (inData,(off_t)EN->offset,0))
	crash ("Can't seek to Experimental notes\n") ;
    for (i = 0 ; i < 4 ; ++i)
	readLine (inData,expLine[i]) ;
    
    /* determine default root name from argv[1]:
     ** I assume this has the format {path}/{name}.alf
     ** Default names will be {name}.1, {name}.2, ..., {name}.MAXCLONES
     */
    lastDot = (size_t)0;
    for (i = strlen(argv[1])-1;i>=0 && argv[1][i] != '/'; i--)
	if (lastDot==0 && argv[1][i] == '.') lastDot = (size_t)i;
    if (lastDot==0) lastDot = strlen(argv[1]);
    lastSlash = (size_t)i;
    
    for (i = 0 ; i < MAXCLONES ; ++i) {
	int guess_name = 0;
	
	readLine (inData,name[i]) ;
	if (0 == strcmp(name[i], "blank") || *name[i] == 0)
	    guess_name = 1;
	else {
	    int j, l = strlen(name[i]);
	    
	    for (j = 0; j < l; j++) {
		if (!(isprint(name[i][j]) &&
		      !isspace(name[i][j]) &&
		      name[i][j] != '/'))
		    break;
	    }
	    if (j != l)
		guess_name = 1;
	}

	if (guess_name) {
	    sprintf(name[i],
		    "%.*s.%d",
		    (int)(lastDot-lastSlash-1),
		    argv[1]+lastSlash+1,
		    (i>10)?i:(i+1)%10);
	}
    }
    for (i = 0 ; i < MAXCLONES ; ++i)
	readLine (inData,note[i]) ;

    for (i = 0 ; i < MAXCLONES ; ++i)
	seqOffset[i] = dataOffset[i] = 0;

    /* gather offset information */
    fseeko (inEnt,(off_t)512,0) ;
    while (readIndexEntry (inEnt,ent)) {
	if (!ent->valid)
	    continue;

	printf ("%s: %d\n",ent->label,ent->offset/512) ;
	if (!strncmp (ent->label,"ALF Sequence data Clone ",24))
	    len = 24 ;
	else if (!strncmp (ent->label,"ALF Processed data Clone ",25))
	    len = 25 ;
	else
	    continue ;
	/* fall through to here if sequence or processed */
	i = atoi (&ent->label[len]) - 1 ;
	
	if (len == 24) {
	    seqOffset[i] = (off_t)ent->offset;
	    seqBlockLen[i] = ent->blockLen;
	    seqDataLen[i] = ent->dataLen;
	} else {
	    dataOffset[i] = (off_t)ent->offset;
	    dataBlockLen[i] = ent->blockLen;
	    dataDataLen[i] = ent->dataLen;
	}
    }

    /* initialise output files for clones */
    for (i = 0 ; i < MAXCLONES ; ++i) {
	if (seqOffset[i]==0 && dataOffset[i]==0) {
	    /* we are missing sequence and/or trace data */
	    printf ("Clone %d: %s - NOT MAKING BECAUSE THERE IS NO TRACE AND SEQUENCE DATA\n",i+1,name[i]);
	} else if (seqOffset[i]==0) {
	    /* we are missing sequence and/or trace data */
	    printf ("Clone %d: %s - NOT MAKING BECAUSE THERE IS NO SEQUENCE DATA\n",i+1,name[i]);
	} else if (dataOffset[i]==0) {
	    /* we are missing sequence and/or trace data */
	    printf ("Clone %d: %s - NOT MAKING BECAUSE THERE IS NO TRACE DATA\n",i+1,name[i]);
	} else  {
	    printf ("Clone %d: %s - %s\n",i+1,name[i],note[i]) ;
	    /* create the file and write the notes */
	    sprintf (fname,"%sALF",name[i]) ;
	    outData[i] = fopen (fname,"w") ;
	    if(!fwrite("ALF ", 4, 1, outData[i]))
		fprintf(stderr, "could not write file: %s\n", fname);
	    if(!fwrite (junk,512-4,1,outData[i]) )
		fprintf(stderr, "could not write file: %s\n", fname);
	    if(!fwrite (junk,512,1,outData[i]) )
		fprintf(stderr, "could not write file: %s\n", fname);
	    len = 0 ;
	    for (j = 0 ; j < 4 ; ++j)
		len += fprintf (outData[i],"%s\n",expLine[j]) ;
	    len += fprintf (outData[i],"%s\n\n\n\n\n\n\n\n\n\n",name[i]) ;
	    len += fprintf (outData[i],"%s\n\n\n\n\n\n\n\n\n\n",note[i]) ;
	    fwrite (junk,512-len,1,outData[i]) ;
	    /* now write the index entry */
	    if (!(outEnt[i] = fopen (fname,"a")))
		crash ("Couldn't open output file %s\n",fname) ;
	    fseeko (outEnt[i],(off_t)512,0) ;
	    EN->offset = 1024 ;
	    EN->dataLen = len ;
	    writeIndexEntry (outEnt[i],EN) ;

	    /*
            ** Copy sequence and trace data
	    */
	    /* trace data */
	    strcpy(ent->label,"ALF Processed data Clone 1");
	    ent->isTraces = 1;
	    ent->offset = (int_4)ftell (outData[i]) ;
	    ent->blockLen = dataBlockLen[i];
	    ent->dataLen = dataDataLen[i];
	    
	    fseeko (inData,dataOffset[i],0) ;
	    len = ent->blockLen/512 ;
	    for (j = 0 ; j < len ; ++j) {
		fread (buf,512,1,inData) ;
		fwrite (buf,512,1,outData[i]) ;
	    }
	    
	    writeIndexEntry (outEnt[i],ent) ;

	    /* sequence data*/
	    strcpy(ent->label,"ALF Sequence data Clone 1");
	    ent->isTraces = 0;
	    ent->offset = (int_4)ftell (outData[i]) ;
	    ent->blockLen = seqBlockLen[i];
	    ent->dataLen = seqDataLen[i];

	    fseeko (inData,seqOffset[i],0) ;
	    len = ent->blockLen/512 ;
	    for (j = 0 ; j < len ; ++j) {
		fread (buf,512,1,inData) ;
		fwrite (buf,512,1,outData[i]) ;
	    }
	    
	    writeIndexEntry (outEnt[i],ent) ;

	}
    }

    return 0;
}
