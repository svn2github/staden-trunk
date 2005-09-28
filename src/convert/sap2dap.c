#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
/*#include <sys/unistd.h>*/
#include "fort.h"

#define l_len 100
#define f_len 200

#define MAXDB 1000

struct _ar_rec {
    int_f idbsiz;
    int_f maxgel;
    int_f idm;
} AR_rec;

struct _rl_rec {
    int_f ngels;
    int_f nconts;
    int_f dum1;
    int_f dum2;
} RL_rec;

struct _tg_rec {
    int_f pos; /* and count */
    int_f len;
    int_f com;
    int_f type;
    int_f next;
} TG_rec;

#define COMMENT_LENGTH 40
struct _cc_rec {
    int_f next;
    char comment[COMMENT_LENGTH];
} CC_rec;

struct _cch_rec{
    int_f next;
    int_f count;
    char comment[COMMENT_LENGTH-sizeof(int_f)];
} CCH_rec;

struct _rd_rec{
    int_f len;
    int_f lcut;
    int_f wlen;
    char type[4];
    char name[12];
} RD_rec;

int main()
{
  char projectName[l_len];
  char versionNumber[l_len];

  fprintf(stdout,"Database conversion program\n");
  fprintf(stdout,"Converts *.RD? file to *.TG? and *.CC? files\n\n");

  fprintf(stdout,"Project name ? ");
  gets(projectName);

  fprintf(stdout,"Version ? ");
  gets(versionNumber);

  if (process(projectName,versionNumber))
    fprintf(stdout,"Error: conversion aborted.\n");
  else
    fprintf(stdout,"Conversion completed.\n");

  return 0;
}

int read_ar(char *AR,int_f *idbsiz)
{

    FILE *AR_fp;

    /*
    ** Check AR file exists
    */
    if ((AR_fp=fopen(AR,"rb"))==NULL) {
	fprintf(stderr,"Cannot open file %s\n",AR);
	return 1;
    }

    /*
    ** Read details from AR file
    */
    fseek(AR_fp,(off_t)((MAXDB-1)*sizeof(AR_rec)),/*SEEK_SET*/0);
    fread(&AR_rec,sizeof(AR_rec),1,AR_fp);

    if (ferror(AR_fp)) {
	fprintf(stderr,"Cannot read file %s\n",AR);
	return 1;
    }

    fclose(AR_fp);

    *idbsiz = AR_rec.idbsiz;
    return 0;
}

int read_rl(char *RL,int_f idbsiz,     int_f * ngels)
{

    FILE *RL_fp;

    /*
    ** Check RL file exists
    */
    if ((RL_fp=fopen(RL,"rb"))==NULL) {
	fprintf(stderr,"Cannot open file %s\n",RL);
	return 1;
    }

    /*
    ** Read details from RL file
    */
    fseek(RL_fp,(off_t)((AR_rec.idbsiz-1)*sizeof(RL_rec)),/*SEEK_SET*/0);
    fread(&RL_rec,sizeof(RL_rec),1,RL_fp);

    if (ferror(RL_fp)) {
	fprintf(stderr,"Cannot read file %s\n",RL);
	return 1;
    }

    fclose(RL_fp);

    *ngels = RL_rec.ngels;
    return 0;
}

void write_tg(FILE *fp,int_f rec, int_f pos, int_f len, int_f com, int_f type, int_f next)
{
    TG_rec.pos = pos;
    TG_rec.len = len;
    TG_rec.com = com;
    TG_rec.type = type;
    TG_rec.next = next;

    fseek(fp,(off_t)((rec-1)*sizeof(TG_rec)),/*SEEK_SET*/0);
    fwrite(&TG_rec,sizeof(TG_rec),1,fp);
}

void write_cc(FILE *fp,int_f rec, int_f next, char *comment)
{
    CC_rec.next = next;
    strncpy(CC_rec.comment,comment,COMMENT_LENGTH);

    fseek(fp,(off_t)((rec-1)*sizeof(CC_rec)),/*SEEK_SET*/0);
    fwrite(&CC_rec,sizeof(CC_rec),1,fp);
}

void write_cc_head(FILE *fp,int_f next, int_f count )
{
    CCH_rec.next = next;
    CCH_rec.count = count;

    fseek(fp,(off_t)0,/*SEEK_SET*/0);
    fwrite(&CCH_rec,sizeof(CC_rec),1,fp);
}

void read_rd(FILE *fp, int_f rec, int_f *len, int_f *lcut, int_f *wlen, char *type,
char *name)
{
    fseek(fp,(off_t)((rec-1)*sizeof(RD_rec)),/*SEEK_SET*/0);
    fread(&RD_rec,sizeof(RD_rec),1,fp);

    *len = RD_rec.len;
    *lcut = RD_rec.lcut;
    *wlen = RD_rec.wlen;
    strncpy(type,RD_rec.type,4);
    strncpy(name,RD_rec.name,12);

}

int process(char *name,char *vers)
{
    char AR[f_len];
    char RD[f_len];
    char CC[f_len];
    char TG[f_len];
    char RL[f_len];

    FILE *RD_fp;
    FILE *CC_fp;
    FILE *TG_fp;

    struct stat statBuff;

    int_f IDBSIZ;
    int_f NGELS;

    if (!*vers) strcpy(vers,"0");

    /* convert bits to upper case */
    {
        char *s;
        for (s=name; *s = islower(*s)?toupper(*s):*s ;s++);
	for (s=vers; *s = islower(*s)?toupper(*s):*s ;s++);
    }

    /*
    ** create file names
    */
    strcpy(RD,name); strcat(RD,".RD"); strncat(RD,vers,1);
    strcpy(TG,name); strcat(TG,".TG"); strncat(TG,vers,1);
    strcpy(CC,name); strcat(CC,".CC"); strncat(CC,vers,1);

    /*
    ** Get AR details
    */
    strcpy(AR,name); strcat(AR,".AR"); strncat(AR,vers,1);
    if (read_ar(AR,&IDBSIZ)) return 1;

    /*
    ** Get RL details
    */
    strcpy(RL,name); strcat(RL,".RL"); strncat(RL,vers,1);
    if (read_rl(RL,IDBSIZ,   &NGELS)) return 1;

    /*
    ** open files TG and CC
    */
    if ( stat(TG,&statBuff) >= 0 ) {
	fprintf(stderr,"%s already exists\n",TG);
	return 1;
    }

    if ( stat(CC,&statBuff) >= 0 ) {
	fprintf(stderr,"%s already exists\n",CC);
	return 1;
    }

    if ((TG_fp=fopen(TG,"wb"))==NULL) {
	 fprintf(stderr,"cannot open %s for writing\n",TG);
	 return 1;
    }

    if ((CC_fp=fopen(CC,"wb"))==NULL) {
	 fprintf(stderr,"cannot open %s for writing\n",CC);
	 fclose(TG_fp);
	 return 1;
    }

    /*
    ** Check RD file exists
    */
    if ((RD_fp=fopen(RD,"rb"))==NULL) {
	/*
	** None:
	** Create anyway
	*/
	write_tg(TG_fp,IDBSIZ/*rec*/,IDBSIZ,0,0,0,0);
	write_cc_head(CC_fp,0/*next*/,1/*count*/);

    } else {
	/*
	** Do the hard graft
	*/
	int i;

	write_tg(TG_fp,IDBSIZ/*rec*/,IDBSIZ,0,0,0,0);
	write_cc_head(CC_fp,0/*next*/,NGELS+1/*count*/);

	for (i=1;i<=NGELS;i++) {
	    int_f len;
	    int_f lcut;
	    int_f wlen;
	    char type[5];
	    char name[13];
	    char comment[COMMENT_LENGTH];

	    read_rd(RD_fp,i,&len,&lcut,&wlen,type,name);
	    type[4] = '\0';
	    name[12] = '\0';

	    sprintf(comment,"%6d%6d%6d%-4s%-18s",len,lcut,wlen,type,name);

	    write_tg(TG_fp,i,0,0,i+1,0,0);
	    write_cc(CC_fp,i+1,0,comment);
	}
    }

    fclose(TG_fp);
    fclose(CC_fp);
    return 0;
}
