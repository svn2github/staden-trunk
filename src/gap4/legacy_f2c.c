/* legacy.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h>
#include "f2c.h"

/* Common Block Declarations */

struct {
    integer points[256];
} shotc_;

#define shotc_1 shotc_

struct {
    integer point1[256];
} iasci1_;

#define iasci1_1 iasci1_

struct {
    integer point2[256];
} iasci2_;

#define iasci2_1 iasci2_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__5 = 5;
static integer c__100 = 100;
static integer c__80 = 80;
static integer c__20 = 20;

/* note to kfs: add extra option to assembly menu "Ignore previous data"? */
/* duplicate dialogue from "normal shotgun assembly" but */
/* add extra argument NOPT to argument list for dbauto. */
/* set to 0 for all existing calls to dbauto and to 1 for new option. */
/*<    >*/
/* Subroutine */ int dbauto_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *maxdb, integer *idbsiz, integer *ngels, 
	integer *nconts, integer *clist, integer *maxgel, char *seq1, char *
	seq2, char *seq3, char *seq4, char *seq5, char *seqc2, char *seqg2, 
	char *seqg3, char *seqc3, char *rnames, integer *maxseq, integer *
	maxglm, integer *savps, integer *savpg, integer *savl, integer *
	maxsav, integer *cends, integer *nends, integer *maxcon, integer *
	idev1, char *namarc, char *nampro, real *percd, integer *iokent, 
	integer *ishow, integer *minmat, integer *maxpg, real *permax, 
	integer *irepsc, integer *iopt, integer *nopt, integer *ansjok, 
	integer *ansfe, integer *iwing, integer *nbad, char *list, integer *
	iok, integer *minovr, ftnlen seq1_len, ftnlen seq2_len, ftnlen 
	seq3_len, ftnlen seq4_len, ftnlen seq5_len, ftnlen seqc2_len, ftnlen 
	seqg2_len, ftnlen seqg3_len, ftnlen seqc3_len, ftnlen rnames_len, 
	ftnlen namarc_len, ftnlen nampro_len, ftnlen list_len)
{
    /* System generated locals */
    integer seqc2_dim1, seqc2_offset, seqg2_dim1, seqg2_offset, i__1, i__2;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, j, n, x;
    static real z__;
    static integer pl[2], kt, pr[2], ilc[2], idm, ltl, ltr;
    extern /* Subroutine */ int ccta_(char *, integer *, ftnlen);
    static integer jobc, jgel;
    static char csen[1];
    static integer lreg, ilcr;
    extern /* Subroutine */ int info_(char *, ftnlen);
    static integer mask, leno, rreg, ilct, ierr, cnum, idim1, idim2, iladd[1];
    extern /* Subroutine */ int swrt1_(char *, char *, ...),
	swrt2_(char *, char *, ...);
    static integer iradd[1], ifail[2];
    extern /* Subroutine */ int swrt5_(char *, char *, ...);
    static integer idim22[2], kfail;
    extern /* Subroutine */ int aline_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, real *, 
	    char *, integer *, real *, integer *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer igelc;
    static char infod[80];
    extern integer gclin_(integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *);
    extern /* Subroutine */ int sindb_(integer *, integer *, char *, char *, 
	    integer *, ftnlen, ftnlen);
    static integer jngel, imatc, igood, joinf, idsav, maxpc, itask, llino[2], 
	    iposc[2], iposg[2], joint[2], idout[2], ioptc, klass, iover, 
	    itype[2], lmost;
    extern /* Subroutine */ int sqrev_(char *, integer *, ftnlen);
    static integer rmost;
    extern /* Subroutine */ int ajoin2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), ajoin3_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), precn1_(char *, char 
	    *, real *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen), dbchek_(integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *), abedin_(integer *, integer *,
	     integer *, integer *, integer *, integer *, char *, integer *, 
	    integer *, char *, integer *, integer *, integer *, integer *, 
	    integer *, ftnlen, ftnlen), swrt2b_(char *, char *, ...);
    extern integer gnread_(char *, ftnlen);
    extern /* Subroutine */ int delcon_(char *, integer *, integer *, integer 
	    *, ftnlen), addtit_(char *, char *, integer *, integer *, ftnlen, 
	    ftnlen), dbautp_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, char *, char *, integer *, integer *, 
	    integer *, char *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, char *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static integer ifcomp;
    extern /* Subroutine */ int aenter_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, char *, integer *, 
	    integer *, integer *, char *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    char *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer lincon[2], ilefts[2], isense[2], itotpc[2];
    static real permis[2];
    static char infoud[80];
    static integer itotpg[2];
    extern integer cmpseq_(integer *, char *, integer *, integer *, integer *,
	     integer *, integer *, char *, char *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer minsli, iempty, jnjoin, maxovr;
    extern /* Subroutine */ int precon_(char *, char *, real *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    ftnlen, ftnlen), erromf_(char *, ftnlen);
    static integer ngelsl, ncontl;
    extern /* Subroutine */ int updout_(void), arrfio_(char *, char *, 
	    integer *, integer *, integer *, ftnlen, ftnlen), aerror_(char *, 
	    char *, integer *, ftnlen, ftnlen), sqcopy_(char *, char *, 
	    integer *, ftnlen, ftnlen), autocn_(char *, integer *, char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    char *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, char *, char *, char *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, real *, integer *, char *, char *, integer *, 
	    integer *, real *, integer *, integer *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), 
	    tolist_(char *, char *, ftnlen, ftnlen), updcon_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    char *, integer *, integer *, integer *, integer *, integer *, 
	    char *, char *, integer *, integer *, integer *, integer *, real *
	    , integer *, integer *, ftnlen, ftnlen, ftnlen), cmplmt_(integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, char *, integer *, integer *, integer *, ftnlen), 
	    sqcomm_(char *, integer *, ftnlen);


/* IOKENT = 0 permit entry, 1 forbid entry */
/* ISHOW  = 1  hide all alignemtns */
/*          2  show passed alignments */
/*          3  show all alignments */
/*          4  show only failed alignments */
/* IREPSC = 0  save alignment scores in file, 1 dont */
/* IOPT   = 1 'Perform normal shotgun assembly' */
/*          2 'Perform shotgun assembly with tagged segments masked' */
/*          3 'Put all sequences in one contig' */
/*          4 'Put all sequences in new contigs' */
/*          5 'Perform shotgun assembly with finished segments masked' */
/* NOPT   = 1 'Perform shotgun assembly but ignore all previous data' */
/* ANSJOK = 0  permit joins, 1 dont */
/* ANSFE  = 1 'Reject failures' */
/*          2 'Enter all readings' */


/*<       INTEGER CLIST(1),CNUM >*/
/*<       INTEGER RELPG(MAXDB),PL(2),PR(2),RMOST >*/
/*<       INTEGER LNGTHG(MAXDB),LNBR(MAXDB),RNBR(MAXDB) >*/
/*<       INTEGER JOINT(2),ITOTPC(2),ITOTPG(2),IDIM22(2),IDOUT(2) >*/
/*<       INTEGER LINCON(2),LLINO(2),ITYPE(2),IFAIL(2) >*/
/*<       INTEGER ILEFTS(2),ILC(2),IPOSC(2),IPOSG(2),ISENSE(2) >*/
/*<       INTEGER LREG,RREG,X,ANSJOK,ANSFE >*/
/*<       INTEGER CENDS(MAXCON),NENDS(MAXCON),ILADD(1),IRADD(1) >*/
/*<       CHARACTER SEQ3(MAXGLM),SEQC2(MAXGLM,2),SEQG2(MAXGLM,2) >*/
/*<       CHARACTER SEQ1(MAXSEQ),SEQ2(MAXGLM),SEQ4(MAXGLM) >*/
/*<       INTEGER SAVPS(MAXSAV),SAVPG(MAXSAV),SAVL(MAXSAV) >*/
/*<       CHARACTER NAMARC*(*),NAMPRO*(*),LIST*(*) >*/
/*<       CHARACTER SEQ5(MAXGLM),SEQG3(MAXGLM),SEQC3(MAXGLM) >*/
/*<       CHARACTER*(*) RNAMES(IDBSIZ) >*/
/*<       CHARACTER CSEN >*/
/*<       REAL PERMIS(2) >*/
/*<       CHARACTER INFOD*80 >*/
/*<       CHARACTER INFOUD*80 >*/
/*<       INTEGER GNREAD,GCLIN,CMPSEQ >*/
/*<       EXTERNAL GNREAD,GCLIN,CMPSEQ >*/
/*<       CALL INFO('Automatic sequence assembler') >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    rnames -= rnames_len;
    --clist;
    --seq1;
    --seqc3;
    --seqg3;
    seqg2_dim1 = *maxglm;
    seqg2_offset = 1 + seqg2_dim1;
    seqg2 -= seqg2_offset;
    seqc2_dim1 = *maxglm;
    seqc2_offset = 1 + seqc2_dim1;
    seqc2 -= seqc2_offset;
    --seq5;
    --seq4;
    --seq3;
    --seq2;
    --savl;
    --savpg;
    --savps;
    --nends;
    --cends;

    /* Function Body */
    info_("Automatic sequence assembler", (ftnlen)28);
/*<       IDM = 5 >*/
    idm = 5;
/*<       MINSLI = 3 >*/
    minsli = 3;
/*<       IFAIL(1) = 0 >*/
    ifail[0] = 0;
/*<       IEMPTY=0 >*/
    iempty = 0;
/*<       MAXPC = MAXPG >*/
    maxpc = *maxpg;
/* note by KFS (22/5/95) - use of CLIST and CNUM by auto_assemble. */
/* auto_assemble is a special case where all contigs are always used */
/* so although CLIST = NULL, CNUM = NCONTS and this fact is catered for in */
/* get_contig_list */
/*<       CNUM = NCONTS >*/
    cnum = *nconts;
/* FIXME */
/*      IF (IOPT.EQ.2) IOPT = 5 */
/*      WRITE(*,*)'permit entry',IOKENT */
/*      WRITE(*,*)'display',ISHOW */
/*      WRITE(*,*)'align',IREPSC */
/*      WRITE(*,*)'option',IOPT */
/*      WRITE(*,*)'joins',ANSJOK */
/*      WRITE(*,*)'failure mode',ANSFE */
/*      WRITE(*,*)'window',IWING */
/*      WRITE(*,*)'dashes',NBAD */
/*      WRITE(*,*)'NGELS',NGELS */
/*      WRITE(*,*)'NOPT=',NOPT */
/*      NOPT = 1 */
/*<       IF((NGELS.LT.1).OR.(NOPT.EQ.1))IEMPTY=1 >*/
    if (*ngels < 1 || *nopt == 1) {
	iempty = 1;
    }
/*<    >*/
    dbchek_(idev1, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &idm, idbsiz, 
	    ngels, nconts, &ierr);
/*<       IF(IERR.GT.1) RETURN >*/
    if (ierr > 1) {
	return 0;
    }
/*<       CALL SINDB(IDEV1,NGELS,RNAMES,NAMARC,1) >*/
    sindb_(idev1, ngels, rnames + rnames_len, namarc, &c__1, rnames_len, 
	    namarc_len);

/* IOKENT = 0 permit entry, 1 forbid entry */
/* ISHOW  = 1  hide all alignemtns */
/*          2  show passed alignments */
/*          3  show all alignments */
/*          4  show only failed alignments */
/* IREPSC = 0  save alignment scores in file, 1 dont */
/* IOPT   = 1 'Perform normal shotgun assembly' */
/*          2 'Perform shotgun assembly with tagged segments masked' */
/*          3 'Put all sequences in one contig' */
/*          4 'Put all sequences in new contigs' */
/*          5 'Perform shotgun assembly with finished segments masked' */
/* ANSJOK = 0  permit joins, 1 dont */
/* ANSFE  = 1 'Reject failures' */
/*          2 'Enter all readings' */
/* MINOVR      The minimum overlap between a reading and a contig. */
/*             If this value is not reached the overlap is reported */
/*             in the count of overlaps but, no alignment is done and */
/*             as far as that match is concerned the reading does not */
/*             overlap ie if it does not overlap anywhere else with */
/*             amount MINOVR the reading will start a new contig. */


/*<       JGEL = 0 >*/
    jgel = 0;
/*<       JNGEL = 0 >*/
    jngel = 0;
/*<       JNJOIN = 0 >*/
    jnjoin = 0;
/*<       JOINF = 0 >*/
    joinf = 0;
/*<       MASK = 0 >*/
    mask = 0;
/*<       IMATC = 0 >*/
    imatc = 0;
/*<       IF((IOPT.EQ.3).OR.(IOPT.EQ.4)) MINMAT = 1 >*/
    if (*iopt == 3 || *iopt == 4) {
	*minmat = 1;
    }
/*<       IF (IOPT.EQ.2) MASK = 3 >*/
    if (*iopt == 2) {
	mask = 3;
    }
/*<       IF (IOPT.EQ.5) MASK = 4 >*/
    if (*iopt == 5) {
	mask = 4;
    }
/*<       IF ((IOPT.EQ.1).OR.(IOPT.EQ.2).OR.(IOPT.EQ.5)) THEN >*/
    if (*iopt == 1 || *iopt == 2 || *iopt == 5) {
/*<         IDIM1=0 >*/
	idim1 = 0;
/*<         MAXOVR=MAXGEL-3*MAX(MAXPC,MAXPG) >*/
	maxovr = *maxgel - max(maxpc,*maxpg) * 3;
/*<         IMATC = 0 >*/
	imatc = 0;
/* get ready for precon */

/* set task (normal consensus+title) */

/*<         ITASK = 5 >*/
	itask = 5;

/* set masking if required */

/*<         IF (MASK.EQ.3) ITASK = ITASK + 32 >*/
	if (mask == 3) {
	    itask += 32;
	}
/*<         IF (MASK.EQ.4) ITASK = 8 + 1 >*/
	if (mask == 4) {
	    itask = 9;
	}
/*<         IF(IEMPTY.EQ.0) THEN >*/
	if (iempty == 0) {
/*<    >*/
	    precon_(seq1 + 1, nampro, percd, idbsiz, &cnum, &clist[1], &itask,
		     idev1, &idim1, maxgel, maxseq, iwing, nbad, iladd, iradd,
		     ifail, (ftnlen)1, nampro_len);
/*<           IF(IFAIL(1).NE.0) THEN >*/
	    if (ifail[0] != 0) {
/*<             CALL ERROMF('Error calculating consensus') >*/
		erromf_("Error calculating consensus", (ftnlen)27);
/*<             GO TO 901 >*/
		goto L901;
/*<           END IF >*/
	    }
/*        IF(IDIM1.GT.0)CALL FMTDB(SEQ1,IDIM1,1,IDIM1,60,30) */
/*        IF (0.EQ.0) RETURN */
/*<         END IF >*/
	}
/*        WRITE(*,*)'INITIAL IDIM1',IDIM1,ITASK */

/* init hashing constants */

/*<       IOPTC = 1 >*/
	ioptc = 1;
/*<    >*/
	idsav = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1], 
		maxsav, seq1 + 1, seq2 + 1, maxseq, maxgel, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/*<       IF (IDSAV.NE.0) THEN >*/
	if (idsav != 0) {
/*<         IOPTC = 6 >*/
	    ioptc = 6;
/*<    >*/
	    idsav = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[
		    1], &idsav, seq1 + 1, seq2 + 1, maxseq, maxgel, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
/*<         RETURN >*/
	    return 0;
/*<       END IF >*/
	}
/*<       END IF >*/
    }

/* set intitial values for contig count and number of readings */
/* just to get thru the first consensus calc */
/*<       NGELSL = NGELS + 2 >*/
    ngelsl = *ngels + 2;
/*<       NCONTL = NCONTS +1 >*/
    ncontl = *nconts + 1;


/*                          MAIN LOOP */


/*< 1     CONTINUE >*/
L1:
/*<       CALL UPDOUT() >*/
    updout_();


/*<       IDIM2=MAXGEL >*/
    idim2 = *maxgel;
/*<       IOK = GNREAD(NAMARC) >*/
    *iok = gnread_(namarc, namarc_len);
/*<       IF(IOK.EQ.1) GO TO 900 >*/
    if (*iok == 1) {
	goto L900;
    }
/*<       IF(IOK.NE.0) GO TO 1 >*/
    if (*iok != 0) {
	goto L1;
    }
/*<       CALL INFO('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>') >*/
    info_(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", (ftnlen)46);
/*      WRITE(*,*)'>>>>>>>>>>>',NAMARC */
/*      WRITE(6,*)'IDIM1',IDIM1,NGELS */
/*         CALL FMTDB1(SEQ1,IDIM1,1,IDIM1,60,6) */
/*      WRITE(*,*)'MAIN' */
/*      WRITE(*,*)(SEQ1(JJJJ),JJJJ=1,IDIM1) */
/*<       JGEL = JGEL + 1 >*/
    ++jgel;
/*      WRITE(INFOD,1006)JGEL */
/* 1006 FORMAT('Processing ',I8,' in batch') */
/* CHECKED */
/*<       CALL SWRT1(INFOD,'Processing %8d in batch%!', JGEL) >*/
    swrt1_(infod, "Processing %8d in batch%!", &jgel, (ftnlen)80, (ftnlen)25);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*       WRITE(INFOD,1007)NAMARC */
/* 1007  FORMAT('File name ',A) */
/* CHECKED */
/*<       CALL SWRT2B(INFOD, 'File name %.*s%!', LEN(NAMARC), NAMARC) >*/
    i__1 = i_len(namarc, namarc_len);
    swrt2b_(infod, "File name %.*s%!", &i__1, namarc, (ftnlen)80, (ftnlen)16, 
	    namarc_len);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/* Added by Simon 23-March-1993 */
/*<       CALL ARRFIO(NAMARC,SEQ2,IDIM2,1,IOK) >*/
    arrfio_(namarc, seq2 + 1, &idim2, &c__1, iok, namarc_len, (ftnlen)1);
/*      CALL OPENRS(IDEV4,NAMARC,IOK,LRECL,2) */
/*<       IF(IOK.NE.0)THEN >*/
    if (*iok != 0) {
/*        IF(INF.EQ.1) RETURN */
/*<          CALL AERROR(LIST,NAMARC,0) >*/
	aerror_(list, namarc, &c__0, list_len, namarc_len);
/*<          GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }
/*       WRITE(INFOD,1800)IDIM2 */
/* 1800  FORMAT('Reading length ',I6) */
/* CHECKED */
/*<       CALL SWRT1(INFOD,'Reading length %6d%!',IDIM2) >*/
    swrt1_(infod, "Reading length %6d%!", &idim2, (ftnlen)80, (ftnlen)20);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);



/*<       IF(IDIM2.LT.MINMAT)THEN >*/
    if (idim2 < *minmat) {
/*<         CALL AERROR(LIST,NAMARC,1) >*/
	aerror_(list, namarc, &c__1, list_len, namarc_len);
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }
/*<       IF((IOPT.EQ.3).OR.(IOPT.EQ.4)) THEN >*/
    if (*iopt == 3 || *iopt == 4) {
/*<    >*/
	dbautp_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq2 + 1, namarc, joint, itype, isense, seqc2 + seqc2_offset, 
		itotpc, &idim2, idout, llino, lincon, ifail, idbsiz, maxdb, 
		idev1, maxgel, &imatc, &iempty, rnames + rnames_len, iopt, (
		ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL AERROR(LIST,NAMARC,3) >*/
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
/*<         ELSE >*/
	} else {
/*<           JNGEL = JNGEL + 1 >*/
	    ++jngel;
/*<         END IF   >*/
	}
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }
/*<       CALL SQCOPY(SEQ2,SEQ3,IDIM2) >*/
    sqcopy_(seq2 + 1, seq3 + 1, &idim2, (ftnlen)1, (ftnlen)1);
/*<       IFCOMP=0 >*/
    ifcomp = 0;
/*<       IMATC=0 >*/
    imatc = 0;
/*<       IFAIL(1)=0 >*/
    ifail[0] = 0;
/*<       IFAIL(2)=0 >*/
    ifail[1] = 0;
/*<       JOBC = 2 >*/
    jobc = 2;
/*<       IF ((NGELSL.LT.NGELS).AND.(NCONTL.LT.NCONTS)) THEN >*/
    if (ngelsl < *ngels && ncontl < *nconts) {
/*<         JOBC = 1 >*/
	jobc = 1;
/*<       ELSE IF (NGELSL.EQ.NGELS) THEN >*/
    } else if (ngelsl == *ngels) {
/*<         JOBC = 0 >*/
	jobc = 0;
/*<       END IF >*/
    }
/*<       NGELSL = NGELS >*/
    ngelsl = *ngels;
/*<       NCONTL = NCONTS >*/
    ncontl = *nconts;
/*<    >*/
    if (iempty == 0) {
	autocn_(seq1 + 1, &idim1, seq2 + 1, &idim2, ilefts, ilc, iposc, iposg,
		 isense, llino, &imatc, &ifcomp, minmat, maxgel, maxglm, seq5 
		+ 1, &savps[1], &savpg[1], &savl[1], maxsav, &cends[1], &
		nends[1], maxcon, seqg2 + seqg2_offset, seqc2 + seqc2_offset, 
		seq4 + 1, idout, idim22, itotpg, itotpc, joint, ifail, itype, 
		&maxpc, maxpg, permax, &minsli, seqg3 + 1, seqc3 + 1, &kfail, 
		&jobc, permis, &leno, ishow, &mask, minovr, (ftnlen)1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
		1, (ftnlen)1);
    }
/*<       IF(IREPSC.EQ.0) THEN >*/
    if (*irepsc == 0) {
/*<         IF(IFCOMP.NE.0) THEN >*/
	if (ifcomp != 0) {
/*<           CALL AERROR(LIST,NAMARC,2) >*/
	    aerror_(list, namarc, &c__2, list_len, namarc_len);
/*<           GO TO 1 >*/
	    goto L1;
/*<         END IF >*/
	}
/*<         IF(IMATC.LE.0) THEN >*/
	if (imatc <= 0) {
/*<           PERMIS(1) = 0. >*/
	    permis[0] = 0.f;
/*<           LENO = 0 >*/
	    leno = 0;
/*<         END IF >*/
	}
/*        WRITE(INFOUD,1022)NAMARC,PERMIS(1),IDIM2,LENO */
/* 1022 FORMAT(A,F5.1,2I6) */
/*<    >*/
	i__1 = i_len(namarc, namarc_len);
	swrt5_(infoud, "%.*s%5.1f%6d%6d%!", &i__1, namarc, permis, &idim2, &
		leno, (ftnlen)80, (ftnlen)17, namarc_len);
/*<         CALL TOLIST(LIST,INFOUD) >*/
	tolist_(list, infoud, list_len, (ftnlen)80);
/*<       END IF >*/
    }
/*<       IF(IOKENT.NE.0) GO TO 1 >*/
    if (*iokent != 0) {
	goto L1;
    }
/*     THIS RETURNS THE FOLLOWING: */
/*     ILEFTS  POSITION IN CONSENSUS OF LEFT END OF MATCHING CONTIGS */
/*     ILC     LENGTHS OF MATCHING CONTIGS */
/*     IPOSC   POSITION OF MATCH RELATIVE TO CONTIG */
/*     IPOSG   POSITION OF MATCH RELATIVE TO NEW GEL */
/*     ISENSE  SENSE OF NEW GEL */
/*     LLINO   LEFT GEL NUMBER IN MATCHING CONTIGS */
/*     IMATC   THE NUMBER OF MATCHING CONTIGS (>2 IS ERROR!) */
/*     IFCOMP  ERROR FLAG FOR COMPARISON (COMPARISON ARRAYS OVERFLOWED) */
/*<       IF(IFCOMP.NE.0) THEN >*/
    if (ifcomp != 0) {
/*<         CALL AERROR(LIST,NAMARC,2) >*/
	aerror_(list, namarc, &c__2, list_len, namarc_len);
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }
/*<       CALL SQCOPY(SEQ3,SEQ2,IDIM2) >*/
    sqcopy_(seq3 + 1, seq2 + 1, &idim2, (ftnlen)1, (ftnlen)1);


/* No overlap below mismatch cutoff */

/*<       IF(IMATC.EQ.0) THEN >*/
    if (imatc == 0) {

/* if masking then count as failure (code 5) unless we are entering all reads */
/* ie if we are masking and there is no match we do not enter unless "enter */
/* all reads (ANSFE=2) is set. */

/*<         IF (MASK.NE.0) THEN >*/
	if (mask != 0) {
/*<           IF (ANSFE.EQ.1) THEN >*/
	    if (*ansfe == 1) {
/*<             CALL AERROR(LIST,NAMARC,5) >*/
		aerror_(list, namarc, &c__5, list_len, namarc_len);
/*<             GO TO 1 >*/
		goto L1;
/*<           END IF >*/
	    }
/*<         END IF >*/
	}

/*                     NO OVERLAP NEW CONTIG */

/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           IF (ANSFE.EQ.1) THEN >*/
	    if (*ansfe == 1) {
/*<             CALL AERROR(LIST,NAMARC,2) >*/
		aerror_(list, namarc, &c__2, list_len, namarc_len);
/*<             GO TO 1 >*/
		goto L1;
/*<           END IF >*/
	    }
/*<           CALL INFO('New reading overlaps poorly: start a new contig') >*/
	    info_("New reading overlaps poorly: start a new contig", (ftnlen)
		    47);
/*<         ELSE >*/
	} else {
/*<           CALL INFO('New reading does not overlap: start a new contig') >*/
	    info_("New reading does not overlap: start a new contig", (ftnlen)
		    48);
/*<         END IF >*/
	}
/*     ITYPE 0 = NO OVERLAP */
/*     ISENSE 1 = SAME SENSE AS ARCHIVE */
/*<         ITYPE(1)=0 >*/
	itype[0] = 0;
/*<         ISENSE(1)=1 >*/
	isense[0] = 1;
/*<         IDOUT(1)=MAXGEL >*/
	idout[0] = *maxgel;
/*<    >*/
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq2 + 1, namarc, &x, itype, isense, seqc2 + (seqc2_dim1 + 1),
		 itotpc, &idim2, idout, llino, lincon, ifail, idbsiz, idev1, 
		maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (ftnlen)1,
		 rnames_len);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL AERROR(LIST,NAMARC,3) >*/
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
/*<           GO TO 1 >*/
	    goto L1;
/*<         END IF >*/
	}
/*<         IEMPTY=0 >*/
	iempty = 0;
/*<         IDIM1=IDIM1+1 >*/
	++idim1;

/* new start */

/*<         LREG = 1 >*/
	lreg = 1;
/*<         RREG = IDIM2 >*/
	rreg = idim2;
/*<         LINCON(1) = IDBSIZ - NCONTS >*/
	lincon[0] = *idbsiz - *nconts;
/*<    >*/
	precn1_(seq1 + 1, nampro, percd, idbsiz, lincon, &lreg, &rreg, &itask,
		 idev1, &idim1, maxgel, maxseq, iwing, nbad, iladd, iradd, 
		ifail, (ftnlen)1, nampro_len);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL ERROMF('Error calculating consensus') >*/
	    erromf_("Error calculating consensus", (ftnlen)27);
/*<           GO TO 900 >*/
	    goto L900;
/*<         END IF >*/
	}
/*         CALL FMTDB1(SEQ1,IDIM1,1,IDIM1,60,1) */
/*         CALL FMTDB(SEQ1,IDIM1,1,IDIM1,60,6) */

/* new end next bit commented out */

/*        IF((IDIM1+19+IDIM2).GT.MAXSEQ)THEN */
/*          WRITE(IDEV,1021)MAXSEQ */
/* 1021      FORMAT(' Database maximum consensus length (',I8,') exceeded') */
/*          GO TO 900 */
/*        END IF */
/*        CALL ADDTIT(SEQ1(IDIM1),NAMPRO,NGELS,IDIM1) */
/*        CALL MSTLKL(SEQ2,IDIM2) */
/*        CALL SQCOPY(SEQ2,SEQ1(IDIM1),IDIM2) */
/*        IDIM1=IDIM1+IDIM2-1 */
/*<         JNGEL = JNGEL + 1 >*/
	++jngel;
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }


/*         OVERLAP SO TRY TO ENTER THE READING */


/*<       DO 100 I=1,IMATC >*/
    i__1 = imatc;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         N=IDBSIZ-NCONTS >*/
	n = *idbsiz - *nconts;
/*<         DO 99 J=N,IDBSIZ-1 >*/
	i__2 = *idbsiz - 1;
	for (j = n; j <= i__2; ++j) {
/*<           IF(LNBR(J).NE.LLINO(I))GO TO 99 >*/
	    if (lnbr[j] != llino[i__ - 1]) {
		goto L99;
	    }
/*<           LINCON(I)=J >*/
	    lincon[i__ - 1] = j;
/*<           GO TO 100 >*/
	    goto L100;
/*< 99      CONTINUE >*/
L99:
	    ;
	}
/*        WRITE(INFOD,10077)LLINO(I) */
/* 10077   FORMAT(' Contig line for contig',I8,' not found!') */
/*<    >*/
	swrt1_(infod, " Contig line for contig%8d not found!%!", &llino[i__ - 
		1], (ftnlen)80, (ftnlen)39);
/*<         CALL ERROMF(INFOD) >*/
	erromf_(infod, (ftnlen)80);
/*<         GO TO 900 >*/
	goto L900;
/*< 100   CONTINUE >*/
L100:
	;
    }

/*<       IF (IMATC.EQ.1) THEN >*/
    if (imatc == 1) {


/*                     SINGLE OVERLAP */



/*        WRITE(INFOD,1014)LLINO(1) */
/* 1014    FORMAT('New reading overlaps contig',I8) */
/* CHECKED */
/*<         CALL SWRT1(INFOD,'New reading overlaps contig%8d%!',LLINO(1)) >*/
	swrt1_(infod, "New reading overlaps contig%8d%!", llino, (ftnlen)80, (
		ftnlen)32);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<         IF(ITOTPG(1).GT.0) CALL CCTA(SEQG2(1,1),IDIM22(1)) >*/
	if (itotpg[0] > 0) {
	    ccta_(seqg2 + (seqg2_dim1 + 1), idim22, (ftnlen)1);
	}
/*      WRITE(*,*)'BEFORE entry' */
/*      WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDIM1) */
/*<    >*/
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (seqg2_dim1 + 1), namarc, joint, itype, isense, seqc2 
		+ (seqc2_dim1 + 1), itotpc, idim22, idout, llino, lincon, 
		ifail, idbsiz, idev1, maxgel, rnames + rnames_len, (ftnlen)1, 
		namarc_len, (ftnlen)1, rnames_len);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL AERROR(LIST,NAMARC,3) >*/
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
/*<           GO TO 1 >*/
	    goto L1;
/*<         END IF >*/
	}
/*<         JNGEL = JNGEL + 1 >*/
	++jngel;
/*      WRITE(*,*)'after entry' */
/*      WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDIM1) */
/*<    >*/
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, ilefts, ilc, lincon, nampro,
		 seq2 + 1, idev1, ifail, maxgel, &idm, percd, &mask, &clist[1]
		, (ftnlen)1, nampro_len, (ftnlen)1);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL ERROMF('Error calculating consensus') >*/
	    erromf_("Error calculating consensus", (ftnlen)27);
/*<           GO TO 900 >*/
	    goto L900;
/*<         END IF >*/
	}
/*<         IF(KFAIL.NE.0) THEN >*/
	if (kfail != 0) {
/*          CALL AERROR(LIST,NAMARC,4) */
/*<         END IF >*/
	}
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }


/*                     DOUBLE OVERLAP */


/*      WRITE(INFOD,1013)LLINO */
/* 1013  FORMAT('Overlap between contigs',I8,' and',I8) */
/* CHECKED */
/*<    >*/
    swrt2_(infod, "Overlap between contigs%8d and%8d%!", llino, &llino[1], (
	    ftnlen)80, (ftnlen)35);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*<       IF(ANSJOK.NE.0) THEN >*/
    if (*ansjok != 0) {

/*  read overlaps 2 contigs but joins are forbidden */
/*  stick in in one of the contigs but do not join */

/*<         IGOOD = 1 >*/
	igood = 1;
/*<    >*/
	info_("Read overlaps 2 contigs: entering it at best site", (ftnlen)49)
		;
/*<         IF(ITOTPG(IGOOD).GT.0) CALL CCTA(SEQG2(1,IGOOD),IDIM22(IGOOD)) >*/
	if (itotpg[igood - 1] > 0) {
	    ccta_(seqg2 + (igood * seqg2_dim1 + 1), &idim22[igood - 1], (
		    ftnlen)1);
	}
/*        WRITE(INFOD,1012)LLINO(IGOOD) */
/*<    >*/
	swrt1_(infod, "Entering the new reading into contig%8d%!", &llino[
		igood - 1], (ftnlen)80, (ftnlen)41);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<    >*/
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (igood * seqg2_dim1 + 1), namarc, &joint[igood - 1], &
		itype[igood - 1], &isense[igood - 1], seqc2 + (igood * 
		seqc2_dim1 + 1), &itotpc[igood - 1], &idim22[igood - 1], &
		idout[igood - 1], &llino[igood - 1], &lincon[igood - 1], &
		ifail[igood - 1], idbsiz, idev1, maxgel, rnames + rnames_len, 
		(ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
/*<         IF(IFAIL(IGOOD).NE.0) THEN >*/
	if (ifail[igood - 1] != 0) {
/*<           CALL AERROR(LIST,NAMARC,3) >*/
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
/*<           GO TO 1 >*/
	    goto L1;
/*<         END IF >*/
	}
/*<         JNGEL = JNGEL + 1 >*/
	++jngel;
/*<    >*/
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, &ilefts[igood - 1], &ilc[
		igood - 1], &lincon[igood - 1], nampro, seq2 + 1, idev1, 
		ifail, maxgel, &idm, percd, &mask, &clist[1], (ftnlen)1, 
		nampro_len, (ftnlen)1);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL ERROMF('Error calculating consensus') >*/
	    erromf_("Error calculating consensus", (ftnlen)27);
/*<           GO TO 900 >*/
	    goto L900;
/*<         END IF >*/
	}
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }



/*<       IF(LLINO(1).EQ.LLINO(2))THEN >*/
    if (llino[0] == llino[1]) {

/*  read overlaps twice in one contig - stick it in the best place */

/*<         IGOOD = 1 >*/
	igood = 1;
/*<    >*/
	info_("Read overlaps twice in one contig: entering it at best site", (
		ftnlen)59);
/*<         IF(ITOTPG(IGOOD).GT.0) CALL CCTA(SEQG2(1,IGOOD),IDIM22(IGOOD)) >*/
	if (itotpg[igood - 1] > 0) {
	    ccta_(seqg2 + (igood * seqg2_dim1 + 1), &idim22[igood - 1], (
		    ftnlen)1);
	}
/*        WRITE(INFOD,1012)LLINO(IGOOD) */
/*<    >*/
	swrt1_(infod, "Entering the new reading into contig%8d%!", &llino[
		igood - 1], (ftnlen)80, (ftnlen)41);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<    >*/
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (igood * seqg2_dim1 + 1), namarc, &joint[igood - 1], &
		itype[igood - 1], &isense[igood - 1], seqc2 + (igood * 
		seqc2_dim1 + 1), &itotpc[igood - 1], &idim22[igood - 1], &
		idout[igood - 1], &llino[igood - 1], &lincon[igood - 1], &
		ifail[igood - 1], idbsiz, idev1, maxgel, rnames + rnames_len, 
		(ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
/*<         IF(IFAIL(IGOOD).NE.0) THEN >*/
	if (ifail[igood - 1] != 0) {
/*<           CALL AERROR(LIST,NAMARC,3) >*/
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
/*<           GO TO 1 >*/
	    goto L1;
/*<         END IF >*/
	}
/*<         JNGEL = JNGEL + 1 >*/
	++jngel;
/*<    >*/
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, &ilefts[igood - 1], &ilc[
		igood - 1], &lincon[igood - 1], nampro, seq2 + 1, idev1, 
		ifail, maxgel, &idm, percd, &mask, &clist[1], (ftnlen)1, 
		nampro_len, (ftnlen)1);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL ERROMF('Error calculating consensus') >*/
	    erromf_("Error calculating consensus", (ftnlen)27);
/*<           GO TO 900 >*/
	    goto L900;
/*<         END IF >*/
	}
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }

/* is the overlap between the contigs too large? */

/*<    >*/
    ajoin3_(&relpg[1], idbsiz, lincon, itype, isense, joint, idim22, &klass, &
	    iover, pl, pr);
/*      WRITE(INFOD,1002)IOVER */
/* 1002  FORMAT('Length of overlap between the contigs',I6) */
/* CHECKED */
/*<    >*/
    swrt1_(infod, "Length of overlap between the contigs%6d%!", &iover, (
	    ftnlen)80, (ftnlen)42);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*<       IF(IOVER.GT.MAXOVR)THEN >*/
    if (iover > maxovr) {
/*<         CALL INFO('Overlap too large: entry only') >*/
	info_("Overlap too large: entry only", (ftnlen)29);

/* cannot align the two contigs, so try to enter the reading into one of them */

/*<         IFAIL(2)=1 >*/
	ifail[1] = 1;
/*<         IGOOD=0 >*/
	igood = 0;
/*<         IF(IFAIL(1).EQ.0)IGOOD=1 >*/
	if (ifail[0] == 0) {
	    igood = 1;
	}
/*<         IF(IFAIL(2).EQ.0)IGOOD=2 >*/
	if (ifail[1] == 0) {
	    igood = 2;
	}
/*<         IF(IGOOD.EQ.0) THEN >*/
	if (igood == 0) {
/*<           CALL AERROR(LIST,NAMARC,2) >*/
	    aerror_(list, namarc, &c__2, list_len, namarc_len);
/*<           JOINF = JOINF + 1 >*/
	    ++joinf;
/*<           GO TO 1 >*/
	    goto L1;
/*<         END IF >*/
	}
/*<         IF(ITOTPG(IGOOD).GT.0) CALL CCTA(SEQG2(1,IGOOD),IDIM22(IGOOD)) >*/
	if (itotpg[igood - 1] > 0) {
	    ccta_(seqg2 + (igood * seqg2_dim1 + 1), &idim22[igood - 1], (
		    ftnlen)1);
	}
/*        WRITE(INFOD,1012)LLINO(IGOOD) */
/*<    >*/
	swrt1_(infod, "Entering the new reading into contig%8d%!", &llino[
		igood - 1], (ftnlen)80, (ftnlen)41);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<    >*/
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (igood * seqg2_dim1 + 1), namarc, &joint[igood - 1], &
		itype[igood - 1], &isense[igood - 1], seqc2 + (igood * 
		seqc2_dim1 + 1), &itotpc[igood - 1], &idim22[igood - 1], &
		idout[igood - 1], &llino[igood - 1], &lincon[igood - 1], &
		ifail[igood - 1], idbsiz, idev1, maxgel, rnames + rnames_len, 
		(ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
/*<         IF(IFAIL(IGOOD).NE.0) THEN >*/
	if (ifail[igood - 1] != 0) {
/*<           CALL AERROR(LIST,NAMARC,3) >*/
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
/*<           JOINF = JOINF + 1 >*/
	    ++joinf;
/*<           GO TO 1 >*/
	    goto L1;
/*<         END IF >*/
	}
/*<         JNGEL = JNGEL + 1 >*/
	++jngel;
/*<    >*/
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, &ilefts[igood - 1], &ilc[
		igood - 1], &lincon[igood - 1], nampro, seq2 + 1, idev1, 
		ifail, maxgel, &idm, percd, &mask, &clist[1], (ftnlen)1, 
		nampro_len, (ftnlen)1);
/*<         IF(IFAIL(1).NE.0) THEN >*/
	if (ifail[0] != 0) {
/*<           CALL ERROMF('Error calculating consensus') >*/
	    erromf_("Error calculating consensus", (ftnlen)27);
/*<           GO TO 900 >*/
	    goto L900;
/*<         END IF >*/
	}
/*        WRITE(INFOD,1020)LLINO */
/*<    >*/
	swrt2_(infod, "Could not join contigs%8d and%8d%!", llino, &llino[1], 
		(ftnlen)80, (ftnlen)34);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/* 1020    FORMAT('Could not join contigs',I8,' and',I8) */
/* 1021    FORMAT('Reading has been entered into contig',I8) */
/*        WRITE(INFOD,1021)LLINO(IGOOD) */
/*<    >*/
	swrt1_(infod, "Reading has been entered into contig%8d%!", &llino[
		igood - 1], (ftnlen)80, (ftnlen)41);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<         JOINF = JOINF + 1 >*/
	++joinf;
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }
/*   WHICH CONTIG IS LEFTMOST? */
/*<       LMOST=1 >*/
    lmost = 1;
/*<       RMOST=2 >*/
    rmost = 2;
/*<       IF(PL(1).GT.PL(2))THEN >*/
    if (pl[0] > pl[1]) {
/*<         LMOST=2 >*/
	lmost = 2;
/*<         RMOST=1 >*/
	rmost = 1;
/*<       END IF >*/
    }
/*   SAVE LENGTH OF RMOST CONTIG FOR DELETION STEP LATER */
/*<       ILCR=ILC(RMOST) >*/
    ilcr = ilc[rmost - 1];
/*<       IF(ITOTPG(LMOST).GT.0) CALL CCTA(SEQG2(1,LMOST),IDIM22(LMOST)) >*/
    if (itotpg[lmost - 1] > 0) {
	ccta_(seqg2 + (lmost * seqg2_dim1 + 1), &idim22[lmost - 1], (ftnlen)1)
		;
    }
/*      WRITE(INFOD,1012)LLINO(LMOST) */
/* CHECKED */
/*<    >*/
    swrt1_(infod, "Entering the new reading into contig%8d%!", &llino[lmost - 
	    1], (ftnlen)80, (ftnlen)41);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/* 1012  FORMAT('Entering the new reading into contig',I8) */
/*<    >*/
    aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, seqg2 + 
	    (lmost * seqg2_dim1 + 1), namarc, &joint[lmost - 1], &itype[lmost 
	    - 1], &isense[lmost - 1], seqc2 + (lmost * seqc2_dim1 + 1), &
	    itotpc[lmost - 1], &idim22[lmost - 1], &idout[lmost - 1], &llino[
	    lmost - 1], &lincon[lmost - 1], &ifail[lmost - 1], idbsiz, idev1, 
	    maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (ftnlen)1, 
	    rnames_len);
/*<       IF(IFAIL(LMOST).NE.0) THEN >*/
    if (ifail[lmost - 1] != 0) {
/*<         CALL AERROR(LIST,NAMARC,3) >*/
	aerror_(list, namarc, &c__3, list_len, namarc_len);
/*<         JOINF = JOINF + 1 >*/
	++joinf;
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }
/*<       JNGEL = JNGEL + 1 >*/
    ++jngel;
/*<    >*/
    updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, nconts, 
	    seq1 + 1, maxseq, &idim1, &ilefts[lmost - 1], &ilc[lmost - 1], &
	    lincon[lmost - 1], nampro, seq2 + 1, idev1, ifail, maxgel, &idm, 
	    percd, &mask, &clist[1], (ftnlen)1, nampro_len, (ftnlen)1);
/*<       IF(IFAIL(1).NE.0) THEN >*/
    if (ifail[0] != 0) {
/*<         CALL ERROMF('Error calculating consensus') >*/
	erromf_("Error calculating consensus", (ftnlen)27);
/*<         GO TO 900 >*/
	goto L900;
/*<       END IF >*/
    }
/*<       IF(ITYPE(LMOST).EQ.1)LLINO(LMOST)=NGELS >*/
    if (itype[lmost - 1] == 1) {
	llino[lmost - 1] = *ngels;
    }
/*<       IF(ILEFTS(LMOST).LT.ILEFTS(RMOST))THEN >*/
    if (ilefts[lmost - 1] < ilefts[rmost - 1]) {
/*<    >*/
	ilefts[rmost - 1] += relpg[lincon[lmost - 1]] - ilc[lmost - 1];
/*<       END IF >*/
    }
/*<       ILC(LMOST) =  RELPG(LINCON(LMOST)) >*/
    ilc[lmost - 1] = relpg[lincon[lmost - 1]];
/*<       DO 500 I=1,2 >*/
    for (i__ = 1; i__ <= 2; ++i__) {
/*<         IF(ISENSE(I).EQ.-1)THEN >*/
	if (isense[i__ - 1] == -1) {
/*<    >*/
	    cmplmt_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    &lincon[i__ - 1], &llino[i__ - 1], seq2 + 1, idbsiz, 
		    idev1, maxgel, (ftnlen)1);
/*<           CALL SQREV(SEQ1(ILEFTS(I)),ILC(I)) >*/
	    sqrev_(seq1 + ilefts[i__ - 1], &ilc[i__ - 1], (ftnlen)1);
/*<           CALL SQCOMM(SEQ1(ILEFTS(I)),ILC(I)) >*/
	    sqcomm_(seq1 + ilefts[i__ - 1], &ilc[i__ - 1], (ftnlen)1);
/*<           KT=IDIM1 >*/
	    kt = idim1;
/*<           CALL ADDTIT(SEQ1((ILEFTS(I)-20)),NAMPRO,LNBR(LINCON(I)),KT) >*/
	    addtit_(seq1 + (ilefts[i__ - 1] - 20), nampro, &lnbr[lincon[i__ - 
		    1]], &kt, (ftnlen)1, nampro_len);
/*<         END IF >*/
	}
/*< 500   CONTINUE >*/
/* L500: */
    }
/*   NEED TO KNOW POSITION OF OVERLAP RELATIVE TO CONTIG, TO CONSENSUS */
/*   WHICH BITS TO SEND TO ALIGNMENT ROUTINES */
/*   SET UP FOR ALINE (NOTE RMOST IS EQUIVALENT TO THE GEL READING AND */
/*   SO IS SLID ALONG THE LMOST CONTIG. THE SECTION SENT TO ALINE MUST */
/*   BE OF LENGTH < MAXGEL-2*MAX(MAXPC,MAXPG) */
/*   IT MUST START AT POSITION 1 IN THE RMOST CONTIG AND EXTEND */
/*<       IPOSC(LMOST)=PL(RMOST)+RELPG(NGELS)-1 >*/
    iposc[lmost - 1] = pl[rmost - 1] + relpg[*ngels] - 1;
/*<    >*/
    ilct = relpg[lincon[lmost - 1]] - relpg[*ngels] - pl[rmost - 1] + 2 + 
	    maxpc;

/* change 5-6-95 line below to line above */

/*      ILCT = RELPG(LINCON(LMOST)) - RELPG(NGELS) - PL(RMOST) + 2 */
/*<       ILC(RMOST)=MIN(ILCT,ILC(RMOST)) >*/
/* Computing MIN */
    i__1 = ilct, i__2 = ilc[rmost - 1];
    ilc[rmost - 1] = min(i__1,i__2);
/*      WRITE(*,*)'ILC(LMOST)',ILC(LMOST) */
/*      WRITE(*,*)'RELPG(LINCON(LMOST))',RELPG(LINCON(LMOST)) */
/*      WRITE(*,*)'RELPG(NGELS)',RELPG(NGELS) */
/*      WRITE(*,*)'PL(RMOST)',PL(RMOST) */
/*      WRITE(*,*)'LENGTH OF OVERLAP',ILC(RMOST) */
/*<       IPOSC(RMOST)=1 >*/
    iposc[rmost - 1] = 1;
/*<       IDOUT(LMOST)=MAXGEL >*/
    idout[lmost - 1] = *maxgel;
/*<       IDOUT(RMOST)=MAXGEL >*/
    idout[rmost - 1] = *maxgel;
/*<       IDSAV=MAXSAV >*/
    idsav = *maxsav;
/*  ON INPUT TO ALINE ILC(RMOST) CONTAINS THE OVERLAP LENGTH */
/*  ON OUTPUT IT CONTAINS THE LENGTH OF THE ALIGNED SECTION (IE INCLUDING */
/*  PADS) */
/*<       CALL INFO('Trying to align the two contigs') >*/
    info_("Trying to align the two contigs", (ftnlen)31);
/*<    >*/
    aline_(seq1 + ilefts[lmost - 1], seq1 + ilefts[rmost - 1], seqc2 + (rmost 
	    * seqc2_dim1 + 1), seqc2 + (lmost * seqc2_dim1 + 1), &savps[1], &
	    savpg[1], &savl[1], &idsav, &ilc[lmost - 1], &ilc[rmost - 1], &
	    idout[lmost - 1], &iposc[lmost - 1], &iposc[rmost - 1], &minsli, &
	    joint[lmost - 1], &itotpc[lmost - 1], &itotpc[rmost - 1], ifail, 
	    itype, &maxpc, &maxpc, permax, seq4 + 1, maxgel, &z__, &leno, 
	    ishow, &mask, &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1);
/* SEQC2(1,LMOST)  NOW CONTAINS THE ALIGNED SECTION OF THE LMOST CONTIG */
/* SEQC2(1,RMOST)  NOW CONTAINS THE ALIGNED SECTION OF THE RMOST CONTIG */
/* ILC(RMOST)  IS NOW THE LENGTH OF ALIGNED SECTION OF THE RMOST CONTIG */
/* IDOUT(LMOST)  IS NOW THE LENGTH OF ALIGNED SECTION OF THE LMOST CONTIG */
/* JOINT(LMOST)  IS THE POSITION OF THE JOIN RLETIVE TO THE LMOST CONTIG */
/* ITYPE IS TYPE OF OVERLAP (-1 = RIGHT END OR INTERNAL, 1 = LEFT END) */
/*  NB SHOULD ALWAYS BE -1 */
/*  IF THIS HAS BEEN DONE OK WE CAN EDIT THE TWO CONTIGS THEN JOIN */
/*<       IF(IFAIL(1).NE.0)THEN >*/
    if (ifail[0] != 0) {
/*<         CALL INFO('Failed to align the two overlapping contigs') >*/
	info_("Failed to align the two overlapping contigs", (ftnlen)43);
/*        CALL AERROR(LIST,NAMARC,4) */
/*<         JOINF = JOINF + 1 >*/
	++joinf;
/*<         GO TO 1 >*/
	goto L1;
/*<       END IF >*/
    }
/*<       IF(ITOTPC(LMOST).GT.0)THEN >*/
    if (itotpc[lmost - 1] > 0) {
/*        WRITE(INFOD,1017)LLINO(LMOST) */
/* CHECKED */
/*<         CALL SWRT1(INFOD,'Editing contig%8d%!',LLINO(LMOST)) >*/
	swrt1_(infod, "Editing contig%8d%!", &llino[lmost - 1], (ftnlen)80, (
		ftnlen)19);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/* 1017    FORMAT('Editing contig',I8) */
/*<    >*/
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq3 + 1, &lincon[lmost - 1], &joint[lmost - 1], seqc2 + (
		lmost * seqc2_dim1 + 1), &itotpc[lmost - 1], &idout[lmost - 1]
		, idbsiz, idev1, maxgel, (ftnlen)1, (ftnlen)1);
/*<       END IF >*/
    }
/*<       JOINT(RMOST)=1 >*/
    joint[rmost - 1] = 1;
/*<       IDOUT(RMOST)=ILC(RMOST) >*/
    idout[rmost - 1] = ilc[rmost - 1];
/*<       IF(ITOTPC(RMOST).GT.0)THEN >*/
    if (itotpc[rmost - 1] > 0) {
/*        WRITE(INFOD,1017)LLINO(RMOST) */
/* CHECKED */
/*<         CALL SWRT1(INFOD,'Editing contig%8d%!',LLINO(RMOST)) >*/
	swrt1_(infod, "Editing contig%8d%!", &llino[rmost - 1], (ftnlen)80, (
		ftnlen)19);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<    >*/
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq3 + 1, &lincon[rmost - 1], &joint[rmost - 1], seqc2 + (
		rmost * seqc2_dim1 + 1), &itotpc[rmost - 1], &idout[rmost - 1]
		, idbsiz, idev1, maxgel, (ftnlen)1, (ftnlen)1);
/*<       END IF >*/
    }
/*<       ILC(RMOST)=ILCR >*/
    ilc[rmost - 1] = ilcr;
/*<       LTL=LNBR(LINCON(LMOST)) >*/
    ltl = lnbr[lincon[lmost - 1]];
/*<       LTR=LNBR(LINCON(RMOST)) >*/
    ltr = lnbr[lincon[rmost - 1]];
/*      WRITE(INFOD,1018)LNBR(LINCON(LMOST)),LNBR(LINCON(RMOST)) */
/* CHECKED */
/*<    >*/
    swrt2_(infod, "Completing the join between contigs%8d and%8d%!", &lnbr[
	    lincon[lmost - 1]], &lnbr[lincon[rmost - 1]], (ftnlen)80, (ftnlen)
	    47);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/* 1018  FORMAT('Completing the join between contigs',I8,' and',I8) */
/*<    >*/
    ajoin2_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    &joint[lmost - 1], &ltl, &ltr, &lincon[lmost - 1], &lincon[rmost 
	    - 1], idev1);
/*<       LLINO(1)=LTL >*/
    llino[0] = ltl;
/*<       IF(ILEFTS(LMOST).GT.ILEFTS(RMOST))THEN >*/
    if (ilefts[lmost - 1] > ilefts[rmost - 1]) {
/*<         CALL DELCON(SEQ1,ILEFTS(LMOST),ILC(LMOST),IDIM1) >*/
	delcon_(seq1 + 1, &ilefts[lmost - 1], &ilc[lmost - 1], &idim1, (
		ftnlen)1);
/*<         CALL DELCON(SEQ1,ILEFTS(RMOST),ILC(RMOST),IDIM1) >*/
	delcon_(seq1 + 1, &ilefts[rmost - 1], &ilc[rmost - 1], &idim1, (
		ftnlen)1);
/*<       END IF >*/
    }
/*<       IF(ILEFTS(RMOST).GE.ILEFTS(LMOST))THEN >*/
    if (ilefts[rmost - 1] >= ilefts[lmost - 1]) {
/*<         CALL DELCON(SEQ1,ILEFTS(RMOST),ILC(RMOST),IDIM1) >*/
	delcon_(seq1 + 1, &ilefts[rmost - 1], &ilc[rmost - 1], &idim1, (
		ftnlen)1);
/*<         CALL DELCON(SEQ1,ILEFTS(LMOST),ILC(LMOST),IDIM1) >*/
	delcon_(seq1 + 1, &ilefts[lmost - 1], &ilc[lmost - 1], &idim1, (
		ftnlen)1);
/*<       END IF >*/
    }
/*<       LREG=1 >*/
    lreg = 1;
/*<       RREG=JOINT(LMOST) >*/
    rreg = joint[lmost - 1];
/*<       IGELC=LLINO(1) >*/
    igelc = llino[0];
/*      JOB = 1 */

/* assume itask set above */

/* also need to add 1 to consensus position */

/*<       IDIM1 = IDIM1 + 1 >*/
    ++idim1;

/* for precon need to find the contig line number after the join */

/*<    >*/
    lincon[lmost - 1] = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], 
	    ngels, nconts, idbsiz, &igelc);
/*<       IF (LINCON(LMOST).EQ.0) THEN >*/
    if (lincon[lmost - 1] == 0) {
/*<         CALL ERROMF('Cannot find contig line! Quitting') >*/
	erromf_("Cannot find contig line! Quitting", (ftnlen)33);
/*<         GO TO 900 >*/
	goto L900;
/*<       END IF >*/
    }
/*<    >*/
    precn1_(seq1 + 1, nampro, percd, idbsiz, &lincon[lmost - 1], &lreg, &rreg,
	     &itask, idev1, &idim1, maxgel, maxseq, iwing, nbad, iladd, iradd,
	     ifail, (ftnlen)1, nampro_len);
/*<       IF(IFAIL(1).NE.0) THEN >*/
    if (ifail[0] != 0) {
/*<         CALL ERROMF('Error calculating consensus') >*/
	erromf_("Error calculating consensus", (ftnlen)27);
/*<         GO TO 900 >*/
	goto L900;
/*<       END IF >*/
    }
/*<       JNJOIN = JNJOIN + 1 >*/
    ++jnjoin;
/*<       IF(KFAIL.NE.0) THEN >*/
    if (kfail != 0) {
/*        CALL AERROR(LIST,NAMARC,4) */
/*        JOINF = JOINF + 1 */
/*<       END IF >*/
    }
/*<       GO TO 1 >*/
    goto L1;
/*< 900   CONTINUE >*/
L900:
/*<       IF ((IOPT.EQ.1).OR.(IOPT.EQ.2).OR.(IOPT.EQ.5)) THEN >*/
    if (*iopt == 1 || *iopt == 2 || *iopt == 5) {
/*<       IOPTC = 6 >*/
	ioptc = 6;
/*<    >*/
	idsav = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1], 
		&idsav, seq1 + 1, seq2 + 1, maxseq, maxgel, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/*<       END IF >*/
    }
/*< 901   CONTINUE >*/
L901:
/*<       CALL INFO('Batch finished') >*/
    info_("Batch finished", (ftnlen)14);
/*      WRITE(INFOD,1030)JGEL */
/* 1030 FORMAT(I8,' sequences processed') */
/* CHECKED */
/*<       CALL SWRT1(INFOD,'%8d sequences processed%!',JGEL) >*/
    swrt1_(infod, "%8d sequences processed%!", &jgel, (ftnlen)80, (ftnlen)25);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*      WRITE(INFOD,1031)JNGEL */
/* 1031 FORMAT(I8,' sequences entered into database') */
/*<       CALL SWRT1(INFOD,'%8d sequences entered into database%!',JNGEL) >*/
    swrt1_(infod, "%8d sequences entered into database%!", &jngel, (ftnlen)80,
	     (ftnlen)37);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*      WRITE(INFOD,1032)JNJOIN */
/* 1032 FORMAT(I8,' joins made') */
/*<       CALL SWRT1(INFOD,'%8d joins made%!',JNJOIN) >*/
    swrt1_(infod, "%8d joins made%!", &jnjoin, (ftnlen)80, (ftnlen)16);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*      WRITE(INFOD,1033)JOINF */
/* 1033 FORMAT(I8,' joins failed') */
/* CHECKED */
/*<       CALL SWRT1(INFOD,'%8d joins failed%!',JOINF) >*/
    swrt1_(infod, "%8d joins failed%!", &joinf, (ftnlen)80, (ftnlen)18);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*<       END >*/
    return 0;
} /* dbauto_ */

/*<    >*/
/* Subroutine */ int dbautp_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *ngels, integer *nconts, char *seq2, char *
	namarc, integer *joint, integer *itype, integer *isense, char *seqc2, 
	integer *itotpc, integer *idim2, integer *idout, integer *llino, 
	integer *lincon, integer *ifail, integer *idbsiz, integer *maxdb, 
	integer *idev1, integer *maxgel, integer *imatc, integer *iempty, 
	char *rnames, integer *iopt, ftnlen seq2_len, ftnlen namarc_len, 
	ftnlen seqc2_len, ftnlen rnames_len)
{
    extern /* Subroutine */ int aenter_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, char *, integer *, 
	    integer *, integer *, char *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    char *, ftnlen, ftnlen, ftnlen, ftnlen);

/*<       INTEGER RELPG(MAXDB) >*/
/*<       INTEGER LNGTHG(MAXDB),LNBR(MAXDB),RNBR(MAXDB) >*/
/*<       CHARACTER SEQ2(MAXGEL),SEQC2(MAXGEL) >*/
/*<       CHARACTER NAMARC*(*) >*/
/*<       CHARACTER*(*) RNAMES(IDBSIZ) >*/
/*  deals with entering all readings into contig 1 (IOPT=3) */
/*  or all readings into new contigs (IOPT=4) */
/*<       IF(IOPT.EQ.3) THEN >*/
    /* Parameter adjustments */
    rnames -= rnames_len;
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --seqc2;
    --seq2;

    /* Function Body */
    if (*iopt == 3) {
/*<         IF(IMATC.EQ.0) THEN >*/
	if (*imatc == 0) {
/*<           ITYPE=0 >*/
	    *itype = 0;
/*<           ISENSE=1 >*/
	    *isense = 1;
/*<           IDOUT=MAXGEL >*/
	    *idout = *maxgel;
/*<    >*/
	    aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    seq2 + 1, namarc, joint, itype, isense, seqc2 + 1, itotpc,
		     idim2, idout, llino, lincon, ifail, idbsiz, idev1, 
		    maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (
		    ftnlen)1, rnames_len);
/*<           IF(IFAIL.NE.0) RETURN >*/
	    if (*ifail != 0) {
		return 0;
	    }
/*<           IEMPTY=0 >*/
	    *iempty = 0;
/*<           IMATC = 1 >*/
	    *imatc = 1;
/*<         ELSE >*/
	} else {
/*<           ITYPE= - 1 >*/
	    *itype = -1;
/*<           ISENSE=1 >*/
	    *isense = 1;
/*<           JOINT = 1 >*/
	    *joint = 1;
/*<           LLINO = NGELS >*/
	    *llino = *ngels;
/*<           LINCON = IDBSIZ - NCONTS >*/
	    *lincon = *idbsiz - *nconts;
/*<           ITOTPC = 0 >*/
	    *itotpc = 0;
/*<           IDOUT=MAXGEL >*/
	    *idout = *maxgel;
/*<    >*/
	    aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    seq2 + 1, namarc, joint, itype, isense, seqc2 + 1, itotpc,
		     idim2, idout, llino, lincon, ifail, idbsiz, idev1, 
		    maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (
		    ftnlen)1, rnames_len);
/*<           IF(IFAIL.NE.0) RETURN >*/
	    if (*ifail != 0) {
		return 0;
	    }
/*<         END IF >*/
	}
/*<       ELSE IF(IOPT.EQ.4) THEN >*/
    } else if (*iopt == 4) {
/*<         ITYPE=0 >*/
	*itype = 0;
/*<         ISENSE=1 >*/
	*isense = 1;
/*<         IDOUT=MAXGEL >*/
	*idout = *maxgel;
/*<    >*/
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq2 + 1, namarc, joint, itype, isense, seqc2 + 1, itotpc, 
		idim2, idout, llino, lincon, ifail, idbsiz, idev1, maxgel, 
		rnames + rnames_len, (ftnlen)1, namarc_len, (ftnlen)1, 
		rnames_len);
/*<         IF(IFAIL.NE.0) RETURN >*/
	if (*ifail != 0) {
	    return 0;
	}
/*<       END IF >*/
    }
/*<       END >*/
    return 0;
} /* dbautp_ */

/*   SUBROUTINE TO ENTER NEW GEL SEQUENCES INTO DATA BASE. */
/*   IT READS IN AN ARCHIVE VERSION AND WRITES OUT A WORKING VERSION. */
/*   IT ALSO SETS UP ANY RELATIONSHIPS WITH OTHER DATA IN THE DATABASE */
/*   BOTH BY POSITION IN A CONTIG AND POINTERS TO LEFT AND RIGHT */
/*   NEIGHBOURS. */
/*<    >*/
/* Subroutine */ int aenter_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *ngels, integer *nconts, char *gel, char *
	namarc, integer *x, integer *itype, integer *isense, char *seqc2, 
	integer *itotpc, integer *idim, integer *idc, integer *ncontc, 
	integer *lincon, integer *ifail, integer *idbsiz, integer *idevr, 
	integer *maxgel, char *rnames, ftnlen gel_len, ftnlen namarc_len, 
	ftnlen seqc2_len, ftnlen rnames_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, n, y, iok;
    extern integer indb_(integer *, char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int info_(char *, ftnlen);
    static integer itmp;
    extern /* Subroutine */ int swrt1_(char *, char *, ...);
    static char namid[40], infod[80];
    extern /* Subroutine */ int sindb_(integer *, integer *, char *, char *, 
	    integer *, ftnlen, ftnlen);
    static integer idevn;
    extern /* Subroutine */ int abedin_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, integer *, integer *, 
	    char *, integer *, integer *, integer *, integer *, integer *, 
	    ftnlen, ftnlen), idline_(char *, char *, ftnlen, ftnlen), erromf_(
	    char *, ftnlen), writec_(integer *, integer *, integer *, integer 
	    *, integer *), writeg_(integer *, integer *, integer *, integer *,
	     integer *, integer *), shiftt_(integer *, integer *, integer *, 
	    integer *), stikit_(integer *, char *, integer *, integer *, char 
	    *, integer *, integer *, integer *, integer *, ftnlen, ftnlen), 
	    writrn_(integer *, integer *, integer *);

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER  RELPG(IDBSIZ),X,Y >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       CHARACTER GEL(MAXGEL),NAMARC*(*) >*/
/*<       CHARACTER SEQC2(IDC) >*/
/*<       CHARACTER NAMID*40,RNAMES(IDBSIZ)*40 >*/
/*<       CHARACTER INFOD*80 >*/
/*<       EXTERNAL INDB >*/
/*      WRITE(*,*)'IN ENTER',NGELS */
/*      WRITE(*,*)'X,ITYPE,ISENSE,IDIM,IDC' */
/*      WRITE(*,*)X,ITYPE,ISENSE,IDIM,IDC */
/*   SET FAIL FLAG */
/*<       IFAIL=0 >*/
    /* Parameter adjustments */
    --seqc2;
    rnames -= 40;
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */
    *ifail = 0;
/*   IS THERE SPACE? */
/*<       IF((IDBSIZ-(NGELS+NCONTS)).GT.2)GO TO 5 >*/
    if (*idbsiz - (*ngels + *nconts) > 2) {
	goto L5;
    }
/*   FULL */
/*<       CALL ERROMF('Database full!') >*/
    erromf_("Database full!", (ftnlen)14);
/*<       IFAIL=7 >*/
    *ifail = 7;
/*<       RETURN >*/
    return 0;
/*< 5     CONTINUE >*/
L5:
/*   NEED TO CHECK TO SEE IF GEL ALREADY IN DB */
/*   LOOK THRU ARC FILE */
/*<       CALL IDLINE(NAMARC, NAMID) >*/
    idline_(namarc, namid, namarc_len, (ftnlen)40);
/*<       J = INDB(NGELS,RNAMES,NAMID) >*/
    j = indb_(ngels, rnames + 40, namid, (ftnlen)40, (ftnlen)40);
/*<       IF (J.NE.0) THEN >*/
    if (j != 0) {
/*   FOUND */
/*        WRITE(INFOD,1013)J */
/* 1013    FORMAT('New reading already in database with number',I8, */
/*     +  ' Entry aborted') */
/* CHECKED */
/*<    >*/
	swrt1_(infod, "New reading already in database with number%8d Entry "
		"aborted%!", &j, (ftnlen)80, (ftnlen)62);
/*<         CALL ERROMF(INFOD) >*/
	erromf_(infod, (ftnlen)80);
/*<         IFAIL=6 >*/
	*ifail = 6;
/*<         RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*   INCREMENT NUMBER OF GELS */
/*<       NGELS=NGELS+1 >*/
    ++(*ngels);

/* set dummy int for idevn */

/*<       IDEVN = 0 >*/
    idevn = 0;
/*<       CALL SINDB(IDEVN,NGELS,RNAMES,NAMID,2) >*/
    sindb_(&idevn, ngels, rnames + 40, namid, &c__2, (ftnlen)40, (ftnlen)40);
/*   SET LENGTH THIS GEL */
/*<       LNGTHG(NGELS)=IDIM*ISENSE >*/
    lngthg[*ngels] = *idim * *isense;
/*      WRITE(INFOD,1003)NGELS */
/*      WRITE(*,1003)NGELS */
/* 1003  FORMAT('This gel reading has been given the number ',I8) */
/* CHECKED */
/*<    >*/
    swrt1_(infod, "This gel reading has been given the number %8d%!", ngels, (
	    ftnlen)80, (ftnlen)48);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*   WRITE NAME OF ARCHIVE TO LIST OF ARCHIVES */
/*   NAMPRO,ARC */
/*      NAMARK=NAMARC(1:16) */
/*      CALL WRITEN(IDEVN,NGELS,NAMARK) */
/* C   WRITE GEL TO WORKING VERSION */
/*      CALL WRITEW(IDEVW,NGELS,GEL,MAXGEL) */
/*      IF(IDEVT.GT.0) CALL ENTRD(IDEVG,IDEVT,IDEVC,NAMARC,NGELS,IOK) */
/* C   CREATE TAGS FOR THIS NASTY */
/*      CALL TAGGEL(NGELS,LNGTHG(NGELS),GEL) */
/*   SET UP RELATIONSHIPS */
/*   DOES THIS GEL OVERLAP? */
/*<       IF(ITYPE.NE.0)GO TO 100 >*/
    if (*itype != 0) {
	goto L100;
    }

/*   DOES NOT OVERLAP SO IT STARTS A CONTIG OF ITS OWN */

/*   SET CONTIG POINTERS AND GENERAL VALUES */
/*   INCREMENT NUMBER OF CONTIGS */
/*<       NCONTS=NCONTS+1 >*/
    ++(*nconts);
/*   POINTER TO THIS CONTIG */
/*<       N=IDBSIZ-NCONTS >*/
    n = *idbsiz - *nconts;
/*   POINTER TO LEFT GEL THIS CONTIG */
/*<       LNBR(N)=NGELS >*/
    lnbr[n] = *ngels;
/*   POINTER TO RIGHT GEL THIS CONTIG */
/*<       RNBR(N)=NGELS >*/
    rnbr[n] = *ngels;
/*   LENGTH OF CONTIG */
/*<       RELPG(N)=IDIM >*/
    relpg[n] = *idim;
/*   WRITE CONTIG DESCRIPTOR */
/*<    >*/
    i__1 = *idbsiz - n;
    writec_(idevr, &i__1, &relpg[n], &lnbr[n], &rnbr[n]);
/*     Setup tags, original positions, conf values, vectors etc */
/*<    >*/
    i__1 = *idbsiz - n;
    stikit_(idevr, namarc, ngels, &lngthg[*ngels], gel + 1, maxgel, &iok, &
	    i__1, &c__1, namarc_len, (ftnlen)1);
/*<       IF (IOK.NE.0) THEN >*/
    if (iok != 0) {
/*<          NCONTS=NCONTS-1 >*/
	--(*nconts);
/*<          NGELS=NGELS-1 >*/
	--(*ngels);
/*<          IFAIL=1 >*/
	*ifail = 1;
/*<          RETURN >*/
	return 0;
/*<       ENDIF >*/
    }

/*   Create gel info */
/*   SET LEFT AND RIGHT POINTERS TO ZERO,RELPG TO 1 */
/*<       LNBR(NGELS)=0 >*/
    lnbr[*ngels] = 0;
/*<       RNBR(NGELS)=0 >*/
    rnbr[*ngels] = 0;
/*<       RELPG(NGELS)=1 >*/
    relpg[*ngels] = 1;
/*   WRITE NEW GEL LINE */
/*<    >*/
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   WRITE DB DESCRIPTOR */
/*<        CALL WRITRN(IDEVR,NGELS,NCONTS) >*/
    writrn_(idevr, ngels, nconts);
/*<       RETURN >*/
    return 0;

/*< 100   CONTINUE >*/
L100:


/*     Shift tags if this new gel adjusts the left end */

/*<       IF (ITYPE.EQ.1) THEN >*/
    if (*itype == 1) {
/*<          CALL SHIFTT(IDEVR, IDBSIZ-LINCON, 1, X-1) >*/
	i__1 = *idbsiz - *lincon;
	i__2 = *x - 1;
	shiftt_(idevr, &i__1, &c__1, &i__2);
/*<          ITMP = 1 >*/
	itmp = 1;
/*<       ELSE >*/
    } else {
/*<          ITMP = X >*/
	itmp = *x;
/*<       ENDIF >*/
    }
/*     Setup tags, original positions, conf values, vectors etc */
/*<    >*/
    i__1 = *idbsiz - *lincon;
    stikit_(idevr, namarc, ngels, &lngthg[*ngels], gel + 1, maxgel, &iok, &
	    i__1, &itmp, namarc_len, (ftnlen)1);
/*<       IF (IOK.NE.0) THEN >*/
    if (iok != 0) {
/*<          NGELS=NGELS-1 >*/
	--(*ngels);
/*<          IFAIL=1 >*/
	*ifail = 1;
/*<          RETURN >*/
	return 0;
/*<       ENDIF >*/
    }

/*   DOES OVERLAP */
/*< 150   CONTINUE >*/
/* L150: */

/*   LEFT END OR RIGHT OVERLAP? */
/*<       IF(ITYPE.EQ.1)GO TO 400 >*/
    if (*itype == 1) {
	goto L400;
    }
/*   RIGHT END OR INTERNAL OVERLAP */

/*< 160   CONTINUE >*/
/* L160: */
/*   NEED TO SEARCH THRU THIS CONTIG TO FIND LEFT AND RIGHT */
/*   NEIGHBOURS FOR THIS NEW GEL */
/*   LINE NUMBER OF LEFT END OF CONTIG */
/*<       N=NCONTC >*/
    n = *ncontc;
/*   LOOK THRU UNTIL CURRENT IS >= THEN IT MUST BE THE PREVIOUS ONE */
/*< 200   CONTINUE >*/
L200:
/*<       IF(RELPG(N).GT.X)GO TO 250 >*/
    if (relpg[n] > *x) {
	goto L250;
    }
/*   IS THIS THE LAST GEL IN CONTIG? */
/*<       IF(RNBR(N).EQ.0)GO TO 350 >*/
    if (rnbr[n] == 0) {
	goto L350;
    }
/*   NO SO LOOK AT NEXT */
/*<       N=RNBR(N) >*/
    n = rnbr[n];
/*<       GO TO 200 >*/
    goto L200;
/*< 250   CONTINUE >*/
L250:
/*   GEL LIES BETWEEN N AND LNBR(N) */
/*   NEED TO EDIT DB HERE */
/*<    >*/
    if (*itotpc > 0) {
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, x, seqc2 + 1, itotpc, idc, idbsiz, idevr, maxgel,
		 (ftnlen)1, (ftnlen)1);
    }


/*   SET POINTERS IN NEW GEL */
/*<       LNBR(NGELS)=LNBR(N) >*/
    lnbr[*ngels] = lnbr[n];
/*<       RNBR(NGELS)=N >*/
    rnbr[*ngels] = n;
/*<       RELPG(NGELS)=X >*/
    relpg[*ngels] = *x;
/*   WRITE NEW GEL LINE */
/*<    >*/
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   SET POINTERS  IN LEFT AND RIGHT NEIGHBOURS */
/*<       K=LNBR(N) >*/
    k = lnbr[n];
/*<       RNBR(K)=NGELS >*/
    rnbr[k] = *ngels;
/*      RNBR(LNBR(N))=NGELS */
/*   WRITE LEFT AND RIGHT NEIGHBOURS */
/*<    >*/
    writeg_(idevr, &k, &relpg[k], &lngthg[k], &lnbr[k], &rnbr[k]);
/*<       LNBR(N)=NGELS >*/
    lnbr[n] = *ngels;
/*<    >*/
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
/*   WRITE NGELS NCONTS */
/*<       CALL WRITRN(IDEVR,NGELS,NCONTS) >*/
    writrn_(idevr, ngels, nconts);
/*   HAVE WE INCREASED LENGTH OF CONTIG? */
/*   ITS LINE NUMBER IS LINCON */
/*   NEED TO UPDATE IDIM IN CASE OF EDITS */
/*<       IDIM=ABS(LNGTHG(NGELS)) >*/
    *idim = (i__1 = lngthg[*ngels], abs(i__1));
/*<       Y=X+IDIM-1 >*/
    y = *x + *idim - 1;
/*<       IF(Y.LE.RELPG(LINCON))RETURN >*/
    if (y <= relpg[*lincon]) {
	return 0;
    }
/*<       RELPG(LINCON)=Y >*/
    relpg[*lincon] = y;
/*   WRITE NEW CONTIG LINE */
/*<    >*/
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
/*<       RETURN >*/
    return 0;
/*< 350   CONTINUE >*/
L350:
/*   MUST BE A RIGHT END OVERLAP */
/*   NEED TO EDIT DB HERE */
/*<    >*/
    if (*itotpc > 0) {
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, x, seqc2 + 1, itotpc, idc, idbsiz, idevr, maxgel,
		 (ftnlen)1, (ftnlen)1);
    }


/*   SET POINTERS FOR NEW GEL */
/*<       LNBR(NGELS)=N >*/
    lnbr[*ngels] = n;
/*<       RNBR(NGELS)=0 >*/
    rnbr[*ngels] = 0;
/*<       RELPG(NGELS)=X >*/
    relpg[*ngels] = *x;
/*   WRITE NEW GEL LINE */
/*<    >*/
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   OLD RIGHT END */
/*<       RNBR(N)=NGELS >*/
    rnbr[n] = *ngels;
/*   WRITE NEW RIGHT LINE */
/*<    >*/
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
/*   RESET RIGHT NAME IN CONTIG */
/*   ITS LINE NUMBER IS LINCON */
/*<       RNBR(LINCON)=NGELS >*/
    rnbr[*lincon] = *ngels;
/*   HAVE WE INCREASED LENGTH OF CONTIG? */
/*   NEED TO UPDATE LENGTH OF GEL IN CASE OF EDITS */
/*<       IDIM=ABS(LNGTHG(NGELS)) >*/
    *idim = (i__1 = lngthg[*ngels], abs(i__1));
/*<       Y=X+IDIM-1 >*/
    y = *x + *idim - 1;
/*<       RELPG(LINCON)=MAX(RELPG(LINCON),Y) >*/
/* Computing MAX */
    i__1 = relpg[*lincon];
    relpg[*lincon] = max(i__1,y);
/*   WRITE HERE */
/*   WRITE CONTIG DESCRIPTOR */
/*<    >*/
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
/*<       CALL WRITRN(IDEVR,NGELS,NCONTS) >*/
    writrn_(idevr, ngels, nconts);
/*<       RETURN >*/
    return 0;

/*< 400   CONTINUE >*/
L400:

/*   ADDING TO LEFT END */
/*< 410   CONTINUE >*/
/* L410: */
/*   NEED TO EDIT DB HERE */
/*<    >*/
    if (*itotpc > 0) {
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, &c__1, seqc2 + 1, itotpc, idc, idbsiz, idevr, 
		maxgel, (ftnlen)1, (ftnlen)1);
    }

/*< 420   CONTINUE >*/
/* L420: */
/*   SET POINTERS IN NEW GEL */
/*<       RELPG(NGELS)=1 >*/
    relpg[*ngels] = 1;
/*<       RNBR(NGELS)=NCONTC >*/
    rnbr[*ngels] = *ncontc;
/*<       LNBR(NGELS)=0 >*/
    lnbr[*ngels] = 0;
/*   WRITE NEW GEL LINE */
/*<    >*/
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   SET POINTERS IN OLD LEFT END */
/*<       LNBR(NCONTC)=NGELS >*/
    lnbr[*ncontc] = *ngels;
/*<       RELPG(NCONTC)=X >*/
    relpg[*ncontc] = *x;
/*   WRITE NEW LEFT END */
/*<    >*/
    writeg_(idevr, ncontc, &relpg[*ncontc], &lngthg[*ncontc], &lnbr[*ncontc], 
	    &rnbr[*ncontc]);
/*   NEW LENGTH OF CONTIG */
/*<       RELPG(LINCON)=RELPG(LINCON)+X-1 >*/
    relpg[*lincon] = relpg[*lincon] + *x - 1;
/*   MAY HAVE JUST ADDED A GEL LONGER THAN CONTIG */
/*<       IDIM=ABS(LNGTHG(NGELS)) >*/
    *idim = (i__1 = lngthg[*ngels], abs(i__1));
/*<       Y=IDIM >*/
    y = *idim;
/*<       IF(Y.GT.RELPG(LINCON))RELPG(LINCON)=Y >*/
    if (y > relpg[*lincon]) {
	relpg[*lincon] = y;
    }
/*   NEW NAME OF LEFT END OF CONTIG */
/*<       LNBR(LINCON)=NGELS >*/
    lnbr[*lincon] = *ngels;
/*   WRITE CONTIG DESCRIPTOR */
/*<    >*/
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
/*<       CALL WRITRN(IDEVR,NGELS,NCONTS) >*/
    writrn_(idevr, ngels, nconts);
/*   NOW GO THRU AND CHANGE ALL RELATIVE POSITIONS */
/*<       N=NCONTC >*/
    n = *ncontc;
/*< 440   CONTINUE >*/
L440:
/*<       IF(RNBR(N).EQ.0)RETURN >*/
    if (rnbr[n] == 0) {
	return 0;
    }
/*<       N=RNBR(N) >*/
    n = rnbr[n];
/*<       RELPG(N)=RELPG(N)+X-1 >*/
    relpg[n] = relpg[n] + *x - 1;
/*   WRITE NEW LINE */
/*<    >*/
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
/*<       GO TO 440 >*/
    goto L440;
/*<       END >*/
} /* aenter_ */

/*      ABEDIN */

/*   ROUTINE TO EDIT THE DB USING A PADDED SEQ */
/*   HAVE AN ARRAY SEQC2 LENGTH IDC OF PADDED SECTION OF CONTIG LINCON */
/*  THE LEFT END OF THE PADDED CONTIG STARTS AT X */
/*   THERE ARE ITOTPC PADS TO MAKE */

/*<    >*/
/* Subroutine */ int abedin_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *ngels, integer *nconts, char *gel, integer *
	lincon, integer *x, char *seqc2, integer *itotpc, integer *idc, 
	integer *idbsiz, integer *idevr, integer *maxgel, ftnlen gel_len, 
	ftnlen seqc2_len)
{
    /* Initialized data */

    static char p[1] = ",";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, iat, ipad, posn, idone;
    extern /* Subroutine */ int padcon_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, ftnlen), erromf_(char 
	    *, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER  RELPG(IDBSIZ),X,POSN >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       CHARACTER SEQC2(IDC),GEL(MAXGEL),P >*/
/*<       SAVE P >*/
/*<       DATA P/','/ >*/
    /* Parameter adjustments */
    --seqc2;
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */

/*   POINT TO CONTIG */
/*<       POSN=X-1 >*/
    posn = *x - 1;
/*   POINT TO SEQC2 */
/*<       IAT=0 >*/
    iat = 0;
/*   COUNT PADS DONE */
/*<       IDONE=0 >*/
    idone = 0;
/*   LOOP FOR ALL SEQC2 */
/*<       DO 100 J=1,IDC >*/
    i__1 = *idc;
    for (j = 1; j <= i__1; ++j) {
/*<       POSN=POSN+1 >*/
	++posn;
/*<       IAT=IAT+1 >*/
	++iat;
/*<       IPAD=0 >*/
	ipad = 0;
/*   IS THIS A PADDING CHAR? */
/*<       IF(SEQC2(IAT).NE.P)GO TO 100 >*/
	if (*(unsigned char *)&seqc2[iat] != *(unsigned char *)&p[0]) {
	    goto L100;
	}
/*< 50    CONTINUE >*/
L50:
/*   COUNT PADS */
/*<       IPAD=IPAD+1 >*/
	++ipad;
/*<       IAT=IAT+1 >*/
	++iat;
/*<       IF(SEQC2(IAT).EQ.P)GO TO 50 >*/
	if (*(unsigned char *)&seqc2[iat] == *(unsigned char *)&p[0]) {
	    goto L50;
	}
/*   END OF THIS STRETCH OF PADS,DO INSERT */
/*   HAVE IPAD INSERTS TO MAKE AT POSN */
/*      WRITE(*,*)'LINCON,POSN,IPAD',LINCON,POSN,IPAD */
/*<    >*/
	padcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, &posn, &ipad, idbsiz, idevr, maxgel, (ftnlen)1);
/*   MOVE POINTER TO CONTIG */
/*<       POSN=POSN+IPAD >*/
	posn += ipad;
/*   COUNT PADS DONE */
/*<       IDONE=IDONE+IPAD >*/
	idone += ipad;
/*   ANY MORE TO DO? */
/*<       IF(IDONE.EQ.ITOTPC)GO TO 101 >*/
	if (idone == *itotpc) {
	    goto L101;
	}
/*< 100   CONTINUE >*/
L100:
	;
    }
/*   ERROR SHOULD HAVE DONE ALL PADS */
/*<       CALL ERROMF('Problem: some pads were not done!') >*/
    erromf_("Problem: some pads were not done!", (ftnlen)33);
/*< 101   CONTINUE >*/
L101:
/*<       END >*/
    return 0;
} /* abedin_ */

/*<       SUBROUTINE ADDTIT(SEQ1,NAMPRO,NGELS,IDIM1) >*/
/* Subroutine */ int addtit_(char *seq1, char *nampro, integer *ngels, 
	integer *idim1, ftnlen seq1_len, ftnlen nampro_len)
{
    extern /* Subroutine */ int cadtit_(char *, char *, integer *, ftnlen, 
	    ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQ1(20),NAMPRO*(*) >*/
/*<       CALL CADTIT(SEQ1, NAMPRO, NGELS) >*/
    /* Parameter adjustments */
    --seq1;

    /* Function Body */
    cadtit_(seq1 + 1, nampro, ngels, (ftnlen)1, nampro_len);
/*<       IDIM1=IDIM1+20 >*/
    *idim1 += 20;
/*<       END >*/
    return 0;
} /* addtit_ */

/*<    >*/
/* Subroutine */ int adism3_(integer *isavps, integer *savpg, integer *cends, 
	integer *nends, integer *idcend, integer *maxcon, integer *ilefts, 
	integer *ilc, integer *iposc, integer *iposg, integer *isense, 
	integer *llino, integer *imatc, integer *istran, integer *nextc, 
	integer *maxc, integer *jj, integer *isavl, integer *lmatch)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, lcl, lcr, savps;
    extern /* Subroutine */ int erromf_(char *, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER ILEFTS(MAXC),ILC(MAXC),IPOSC(MAXC),IPOSG(MAXC) >*/
/*<       INTEGER ISENSE(MAXC),LLINO(MAXC) >*/
/*<       INTEGER SAVPS,SAVPG,CENDS(MAXCON) >*/
/*<       INTEGER NENDS(MAXCON) >*/
/*<       SAVPS=ISAVPS-19 >*/
    /* Parameter adjustments */
    --nends;
    --cends;
    --llino;
    --isense;
    --iposg;
    --iposc;
    --ilc;
    --ilefts;

    /* Function Body */
    savps = *isavps - 19;
/*<       JJ=1 >*/
    *jj = 1;

/* we have a match isavps, isavpg, isavl (pos in consensus, gel, length) */
/* 1. which contig is it in? --> JJ */
/* 2. save pos of contig, pos in contig, contig length, contig number, sense */

/*      WRITE(*,*)'ENTER ADISM3, IMATC',IMATC */
/*<       DO 5 J=2,IDCEND >*/
    i__1 = *idcend;
    for (j = 2; j <= i__1; ++j) {
/*<         IF(SAVPS.GT.CENDS(J))GO TO 5 >*/
	if (savps > cends[j]) {
	    goto L5;
	}
/*<         JJ=J-1 >*/
	*jj = j - 1;
/*<         GO TO 6 >*/
	goto L6;
/*< 5     CONTINUE >*/
L5:
	;
    }
/*<       JJ=IDCEND >*/
    *jj = *idcend;
/*< 6     CONTINUE >*/
L6:
/*<       SAVPS=SAVPS-1 >*/
    --savps;
/*<       LCL=SAVPS-CENDS(JJ) >*/
    lcl = savps - cends[*jj];
/*<       LCR=CENDS(JJ+1)-ISAVPS-1 >*/
    lcr = cends[*jj + 1] - *isavps - 1;
/*<       NEXTC=CENDS(JJ+1)+20 >*/
    *nextc = cends[*jj + 1] + 20;
/*<       IF(IMATC.LE.MAXC) THEN >*/
    if (*imatc <= *maxc) {
/*<         ILEFTS(IMATC)=CENDS(JJ)+20 >*/
	ilefts[*imatc] = cends[*jj] + 20;
/*<         ILC(IMATC)=LCL+LCR+1 >*/
	ilc[*imatc] = lcl + lcr + 1;
/*<         IPOSC(IMATC)=LCL+1 >*/
	iposc[*imatc] = lcl + 1;
/*<         IPOSG(IMATC)=SAVPG >*/
	iposg[*imatc] = *savpg;
/*<         LLINO(IMATC)=NENDS(JJ) >*/
	llino[*imatc] = nends[*jj];
/*<         ISENSE(IMATC)=1 >*/
	isense[*imatc] = 1;
/*<         IF(ISTRAN.EQ.2)ISENSE(IMATC)=-1 >*/
	if (*istran == 2) {
	    isense[*imatc] = -1;
	}
/*<         LMATCH = ISAVL >*/
	*lmatch = *isavl;
/*        WRITE(INFOD,1000)LLINO(IMATC),IPOSC(IMATC),ISTRAN, */
/*     +  IPOSG(IMATC) */
/* 1000   FORMAT */
/*     +  ('Contig',I8,' position',I8,' matches strand',I2, */
/*     +  ' at position',I8) */
/*        CALL INFO(INFOD) */
/*<       ELSE >*/
    } else {
/*<         CALL ERROMF('Warning: too many overlaps') >*/
	erromf_("Warning: too many overlaps", (ftnlen)26);
/*<       END IF >*/
    }
/*<       END >*/
    return 0;
} /* adism3_ */

/*<    >*/
/* Subroutine */ int adism4_(integer *idim, integer *idimg, integer *savps, 
	integer *savpg, integer *savl, integer *idsav, integer *cends, 
	integer *nends, integer *idcend, integer *maxcon, integer *ilefts, 
	integer *ilc, integer *iposc, integer *iposg, integer *isense, 
	integer *llino, integer *imatc, integer *istran, integer *maxc)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, lend, lastc, nextc;
    extern /* Subroutine */ int bub3as_(integer *, integer *, integer *, 
	    integer *), adism3_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer lmatch;

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER ILEFTS(MAXC),ILC(MAXC),IPOSC(MAXC),IPOSG(MAXC) >*/
/*<       INTEGER ISENSE(MAXC),LLINO(MAXC) >*/
/*<       INTEGER CENDS(MAXCON) >*/
/*<       INTEGER NENDS(MAXCON) >*/
/*<       INTEGER SAVPS(IDSAV),SAVPG(IDSAV),SAVL(IDSAV) >*/
/*      CHARACTER INFOD*80 */
/*      WRITE(*,*)'ENTER ADISM4    , IMATC',IMATC */
/*<       NEXTC=IDIM+1 >*/
    /* Parameter adjustments */
    --savl;
    --savpg;
    --savps;
    --nends;
    --cends;
    --llino;
    --isense;
    --iposg;
    --iposc;
    --ilc;
    --ilefts;

    /* Function Body */
    nextc = *idim + 1;

/* sort on position in consensus */

/*<       CALL BUB3AS(SAVPS,SAVPG,SAVL,IDSAV) >*/
    bub3as_(&savps[1], &savpg[1], &savl[1], idsav);
/*      DO 123 II = 1,IDSAV */
/*        WRITE(*,*)II,SAVPS(II),SAVPG(II),SAVL(II) */
/* 123    CONTINUE */
/*<         IMATC=IMATC+1 >*/
    ++(*imatc);
/*      WRITE(*,*)'IN ADISM4, UPDATED IMATC',IMATC */

/* get the contig info for the first match */

/*        WRITE(*,*)'sav1',SAVL(1) */
/* we have a match savps, savpg, savl (pos in consensus, gel, length) */
/* 2. save pos of contig, pos in contig, contig length, contig number, sense */
/*<    >*/
    adism3_(&savps[1], &savpg[1], &cends[1], &nends[1], idcend, maxcon, &
	    ilefts[1], &ilc[1], &iposc[1], &iposg[1], &isense[1], &llino[1], 
	    imatc, istran, &nextc, maxc, &lastc, &savl[1], &lmatch);

/* now decide when a match is with a new contig and get the relevant info. */
/* Decide its the same overlap if it is covered by the previous gel position */
/* If we want to record the longest match for each overlap we should test */
/* here to see if overlapping ones are longer than the one weve recorded */

/*<       LEND=IDIMG-SAVPG(1)+SAVPS(1) >*/
    lend = *idimg - savpg[1] + savps[1];
/*<       DO 10 I=2,IDSAV >*/
    i__1 = *idsav;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*         WRITE(*,*)SAVPS(I),SAVPG(I) */
/*         WRITE(*,*)'SAVPS(I)-SAVPG(I)',SAVPS(I)-SAVPG(I) */
/*        WRITE(*,*)'savl(I),lend',SAVL(I),LEND */
/*<         IF((SAVPS(I).LT.LEND).AND.(SAVPS(I).LT.NEXTC)) THEN >*/
	if (savps[i__] < lend && savps[i__] < nextc) {

/* test here if this match is longer */

/*<           IF(SAVL(I).GT.LMATCH) THEN >*/
	    if (savl[i__] > lmatch) {

/* next test added 22-11-94 because the trap in adism3 is insufficient */

/*<             IF(IMATC.LE.MAXC) THEN >*/
		if (*imatc <= *maxc) {
/*<                IPOSC(IMATC) = SAVPS(I) - CENDS(LASTC) + 1 - 20 >*/
		    iposc[*imatc] = savps[i__] - cends[lastc] - 19;
/*<                IPOSG(IMATC) = SAVPG(I) >*/
		    iposg[*imatc] = savpg[i__];
/*<                LMATCH = SAVL(I) >*/
		    lmatch = savl[i__];
/*              WRITE(*,*)'new best g,c,l',IPOSG(IMATC),IPOSC(IMATC),LMATCH */
/*<             END IF >*/
		}
/*<           END IF >*/
	    }
/*<           GO TO 10 >*/
	    goto L10;
/*<         END IF >*/
	}
/*         WRITE(*,*)'2SAVPS(I)-SAVPG(I)',SAVPS(I)-SAVPG(I) */
/*         WRITE(*,*)IPOSC(IMATC),IPOSG(IMATC),SAVPS(I),SAVPG(I) */
/*<         IMATC=IMATC+1 >*/
	++(*imatc);
/*      WRITE(*,*)'IN ADISM4, UPDATED AGAIN IMATC',IMATC */
/* we have a match savps, savpg, savl (pos in consensus, gel, length) */
/* 2. save pos of contig, pos in contig, contig length, contig number, sense */
/*<    >*/
	adism3_(&savps[i__], &savpg[i__], &cends[1], &nends[1], idcend, 
		maxcon, &ilefts[1], &ilc[1], &iposc[1], &iposg[1], &isense[1],
		 &llino[1], imatc, istran, &nextc, maxc, &lastc, &savl[i__], &
		lmatch);
/*<         LEND=IDIMG-SAVPG(I)+SAVPS(I) >*/
	lend = *idimg - savpg[i__] + savps[i__];
/*        RSTART = SAVPS(I) - SAVPG(I) */
/*< 10    CONTINUE >*/
L10:
	;
    }
/*<       IMATC = MIN(IMATC,MAXC) >*/
    *imatc = min(*imatc,*maxc);
/*      WRITE(*,*)'IN ADISM4, LAST IMATC',IMATC */
/*<       END >*/
    return 0;
} /* adism4_ */

/* C    AJOIN2 */
/* C   COMPLETES JOIN AND RETURNS LENGTH OF NEW CONTIG IN LLINOR */
/*<    >*/
/* Subroutine */ int ajoin2_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *ngels, integer *nconts, integer *idbsiz, 
	integer *relx, integer *llinol, integer *llinor, integer *lnconl, 
	integer *lnconr, integer *idevr)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer n;
    static real x;
    extern /* Subroutine */ int merge_(integer *, integer *, integer *, 
	    integer *, integer *, integer *), remcnl_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), mrgtag_(integer *, integer *, integer *, integer *), 
	    writec_(integer *, integer *, integer *, integer *, integer *), 
	    writeg_(integer *, integer *, integer *, integer *, integer *, 
	    integer *), mrgnot_(integer *, integer *, integer *);

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER RELPG(IDBSIZ) >*/
/*<       INTEGER LNBR(IDBSIZ),RNBR(IDBSIZ),LNGTHG(IDBSIZ) >*/
/*<       INTEGER RELX >*/
/*   RELX IS THE POSITION OF THE JOINT */
/*   LLINOL IS THE LEFT GEL NUMBER OF THE LEFT CONTIG */
/*   LLINOR IS THE LEFT GEL OF THE RIGHT CONTIG */
/*   LNCONL IS THE LEFT CONTIG LINE NUMBER */
/*   LNCONR IS THE RIGHT CONTIG LINE NUMBER */

/*   ADJUST ALL RELATIVE POSITIONS IN RIGHT CONTIG */
/*<       N=LLINOR >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    n = *llinor;
/*<       RELPG(N)=RELX >*/
    relpg[n] = *relx;
/*< 50    CONTINUE >*/
L50:
/*<       IF(RNBR(N).EQ.0)GO TO 60 >*/
    if (rnbr[n] == 0) {
	goto L60;
    }
/*<       N=RNBR(N) >*/
    n = rnbr[n];
/*<       RELPG(N)=RELPG(N)+RELX-1 >*/
    relpg[n] = relpg[n] + *relx - 1;
/*<       GO TO 50 >*/
    goto L50;
/*< 60    CONTINUE >*/
L60:

/*   FIX UP NEW GEL LINE FOR OLD LEFT OF RIGHT CONTIG */
/*<       LNBR(LLINOR)=RNBR(LNCONL) >*/
    lnbr[*llinor] = rnbr[*lnconl];
/*   FIX UP RIGHT GEL OF LEFT CONTIG */
/*<       N=RNBR(LNCONL) >*/
    n = rnbr[*lnconl];
/*<       RNBR(N)=LLINOR >*/
    rnbr[n] = *llinor;
/*   MERGE WILL SORT OUT THE CORRECT NEIGHBOURS */

/*<       CALL MERGE(RELPG,LNGTHG,LNBR,RNBR,LNCONL,IDBSIZ) >*/
    merge_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], lnconl, idbsiz);
/*   MERGE DOES NOT WRITE TO DISK */
/*<       N=LNBR(LNCONL) >*/
    n = lnbr[*lnconl];
/*< 65    CONTINUE >*/
L65:
/*      WRITE(IDEVR,REC=N)RELPG(N),LNGTHG(N),LNBR(N),RNBR(N) */
/*<       CALL WRITEG(IDEVR,N,RELPG(N),LNGTHG(N),LNBR(N),RNBR(N)) >*/
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
/*<       N=RNBR(N) >*/
    n = rnbr[n];
/*<       IF(N.NE.0)GO TO 65 >*/
    if (n != 0) {
	goto L65;
    }
/*   Merge annotation lists. */
/*<       CALL MRGTAG(IDEVR, IDBSIZ-LNCONL, IDBSIZ-LNCONR, RELX-1) >*/
    i__1 = *idbsiz - *lnconl;
    i__2 = *idbsiz - *lnconr;
    i__3 = *relx - 1;
    mrgtag_(idevr, &i__1, &i__2, &i__3);
/*<       CALL MRGNOT(IDEVR, IDBSIZ-LNCONR, IDBSIZ-LNCONL) >*/
    i__1 = *idbsiz - *lnconr;
    i__2 = *idbsiz - *lnconl;
    mrgnot_(idevr, &i__1, &i__2);
/*   CONTIG LINES */
/*<       X=RELPG(LNCONR)+RELX-1 >*/
    x = (real) (relpg[*lnconr] + *relx - 1);
/*   LENGTH MAY NOT HAVE INCREASED! */
/*<       IF(X.GT.RELPG(LNCONL))RELPG(LNCONL)=X >*/
    if (x > (real) relpg[*lnconl]) {
	relpg[*lnconl] = x;
    }
/*   SAVE LENGTH OF NEW CONTIG */
/*<       RELX=RELPG(LNCONL) >*/
    *relx = relpg[*lnconl];
/*      WRITE(IDEVR,REC=LNCONL)RELPG(LNCONL),LNGTHG(LNCONL),LNBR(LNCONL), */
/*     1RNBR(LNCONL) */
/*<    >*/
    i__1 = *idbsiz - *lnconl;
    writec_(idevr, &i__1, &relpg[*lnconl], &lnbr[*lnconl], &rnbr[*lnconl]);
/*   Now remove the old contig. We must use the C routine for this so that */
/*   it can update the contig order, tag lists, etc. */
/*<    >*/
    remcnl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    lnconr, idevr);
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* ajoin2_ */

/*     SUBROUTINE AJOIN3 */
/*<    >*/
/* Subroutine */ int ajoin3_(integer *relpg, integer *idbsiz, integer *lincon,
	 integer *itype, integer *isense, integer *joint, integer *idim22, 
	integer *klass, integer *iover, integer *pl, integer *pr)
{
    static integer i__;

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER RELPG(IDBSIZ),LINCON(2),IDIM22(2) >*/
/*<       INTEGER ITYPE(2),ISENSE(2),JOINT(2),PL(2),PR(2) >*/

/*   CALC POSITIONS OF CONTIGS RELATIVE TO FIXED GEL */
/*<       DO 20 I=1,2 >*/
    /* Parameter adjustments */
    --relpg;
    --lincon;
    --itype;
    --isense;
    --joint;
    --idim22;
    --pl;
    --pr;

    /* Function Body */
    for (i__ = 1; i__ <= 2; ++i__) {
/*   R+ */
/*<       IF((ITYPE(I).NE.-1).OR.(ISENSE(I).NE.1))GO TO 11 >*/
	if (itype[i__] != -1 || isense[i__] != 1) {
	    goto L11;
	}
/*<       PL(I)=-1*JOINT(I)+2 >*/
	pl[i__] = -joint[i__] + 2;
/*<       PR(I)=PL(I)+RELPG(LINCON(I))-1 >*/
	pr[i__] = pl[i__] + relpg[lincon[i__]] - 1;
/*<       GO TO 20 >*/
	goto L20;
/*   L+ */
/*< 11    CONTINUE >*/
L11:
/*<       IF((ITYPE(I).NE.1).OR.(ISENSE(I).NE.1))GO TO 12 >*/
	if (itype[i__] != 1 || isense[i__] != 1) {
	    goto L12;
	}
/*<       PL(I)=JOINT(I) >*/
	pl[i__] = joint[i__];
/*<       PR(I)=PL(I)+RELPG(LINCON(I))-1 >*/
	pr[i__] = pl[i__] + relpg[lincon[i__]] - 1;
/*<       GO TO 20 >*/
	goto L20;
/*   R- */
/*< 12    CONTINUE >*/
L12:
/*<       IF((ITYPE(I).NE.-1).OR.(ISENSE(I).NE.-1))GO TO 13 >*/
	if (itype[i__] != -1 || isense[i__] != -1) {
	    goto L13;
	}
/*<       PR(I)=JOINT(I)+IDIM22(I)-1 >*/
	pr[i__] = joint[i__] + idim22[i__] - 1;
/*<       PL(I)=PR(I)-RELPG(LINCON(I))+1 >*/
	pl[i__] = pr[i__] - relpg[lincon[i__]] + 1;
/*<       GO TO 20 >*/
	goto L20;
/*   L- */
/*< 13    CONTINUE >*/
L13:
/*<       PR(I)=IDIM22(I)-JOINT(I)+1 >*/
	pr[i__] = idim22[i__] - joint[i__] + 1;
/*<       PL(I)=PR(I)-RELPG(LINCON(I))+1 >*/
	pl[i__] = pr[i__] - relpg[lincon[i__]] + 1;
/*< 20    CONTINUE >*/
L20:
	;
    }
/*  LENGTH OF OVERLAP */
/*<       IOVER=MIN(PR(1),PR(2))-MAX(PL(1),PL(2))+1 >*/
    *iover = min(pr[1],pr[2]) - max(pl[1],pl[2]) + 1;

/*  CLASS NUMBER 1-16 */
/*<       KLASS=1 >*/
    *klass = 1;
/*<       IF(ITYPE(1).EQ.1)KLASS=KLASS+8 >*/
    if (itype[1] == 1) {
	*klass += 8;
    }
/*<       IF(ISENSE(1).EQ.-1)KLASS=KLASS+4 >*/
    if (isense[1] == -1) {
	*klass += 4;
    }
/*<       IF(ITYPE(2).EQ.1)KLASS=KLASS+2 >*/
    if (itype[2] == 1) {
	*klass += 2;
    }
/*<       IF(ISENSE(2).EQ.-1)KLASS=KLASS+1 >*/
    if (isense[2] == -1) {
	++(*klass);
    }
/*<       END >*/
    return 0;
} /* ajoin3_ */

/*      ALINE */

/*    ROUTINE TO LINE UP 2 SEQS. */
/*   IT SLIDES,REMOVES OVERLAPPING MATCHES, */
/*   SORTS MATCHES INTO ASCENDING ORDER, THEN DOES DOES A TOPOLOGICAL */
/*   CHECK, AND THEN PRODUCES 2 LINED UP SEQS WITH PADDING CHARS */
/*   VARIABLES */
/*       SEQ1 CONSENSUS */
/*       SEQ2 GEL ORIGINAL IN CORRECT ORIENTATION */
/*       SEQG2 ALIGNED GEL */
/*       SEQC2 ALIGNED CONSENSUS */
/*       SEQ3 SAVED GEL RAW DATA */
/*       ISAV1,2,3 STORE MATCHES AND POSITIONS */
/*       IDSAV NUMBER ISAV'S */
/*       IDC LENGTH OF INPUT SEQ1 */
/*       IDIM2 LENGTH OF INPUT SEQ2 */
/*       IDOUT LENGTH OF OUTPUT ALIGNED SEQ1 */
/*       IDIM2 LENGTH OF SEQ2 ON OUTPUT AFTER ALIGNMENT */
/*       MINSLI MIN MATCH FOR SLIDING */
/*       IFAIL FLAG TO SHOW IF ALIGNMENT FAILED DUE TO TOO */
/*   MANY MISMATCHES OR TOPOLIGICAL CHECK OR TOO MANY OR TOO MANY */
/*   PADDING CHARS. 1=FAIL,0=PASS */

/*<    >*/
/* Subroutine */ int aline_(char *seq1, char *seq2, char *seqg2, char *seqc2, 
	integer *isav1, integer *isav2, integer *isav3, integer *idsav, 
	integer *idc, integer *idim2, integer *idout, integer *ic1, integer *
	ig1, integer *minsli, integer *joint, integer *itotpc, integer *
	itotpg, integer *ifail, integer *itype, integer *maxpc, integer *
	maxpg, real *permax, char *seq3, integer *maxgel, real *percm, 
	integer *leno, integer *ishow, integer *mask, integer *jrorc, ftnlen 
	seq1_len, ftnlen seq2_len, ftnlen seqg2_len, ftnlen seqc2_len, ftnlen 
	seq3_len)
{
    static integer ipp;
    extern /* Subroutine */ int maskc_(char *, integer *, integer *, ftnlen);
    static integer idim2i;
    extern /* Subroutine */ int bub3as_(integer *, integer *, integer *, 
	    integer *), dalign_(char *, char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, real *, integer *, integer *, 
	    real *, integer *, integer *, integer *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), tpchek_(integer *, integer *, integer *, 
	    integer *), upchek_(integer *, integer *, integer *, integer *), 
	    slides_(char *, integer *, char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, ftnlen, ftnlen), lineup_(char *, char *, char *, 
	    char *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), removl_(
	    integer *, integer *, integer *, integer *);
    static integer minslt;
    extern /* Subroutine */ int mstlkl_(char *, integer *, ftnlen), sqcopy_(
	    char *, char *, integer *, ftnlen, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQ1(IDC),SEQ2(IDIM2),SEQG2(IDOUT),SEQC2(IDOUT) >*/
/*<       CHARACTER SEQ3(MAXGEL) >*/
/*<       INTEGER ISAV1(IDSAV),ISAV2(IDSAV),ISAV3(IDSAV) >*/
/*<       MINSLT=MINSLI >*/
    /* Parameter adjustments */
    --isav3;
    --isav2;
    --isav1;
    --seq1;
    --seq2;
    --seqc2;
    --seqg2;
    --seq3;

    /* Function Body */
    minslt = *minsli;
/*<       IDIM2I = IDIM2 >*/
    idim2i = *idim2;

/* need to unmask both(for contig joins) sequences */

/*        CALL FMTDB(SEQ2,IDIM2,1,IDIM2,60,6) */
/*<       IF (MASK.NE.0) THEN >*/
    if (*mask != 0) {
/*        WRITE(*,*)'SEQ1 B' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
/*<         CALL MASKC(SEQ1,IDC,2) >*/
	maskc_(seq1 + 1, idc, &c__2, (ftnlen)1);
/*        WRITE(*,*)'SEQ1 A' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
/*        WRITE(*,*)'SEQ2 B' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2) */
/*<         CALL MASKC(SEQ2,IDIM2,2) >*/
	maskc_(seq2 + 1, idim2, &c__2, (ftnlen)1);
/*        WRITE(*,*)'SEQ2 A' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2) */
/*<       END IF >*/
    }
/*        CALL FMTDB(SEQ2,IDIM2,1,IDIM2,60,6) */
/*   SAVE SEQ2 */
/*<       CALL SQCOPY(SEQ2,SEQ3,IDIM2) >*/
    sqcopy_(seq2 + 1, seq3 + 1, idim2, (ftnlen)1, (ftnlen)1);
/*<       CALL MSTLKL(SEQ3,IDIM2) >*/
    mstlkl_(seq3 + 1, idim2, (ftnlen)1);
/*        CALL FMTDB(SEQ3,IDIM2,1,IDIM2,60,6) */
/*<       IFAIL=1 >*/
    *ifail = 1;
/*   FIND MATCHES */
/*<       IPP=IDSAV >*/
    ipp = *idsav;
/*      WRITE(*,*)'IC1,IG1',IC1,IG1,MAXPG,MAXPC,MINSLT */
/*<    >*/
    slides_(seq1 + 1, idc, seq3 + 1, idim2, ic1, ig1, maxpg, maxpc, &minslt, &
	    isav1[1], &isav2[1], &isav3[1], &ipp, (ftnlen)1, (ftnlen)1);
/*      WRITE(*,*)'IPP',IPP,IDSAV */
/*<       IF(IPP.GT.IDSAV) GO TO 50 >*/
    if (ipp > *idsav) {
	goto L50;
    }
/*<       IF(IPP.LT.1) GO TO 50 >*/
    if (ipp < 1) {
	goto L50;
    }
/*<       CALL REMOVL(ISAV2,ISAV3,ISAV1,IPP) >*/
    removl_(&isav2[1], &isav3[1], &isav1[1], &ipp);
/*      WRITE(*,*)'IPP',IPP,IDSAV */
/*<       CALL BUB3AS(ISAV2,ISAV3,ISAV1,IPP) >*/
    bub3as_(&isav2[1], &isav3[1], &isav1[1], &ipp);
/*   DO TOPOLOGICAL CHECK */
/*<       CALL TPCHEK(ISAV2,ISAV3,ISAV1,IPP) >*/
    tpchek_(&isav2[1], &isav3[1], &isav1[1], &ipp);

/* added next routine 27-2-93 */

/*      WRITE(*,*)'IPP',IPP,IDSAV */
/*<       CALL UPCHEK(ISAV2,ISAV3,ISAV1,IPP) >*/
    upchek_(&isav2[1], &isav3[1], &isav1[1], &ipp);
/*      WRITE(*,*)'IPP',IPP,IDSAV */
/*<    >*/
    lineup_(seq2 + 1, seq1 + 1, seqg2 + 1, seqc2 + 1, idc, idim2, idout, &
	    isav3[1], &isav2[1], &isav1[1], &ipp, itotpc, itotpg, joint, 
	    itype, maxgel, ifail, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/*       WRITE(*,*)'ITOTPC,ITOTPG',ITOTPC,ITOTPG,IFAIL */
/*      IF(ITOTPC.GT.MAXPC)IFAIL=1 */
/*      IF(ITOTPG.GT.MAXPG)IFAIL=1 */
/*<       IF(IFAIL.NE.0) GO TO 50 >*/
    if (*ifail != 0) {
	goto L50;
    }
/*   IDIM2 IS NOW LENGTH OF ALIGNED GEL */
/*<    >*/
    dalign_(seqc2 + 1, seqg2 + 1, seq3 + 1, maxgel, idout, idim2, joint, 
	    itype, percm, ifail, leno, permax, ishow, maxpg, maxpc, itotpg, 
	    itotpc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/*<  50   CONTINUE >*/
L50:

/* need to remask both(for contig joins) sequences */

/*<       IF (MASK.NE.0) THEN >*/
    if (*mask != 0) {
/*        WRITE(*,*)'SEQ1 B1' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
/*<         CALL MASKC(SEQ1,IDC,3) >*/
	maskc_(seq1 + 1, idc, &c__3, (ftnlen)1);
/*        WRITE(*,*)'SEQ1 A1' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */

/* only mask consensus data which is only lowercase where marked */

/*        WRITE(*,*)'SEQ2 B1' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2I) */
/*<         IF (JRORC.EQ.1) CALL MASKC(SEQ2,IDIM2I,3) >*/
	if (*jrorc == 1) {
	    maskc_(seq2 + 1, &idim2i, &c__3, (ftnlen)1);
	}
/*        WRITE(*,*)'SEQ2 A1' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2I) */
/*<       END IF >*/
    }
/*<       END >*/
    return 0;
} /* aline_ */

/*<    >*/
/* Subroutine */ int autocn_(char *seq1, integer *idim, char *gel, integer *
	idimg, integer *ilefts, integer *ilc, integer *iposc, integer *iposg, 
	integer *isense, integer *llino, integer *imatc, integer *ifcomp, 
	integer *minmat, integer *maxgel, integer *maxglm, char *gelcop, 
	integer *savps, integer *savpg, integer *savl, integer *maxsav, 
	integer *cends, integer *nends, integer *maxcon, char *seqg2, char *
	seqc2, char *seq4, integer *idout, integer *idim22, integer *itotpg, 
	integer *itotpc, integer *joint, integer *ifail, integer *itype, 
	integer *maxpc, integer *maxpg, real *permax, integer *minsli, char *
	seqg3, char *seqc3, integer *kfail, integer *jobc, real *permis, 
	integer *leno, integer *ishow, integer *mask, integer *minovr, ftnlen 
	seq1_len, ftnlen gel_len, ftnlen gelcop_len, ftnlen seqg2_len, ftnlen 
	seqc2_len, ftnlen seq4_len, ftnlen seqg3_len, ftnlen seqc3_len)
{
    /* System generated locals */
    integer seqg2_dim1, seqg2_offset, seqc2_dim1, seqc2_offset, i__1, i__2, 
	    i__3, i__4, i__5;

    /* Local variables */
    static integer i__, jlc[100];
    static char csen[1];
    extern /* Subroutine */ int info_(char *, ftnlen),
	swrt1_(char *, char *, ...), swrt4_(char *, char *, ...);
    static integer jfail, jdim22;
    extern /* Subroutine */ int aline_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, real *, 
	    char *, integer *, real *, integer *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jmatc;
    static char infod[80];
    static integer idsav, jposc[100], jposg[100], ioptc;
    extern /* Subroutine */ int sqcom_(char *, integer *, ftnlen);
    static integer jdout;
    extern /* Subroutine */ int busyf_(void);
    static integer jtype, start;
    static real perms;
    extern /* Subroutine */ int copym_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, char *, char *, char *, char *, real *, real *, ftnlen,
	     ftnlen, ftnlen, ftnlen), sqrev_(char *, integer *, ftnlen), 
	    adism4_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer idcend;
    extern /* Subroutine */ int fndcon_(char *, integer *, integer *, integer 
	    *, integer *, integer *, ftnlen);
    static integer jlefts[100], jsense[100], jllino[100];
    extern integer cmpseq_(integer *, char *, integer *, integer *, integer *,
	     integer *, integer *, char *, char *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer istran, ksense;
    extern /* Subroutine */ int mstlkl_(char *, integer *, ftnlen);
    static integer lenovr, jjoint, jtotpc, jtotpg;
    extern /* Subroutine */ int sqcopy_(char *, char *, integer *, ftnlen, 
	    ftnlen);

/*   AUTHOR: RODGER STADEN */
/*   changed 29-11-90 to make first in list of alignments the best */
/*<       INTEGER ILEFTS(2),ILC(2),IPOSC(2),IPOSG(2),ISENSE(2),LLINO(2) >*/
/*<       INTEGER SAVPS(MAXSAV) >*/
/*<       INTEGER SAVPG(MAXSAV),SAVL(MAXSAV) >*/
/*<       CHARACTER GELCOP(MAXGLM) >*/
/*<       INTEGER CENDS(MAXCON),NENDS(MAXCON) >*/
/*<       CHARACTER SEQ1(IDIM),GEL(MAXGLM) >*/

/*<       CHARACTER SEQG2(MAXGLM,2),SEQC2(MAXGLM,2),SEQ4(MAXGLM) >*/
/*<       INTEGER IDOUT(2),IDIM22(2),ITOTPG(2),ITOTPC(2),JOINT(2) >*/
/*<       INTEGER IFAIL(2),ITYPE(2) >*/
/*<       PARAMETER (MAXC = 100) >*/
/*<       CHARACTER SEQG3(MAXGLM),SEQC3(MAXGLM) >*/
/*<       INTEGER JLEFTS(MAXC),JLC(MAXC),JPOSC(MAXC),JPOSG(MAXC) >*/
/*<       INTEGER JSENSE(MAXC),JLLINO(MAXC),START >*/
/*<       REAL PERMIS(2) >*/
/*<       CHARACTER CSEN >*/
/*<       CHARACTER INFOD*80 >*/
/*<       INTEGER CMPSEQ >*/
/*<       EXTERNAL CMPSEQ >*/
/*<       CSEN = 'f' >*/
    /* Parameter adjustments */
    --seq1;
    --ilefts;
    --ilc;
    --iposc;
    --iposg;
    --isense;
    --llino;
    --seqc3;
    --seqg3;
    --seq4;
    seqc2_dim1 = *maxglm;
    seqc2_offset = 1 + seqc2_dim1;
    seqc2 -= seqc2_offset;
    seqg2_dim1 = *maxglm;
    seqg2_offset = 1 + seqg2_dim1;
    seqg2 -= seqg2_offset;
    --gelcop;
    --gel;
    --savl;
    --savpg;
    --savps;
    --nends;
    --cends;
    --idout;
    --idim22;
    --itotpg;
    --itotpc;
    --joint;
    --ifail;
    --itype;
    --permis;

    /* Function Body */
    *(unsigned char *)csen = 'f';

/* jobc tells how to update the hash tables: */
/* 0 means dont do anything because the consensus hasnt changed */
/* 1 means add the last contig because a new one has been stuck on the end */
/* 2 means do the whole consensus */

/*<       IFAIL(1) = 1 >*/
    ifail[1] = 1;
/*<       IFAIL(2) = 1 >*/
    ifail[2] = 1;
/*<       KFAIL = 0 >*/
    *kfail = 0;
/*  23-8-90 Need to deal with failures in a better way. Problem is */
/*          case where overlaps are found but fail to align. In future */
/*          signal them with new variable KFAIL which will be nonzero */
/*          if any alignment fails. */
/*  29-11-90 Changed sorting of overlaps so that the best is first in the */
/*           list returned to caller. */
/*      WRITE(*,*)'MINMAT,ITOTPG,ITOTPC',MINMAT,ITOTPG,ITOTPC */
/*      WRITE(*,*)'ISHOW,MASK,MINSLI',ISHOW,MASK,MINSLI */
/*      WRITE(*,*)'MAXPG,MAXPC,PERMAX',MAXPG,MAXPC,PERMAX */
/*   SAVE GEL */
/*<       CALL SQCOPY(GEL,GELCOP,IDIMG) >*/
    sqcopy_(gel + 1, gelcop + 1, idimg, (ftnlen)1, (ftnlen)1);
/*  COUNT NUMBER OF CONTIGS THAT MATCH */
/*<       IMATC=0 >*/
    *imatc = 0;
/*<       IDCEND=MAXCON >*/
    idcend = *maxcon;
/*<       CALL BUSYF() >*/
    busyf_();
/*      WRITE(*,*)'IDIM',IDIM,IDCEND */
/*<       CALL FNDCON(SEQ1,IDIM,CENDS,NENDS,IDCEND,MAXCON) >*/
    fndcon_(seq1 + 1, idim, &cends[1], &nends[1], &idcend, maxcon, (ftnlen)1);
/*<       IF (JOBC.NE.0) THEN >*/
    if (*jobc != 0) {
/*<         START = 1 >*/
	start = 1;
/*<         IF(JOBC.EQ.1) START = CENDS(IDCEND) >*/
	if (*jobc == 1) {
	    start = cends[idcend];
	}
/*<         IOPTC = 2 >*/
	ioptc = 2;
/*<    >*/
	*ifcomp = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1]
		, &idsav, seq1 + 1, gel + 1, idim, idimg, (ftnlen)1, (ftnlen)
		1, (ftnlen)1);
/*<         IF (IFCOMP.NE.0) RETURN >*/
	if (*ifcomp != 0) {
	    return 0;
	}
/*<       END IF >*/
    }
/*< 1     CONTINUE >*/
/* L1: */
/*<       ISTRAN=1 >*/
    istran = 1;
/*< 2     CONTINUE >*/
L2:
/*<       CALL MSTLKL(GEL,IDIMG) >*/
    mstlkl_(gel + 1, idimg, (ftnlen)1);
/*<       IDSAV=MAXSAV >*/
    idsav = *maxsav;
/*<       IOPTC = 3 >*/
    ioptc = 3;
/*<    >*/
    *ifcomp = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1], &
	    idsav, seq1 + 1, gel + 1, idim, idimg, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1);
/*<       IF (IFCOMP.LT.0) RETURN >*/
    if (*ifcomp < 0) {
	return 0;
    }
/*<       IDSAV = IFCOMP >*/
    idsav = *ifcomp;
/*<       IF(IDSAV.NE.0)THEN >*/
    if (idsav != 0) {
/*<    >*/
	adism4_(idim, idimg, &savps[1], &savpg[1], &savl[1], &idsav, &cends[1]
		, &nends[1], &idcend, maxcon, jlefts, jlc, jposc, jposg, 
		jsense, jllino, imatc, &istran, &c__100);
/*<       END IF >*/
    }
/*<       ISTRAN=ISTRAN+1 >*/
    ++istran;
/*<       IF(ISTRAN.EQ.2) THEN >*/
    if (istran == 2) {
/*<         CALL SQCOPY(GELCOP,GEL,IDIMG) >*/
	sqcopy_(gelcop + 1, gel + 1, idimg, (ftnlen)1, (ftnlen)1);
/*<         CALL SQREV(GEL,IDIMG) >*/
	sqrev_(gel + 1, idimg, (ftnlen)1);
/*<         CALL SQCOM(GEL,IDIMG) >*/
	sqcom_(gel + 1, idimg, (ftnlen)1);
/*<         GO TO 2 >*/
	goto L2;
/*<       END IF >*/
    }
/*<       CALL SQCOPY(GELCOP,GEL,IDIMG) >*/
    sqcopy_(gelcop + 1, gel + 1, idimg, (ftnlen)1, (ftnlen)1);
/*<       KSENSE = 0 >*/
    ksense = 0;
/*      WRITE(INFOD,1000)IMATC */
/* 1000 FORMAT('Total matches found',I6) */
/* CHECKED */
/*<       CALL SWRT1(INFOD, 'Total matches found%6d%!', IMATC) >*/
    swrt1_(infod, "Total matches found%6d%!", imatc, (ftnlen)80, (ftnlen)24);
/*<       CALL INFO(INFOD) >*/
    info_(infod, (ftnlen)80);
/*<       IF(IMATC.EQ.0) THEN >*/
    if (*imatc == 0) {
/*<         IFAIL(1) = 0 >*/
	ifail[1] = 0;
/*<         IFCOMP = 0 >*/
	*ifcomp = 0;
/*<         RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       DO 99 I = 1,IMATC >*/
    i__1 = *imatc;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        WRITE(INFOD,1002)JLLINO(I),JPOSC(I),JSENSE(I), */
/*     +  JPOSG(I) */
/* 1002   FORMAT */
/*     +  ('Contig',I8,' position',I8,' matches strand ',I2, */
/*     +  ' at position',I8) */
/* CHECKED */
/*<    >*/
	swrt4_(infod, "Contig%8d position%8d matches strand %2d at position%"
		"8d%!", &jllino[i__ - 1], &jposc[i__ - 1], &jsense[i__ - 1], &
		jposg[i__ - 1], (ftnlen)80, (ftnlen)57);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<  99   CONTINUE >*/
/* L99: */
    }
/*<       JMATC = 0 >*/
    jmatc = 0;
/*<       DO 100 I = 1,IMATC >*/
    i__1 = *imatc;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*  3-10-95 New idea! have minimum overlap before allowing entry. */
/*  Simplest place to apply is here (though not where it should */
/*  be if we had started afresh) */
/*<    >*/
/* Computing MIN */
	i__2 = jposg[i__ - 1], i__3 = jposc[i__ - 1];
/* Computing MIN */
	i__4 = *idimg - jposg[i__ - 1], i__5 = jlc[i__ - 1] - jposc[i__ - 1];
	lenovr = min(i__2,i__3) + min(i__4,i__5);
/*<        IF (LENOVR.LT.MINOVR) THEN >*/
	if (lenovr < *minovr) {
/*         WRITE(*,*)'SHORT OVERLAP', */
/*     +   LENOVR,JPOSG(I),IDIMG,JPOSC(I),JLC(I) */
/*<          GO TO 100 >*/
	    goto L100;
/*<        END IF >*/
	}



/*         WRITE(*,*)'*******LONG OVERLAP', */
/*     +   LENOVR,JPOSG(I),IDIMG,JPOSC(I),JLC(I) */
/*<         IF(JSENSE(I).EQ.-1) THEN >*/
	if (jsense[i__ - 1] == -1) {
/*<           IF(KSENSE.EQ.0) THEN  >*/
	    if (ksense == 0) {
/*<             CALL SQREV(GEL,IDIMG) >*/
		sqrev_(gel + 1, idimg, (ftnlen)1);
/*<             CALL SQCOM(GEL,IDIMG) >*/
		sqcom_(gel + 1, idimg, (ftnlen)1);
/*<             KSENSE = 1 >*/
		ksense = 1;
/*<           END IF >*/
	    }
/*<         END IF >*/
	}
/*<         JDIM22 = IDIMG >*/
	jdim22 = *idimg;
/*<         JDOUT = MAXGEL >*/
	jdout = *maxgel;
/*<         IDSAV = MAXSAV >*/
	idsav = *maxsav;
/*        WRITE(INFOD,1001)JLLINO(I) */
/* 1001   FORMAT('Trying to align with contig ',I8) */
/* CHECKED */
/*<         CALL SWRT1(INFOD,'Trying to align with contig %8d%!',JLLINO(I)) >*/
	swrt1_(infod, "Trying to align with contig %8d%!", &jllino[i__ - 1], (
		ftnlen)80, (ftnlen)33);
/*<         CALL INFO(INFOD) >*/
	info_(infod, (ftnlen)80);
/*<    >*/
	aline_(seq1 + jlefts[i__ - 1], gel + 1, seqg3 + 1, seqc3 + 1, &savps[
		1], &savpg[1], &savl[1], &idsav, &jlc[i__ - 1], &jdim22, &
		jdout, &jposc[i__ - 1], &jposg[i__ - 1], minsli, &jjoint, &
		jtotpc, &jtotpg, &jfail, &jtype, maxpc, maxpg, permax, seq4 + 
		1, maxgel, &perms, leno, ishow, mask, &c__0, (ftnlen)1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/*<         IF(JFAIL.EQ.0) THEN >*/
	if (jfail == 0) {
/*<           JMATC = JMATC + 1 >*/
	    ++jmatc;
/*<           IF(JMATC.EQ.1) THEN >*/
	    if (jmatc == 1) {
/*    Save in elements 1 */
/*<    >*/
		copym_(&jlefts[i__ - 1], &ilefts[1], &jlc[i__ - 1], &ilc[1], &
			jposc[i__ - 1], &iposc[1], &jsense[i__ - 1], &isense[
			1], &jllino[i__ - 1], &llino[1], &jjoint, &joint[1], &
			jtotpc, &itotpc[1], &jtotpg, &itotpg[1], &jtype, &
			itype[1], &jdout, &idout[1], &jdim22, &idim22[1], 
			seqg3 + 1, seqg2 + (seqg2_dim1 + 1), seqc3 + 1, seqc2 
			+ (seqc2_dim1 + 1), &perms, &permis[1], (ftnlen)1, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
/*<             IFAIL(1) = 0 >*/
		ifail[1] = 0;
/*<           ELSE IF(JMATC.EQ.2) THEN >*/
	    } else if (jmatc == 2) {
/*<             IF(PERMS.LT.PERMIS(1)) THEN >*/
		if (perms < permis[1]) {
/*    Better match so save in elements 1, so copy 1 to 2 first */
/*<    >*/
		    copym_(&ilefts[1], &ilefts[2], &ilc[1], &ilc[2], &iposc[1]
			    , &iposc[2], &isense[1], &isense[2], &llino[1], &
			    llino[2], &joint[1], &joint[2], &itotpc[1], &
			    itotpc[2], &itotpg[1], &itotpg[2], &itype[1], &
			    itype[2], &idout[1], &idout[2], &idim22[1], &
			    idim22[2], seqg2 + (seqg2_dim1 + 1), seqg2 + ((
			    seqg2_dim1 << 1) + 1), seqc2 + (seqc2_dim1 + 1), 
			    seqc2 + ((seqc2_dim1 << 1) + 1), &permis[1], &
			    permis[2], (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			    ftnlen)1);
/*<                 IFAIL(2) = 0 >*/
		    ifail[2] = 0;
/*    Now save in 1 */
/*<    >*/
		    copym_(&jlefts[i__ - 1], &ilefts[1], &jlc[i__ - 1], &ilc[
			    1], &jposc[i__ - 1], &iposc[1], &jsense[i__ - 1], 
			    &isense[1], &jllino[i__ - 1], &llino[1], &jjoint, 
			    &joint[1], &jtotpc, &itotpc[1], &jtotpg, &itotpg[
			    1], &jtype, &itype[1], &jdout, &idout[1], &jdim22,
			     &idim22[1], seqg3 + 1, seqg2 + (seqg2_dim1 + 1), 
			    seqc3 + 1, seqc2 + (seqc2_dim1 + 1), &perms, &
			    permis[1], (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			    ftnlen)1);
/*<             ELSE >*/
		} else {
/*    Save in element 2 */
/*<    >*/
		    copym_(&jlefts[i__ - 1], &ilefts[2], &jlc[i__ - 1], &ilc[
			    2], &jposc[i__ - 1], &iposc[2], &jsense[i__ - 1], 
			    &isense[2], &jllino[i__ - 1], &llino[2], &jjoint, 
			    &joint[2], &jtotpc, &itotpc[2], &jtotpg, &itotpg[
			    2], &jtype, &itype[2], &jdout, &idout[2], &jdim22,
			     &idim22[2], seqg3 + 1, seqg2 + ((seqg2_dim1 << 1)
			     + 1), seqc3 + 1, seqc2 + ((seqc2_dim1 << 1) + 1),
			     &perms, &permis[2], (ftnlen)1, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
/*<               IFAIL(2) = 0 >*/
		    ifail[2] = 0;
/*<             END IF >*/
		}
/*<           ELSE >*/
	    } else {
/*<             IF(PERMS.LT.PERMIS(1)) THEN >*/
		if (perms < permis[1]) {
/*    Better match so save in elements 1, so copy 1 to 2 first */
/*<    >*/
		    copym_(&ilefts[1], &ilefts[2], &ilc[1], &ilc[2], &iposc[1]
			    , &iposc[2], &isense[1], &isense[2], &llino[1], &
			    llino[2], &joint[1], &joint[2], &itotpc[1], &
			    itotpc[2], &itotpg[1], &itotpg[2], &itype[1], &
			    itype[2], &idout[1], &idout[2], &idim22[1], &
			    idim22[2], seqg2 + (seqg2_dim1 + 1), seqg2 + ((
			    seqg2_dim1 << 1) + 1), seqc2 + (seqc2_dim1 + 1), 
			    seqc2 + ((seqc2_dim1 << 1) + 1), &permis[1], &
			    permis[2], (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			    ftnlen)1);
/*<                 IFAIL(2) = 0 >*/
		    ifail[2] = 0;
/*    Now save in 1 */
/*<    >*/
		    copym_(&jlefts[i__ - 1], &ilefts[1], &jlc[i__ - 1], &ilc[
			    1], &jposc[i__ - 1], &iposc[1], &jsense[i__ - 1], 
			    &isense[1], &jllino[i__ - 1], &llino[1], &jjoint, 
			    &joint[1], &jtotpc, &itotpc[1], &jtotpg, &itotpg[
			    1], &jtype, &itype[1], &jdout, &idout[1], &jdim22,
			     &idim22[1], seqg3 + 1, seqg2 + (seqg2_dim1 + 1), 
			    seqc3 + 1, seqc2 + (seqc2_dim1 + 1), &perms, &
			    permis[1], (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			    ftnlen)1);
/*<             ELSE IF(PERMS.LT.PERMIS(2)) THEN >*/
		} else if (perms < permis[2]) {
/*    Save in element 2 */
/*<    >*/
		    copym_(&jlefts[i__ - 1], &ilefts[2], &jlc[i__ - 1], &ilc[
			    2], &jposc[i__ - 1], &iposc[2], &jsense[i__ - 1], 
			    &isense[2], &jllino[i__ - 1], &llino[2], &jjoint, 
			    &joint[2], &jtotpc, &itotpc[2], &jtotpg, &itotpg[
			    2], &jtype, &itype[2], &jdout, &idout[2], &jdim22,
			     &idim22[2], seqg3 + 1, seqg2 + ((seqg2_dim1 << 1)
			     + 1), seqc3 + 1, seqc2 + ((seqc2_dim1 << 1) + 1),
			     &perms, &permis[2], (ftnlen)1, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
/*<             END IF >*/
		}
/*<           END IF >*/
	    }
/*<         ELSE >*/
	} else {
/*<           KFAIL = 1 >*/
	    *kfail = 1;
/*<         END IF >*/
	}
/*< 100   CONTINUE >*/
L100:
	;
    }
/*<       IMATC = MIN(2,JMATC) >*/
    *imatc = min(2,jmatc);
/*<       IFCOMP = 0 >*/
    *ifcomp = 0;
/*<       END >*/
    return 0;
} /* autocn_ */

/*     BUBBL3 */
/*   SUBROUTINE TO SORT INTEGER ARRAY (LIST) INTO ASCENDING  ORDER */

/*<       SUBROUTINE CCTA(SEQ,ID) >*/
/* Subroutine */ int ccta_(char *seq, integer *id, ftnlen seq_len)
{
    /* Initialized data */

    static char com[1] = ",";
    static char as[1] = "*";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*<       CHARACTER SEQ(ID),COM,AS >*/
/*<       SAVE COM,AS >*/
/*<       DATA COM/','/,AS/'*'/ >*/
    /* Parameter adjustments */
    --seq;

    /* Function Body */
/*<       DO 10 I = 1,ID >*/
    i__1 = *id;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         IF(SEQ(I).EQ.COM) SEQ(I) = AS >*/
	if (*(unsigned char *)&seq[i__] == *(unsigned char *)&com[0]) {
	    *(unsigned char *)&seq[i__] = *(unsigned char *)&as[0];
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       END >*/
    return 0;
} /* ccta_ */

/*<    >*/
/* Subroutine */ int copym_(integer *jlefts, integer *ilefts, integer *jlc, 
	integer *ilc, integer *jposc, integer *iposc, integer *jsense, 
	integer *isense, integer *jllino, integer *llino, integer *jjoint, 
	integer *joint, integer *jtotpc, integer *itotpc, integer *jtotpg, 
	integer *itotpg, integer *jtype, integer *itype, integer *jdout, 
	integer *idout, integer *jdim22, integer *idim22, char *seqg3, char *
	seqg2, char *seqc3, char *seqc2, real *perms, real *permis, ftnlen 
	seqg3_len, ftnlen seqg2_len, ftnlen seqc3_len, ftnlen seqc2_len)
{
    extern /* Subroutine */ int sqcopy_(char *, char *, integer *, ftnlen, 
	    ftnlen);

/*<       CHARACTER SEQG3(JDIM22),SEQG2(JDIM22),SEQC3(JDOUT),SEQC2(JDOUT) >*/
/*<       ILEFTS = JLEFTS >*/
    /* Parameter adjustments */
    --seqc2;
    --seqc3;
    --seqg2;
    --seqg3;

    /* Function Body */
    *ilefts = *jlefts;
/*<       ILC = JLC >*/
    *ilc = *jlc;
/*<       IPOSC = JPOSC >*/
    *iposc = *jposc;
/*<       ISENSE = JSENSE >*/
    *isense = *jsense;
/*<       LLINO = JLLINO >*/
    *llino = *jllino;
/*<       JOINT = JJOINT >*/
    *joint = *jjoint;
/*<       ITOTPC = JTOTPC >*/
    *itotpc = *jtotpc;
/*<       ITOTPG = JTOTPG >*/
    *itotpg = *jtotpg;
/*<       ITYPE = JTYPE >*/
    *itype = *jtype;
/*<       IDOUT = JDOUT >*/
    *idout = *jdout;
/*<       IDIM22 = JDIM22 >*/
    *idim22 = *jdim22;
/*<       CALL SQCOPY(SEQG3,SEQG2,JDIM22) >*/
    sqcopy_(seqg3 + 1, seqg2 + 1, jdim22, (ftnlen)1, (ftnlen)1);
/*<       CALL SQCOPY(SEQC3,SEQC2,JDOUT) >*/
    sqcopy_(seqc3 + 1, seqc2 + 1, jdout, (ftnlen)1, (ftnlen)1);
/*<       PERMIS = PERMS >*/
    *permis = *perms;
/*<       END >*/
    return 0;
} /* copym_ */

/*     SUBROUTINE DALIGN */

/*   COUNTS MISMATCHES AND DISPLAYS OVERLAP. */
/*<    >*/
/* Subroutine */ int dalign_(char *seqc2, char *seqg2, char *seq3, integer *
	maxgel, integer *idout, integer *idim2, integer *joint, integer *
	itype, real *x, integer *ifail, integer *lo, real *permax, integer *
	ishow, integer *maxpg, integer *maxpc, integer *itotpg, integer *
	itotpc, ftnlen seqc2_len, ftnlen seqg2_len, ftnlen seq3_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static real y;
    static integer kc, lg;
    extern /* Subroutine */ int info_(char *, ftnlen);
    static char name1[15], name2[15];
    extern /* Subroutine */ int swrt0_(char *, char *, ...);
    static integer iendc;
    extern /* Subroutine */ int swrt3_(char *, char *, ...);
    static integer iendg;
    static char infod[80];
    extern integer forta_(char *, char *, integer *, char *, char *, integer *
	    , integer *, integer *, char *, integer *, ftnlen, ftnlen, ftnlen,
	     ftnlen, ftnlen);
    extern /* Subroutine */ int erromf_(char *, ftnlen);
    extern integer ctonum_(char *, ftnlen);
    extern /* Subroutine */ int mstlkl_(char *, integer *, ftnlen), sqcopy_(
	    char *, char *, integer *, ftnlen, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQC2(MAXGEL),SEQG2(MAXGEL),SEQ3(MAXGEL) >*/
/*      CHARACTER PAD,DASH */
/*<       CHARACTER INFOD*80,NAME1*15,NAME2*15 >*/
/*<       INTEGER CTONUM, FORTA >*/
/*<       EXTERNAL CTONUM, FORTA >*/
/*      SAVE PAD,DASH */
/*      DATA PAD,DASH/',','-'/ */
/*<       IENDG=1 >*/
    /* Parameter adjustments */
    --seq3;
    --seqg2;
    --seqc2;

    /* Function Body */
    iendg = 1;
/*<       IENDC=JOINT >*/
    iendc = *joint;
/*   ONLY LOOK AT OVERLAP WHICH IS FROM JOINT FOR LEFT TYPE JOIN */
/*<       IF(ITYPE.EQ.1)THEN >*/
    if (*itype == 1) {
/*<         IENDG=JOINT >*/
	iendg = *joint;
/*<         IENDC=1 >*/
	iendc = 1;
/*<       END IF >*/
    }
/*< 100   CONTINUE >*/
/* L100: */
/*   LENGTH OF OVERLAP? */
/*<       LG=IDIM2-IENDG+1 >*/
    lg = *idim2 - iendg + 1;
/*<       LO=MIN(IDOUT,LG) >*/
    *lo = min(*idout,lg);
/*   SAVE RAW DATA */
/*<       CALL SQCOPY(SEQG2,SEQ3,IDIM2) >*/
    sqcopy_(seqg2 + 1, seq3 + 1, idim2, (ftnlen)1, (ftnlen)1);
/*<       CALL MSTLKL(SEQ3,IDIM2) >*/
    mstlkl_(seq3 + 1, idim2, (ftnlen)1);
/*<       X=FLOAT(LO) >*/
    *x = (real) (*lo);
/*<       Y=X >*/
    y = *x;
/*<       K=IENDG+LO-1 >*/
    k = iendg + *lo - 1;
/*   POINT TO CONSENSUS */
/*<       J=0 >*/
    j = 0;
/*   CHECK FOR OVERFLOW */
/*<       IF(K.GT.MAXGEL)THEN >*/
    if (k > *maxgel) {
/*<         CALL ERROMF('DALIGN: matching region too long') >*/
	erromf_("DALIGN: matching region too long", (ftnlen)32);
/*<         IFAIL = 1 >*/
	*ifail = 1;
/*<         RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       DO 200 I=IENDG,K >*/
    i__1 = k;
    for (i__ = iendg; i__ <= i__1; ++i__) {
/*<         J=J+1 >*/
	++j;
/*<         IF(CTONUM(SEQC2(J)).EQ.CTONUM(SEQ3(I))) THEN >*/
	if (ctonum_(seqc2 + j, (ftnlen)1) == ctonum_(seq3 + i__, (ftnlen)1)) {
/*<           IF(CTONUM(SEQC2(J)).LT.5) GO TO 200 >*/
	    if (ctonum_(seqc2 + j, (ftnlen)1) < 5) {
		goto L200;
	    }
/*         SO FAR THEY ARE = AND ACGT, WHAT IS LEFT IS = AND 5 */
/*<    >*/
	    if ((*(unsigned char *)&seqc2[j] == '*' || *(unsigned char *)&
		    seqc2[j] == ',') && (*(unsigned char *)&seq3[i__] == '*' 
		    || *(unsigned char *)&seq3[i__] == ',')) {
		goto L200;
	    }
/*<         END IF >*/
	}
/*<         X=X-1. >*/
	*x += -1.f;
/*< 200   CONTINUE >*/
L200:
	;
    }
/*<       X=(Y-X)*100./Y >*/
    *x = (y - *x) * 100.f / y;
/*<       IFAIL=0 >*/
    *ifail = 0;
/*<       IF (X.GT.PERMAX) IFAIL = 1 >*/
    if (*x > *permax) {
	*ifail = 1;
    }
/*<       IF(ITOTPC.GT.MAXPC) IFAIL = 1 >*/
    if (*itotpc > *maxpc) {
	*ifail = 1;
    }
/*<       IF(ITOTPG.GT.MAXPG) IFAIL = 1 >*/
    if (*itotpg > *maxpg) {
	*ifail = 1;
    }

/* ISHOW 1 hide all alignments */
/*       2 show passes */
/*       3 show all alignments */
/*       4 show failures only */

/*          WRITE(*,*)X,ITOTPC,ITOTPG */
/*<       IF (ISHOW.EQ.1) THEN >*/
    if (*ishow == 1) {
/*<         IF(IFAIL.EQ.0) THEN >*/
	if (*ifail == 0) {
/*          WRITE(INFOD,1052)X,ITOTPC,ITOTPG */
/* 1052      FORMAT('Percent mismatch ',F4.1,', pads in contig',I3, */
/*     +    ', pads in gel',I3) */
/* CHECKED */
/*<    >*/
	    swrt3_(infod, "Percent mismatch %4.1f, pads in contig%3d, pads i"
		    "n gel%3d%!", x, itotpc, itotpg, (ftnlen)80, (ftnlen)59);
/*<           CALL INFO(INFOD) >*/
	    info_(infod, (ftnlen)80);
/*<         END IF >*/
	}
/*<         RETURN >*/
	return 0;
/*<       ELSE IF(ISHOW.EQ.2) THEN >*/
    } else if (*ishow == 2) {
/*<         IF (IFAIL.NE.0) RETURN >*/
	if (*ifail != 0) {
	    return 0;
	}
/*<       ELSE IF(ISHOW.EQ.4) THEN >*/
    } else if (*ishow == 4) {
/*<         IF (IFAIL.EQ.0) THEN >*/
	if (*ifail == 0) {
/*          WRITE(INFOD,1052)X,ITOTPC,ITOTPG */
/* CHECKED */
/*<    >*/
	    swrt3_(infod, "Percent mismatch %4.1f, pads in contig%3d, pads i"
		    "n gel%3d%!", x, itotpc, itotpg, (ftnlen)80, (ftnlen)59);
/*<           CALL INFO(INFOD) >*/
	    info_(infod, (ftnlen)80);
/*<           RETURN >*/
	    return 0;
/*<         END IF >*/
	}
/*<       END IF >*/
    }
/*      WRITE(INFOD,1052)X,ITOTPC,ITOTPG */
/* CHECKED */
/*<    >*/
    swrt3_(infod, "Percent mismatch %4.1f, pads in contig%3d, pads in gel%3d"
	    "%!", x, itotpc, itotpg, (ftnlen)80, (ftnlen)59);
/*      WRITE(NAME2,1000)'     Consensus' */
/*      WRITE(NAME1,1000)'       Reading' */
/* 1000  FORMAT(A) */
/* CHECKED */
/*<       CALL SWRT0(NAME2, '    Consensus %!') >*/
    swrt0_(name2, "    Consensus %!", (ftnlen)15, (ftnlen)16);
/*<       CALL SWRT0(NAME1, '      Reading %!') >*/
    swrt0_(name1, "      Reading %!", (ftnlen)15, (ftnlen)16);
/*<    >*/
    i__1 = i_len(name1, (ftnlen)15);
    kc = forta_(seqc2 + 1, seqg2 + iendg, lo, name2, name1, &i__1, &iendc, &
	    iendg, infod, &c__80, (ftnlen)1, (ftnlen)1, (ftnlen)15, (ftnlen)
	    15, (ftnlen)80);
/*<       END >*/
    return 0;
} /* dalign_ */

/*     DELCON */

/*   DELETES CONTIG FROM CONSENSUS SEQUENCE */
/*<       SUBROUTINE DELCON(SEQ1,ILEFT,ILC,IDIM1) >*/
/* Subroutine */ int delcon_(char *seq1, integer *ileft, integer *ilc, 
	integer *idim1, ftnlen seq1_len)
{
    static integer i1, i2, id;
    extern /* Subroutine */ int sqcopy_(char *, char *, integer *, ftnlen, 
	    ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQ1(IDIM1) >*/
/*   FIRST CHAR TO REPLACE */
/*<       I1=ILEFT-20 >*/
    /* Parameter adjustments */
    --seq1;

    /* Function Body */
    i1 = *ileft - 20;
/*   FIRST CHAR TO MOVE */
/*<       I2=ILEFT+ILC >*/
    i2 = *ileft + *ilc;
/*   IS THIS RIGHTMOST CONTIG ANYWAY? */
/*<       IF(I2.GT.IDIM1)GO TO 10 >*/
    if (i2 > *idim1) {
	goto L10;
    }
/*   NUMBER TO MOVE */
/*<       ID=IDIM1-I2+1 >*/
    id = *idim1 - i2 + 1;
/*   MOVE */
/*<       CALL SQCOPY(SEQ1(I2),SEQ1(I1),ID) >*/
    sqcopy_(seq1 + i2, seq1 + i1, &id, (ftnlen)1, (ftnlen)1);
/*   RESET LENGTH */
/*<       IDIM1=I1+ID-1 >*/
    *idim1 = i1 + id - 1;
/*<       RETURN >*/
    return 0;
/*< 10    CONTINUE >*/
L10:
/*   RIGHTMOST CONTIG SO DONT MOVE */
/*<       IDIM1=I1-1 >*/
    *idim1 = i1 - 1;

/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* delcon_ */

/*     LINEUP */

/*   TAKES 2 SEQS SET OF MATCHES AND PRODUCES LINED UP SEQS */
/*   FINDS IF WE HAVE A LEFT OVERLAP */
/*   RETURNS POSITION OF JOINT. THIS IS RELATIVE TO THE CONTIG */
/*   FOR MOST MATCHES BUT I RELATIVE TO THE GEL FOR A LEFT OVERLAP */
/*<    >*/
/* Subroutine */ int lineup_(char *seqg, char *seqc, char *seqg2, char *seqc2,
	 integer *idc, integer *idg, integer *idout, integer *matg, integer *
	matc, integer *matl, integer *ip, integer *itotpc, integer *itotpg, 
	integer *joint, integer *itype, integer *maxgel, integer *ifail, 
	ftnlen seqg_len, ftnlen seqc_len, ftnlen seqg2_len, ftnlen seqc2_len)
{
    /* Initialized data */

    static char pad[1] = ",";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, l, m, l5, ic1, ic2, lc1, ig1, ig2, lc2, lg1, lg2, 
	    is1, is2, nmtch;
    extern /* Subroutine */ int padcop_(char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), erromf_(char *, 
	    ftnlen), sqcopy_(char *, char *, integer *, ftnlen, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQG(IDG),SEQC(IDC),SEQG2(IDOUT),SEQC2(IDOUT),PAD >*/
/*<       INTEGER MATG(IP),MATC(IP),MATL(IP) >*/
/*<       SAVE PAD >*/
/*<       DATA PAD/','/ >*/
    /* Parameter adjustments */
    --seqc;
    --seqg;
    --seqc2;
    --seqg2;
    --matl;
    --matc;
    --matg;

    /* Function Body */
/*<       IFAIL=0 >*/
    *ifail = 0;
/*   ZERO PADDING CHARS IN CONTIG (GEL DONE AT END BY DIFFERENCE */
/*   IN INPUT AND OUTPUT LENGTHS) */
/*<       ITOTPC=0 >*/
    *itotpc = 0;
/*   FILL OUTPUT WITH PADDING */
/*<       DO 10 I=1,IDOUT >*/
    i__1 = *idout;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         SEQG2(I)=PAD >*/
	*(unsigned char *)&seqg2[i__] = *(unsigned char *)&pad[0];
/*<         SEQC2(I)=PAD >*/
	*(unsigned char *)&seqc2[i__] = *(unsigned char *)&pad[0];
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       NMTCH=0 >*/
    nmtch = 0;
/*   SET INITIAL POINTERS TO OUTPUT */
/*   CONSENSUS */
/*<       IS1=1 >*/
    is1 = 1;
/*   GEL */
/*<       IS2=1 >*/
    is2 = 1;
/*   FIND DISTANCE FROM LEFT MATCH IN GEL TO LEFT OF GEL */
/*<       IG2=MATG(1)-1 >*/
    ig2 = matg[1] - 1;
/*<       IF(IG2.EQ.0)THEN >*/
    if (ig2 == 0) {
/*       THE LEFT END OF THE GEL MATCHES SO THIS IS NOT A LEFT OVERLAP */
/*       SET TYPE */
/*<         ITYPE=-1 >*/
	*itype = -1;
/*       SET JOINT */
/*<         JOINT=MATC(1) >*/
	*joint = matc[1];
/*       SKIP NEXT SECTION */
/*<         GO TO 50 >*/
	goto L50;
/*<       END IF >*/
    }
/*   FIND DISTANCE FROM LEFT MATCH IN CONTIG TO LEFT OF CONTIG */
/*<       IC2=MATC(1)-1 >*/
    ic2 = matc[1] - 1;
/*   GET DISTANCE FROM FIRST MATCH IN CONTIG TO FIRST MATCH IN GEL. */
/*   IF THIS DISTANCE <0 THEN WE HAVE A LEFT OVERLAP */
/*<       IC1=IC2-IG2+1 >*/
    ic1 = ic2 - ig2 + 1;
/*<       IF(IC1.GT.0)THEN >*/
    if (ic1 > 0) {
/*       THIS IS NOT A LEFT OVERLAP */
/*       SET TYPE */
/*<         ITYPE=-1 >*/
	*itype = -1;
/*       SET LEFT END */
/*<         JOINT=IC1 >*/
	*joint = ic1;
/*       COPY THE GEL UPTO THE FIRST MATCH, INTO THE OUTPUT ARRAY */
/*       CHECK FOR OVERFLOW */
/*<         IF(IG2.GT.MAXGEL)GO TO 700 >*/
	if (ig2 > *maxgel) {
	    goto L700;
	}
/*<         CALL SQCOPY(SEQG(1),SEQG2(1),IG2) >*/
	sqcopy_(seqg + 1, seqg2 + 1, &ig2, (ftnlen)1, (ftnlen)1);
/*       COPY THE CONTIG FOR THE SAME REGION */
/*<         IF(IG2.GT.MAXGEL)GO TO 700 >*/
	if (ig2 > *maxgel) {
	    goto L700;
	}
/*<         CALL SQCOPY(SEQC(IC1),SEQC2(1),IG2) >*/
	sqcopy_(seqc + ic1, seqc2 + 1, &ig2, (ftnlen)1, (ftnlen)1);
/*<         IS1=IS1+IG2 >*/
	is1 += ig2;
/*<         IS2=IS2+IG2 >*/
	is2 += ig2;
/*<         GO TO 50 >*/
	goto L50;
/*<       END IF >*/
    }
/*   MUST BE LEFT END OVERLAP */
/*   SET TYPE */
/*<       ITYPE=1 >*/
    *itype = 1;
/*   SET POSITION OF JOINT RELATIVE TO GEL */
/*<       JOINT=ABS(IC1)+2 >*/
    *joint = abs(ic1) + 2;
/*   COPY OVER THE GEL UPTO THE JOINT */
/*   CHECK FOR OVERFLOW */
/*<       IF(IG2.GT.MAXGEL)GO TO 700 >*/
    if (ig2 > *maxgel) {
	goto L700;
    }
/*<       CALL SQCOPY(SEQG(1),SEQG2(1),IG2) >*/
    sqcopy_(seqg + 1, seqg2 + 1, &ig2, (ftnlen)1, (ftnlen)1);
/*<       IS2=IS2+IG2 >*/
    is2 += ig2;
/*   WE MAY ALSO HAVE MISMATCHING */
/*   DATA AT THE JOIN SO DEAL WITH THAT NOW */
/*   IF IC2 >0 THE LEFT END OF THE CONTIG MATCHES THE GEL BUT OTHERWISE */
/*   WE HAVE SOME MISMATCHED DATA TO DEAL WITH - WE NEED TO TRANSFER */
/*   THE MISMATCHED REGION OF THE CONTIG TO THE OUTPUT ARRAY */
/*<       IF(IC2.GT.0)THEN >*/
    if (ic2 > 0) {
/*<         IF(IC2.GT.MAXGEL)GO TO 700 >*/
	if (ic2 > *maxgel) {
	    goto L700;
	}
/*<         CALL SQCOPY(SEQC(1),SEQC2(1),IC2) >*/
	sqcopy_(seqc + 1, seqc2 + 1, &ic2, (ftnlen)1, (ftnlen)1);
/*<         IS1=IS1+IC2 >*/
	is1 += ic2;
/*<       END IF >*/
    }
/*   WHEN WE GET HERE WE HAVE SORTED OUT THE LEFT ENDS FOR LEFT OVERLAP */
/*   AND MISMATCHED LEFT ENDS, WE NOW DEAL WITH THE REST OF THE SEQUENCE */
/*   STARTING WITH THE FIRST BLOCK OF IDENTITY */

/* IG1 POSITION IN INPUT GEL */
/* IS2 POSITION IN OUTPUT GEL */
/* IC1 POSITION IN INPUT CONTIG */
/* IS1 POSITION IN OUTPUT CONTIG */
/* LG1 POSITION OF END OF CURRENT MATCH IN OUTPUT GEL */
/* LC1 POSITION OF END OF CURRENT MATCH IN OUTPUT CONTIG */
/* LG2 DISTANCE FROM CURRENT MATCH IN INPUT GEL TO NEXT MATCH */
/* LC2 DISTANCE FROM CURRENT MATCH IN INPUT CONTIG TO NEXT MATCH */

/*< 50    CONTINUE >*/
L50:
/*   POINT TO NEXT MATCH */
/*<       NMTCH=NMTCH+1 >*/
    ++nmtch;
/*   COPY NEXT MATCH */
/*<       IG1=MATG(NMTCH) >*/
    ig1 = matg[nmtch];
/*<       IC1=MATC(NMTCH) >*/
    ic1 = matc[nmtch];
/*<       L=MATL(NMTCH) >*/
    l = matl[nmtch];
/*   CHECK FOR OVERFLOW */
/*<       IF(IS2+L-1.GT.MAXGEL)GO TO 700 >*/
    if (is2 + l - 1 > *maxgel) {
	goto L700;
    }
/*<       CALL SQCOPY(SEQG(IG1),SEQG2(IS2),L) >*/
    sqcopy_(seqg + ig1, seqg2 + is2, &l, (ftnlen)1, (ftnlen)1);
/*   CHECK FOR OVERFLOW */
/*<       IF(IS1+L-1.GT.MAXGEL)GO TO 700 >*/
    if (is1 + l - 1 > *maxgel) {
	goto L700;
    }
/*<       CALL SQCOPY(SEQC(IC1),SEQC2(IS1),L) >*/
    sqcopy_(seqc + ic1, seqc2 + is1, &l, (ftnlen)1, (ftnlen)1);
/*   POINT TO NEXT OUTPUT POSITIONS */
/*<       IS1=IS1+L >*/
    is1 += l;
/*<       IS2=IS2+L >*/
    is2 += l;
/*   END OF CURRENT MATCH */
/*<       LG1=IG1+L >*/
    lg1 = ig1 + l;
/*<       LC1=IC1+L >*/
    lc1 = ic1 + l;
/*   ANY MORE MATCHES */
/*<       IF(NMTCH.EQ.IP)GO TO 500 >*/
    if (nmtch == *ip) {
	goto L500;
    }
/*<       K=NMTCH+1 >*/
    k = nmtch + 1;
/*<       LG2=MATG(K)-LG1 >*/
    lg2 = matg[k] - lg1;
/*<       LC2=MATC(K)-LC1 >*/
    lc2 = matc[k] - lc1;
/*   ANY DIFFERENCE IN LENGTH? IF SO WE HAVE TO PAD SO THEY BECOME THE SAME */
/*<       L5=ABS(LG2-LC2) >*/
    l5 = (i__1 = lg2 - lc2, abs(i__1));
/*   COUNT PADDING CHARS IN CONTIG */
/*<       IF(LG2.GT.LC2)ITOTPC=ITOTPC+L5 >*/
    if (lg2 > lc2) {
	*itotpc += l5;
    }
/*   IF DIFFERENCE INCREMENT SHORTER */
/*<       IF(LG2.GT.LC2)IS1=IS1+L5 >*/
    if (lg2 > lc2) {
	is1 += l5;
    }
/*   IF GEL NEEDS PADDING TRY TO PUT PADS NEXT TO DOUBLE CODES */
/*<    >*/
    if (lc2 > lg2) {
	padcop_(seqg + 1, seqg2 + 1, &lg1, &matg[k], &l5, &is2, &lg2, maxgel, 
		ifail, seqc + 1, idc, &lc1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    }
/*   CHECK FOR OVERFLOW */
/*<       IF(IFAIL.EQ.1)GO TO 700 >*/
    if (*ifail == 1) {
	goto L700;
    }
/*   NOW COPY MISSMATCHED REGION */
/*   CHECK FOR OVERFLOW */
/*<       IF(IS2+LG2-1.GT.MAXGEL)GO TO 700 >*/
    if (is2 + lg2 - 1 > *maxgel) {
	goto L700;
    }
/*<       IF(LG2.GT.0)CALL SQCOPY(SEQG(LG1),SEQG2(IS2),LG2) >*/
    if (lg2 > 0) {
	sqcopy_(seqg + lg1, seqg2 + is2, &lg2, (ftnlen)1, (ftnlen)1);
    }
/*   CHECK FOR OVERFLOW */
/*<       IF(IS1+LC2-1.GT.MAXGEL)GO TO 700 >*/
    if (is1 + lc2 - 1 > *maxgel) {
	goto L700;
    }
/*<       IF(LC2.GT.0)CALL SQCOPY(SEQC(LC1),SEQC2(IS1),LC2) >*/
    if (lc2 > 0) {
	sqcopy_(seqc + lc1, seqc2 + is1, &lc2, (ftnlen)1, (ftnlen)1);
    }
/*   POINT TO NEXT OUTPUT POSITIONS */
/*<       IS1=IS1+LC2 >*/
    is1 += lc2;
/*<       IS2=IS2+LG2 >*/
    is2 += lg2;
/*   GET NEXT MATCH */
/*<       GO TO 50 >*/
    goto L50;
/*< 500   CONTINUE >*/
L500:

/*   FINISH RIGHT ENDS */
/*   ONLY COPY TO END OF GEL IN GEL AND TO THE SAME RELATIVE POSITION */
/*   IN THE CONTIG FOR DISPLAY PURPOSES AND FOR COUNTING MISMATCH */
/*   CURRENT ENDS AT LG1,LC1 */
/*   HOW FAR TO END OF GEL? */
/*   SET M */
/*<       M=0 >*/
    m = 0;
/*<       L=IDG-LG1+1 >*/
    l = *idg - lg1 + 1;
/*<       IF(L.LT.1)GO TO 600 >*/
    if (l < 1) {
	goto L600;
    }
/*   CHECK FOR OVERFLOW */
/*<       IF(IS2+L-1.GT.MAXGEL)GO TO 700 >*/
    if (is2 + l - 1 > *maxgel) {
	goto L700;
    }
/*<       CALL SQCOPY(SEQG(LG1),SEQG2(IS2),L) >*/
    sqcopy_(seqg + lg1, seqg2 + is2, &l, (ftnlen)1, (ftnlen)1);
/*   NEED TO COPY TO END OF GEL IN CONTIG FOR DISPLAY */
/*   POINT TO POSN IN CONTIG LEVEL WITH END OF GEL */
/*<       M=LC1+L-1 >*/
    m = lc1 + l - 1;
/*   IS THIS OVER END OF CONTIG? */
/*<       IF(M.GT.IDC)M=IDC >*/
    if (m > *idc) {
	m = *idc;
    }
/*   NUMBER TO COPY */
/*<       M=M-LC1+1 >*/
    m = m - lc1 + 1;
/*   CHECK FOR OVERFLOW */
/*<       IF(IS1+M-1.GT.MAXGEL)GO TO 700 >*/
    if (is1 + m - 1 > *maxgel) {
	goto L700;
    }
/*<       IF(M.GT.0)CALL SQCOPY(SEQC(LC1),SEQC2(IS1),M) >*/
    if (m > 0) {
	sqcopy_(seqc + lc1, seqc2 + is1, &m, (ftnlen)1, (ftnlen)1);
    }
/*< 600   CONTINUE >*/
L600:
/*   COUNT PADDING IN GEL */
/*<       ITOTPG=IS2+L-1-IDG >*/
    *itotpg = is2 + l - 1 - *idg;
/*   SET NEW LENGTHS FOR RETURN TO CALLING ROUTINE */
/*<       IDOUT=IS1+M-1 >*/
    *idout = is1 + m - 1;
/*<       IDG=IS2+L-1 >*/
    *idg = is2 + l - 1;
/*<       IFAIL=0 >*/
    *ifail = 0;
/*<       RETURN >*/
    return 0;
/*< 700   CONTINUE >*/
L700:
/*<    >*/
    erromf_("Matching region too long in lineup: alignment aborted", (ftnlen)
	    53);
/*<       IFAIL=1 >*/
    *ifail = 1;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* lineup_ */

/*     MERGE */

/*   ROUTINE SENT CONTIG WHOSE GELS MAY BE OUT OF ORDER */
/*   REORDERS GELS ON POSITION OF LEFT ENDS AND SETS LEFT */
/*   GEL NUMBER FOR THE REORDERED CONTIG */

/*<       SUBROUTINE MERGE(RELPG,LNGTHG,LNBR,RNBR,LINCON,IDBSIZ) >*/
/* Subroutine */ int merge_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *lincon, integer *idbsiz)
{
    static integer m, n, i1, i2, nr;

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER RELPG(IDBSIZ) >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/

/*   START AT LEFT END */
/*<       N=LNBR(LINCON) >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    n = lnbr[*lincon];
/*<       GO TO 22 >*/
    goto L22;
/*< 21    CONTINUE >*/
L21:
/*   SET POINTER TO NEXT GEL TO RIGHT IN LIST */
/*<       N=NR >*/
    n = nr;
/*<       IF(I1.GT.0)N=I2 >*/
    if (i1 > 0) {
	n = i2;
    }
/*< 22    CONTINUE >*/
L22:
/*   SET POINTER TO NEXT GEL TO RIGHT */
/*<       NR=RNBR(N) >*/
    nr = rnbr[n];
/*<       IF(NR.EQ.0)GO TO 30 >*/
    if (nr == 0) {
	goto L30;
    }
/*   HAVENT REACHED END YET */
/*<       I1=0 >*/
    i1 = 0;
/*< 23    CONTINUE >*/
L23:
/*   ARE THESE 2 IN CORRECT ORDER IE N<=NR ? */
/*<       IF(RELPG(N).LE.RELPG(NR))GO TO 21 >*/
    if (relpg[n] <= relpg[nr]) {
	goto L21;
    }
/*   NOT IN ORDER SO CHAIN LEFT UNTIL CORRECTLY POSITIONED */
/*   THEN COME BACK TO THIS POINT AND CONTINUE */
/*   IF FIRST MOVE SAVE POSITION */
/*<       IF(I1.EQ.0)I2=N >*/
    if (i1 == 0) {
	i2 = n;
    }
/*<       I1=1 >*/
    i1 = 1;
/*   EXCHANGE NEIGHBOURS */
/*<       M=RNBR(NR) >*/
    m = rnbr[nr];
/*<       IF(M.NE.0)LNBR(M)=N >*/
    if (m != 0) {
	lnbr[m] = n;
    }
/*<       M=LNBR(N) >*/
    m = lnbr[n];
/*<       IF(M.NE.0)RNBR(M)=NR >*/
    if (m != 0) {
	rnbr[m] = nr;
    }
/*<       RNBR(N)=RNBR(NR) >*/
    rnbr[n] = rnbr[nr];
/*<       RNBR(NR)=N >*/
    rnbr[nr] = n;
/*<       LNBR(NR)=LNBR(N) >*/
    lnbr[nr] = lnbr[n];
/*<       LNBR(N)=NR >*/
    lnbr[n] = nr;
/*   CHAIN BACK THRU LIST */
/*<       N=LNBR(NR) >*/
    n = lnbr[nr];
/*<       IF(N.EQ.0)GO TO 21 >*/
    if (n == 0) {
	goto L21;
    }
/*   END NOT REACHED */
/*<       GO TO 23 >*/
    goto L23;
/*< 30    CONTINUE >*/
L30:
/*  ALL DONE POINTER AT RIGHT GEL */
/*<       RNBR(LINCON)=N >*/
    rnbr[*lincon] = n;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* merge_ */

/*<       SUBROUTINE REMOVL(MATC,MATG,MATL,IP) >*/
/* Subroutine */ int removl_(integer *matc, integer *matg, integer *matl, 
	integer *ip)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k1, k2, k3, k4, k5, k6, ipp, idelt, nmtch;
    extern /* Subroutine */ int bubbl3_(integer *, integer *, integer *, 
	    integer *);

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER MATC(IP),MATG(IP),MATL(IP) >*/

/*   SET POINTER TO FIRST MATCH */
/*<       NMTCH=0 >*/
    /* Parameter adjustments */
    --matl;
    --matg;
    --matc;

    /* Function Body */
    nmtch = 0;
/*< 10    CONTINUE >*/
L10:
/*   POINT TO NEXT MATCH */
/*<       NMTCH=NMTCH+1 >*/
    ++nmtch;
/*   SORT MATCHES ON LENGTH */
/*<       IPP=IP-NMTCH+1 >*/
    ipp = *ip - nmtch + 1;
/*<       CALL BUBBL3(MATL(NMTCH),MATG(NMTCH),MATC(NMTCH),IPP) >*/
    bubbl3_(&matl[nmtch], &matg[nmtch], &matc[nmtch], &ipp);
/*   LOOK FOR END OF POSITIVES */
/*<       DO 20 I=NMTCH,IP >*/
    i__1 = *ip;
    for (i__ = nmtch; i__ <= i__1; ++i__) {
/*<       J=I >*/
	j = i__;
/*< 20    IF(MATL(I).LT.1)GO TO 30 >*/
/* L20: */
	if (matl[i__] < 1) {
	    goto L30;
	}
    }
/*<       J=J+1 >*/
    ++j;
/*< 30    CONTINUE >*/
L30:
/*<       IP=J-1 >*/
    *ip = j - 1;
/*   END OF POSITIVES AT IP */
/*<       IF(NMTCH.GE.IP)RETURN >*/
    if (nmtch >= *ip) {
	return 0;
    }
/*<       K1=MATC(NMTCH) >*/
    k1 = matc[nmtch];
/*<       K2=K1+MATL(NMTCH)-1 >*/
    k2 = k1 + matl[nmtch] - 1;
/*<       K3=MATG(NMTCH) >*/
    k3 = matg[nmtch];
/*<       K4=K3+MATL(NMTCH)-1 >*/
    k4 = k3 + matl[nmtch] - 1;
/*   POINT TO FIRST MATCH TO TEST */
/*<       K6=NMTCH+1 >*/
    k6 = nmtch + 1;
/*<       DO 200 I=K6,IP >*/
    i__1 = *ip;
    for (i__ = k6; i__ <= i__1; ++i__) {
/*   DO CONSENSUS FIRST */
/*   OVERLAP? */
/*<       IF(MATC(I).GT.K2)GO TO 100 >*/
	if (matc[i__] > k2) {
	    goto L100;
	}
/*<       K5=MATC(I)+MATL(I)-1 >*/
	k5 = matc[i__] + matl[i__] - 1;
/*<       IF(K5.LT.K1)GO TO 100 >*/
	if (k5 < k1) {
	    goto L100;
	}
/*   DOES OVERLAP */
/*   WHICH END */
/*<       IF(K5.LE.K2)GO TO 80 >*/
	if (k5 <= k2) {
	    goto L80;
	}
/*   LENGTH TO REDUCE MATCH BY IS IDELT */
/*<       IDELT=K2-MATC(I)+1 >*/
	idelt = k2 - matc[i__] + 1;
/*   NEW LENGTH */
/*<       MATL(I)=MATL(I)-IDELT >*/
	matl[i__] -= idelt;
/*  MOVE LEFT ENDS */
/*<       MATC(I)=MATC(I)+IDELT >*/
	matc[i__] += idelt;
/*<       MATG(I)=MATG(I)+IDELT >*/
	matg[i__] += idelt;
/*<       GO TO 100 >*/
	goto L100;
/*< 80    CONTINUE >*/
L80:
/*   LENGTH */
/*<       MATL(I)=K1-MATC(I) >*/
	matl[i__] = k1 - matc[i__];
/*< 100   CONTINUE >*/
L100:
/*   NOW LOOK FOR OVERLAPS WITH GEL */
/*   OVERLAP? */
/*<       IF(MATG(I).GT.K4)GO TO 200 >*/
	if (matg[i__] > k4) {
	    goto L200;
	}
/*<       K5=MATG(I)+MATL(I)-1 >*/
	k5 = matg[i__] + matl[i__] - 1;
/*<       IF(K5.LT.K3)GO TO 200 >*/
	if (k5 < k3) {
	    goto L200;
	}
/*   DOES OVERLAP */
/*   WHICH END? */
/*<       IF(K5.LE.K4)GO TO 180 >*/
	if (k5 <= k4) {
	    goto L180;
	}
/*   LENGTH TO REDUCE MATCH BY IS IDELT */
/*<       IDELT=K4-MATG(I)+1 >*/
	idelt = k4 - matg[i__] + 1;
/*   NEW LENGTH */
/*<       MATL(I)=MATL(I)-IDELT >*/
	matl[i__] -= idelt;
/*   MOVE LEFT ENDS */
/*<       MATC(I)=MATC(I)+IDELT >*/
	matc[i__] += idelt;
/*<       MATG(I)=MATG(I)+IDELT >*/
	matg[i__] += idelt;
/*<       GO TO 200 >*/
	goto L200;
/*< 180   CONTINUE >*/
L180:
/*   LENGTH */
/*<       MATL(I)=K3-MATG(I) >*/
	matl[i__] = k3 - matg[i__];
/*< 200   CONTINUE >*/
L200:
	;
    }
/*<       GO TO 10 >*/
    goto L10;
/*<       END >*/
} /* removl_ */

/*<       SUBROUTINE TPCHEK(PC,PG,L,N) >*/
/* Subroutine */ int tpchek_(integer *pc, integer *pg, integer *l, integer *n)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j1, k1;
    extern /* Subroutine */ int ml_(integer *, integer *, integer *, integer *
	    , integer *);

/*<       INTEGER PC(N),PG(N),L(N) >*/
/*     AUTHOR RODGER STADEN */
/*     IF OVERLAPPING BLOCKS ARE FOUND REMOVE THE SHORTER ONE */
/*     THEN REMOVE LARGE GAPS AT ENDS (THOSE AS LARGE AS THE END BLOCK) */
/*<       K1 = 2 >*/
    /* Parameter adjustments */
    --l;
    --pg;
    --pc;

    /* Function Body */
    k1 = 2;
/*< 1     CONTINUE >*/
L1:
/*<       DO 10 I = K1,N >*/
    i__1 = *n;
    for (i__ = k1; i__ <= i__1; ++i__) {
/*<         J1 = I >*/
	j1 = i__;
/*<         IF(PC(I).LE.PC(I-1)) GO TO 20 >*/
	if (pc[i__] <= pc[i__ - 1]) {
	    goto L20;
	}
/*<         IF(PG(I).LE.PG(I-1)) GO TO 20 >*/
	if (pg[i__] <= pg[i__ - 1]) {
	    goto L20;
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*     REMOVE LARGE GAPS FROM ENDS */
/*     THIS RULE OF THUMB COULD BE CHANGED TO USE A DIFFERENCE */
/*     BETWEEN THE NUMBERS OF MISMATCHING CHARACTERS */
/*<       IF(N.GT.1) THEN >*/
    if (*n > 1) {
/*<         K1 = PC(2) - PC(1) - L(1)  >*/
	k1 = pc[2] - pc[1] - l[1];
/*<         J1 = PG(2) - PG(1) - L(1) >*/
	j1 = pg[2] - pg[1] - l[1];
/*<         IF(MAX(K1,J1).GT.L(1)) THEN >*/
	if (max(k1,j1) > l[1]) {
/*<           CALL ML(PC,PG,L,N,1) >*/
	    ml_(&pc[1], &pg[1], &l[1], n, &c__1);
/*<           N = N - 1 >*/
	    --(*n);
/*<         END IF >*/
	}
/*<         IF(N.GT.1) THEN >*/
	if (*n > 1) {
/*<           K1 = PC(N) - PC(N-1) - L(N-1) >*/
	    k1 = pc[*n] - pc[*n - 1] - l[*n - 1];
/*<           J1 = PG(N) - PG(N-1) - L(N-1) >*/
	    j1 = pg[*n] - pg[*n - 1] - l[*n - 1];
/*<           IF(MAX(K1,J1).GT.L(N)) THEN >*/
	    if (max(k1,j1) > l[*n]) {
/*<             CALL ML(PC,PG,L,N,N) >*/
		ml_(&pc[1], &pg[1], &l[1], n, n);
/*<             N = N - 1 >*/
		--(*n);
/*<           END IF >*/
	    }
/*<         END IF >*/
	}
/*<       END IF >*/
    }
/*<       RETURN >*/
    return 0;
/*< 20    CONTINUE >*/
L20:
/*<       IF(L(J1-1).GT.L(J1)) THEN >*/
    if (l[j1 - 1] > l[j1]) {
/*<         CALL ML(PC,PG,L,N,J1) >*/
	ml_(&pc[1], &pg[1], &l[1], n, &j1);
/*<       ELSE >*/
    } else {
/*<         CALL ML(PC,PG,L,N,J1-1) >*/
	i__1 = j1 - 1;
	ml_(&pc[1], &pg[1], &l[1], n, &i__1);
/*<       END IF >*/
    }
/*  Until 25-11-90 next line was k1=j1 but this does not deal with all */
/*  cases: when a line is deleted we must compare it with the previous */
/*  one before dealing with the rest, because it could be left of that */
/*   one as well! */
/*<       K1 = MAX(2,J1-1) >*/
/* Computing MAX */
    i__1 = 2, i__2 = j1 - 1;
    k1 = max(i__1,i__2);
/*<       N = N - 1 >*/
    --(*n);
/*<       GO TO 1 >*/
    goto L1;
/*<       END >*/
} /* tpchek_ */

/*<    >*/
/* Subroutine */ int updcon_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *idbsiz, integer *ngels, integer *nconts, char 
	*seq, integer *maxseq, integer *idim1, integer *cstart, integer *
	cleno, integer *lincon, char *nampro, char *seq2, integer *idevr, 
	integer *ifail, integer *maxgel, integer *idm, real *percd, integer *
	mask, integer *clist, ftnlen seq_len, ftnlen nampro_len, ftnlen 
	seq2_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer l, b1, s1, ld, nbad, lreg, rreg, iladd[1], iradd[1], igelc;
    extern integer chnrp_(integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *);
    static integer itask, iwing;
    extern /* Subroutine */ int precn1_(char *, char *, real *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, ftnlen, ftnlen), makhca_(char *, integer *, integer *, 
	    integer *, integer *, integer *, ftnlen), addtit_(char *, char *, 
	    integer *, integer *, ftnlen, ftnlen), erromf_(char *, ftnlen);

/*<       INTEGER RELPG(IDBSIZ),LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       INTEGER CSTART,CLENO,S1,B1,RREG >*/
/*<       CHARACTER SEQ(MAXSEQ),SEQ2(MAXGEL) >*/
/*<       CHARACTER NAMPRO*(*) >*/
/*<       INTEGER ILADD(1),IRADD(1) >*/
/*<       INTEGER CLIST(1) >*/
/*<       INTEGER CHNRP >*/
/*<       EXTERNAL CHNRP >*/
/* cstart consensus start point (before new reading) */
/* cleno consensus length (before new reading) */
/* lincon element number of contig */
/* s1 number of first reading to shift */
/* b1 number of first base to shift (in overall consensus positioning) */

/* there are 2 tasks: 1. make space for the new and altered region */
/*                    2. calculate the new consensus and put it in the space */
/* we do not have to make space if: */
/* a. we are dealing with the last contig in the consensus and there are no */
/*    readings starting to the right of the new data */
/* b. the contig has not been padded */

/* New code to update the consensus only for the region affected by the */
/* new reading. Find the next reading to the right of the new one, which */
/* the new one does not overlap (might not be one!). Make a consensus from */
/* start of new reading to here. Prior to this make space for it by moving */
/* the consensus right (only if the contig is longer (padding or extra data */
/* at its ends). Let s1 be the first reading to shift. We shift from its */
/* left end to the end of the contig - where is this in the overall consensus? */
/* The distance of the left end of s1 to the right end of the contig is */
/* unchanged. This means that the new relpg(s1) is the same distance from */
/* the right end of the old consensus as the old relpg(s1) was from the right */
/* end of the old consensus. So from this we can calculate the position of the */
/* the first base to move. */
/* Let L be the position in the overall consensus of  the last base in this contig */
/*            L = cstart - cleno - 1 */
/* Let D = distance to end of contig */
/*            D = RELPG(LINCON) - relpg(s1) + 1. */
/* First base to shift B1 = L - D + 1 */
/* Last base to shift is idim1 */
/* Distance to move to right is relpg(lincon) - cleno ie the number of extra bases */
/* make consensus from relpg(ngels) to relpg(s1) - 1 */
/* put it at cstart + relpg(ngels) - 1 */

/* Potential problems: */
/* 1) reading at right end of contig */
/* the search for the first nonoverlapping read to the right will return 0 */
/* shift al the next contig: ie cstart + cleno onwards */
/* make consensus from relpg(ngels) to end of contig */
/* put it at cstart + relpg(ngels) -1 */

/* 2) reading at left end of contig */
/* shift whole contig ie cstart - 20 */
/* add new title */
/* shift consensus relpg(lincon) - cleno to the right */

/* 3) new reading contains contig - cases 1 and 2 combined */
/* the search for the first nonoverlapping read to the right will return 0 */
/* shift whole of next contig and make consensus from relpg(ngels) to end of */
/* contig. */

/* 4) Might not be a next contig to shift */

/* get number of first reading to shift */

/*<    >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --seq;
    --seq2;
    --clist;

    /* Function Body */
    i__2 = relpg[*ngels] + (i__1 = lngthg[*ngels], abs(i__1)) - 1;
    s1 = chnrp_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
	    nconts, &i__2);
/*      WRITE(*,*)'S1',S1 */

/* is the altered region longer than the original: only then do we need to shift */

/*           WRITE(*,*)'IDIM1',IDIM1 */
/*           WRITE(*,*)'RELPG(LINCON)',RELPG(LINCON) */
/*           WRITE(*,*)'CSTART,CLENO',CSTART,CLENO */
/*<       IF (RELPG(LINCON) - CLENO.GT.0) THEN >*/
    if (relpg[*lincon] - *cleno > 0) {

/* it is longer so we probably need to shift */

/*<         IF (S1.EQ.0) THEN >*/
	if (s1 == 0) {

/* no readings start to the right of the new data */

/*<           IF (CSTART+CLENO-1.LT.IDIM1) THEN >*/
	    if (*cstart + *cleno - 1 < *idim1) {

/* there are other contigs to the right */

/*           WRITE(*,*)'CSTART,CLENO',CSTART,CLENO */
/*<             B1 = CSTART + CLENO >*/
		b1 = *cstart + *cleno;
/*            WRITE(*,*)'B1',B1 */
/*<             CALL MAKHCA(SEQ,MAXSEQ,B1,RELPG(LINCON)-CLENO,IDIM1,IFAIL) >*/
		i__1 = relpg[*lincon] - *cleno;
		makhca_(seq + 1, maxseq, &b1, &i__1, idim1, ifail, (ftnlen)1);
/*<             IF(IFAIL.NE.0) THEN >*/
		if (*ifail != 0) {
/*<               CALL ERROMF('Error: consensus too long') >*/
		    erromf_("Error: consensus too long", (ftnlen)25);
/*<               RETURN >*/
		    return 0;
/*<             END IF >*/
		}
/*<           ELSE >*/
	    } else {

/* there are no contigs to the right and no readings start to the right of */
/* the new one so nothing to shift */

/*<           END IF >*/
	    }
/*<         ELSE >*/
	} else {

/* there are readings starting to the right of the new one */

/* shift from start of next reading to right */

/*<            L = CSTART + CLENO - 1 >*/
	    l = *cstart + *cleno - 1;
/*           WRITE(*,*)'CSTART,CLENO,L',CSTART,CLENO,L */
/*<            LD = RELPG(LINCON) - RELPG(S1) + 1 >*/
	    ld = relpg[*lincon] - relpg[s1] + 1;
/*           WRITE(*,*)'LD',LD */
/*<            B1 = L - LD + 1 >*/
	    b1 = l - ld + 1;
/*            WRITE(*,*)'B1',B1 */
/*<            CALL MAKHCA(SEQ,MAXSEQ,B1,RELPG(LINCON)-CLENO,IDIM1,IFAIL) >*/
	    i__1 = relpg[*lincon] - *cleno;
	    makhca_(seq + 1, maxseq, &b1, &i__1, idim1, ifail, (ftnlen)1);
/*<            IF(IFAIL.NE.0) THEN >*/
	    if (*ifail != 0) {
/*<              CALL ERROMF('Error: consensus too long') >*/
		erromf_("Error: consensus too long", (ftnlen)25);
/*<              RETURN >*/
		return 0;
/*<            END IF >*/
	    }
/*<         END IF >*/
	}
/*<       END IF >*/
    }

/* now make new consensus (where do we put it,  do we need */
/* to give it a header, and what region do we make it for ? */
/* in the simplest case make it for relpg(ngels) to relpg(s1) -1 */
/* if s1=0 make it for relpg(ngels) to end of contig (relpg(lincon)) */
/* we give it a header if it is at the left end of the contig ie lnbr(ngels)=0 */

/* we always start at the left end of the new reading */

/*<       LREG = RELPG(NGELS) >*/
    lreg = relpg[*ngels];

/* we end at the next reading to the right or the end of the contig */

/*<       IF (S1.NE.0) THEN >*/
    if (s1 != 0) {
/*<         RREG = RELPG(S1) - 1 >*/
	rreg = relpg[s1] - 1;
/*<       ELSE >*/
    } else {
/*<         RREG = RELPG(LINCON) >*/
	rreg = relpg[*lincon];
/*<       END IF >*/
    }

/* where do we put the new consensus ? */

/*<       B1 = CSTART + RELPG(NGELS) - 1 >*/
    b1 = *cstart + relpg[*ngels] - 1;
/*      WRITE(*,*)'LREG,RREG',LREG,RREG */
/*            WRITE(*,*)'B1',B1 */

/* do we need to add a title */

/*<       IF (LNBR(NGELS).EQ.0) THEN >*/
    if (lnbr[*ngels] == 0) {
/*<         B1 = CSTART - 20 >*/
	b1 = *cstart - 20;
/*        WRITE(*,*)'ADD NEW TIT AT',B1 */
/*<         CALL ADDTIT(SEQ(B1),NAMPRO,NGELS,B1) >*/
	addtit_(seq + b1, nampro, ngels, &b1, (ftnlen)1, nampro_len);
/*<       END IF >*/
    }
/*<       IGELC = LNBR(LINCON) >*/
    igelc = lnbr[*lincon];

/* note aconsn will chain along until it find the first useful reading!! */

/*      JOB = 2 */

/* set dummy values for precon (and iladd,iradd above) */

/*<       NBAD = 0 >*/
    nbad = 0;
/*<       IWING = 0 >*/
    iwing = 0;

/* set task (normal consensus) */

/*<       ITASK = 4 >*/
    itask = 4;

/* add masking if required */

/*<       IF (MASK.EQ.3) ITASK = ITASK + 32 >*/
    if (*mask == 3) {
	itask += 32;
    }
/*<       IF (MASK.EQ.4) ITASK = 8 >*/
    if (*mask == 4) {
	itask = 8;
    }
/*      WRITE(*,*)'BEFORE' */
/*      WRITE(*,*)(SEQ(JJJ),JJJ=1,IDIM1+RELPG(LINCON)-CLENO) */
/*         CALL FMTDB1(SEQ,IDIM1,1,IDIM1,60,6) */
/*      WRITE(*,*)'NOCONT,LREG,RREG,ITASK,B1' */
/*      WRITE(*,*)NOCONT,LREG,RREG,ITASK,B1 */
/*<    >*/
    precn1_(seq + 1, nampro, percd, idbsiz, lincon, &lreg, &rreg, &itask, 
	    idevr, &b1, maxgel, maxseq, &iwing, &nbad, iladd, iradd, ifail, (
	    ftnlen)1, nampro_len);
/*<       IF(IFAIL.NE.0) THEN >*/
    if (*ifail != 0) {
/*<         CALL ERROMF('Error calculating consensus') >*/
	erromf_("Error calculating consensus", (ftnlen)27);
/*<         RETURN >*/
	return 0;
/*<       END IF >*/
    }

/* before we leave we must make the overall consensus length correct */
/*  so add on the extra length (if any) which is the new length - old length */

/*      WRITE(*,*)'OLD IDIM1',IDIM1 */
/*<       IDIM1 = IDIM1 + RELPG(LINCON) - CLENO >*/
    *idim1 = *idim1 + relpg[*lincon] - *cleno;
/*      WRITE(*,*)'after NEW IDIM1/2',IDIM1 */
/*      WRITE(*,*)(SEQ(JJJ),JJJ=1,IDIM1) */
/*         CALL FMTDB1(SEQ,IDIM1,1,IDIM1,60,6) */
/*<       END >*/
    return 0;
} /* updcon_ */

/*<       SUBROUTINE MAKHCA(STRING,MAXAR,FROM,HSIZE,ASIZE,IFAIL) >*/
/* Subroutine */ int makhca_(char *string, integer *maxar, integer *from, 
	integer *hsize, integer *asize, integer *ifail, ftnlen string_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;

/*<       CHARACTER STRING(MAXAR) >*/
/*<       INTEGER FROM,HSIZE,ASIZE >*/

/* make a hole of size hsize in character array size asize */

/*<       J = ASIZE + HSIZE >*/
    /* Parameter adjustments */
    --string;

    /* Function Body */
    j = *asize + *hsize;
/*<       IF (J.GT.MAXAR) THEN >*/
    if (j > *maxar) {
/*<         IFAIL = 1 >*/
	*ifail = 1;
/*<         RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       DO 10 I=ASIZE,FROM,-1 >*/
    i__1 = *from;
    for (i__ = *asize; i__ >= i__1; --i__) {
/*<         STRING(J) = STRING(I) >*/
	*(unsigned char *)&string[j] = *(unsigned char *)&string[i__];
/*<         J = J - 1 >*/
	--j;
/*<  10     CONTINUE >*/
/* L10: */
    }
/*<       IFAIL = 0 >*/
    *ifail = 0;
/*<       END >*/
    return 0;
} /* makhca_ */

/*<    >*/
integer chnrp_(integer *relpg, integer *lngthg, integer *lnbr, integer *rnbr, 
	integer *idbsiz, integer *lgel, integer *ncont, integer *lreg)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;

/*<       INTEGER RELPG(IDBSIZ),LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/

/* find first reading starting past lreg (0=none found) */

/*<       I = LGEL >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *lgel;
/*<       CHNRP = 0 >*/
    ret_val = 0;
/*< 10    CONTINUE >*/
L10:
/*<       IF(I.NE.0) THEN >*/
    if (i__ != 0) {
/*<         IF(RELPG(I).LE.LREG) THEN >*/
	if (relpg[i__] <= *lreg) {
/*<           I = RNBR(I) >*/
	    i__ = rnbr[i__];
/*<           GO TO 10 >*/
	    goto L10;
/*<         END IF >*/
	}
/*<         CHNRP = I >*/
	ret_val = i__;
/*<         RETURN >*/
	return ret_val;
/*<       END IF >*/
    }
/*<       END >*/
    return ret_val;
} /* chnrp_ */

/*<    >*/
integer chnrp1_(integer *relpg, integer *lngthg, integer *lnbr, integer *rnbr,
	 integer *idbsiz, integer *lgel, integer *lreg)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;

/*<       INTEGER RELPG(IDBSIZ),LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/

/* find first reading with data covering or past lreg (0=none found) */

/*<       I = LGEL >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *lgel;
/*<       CHNRP1 = 0 >*/
    ret_val = 0;
/*< 10    CONTINUE >*/
L10:
/*<       IF(I.NE.0) THEN >*/
    if (i__ != 0) {
/*<         IF(RELPG(I)+ABS(LNGTHG(I))-1.LT.LREG) THEN >*/
	if (relpg[i__] + (i__1 = lngthg[i__], abs(i__1)) - 1 < *lreg) {
/*<           I = RNBR(I) >*/
	    i__ = rnbr[i__];
/*<           GO TO 10 >*/
	    goto L10;
/*<         END IF >*/
	}
/*<         CHNRP1 = I >*/
	ret_val = i__;
/*<         RETURN >*/
	return ret_val;
/*<       END IF >*/
    }
/*<       END >*/
    return ret_val;
} /* chnrp1_ */

/*<       SUBROUTINE AERROR(LIST,NAME,IERR) >*/
/* Subroutine */ int aerror_(char *list, char *name__, integer *ierr, ftnlen 
	list_len, ftnlen name_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer j, l;
    extern /* Subroutine */ int info_(char *, ftnlen);
    static char infod[60];
    extern /* Subroutine */ int swrt2b_(char *, char *, ...),
	swrt3b_(char *, char *, ...), erromf_(char *, ftnlen);
    static char errmsg[333];
    extern /* Subroutine */ int tolist_(char *, char *, ftnlen, ftnlen);

/*<       CHARACTER LIST*(*),NAME*(*) >*/
/*<       CHARACTER INFOD*60 >*/
/*<       CHARACTER ERRMSG*333 >*/

/* handle errors for assembly */

/* errors are: */
/* 0 file not found or file is of invalid format */
/* 1 read too short */
/* 2 failed to align and not entered */
/* 3 failed on entry */
/* 4 failed to align but entered */
/* 5 no match found during masked assembly */
/*<       L=1 >*/
    l = 1;
/*<       DO 5 J=1,LEN(NAME) >*/
    i__1 = i_len(name__, name_len);
    for (j = 1; j <= i__1; ++j) {
/*<          L=J >*/
	l = j;
/*<          IF (NAME(J:J).EQ.' ') THEN >*/
	if (*(unsigned char *)&name__[j - 1] == ' ') {
/*<             GO TO 6 >*/
	    goto L6;
/*<          END IF >*/
	}
/*<  5    CONTINUE >*/
/* L5: */
    }
/*<  6    CONTINUE >*/
L6:
/*      WRITE(INFOD,1000)NAME(1:L),IERR */
/* 1000 FORMAT(A,I2) */
/* CHECKED */
/*<       CALL SWRT3B(INFOD,'%.*s%2d%!',LEN(NAME(1:L)),NAME(1:L),IERR) >*/
    i__1 = i_len(name__, l);
    swrt3b_(infod, "%.*s%2d%!", &i__1, name__, ierr, (ftnlen)60, (ftnlen)9, l)
	    ;
/*      WRITE(ERRMSG,1010)'Failed file ',NAME(1:L), */
/*     +     'written to error file' */
/* 1010 FORMAT(A,A,A) */
/* CHECKED */
/*<    >*/
    i__1 = i_len(name__, l);
    swrt2b_(errmsg, "Failed file %.*swritten to error file%!", &i__1, name__, 
	    (ftnlen)333, (ftnlen)39, l);
/*<       CALL ERROMF(ERRMSG) >*/
    erromf_(errmsg, (ftnlen)333);
/*<       CALL TOLIST(LIST,INFOD) >*/
    tolist_(list, infod, list_len, (ftnlen)60);
/*<       CALL INFO(ERRMSG) >*/
    info_(errmsg, (ftnlen)333);
/*<       END >*/
    return 0;
} /* aerror_ */

/*<       SUBROUTINE SINDB(IDEVN,NGELS,RNAMES,NAME,JOB) >*/
/* Subroutine */ int sindb_(integer *idevn, integer *ngels, char *rnames, 
	char *name__, integer *job, ftnlen rnames_len, ftnlen name_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int readn_(integer *, integer *, char *, ftnlen);

/*<       CHARACTER*(*) RNAMES(NGELS),NAME*(*) >*/
/*<       IF (JOB.EQ.1) THEN >*/
    /* Parameter adjustments */
    rnames -= rnames_len;

    /* Function Body */
    if (*job == 1) {
/*<         DO 10 J=1,NGELS >*/
	i__1 = *ngels;
	for (j = 1; j <= i__1; ++j) {
/*<           CALL READN(IDEVN,J,RNAMES(J)) >*/
	    readn_(idevn, &j, rnames + j * rnames_len, rnames_len);
/*          WRITE(*,*)'INITIALISING RNAMES ',RNAMES(J) */
/*< 10      CONTINUE >*/
/* L10: */
	}
/*<         RETURN >*/
	return 0;
/*<       ELSE IF (JOB.EQ.2) THEN >*/
    } else if (*job == 2) {
/*<         RNAMES(NGELS) = NAME >*/
	s_copy(rnames + *ngels * rnames_len, name__, rnames_len, name_len);
/*          WRITE(*,*)' ADDING TO RNAMES ',RNAMES(NGELS) */
/*<       END IF >*/
    }
/*<       END >*/
    return 0;
} /* sindb_ */

/*<       INTEGER FUNCTION INDB(NGELS,RNAMES,NAME) >*/
integer indb_(integer *ngels, char *rnames, char *name__, ftnlen rnames_len, 
	ftnlen name_len)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer j;

/*<       CHARACTER RNAMES(NGELS)*40,NAME*(*) >*/
/*<       DO 10 J=1,NGELS >*/
    /* Parameter adjustments */
    rnames -= 40;

    /* Function Body */
    i__1 = *ngels;
    for (j = 1; j <= i__1; ++j) {
/*        WRITE(*,*)'CHECKING RNAMES ',NAME,' ',RNAMES(J) */
/*<         IF(NAME.EQ.RNAMES(J)) THEN >*/
	if (s_cmp(name__, rnames + j * 40, name_len, (ftnlen)40) == 0) {
/*<           INDB = J >*/
	    ret_val = j;
/*<           RETURN >*/
	    return ret_val;
/*<         END IF >*/
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       INDB = 0 >*/
    ret_val = 0;
/*<       END >*/
    return ret_val;
} /* indb_ */

/*<    >*/
/* Subroutine */ int slides_(char *seq1, integer *idc, char *seq2, integer *
	idim2, integer *ms1, integer *ms2, integer *maxpg, integer *maxpc, 
	integer *minsli, integer *matl, integer *matc, integer *matg, integer 
	*ip, ftnlen seq1_len, ftnlen seq2_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, n, p1, p2, ip1, p1s;
    extern /* Subroutine */ int savit_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    extern integer ctonum_(char *, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQ1(IDC),SEQ2(IDIM2) >*/
/*<       INTEGER MATL(IP),MATC(IP),MATG(IP),P1S,P1,P2 >*/
/*<       INTEGER CTONUM >*/
/*<       EXTERNAL CTONUM >*/
/*<       IP1 = IP >*/
    /* Parameter adjustments */
    --seq1;
    --seq2;
    --matg;
    --matc;
    --matl;

    /* Function Body */
    ip1 = *ip;
/*<       IP = 0 >*/
    *ip = 0;
/*   LEFT END S2 RELATIVE S1 - MAX PADS -2 READY FOR LOOP */
/*<       P1S = MS1 - MS2 - MAXPC - 1 >*/
    p1s = *ms1 - *ms2 - *maxpc - 1;
/*   TRY NSLIDE START POSNS FOR SEQ2 */
/*      WRITE(*,*)'IDC,IDIM2',IDC,IDIM2 */
/*      WRITE(*,*)'SEQ1' */
/*      WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
/*      WRITE(*,*)'SEQ2' */
/*      WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2) */
/*<       DO 100 I=1,MAXPG+MAXPC+1 >*/
    i__1 = *maxpg + *maxpc + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       POINT TO SEQ1 START */
/*<         P1S = P1S + 1 >*/
	++p1s;
/*       POINT TO CURRENT SEQ1 POSN */
/*<         P1 = P1S >*/
	p1 = p1s;
/*<         N = 0 >*/
	n = 0;
/*       COMPARE WHOLE LENGTH OF SEQ2 (IF P1 WITHIN RANGE) */
/*<         DO 50 J=1,IDIM2 >*/
	i__2 = *idim2;
	for (j = 1; j <= i__2; ++j) {
/*<           P2 = J >*/
	    p2 = j;
/*<           P1 = P1 + 1 >*/
	    ++p1;
/*<           IF(P1.LT.1)GO TO 50 >*/
	    if (p1 < 1) {
		goto L50;
	    }
/*         OFF RIGHT END? IF SO MAY HAVE BEEN A MATCH */
/*<           IF(P1.GT.IDC)GO TO 40 >*/
	    if (p1 > *idc) {
		goto L40;
	    }
/*<           IF(CTONUM(SEQ1(P1)).EQ.CTONUM(SEQ2(P2)))GO TO 45 >*/
	    if (ctonum_(seq1 + p1, (ftnlen)1) == ctonum_(seq2 + p2, (ftnlen)1)
		    ) {
		goto L45;
	    }
/*< 40        CONTINUE >*/
L40:
/*<           IF(N.GE.MINSLI)CALL SAVIT(N,P1,P2,IP,MATL,MATC,MATG,IP1) >*/
	    if (n >= *minsli) {
		savit_(&n, &p1, &p2, ip, &matl[1], &matc[1], &matg[1], &ip1);
	    }
/*<           N = 0 >*/
	    n = 0;
/*<           GO TO 50 >*/
	    goto L50;
/*< 45        CONTINUE >*/
L45:
/*<           N = N + 1 >*/
	    ++n;
/*< 50      CONTINUE >*/
L50:
	    ;
	}
/*       GOOD SCORE AT END? NEED TO INCREMENT POINTERS FOR SAVIT */
/*<         P1 = P1 + 1 >*/
	++p1;
/*<         P2 = P2 + 1 >*/
	++p2;
/*<         IF(N.GE.MINSLI)CALL SAVIT(N,P1,P2,IP,MATL,MATC,MATG,IP1) >*/
	if (n >= *minsli) {
	    savit_(&n, &p1, &p2, ip, &matl[1], &matc[1], &matg[1], &ip1);
	}
/*< 100   CONTINUE >*/
/* L100: */
    }
/*<       END >*/
    return 0;
} /* slides_ */

/*<    >*/
/* Subroutine */ int padcon_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *ngels, integer *nconts, char *gel, integer *
	lincon, integer *posn, integer *nc, integer *idbsiz, integer *idevr, 
	integer *maxgel, ftnlen gel_len)
{
    /* Initialized data */

    static char pad[1] = "*";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_sign(integer *, integer *);

    /* Local variables */
    static integer i__, k, x;
    extern /* Subroutine */ int readw_(integer *, integer *, char *, integer *
	    , ftnlen);
    static integer llino;
    extern /* Subroutine */ int insbas_(integer *, integer *, integer *, char 
	    *, ftnlen), writec_(integer *, integer *, integer *, integer *, 
	    integer *), writeg_(integer *, integer *, integer *, integer *, 
	    integer *, integer *), shiftt_(integer *, integer *, integer *, 
	    integer *);

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER RELPG(IDBSIZ),POSN,X >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       CHARACTER GEL(MAXGEL) >*/
/*<       CHARACTER PAD >*/
/*<       SAVE PAD >*/
/*<       DATA PAD/'*'/ >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */
/*   NOW FIND FIRST CHAR THAT OVERLAPS REGION */
/*<       LLINO=LNBR(LINCON) >*/
    llino = lnbr[(0 + (0 + (*lincon << 2))) / 4];
/*< 30    CONTINUE >*/
L30:
/*<       X=RELPG(LLINO)+ABS(LNGTHG(LLINO))-1 >*/
    x = relpg[llino] + (i__1 = lngthg[llino], abs(i__1)) - 1;
/*<       IF(X.GE.POSN)GO TO 40 >*/
    if (x >= *posn) {
	goto L40;
    }
/*   NOT IN REGION */
/*<       LLINO=RNBR(LLINO) >*/
    llino = rnbr[llino];
/*<       GO TO 30 >*/
    goto L30;
/*< 40    CONTINUE >*/
L40:
/*   NOW GET THIS GEL FROM DISK */
/*<       CALL READW(IDEVR,LLINO,GEL,MAXGEL) >*/
    readw_(idevr, &llino, gel + 1, maxgel, (ftnlen)1);
/*   CALC POSN IN THIS GEL TO EDIT */
/*<       X=POSN-RELPG(LLINO)+1 >*/
    x = *posn - relpg[llino] + 1;
/*<       K=X >*/
    k = x;
/*<       LNGTHG(LLINO)=LNGTHG(LLINO)+SIGN(NC,LNGTHG(LLINO)) >*/
    lngthg[llino] += i_sign(nc, &lngthg[llino]);
/*<    >*/
    if ((i__1 = lngthg[llino], abs(i__1)) > *maxgel) {
	lngthg[llino] = i_sign(maxgel, &lngthg[llino]);
    }
/*   INSBAS TAKES CARE OF OVERFLOW OF READ IN ITS OWN WAY */
/*<       DO 61 I=1,NC >*/
    i__1 = *nc;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        WRITE(*,*)'INSERT TO ',LLINO,' AT ',K */
/*<         CALL INSBAS(IDEVR,LLINO,K,PAD) >*/
	insbas_(idevr, &llino, &k, pad, (ftnlen)1);
/*<  61   CONTINUE >*/
/* L61: */
    }
/*   WRITE NEW LINE */
/*<    >*/
    writeg_(idevr, &llino, &relpg[llino], &lngthg[llino], &lnbr[llino], &rnbr[
	    llino]);
/*< 65    CONTINUE >*/
L65:
/*   NOW GET NEXT GEL */
/*<       LLINO=RNBR(LLINO) >*/
    llino = rnbr[llino];
/*   LAST GEL? */
/*<       IF(LLINO.EQ.0)GO TO 70 >*/
    if (llino == 0) {
	goto L70;
    }
/*   DOES IT HAVE DATA IN REGION? */
/*   IE DO RELPG  AND RELPG+LNGTHG-1 LIE EITHER SIDE OF POSN? */
/*<       IF(RELPG(LLINO).GE.POSN)GO TO 70 >*/
    if (relpg[llino] >= *posn) {
	goto L70;
    }
/*      WRITE(*,*)LLINO,RELPG(LLINO),POSN */
/*<       X=RELPG(LLINO)+ABS(LNGTHG(LLINO))-1 >*/
    x = relpg[llino] + (i__1 = lngthg[llino], abs(i__1)) - 1;
/*<       IF(X.LT.POSN)GO TO 65 >*/
    if (x < *posn) {
	goto L65;
    }
/*  WITHIN */
/*<       GO TO 40 >*/
    goto L40;
/*< 70    CONTINUE >*/
L70:
/*      WRITE(*,*)'NOW MOVING' */
/*   INSERTS FINISHED SO NEED TO INCREMENT ALL THOSE GELS TO RIGHT */
/*<       LLINO=LNBR(LINCON) >*/
    llino = lnbr[*lincon];
/*< 75    CONTINUE >*/
L75:
/*<       IF(RELPG(LLINO).GE.POSN)GO TO 80 >*/
    if (relpg[llino] >= *posn) {
	goto L80;
    }
/*< 76    CONTINUE >*/
L76:
/*      WRITE(*,*)'SKIPPING ',LLINO */
/*<       LLINO=RNBR(LLINO) >*/
    llino = rnbr[llino];
/*<       IF(LLINO.EQ.0)GO TO 90 >*/
    if (llino == 0) {
	goto L90;
    }
/*<       GO TO 75 >*/
    goto L75;
/*< 80    CONTINUE >*/
L80:
/*      WRITE(*,*)'ADDING TO ',LLINO */
/*<       RELPG(LLINO)=RELPG(LLINO)+NC >*/
    relpg[llino] += *nc;
/*   WRITE NEW LINE */
/*<    >*/
    writeg_(idevr, &llino, &relpg[llino], &lngthg[llino], &lnbr[llino], &rnbr[
	    llino]);
/*<       GO TO 76 >*/
    goto L76;
/*< 90    CONTINUE >*/
L90:
/*   NEED TO INCREMENT CONTIG LINE */
/*<       RELPG(LINCON)=RELPG(LINCON)+NC >*/
    relpg[*lincon] += *nc;
/*<    >*/
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
/*   Now move tags along on the consensus */
/*<       CALL SHIFTT(IDEVR, IDBSIZ-LINCON, POSN, NC) >*/
    i__1 = *idbsiz - *lincon;
    shiftt_(idevr, &i__1, posn, nc);
/*<       END >*/
    return 0;
} /* padcon_ */

/*<       SUBROUTINE UPCHEK(PC,PG,L,N) >*/
/* Subroutine */ int upchek_(integer *pc, integer *pg, integer *l, integer *n)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j1, k1, dc, dg;
    extern /* Subroutine */ int ml_(integer *, integer *, integer *, integer *
	    , integer *);

/*<       INTEGER PC(N),PG(N),L(N),DC,DG >*/
/*     AUTHOR RODGER STADEN */

/* only allow gaps that are shorter than the next block of identity */

/*<       K1 = 2 >*/
    /* Parameter adjustments */
    --l;
    --pg;
    --pc;

    /* Function Body */
    k1 = 2;
/*< 1     CONTINUE >*/
L1:
/*<       DO 10 I = K1,N >*/
    i__1 = *n;
    for (i__ = k1; i__ <= i__1; ++i__) {
/*<         J1 = I >*/
	j1 = i__;
/*<         DC = PC(I) - PC(I-1) - L(I-1) >*/
	dc = pc[i__] - pc[i__ - 1] - l[i__ - 1];
/*<         DG = PG(I) - PG(I-1) - L(I-1) >*/
	dg = pg[i__] - pg[i__ - 1] - l[i__ - 1];
/*<         IF(ABS(DC-DG).GE.L(I)) GO TO 20 >*/
	if ((i__2 = dc - dg, abs(i__2)) >= l[i__]) {
	    goto L20;
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;
/*< 20    CONTINUE >*/
L20:
/*      WRITE(*,*)'REMOVING!!' */
/*<         CALL ML(PC,PG,L,N,J1) >*/
    ml_(&pc[1], &pg[1], &l[1], n, &j1);
/*      IF(L(J1-1).GT.L(J1)) THEN */
/*        CALL ML(PC,PG,L,N,J1) */
/*      ELSE */
/*        CALL ML(PC,PG,L,N,J1-1) */
/*      END IF */
/*<       K1 = MAX(2,J1-1) >*/
/* Computing MAX */
    i__1 = 2, i__2 = j1 - 1;
    k1 = max(i__1,i__2);
/*<       N = N - 1 >*/
    --(*n);
/*<       GO TO 1 >*/
    goto L1;
/*<       END >*/
} /* upchek_ */

/*<    >*/
/* Subroutine */ int padcop_(char *seqg, char *seqg2, integer *lg1, integer *
	mg, integer *l5, integer *is2, integer *lg2, integer *maxgel, integer 
	*ifail, char *seqc, integer *idc, integer *ic1, ftnlen seqg_len, 
	ftnlen seqg2_len, ftnlen seqc_len)
{
    /* Initialized data */

    static char dubbl[1*4] = "D" "B" "V" "H";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, m, jc1, mgm1;
    extern /* Subroutine */ int info_(char *, ftnlen);
    static integer idone, maxreq;

/*   AUTHOR: RODGER STADEN */
/*<       PARAMETER (NDUBL = 4) >*/
/*<       CHARACTER SEQG(MAXGEL),SEQG2(MAXGEL),DUBBL(NDUBL),SEQC(IDC) >*/
/*<       SAVE DUBBL >*/
/*<       DATA DUBBL/'D','B','V','H'/ >*/
    /* Parameter adjustments */
    --seqg2;
    --seqg;
    --seqc;

    /* Function Body */
/*<       JC1 = IC1 >*/
    jc1 = *ic1;
/* Make seqg2 from seqg placing L5 padding chars before position MG */
/* which is the start of the next block of identity. Try to put the */
/* padding either in line with consensus pads, or next to double */
/* codes. The positions in seqg are LG1 to MG-1. seqg2 needs to be long */
/* enough to be extended from IS2 to IS2 + L5 -1 + MGM1-LG1 +1 */
/* ie we add L5 pads, plus the chars between and including  LG1 and MGM1 */
/*<       IDONE=0 >*/
    idone = 0;
/*   POINT TO END OF MISMATCH */
/*<       MGM1=MG-1 >*/
    mgm1 = *mg - 1;
/*   MAY BE NO CHARS TO COPY */
/*<       IF(MGM1.LT.LG1)GO TO 111 >*/
    if (mgm1 < *lg1) {
	goto L111;
    }
/*  Next check added 26-2-91 */
/*<       MAXREQ = IS2 + L5 - 1 + MGM1 - LG1 + 1 >*/
    maxreq = *is2 + *l5 - 1 + mgm1 - *lg1 + 1;
/*<       IF((MGM1.GT.MAXGEL).OR.(MAXREQ.GT.MAXGEL)) THEN >*/
    if (mgm1 > *maxgel || maxreq > *maxgel) {
/*<    >*/
	info_("Matching region too large in padcop: alignment aborted", (
		ftnlen)54);
/*<         IFAIL=1 >*/
	*ifail = 1;
/*<         RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       DO 110 J=LG1,MGM1 >*/
    i__1 = mgm1;
    for (j = *lg1; j <= i__1; ++j) {
/*<         IF(IDONE.LT.L5) THEN >*/
	if (idone < *l5) {
/*<           IF((JC1.GT.0).AND.(JC1.LT.IDC)) THEN >*/
	    if (jc1 > 0 && jc1 < *idc) {
/*<           IF(SEQC(JC1).EQ.'*') THEN >*/
		if (*(unsigned char *)&seqc[jc1] == '*') {
/*<             IS2 = IS2 + 1 >*/
		    ++(*is2);
/*<             JC1 = JC1 + 1 >*/
		    ++jc1;
/*<             IDONE = IDONE + 1 >*/
		    ++idone;
/*<             GO TO 109 >*/
		    goto L109;
/*<           END IF >*/
		}
/*<           END IF >*/
	    }
/*<           DO 108 M=1,NDUBL >*/
	    for (m = 1; m <= 4; ++m) {
/*<             IF(SEQG(J).EQ.DUBBL(M)) THEN >*/
		if (*(unsigned char *)&seqg[j] == *(unsigned char *)&dubbl[m 
			- 1]) {
/*<               IS2 = IS2 + 1 >*/
		    ++(*is2);
/*<               JC1 = JC1 + 1 >*/
		    ++jc1;
/*<               IDONE = IDONE + 1 >*/
		    ++idone;
/*<               GO TO 109 >*/
		    goto L109;
/*<             END IF >*/
		}
/*< 108       CONTINUE >*/
/* L108: */
	    }
/*< 109       CONTINUE >*/
L109:
/*<         END IF >*/
	    ;
	}
/*<         SEQG2(IS2) = SEQG(J) >*/
	*(unsigned char *)&seqg2[*is2] = *(unsigned char *)&seqg[j];
/*<         IS2 = IS2 + 1 >*/
	++(*is2);
/*<         JC1 = JC1 + 1 >*/
	++jc1;
/*< 110   CONTINUE >*/
/* L110: */
    }
/*< 111   CONTINUE >*/
L111:
/*   ALL CHARS COPIED. ENOUGH PADDING? */
/*<       IF(IDONE.LT.L5)IS2=IS2+L5-IDONE >*/
    if (idone < *l5) {
	*is2 = *is2 + *l5 - idone;
    }
/*   IS2 SHOULD NOW BE POINTING AT NEXT CHAR */
/*   ZERO LG2 TO SHOW CALLING ROUTINE COPYING DONE */
/*<       LG2=0 >*/
    *lg2 = 0;
/*<       IFAIL=0 >*/
    *ifail = 0;
/*<       END >*/
    return 0;
} /* padcop_ */

/*<       SUBROUTINE ML(PC,PG,L,N,J) >*/
/* Subroutine */ int ml_(integer *pc, integer *pg, integer *l, integer *n, 
	integer *j)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*<       INTEGER PC(N),PG(N),L(N) >*/
/*<       DO 10 I = J,N-1 >*/
    /* Parameter adjustments */
    --l;
    --pg;
    --pc;

    /* Function Body */
    i__1 = *n - 1;
    for (i__ = *j; i__ <= i__1; ++i__) {
/*<         PC(I) = PC(I+1) >*/
	pc[i__] = pc[i__ + 1];
/*<         PG(I) = PG(I+1) >*/
	pg[i__] = pg[i__ + 1];
/*<         L(I) = L(I+1) >*/
	l[i__] = l[i__ + 1];
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       END >*/
    return 0;
} /* ml_ */

/*<       SUBROUTINE MSTLKL(SEQ,IDIM) >*/
/* Subroutine */ int mstlkl_(char *seq, integer *idim, ftnlen seq_len)
{
    /* System generated locals */
    integer i__1;
    char ch__1[1];

    /* Local variables */
    static integer i__, j, k;
    extern /* Character */ VOID charsu_(char *, ftnlen, integer *);
    extern integer indexs_(char *, integer *, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQ(IDIM) >*/
/*<       CHARACTER CHARSU >*/
/*<       EXTERNAL CHARSU,INDEXS >*/
/*<       DO 100 I=1,IDIM >*/
    /* Parameter adjustments */
    --seq;

    /* Function Body */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         J = INDEXS(SEQ(I),K) >*/
	j = indexs_(seq + i__, &k, (ftnlen)1);
/*<         SEQ(I) = CHARSU(J) >*/
	charsu_(ch__1, (ftnlen)1, &j);
	*(unsigned char *)&seq[i__] = *(unsigned char *)&ch__1[0];
/*< 100   CONTINUE >*/
/* L100: */
    }
/*<       END >*/
    return 0;
} /* mstlkl_ */

/*<       INTEGER FUNCTION INDEXS(C,S) >*/
integer indexs_(char *c__, integer *s, ftnlen c_len)
{
    /* Initialized data */

    static integer ind[29] = { 1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,6,6,6,6,6,6,1,
	    2,3,4,5,5,6 };
    static integer scores[29] = { 100,100,100,100,75,75,75,75,100,100,100,100,
	    100,100,100,100,10,10,10,10,10,10,100,100,100,100,100,100,10 };

    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;

/*<       PARAMETER (IDM = 29) >*/
/*<       CHARACTER C >*/
/*<       INTEGER POINTS(0:255),SCORES(IDM),IND(IDM),S >*/
/*<       COMMON /SHOTC/POINTS >*/
/*<       SAVE /SHOTC/ >*/
/*<       SAVE SCORES,IND >*/
/*<    >*/
/*      DATA DUP/'CTAG1234DVBHKLMNRY5678ctag*,-'/ */
/*  changed 28-7-91 to give 10 to old zeroes and 100 to lowercase */
/*<    >*/
/*<       I = ICHAR(C) >*/
    i__ = *(unsigned char *)c__;
/*<       I = POINTS(I) >*/
    i__ = shotc_1.points[i__];
/*<       S = SCORES(I) >*/
    *s = scores[i__ - 1];
/*<       INDEXS = IND(I) >*/
    ret_val = ind[i__ - 1];
/*<       END >*/
    return ret_val;
} /* indexs_ */

/*<       CHARACTER*1 FUNCTION CHARSU(I) >*/
/* Character */ VOID charsu_(char *ret_val, ftnlen ret_val_len, integer *i__)
{
    /* Initialized data */

    static char c__[6] = "CTAG*-";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

/*<       CHARACTER C*6 >*/
/*<       SAVE C >*/
/*<       DATA C/'CTAG*-'/ >*/
/*<       CHARSU = C(I:I) >*/
    s_copy(ret_val, c__ + (0 + (0 + (*i__ - 1))), (ftnlen)1, *i__ - (*i__ - 1)
	    );
/*<       END >*/
} /* charsu_ */

/*     SAVIT */

/*<       SUBROUTINE SAVIT(N,J,K,IP,S1,S2,S3,IP1) >*/
/* Subroutine */ int savit_(integer *n, integer *j, integer *k, integer *ip, 
	integer *s1, integer *s2, integer *s3, integer *ip1)
{
/*   AUTHOR: RODGER STADEN */
/*<       INTEGER S1(IP1),S2(IP1),S3(IP1) >*/

/*<       IP=IP+1 >*/
    /* Parameter adjustments */
    --s3;
    --s2;
    --s1;

    /* Function Body */
    ++(*ip);
/*   TEST FOR OVERFLOW */
/*<       IF(IP.GT.IP1)RETURN >*/
    if (*ip > *ip1) {
	return 0;
    }
/*<       S1(IP)=N >*/
    s1[*ip] = *n;
/*<       S2(IP)=J-N >*/
    s2[*ip] = *j - *n;
/*<       S3(IP)=K-N >*/
    s3[*ip] = *k - *n;

/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* savit_ */

/*<    >*/
integer clen_(integer *relpg, integer *lngthg, integer *lnbr, integer *rnbr, 
	integer *ngels, integer *nconts, integer *idbsiz, integer *iin)
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, len;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LEFT GEL NUMBER OR ZERO FOR ERROR */
/*<       INTEGER RELPG(IDBSIZ) >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       I = IIN >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *iin;
/*<       CLEN= 0 >*/
    ret_val = 0;
/*<       LEN = 0 >*/
    len = 0;
/*< 10    CONTINUE >*/
L10:
/*<       IF(I.NE.0)THEN >*/
    if (i__ != 0) {
/*<         LEN = MAX(LEN,(RELPG(I) + ABS(LNGTHG(I)) - 1)) >*/
/* Computing MAX */
	i__2 = len, i__3 = relpg[i__] + (i__1 = lngthg[i__], abs(i__1)) - 1;
	len = max(i__2,i__3);
/*<         I = RNBR(I) >*/
	i__ = rnbr[i__];
/*<         IF(I.EQ.IIN)RETURN >*/
	if (i__ == *iin) {
	    return ret_val;
	}
/*<         GO TO 10 >*/
	goto L10;
/*<       END IF >*/
    }
/*<       CLEN = LEN >*/
    ret_val = len;
/*<       END >*/
    return ret_val;
} /* clen_ */

/*<    >*/
/* Subroutine */ int gllino_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *idbsiz, integer *nconts, integer *llino, 
	integer *lincon)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n, mxt;

/*<       INTEGER RELPG(IDBSIZ),LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/

/* routine to get the left gel number and contig line of the longest contig */

/*<       LLINO = 0 >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    *llino = 0;
/*<       LINCON = 0 >*/
    *lincon = 0;
/*<       N=IDBSIZ-NCONTS >*/
    n = *idbsiz - *nconts;
/*<       MXT = 0 >*/
    mxt = 0;
/*<       DO 4 I=N,IDBSIZ-1 >*/
    i__1 = *idbsiz - 1;
    for (i__ = n; i__ <= i__1; ++i__) {
/*<         IF(RELPG(I).GT.MXT) THEN >*/
	if (relpg[i__] > mxt) {
/*<           MXT = RELPG(I) >*/
	    mxt = relpg[i__];
/*<           LLINO = LNBR(I) >*/
	    *llino = lnbr[i__];
/*<           LINCON = I >*/
	    *lincon = i__;
/*<         END IF >*/
	}
/*< 4     CONTINUE >*/
/* L4: */
    }
/*<       END >*/
    return 0;
} /* gllino_ */

/*<       INTEGER FUNCTION CLINNO(LNBR,IDBSIZ,NCONTS,IIN) >*/
integer clinno_(integer *lnbr, integer *idbsiz, integer *nconts, integer *iin)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer j, n;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LINE NUMBER OR ZERO FOR ERROR */
/*<       INTEGER LNBR(IDBSIZ) >*/
/*<       CLINNO = 0 >*/
    /* Parameter adjustments */
    --lnbr;

    /* Function Body */
    ret_val = 0;
/*<       N=IDBSIZ-NCONTS >*/
    n = *idbsiz - *nconts;
/*<       DO 10 J=N,IDBSIZ-1 >*/
    i__1 = *idbsiz - 1;
    for (j = n; j <= i__1; ++j) {
/*<         IF(LNBR(J).EQ.IIN) THEN >*/
	if (lnbr[j] == *iin) {
/*<           CLINNO = J >*/
	    ret_val = j;
/*<           RETURN >*/
	    return ret_val;
/*<         END IF >*/
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       END >*/
    return ret_val;
} /* clinno_ */

/*<       SUBROUTINE SHFTLA(STRING,MAXAR,FROMS,TO,FROME) >*/
/* Subroutine */ int shftla_(char *string, integer *maxar, integer *froms, 
	integer *to, integer *frome, ftnlen string_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;

/*<       CHARACTER STRING(MAXAR) >*/
/*<       INTEGER FROMS,TO,FROME >*/

/* shift an array left from froms to to */

/*<       J = TO >*/
    /* Parameter adjustments */
    --string;

    /* Function Body */
    j = *to;
/*<       DO 10 I=FROMS,FROME >*/
    i__1 = *frome;
    for (i__ = *froms; i__ <= i__1; ++i__) {
/*<         STRING(J) = STRING(I) >*/
	*(unsigned char *)&string[j] = *(unsigned char *)&string[i__];
/*<         J = J + 1 >*/
	++j;
/*<  10   CONTINUE >*/
/* L10: */
    }
/*<       END >*/
    return 0;
} /* shftla_ */

/*<    >*/
integer randc_(integer *relpg, integer *lngthg, integer *lnbr, integer *rnbr, 
	integer *ngels, integer *nconts, integer *idbsiz, integer *iin, 
	integer *lincon, integer *llino)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;
    extern integer gclin_(integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *), chainl_(integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *)
	    ;

/*<       INTEGER RELPG(IDBSIZ),LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       INTEGER CHAINL,GCLIN >*/
/*<       EXTERNAL CHAINL,GCLIN >*/

/* return reading and contig number for contig containing IIN */
/* -1 = ERROR, 0 = OK */

/*<       I = CHAINL(RELPG,LNGTHG,LNBR,RNBR,NGELS,NCONTS,IDBSIZ,IIN) >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = chainl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, iin);
/*<       IF (I.EQ.0) THEN >*/
    if (i__ == 0) {
/*<         RANDC = -1 >*/
	ret_val = -1;
/*<         RETURN >*/
	return ret_val;
/*<       END IF >*/
    }
/*<       LLINO = I >*/
    *llino = i__;
/*<       I = GCLIN(RELPG,LNGTHG,LNBR,RNBR,NGELS,NCONTS,IDBSIZ,LLINO) >*/
    i__ = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, llino);
/*<       IF (I.EQ.0) THEN >*/
    if (i__ == 0) {
/*<         RANDC = -2 >*/
	ret_val = -2;
/*<         RETURN >*/
	return ret_val;
/*<       END IF >*/
    }
/*<       LINCON = I >*/
    *lincon = i__;
/*<       RANDC = 0 >*/
    ret_val = 0;
/*<       END >*/
    return ret_val;
} /* randc_ */

/*<       SUBROUTINE INITS >*/
/* Subroutine */ int inits_(void)
{
    /* Initialized data */

    static char dup[29] = "CTAG1234DVBHKLMNRY5678ctag*,-";

    static integer i__, j;

/*  AUTHOR RODGER STADEN */
/*<       INTEGER POINTS(0:255) >*/
/*<       PARAMETER (IDM = 29) >*/
/*<       CHARACTER DUP*29 >*/
/*<       COMMON /SHOTC/POINTS >*/
/*<       SAVE /SHOTC/ >*/
/*<       DATA DUP/'CTAG1234DVBHKLMNRY5678ctag*,-'/ >*/
/*  ICHAR RETURNS THE COLLATING SEQUENCE NUMBER */
/*  I WANT 1-4 FOR ACGT */
/*                 acgt */
/*                 1234 */
/*                 BDHV */
/*                 KLMN */
/*      5 FOR      * */
/*      6 FOR      5678- AND ELSE */
/*  THE ACTUAL VALUE RETURNED BY ICHAR IS NOT PORTABLE */
/*  SO I NEED TO INITIALIZE POINTR SO THAT THE CORRECT */
/*  ELEMENTS CONTAIN VALUES 1 - 6 */

/*<         DO 30 I = 0,255 >*/
    for (i__ = 0; i__ <= 255; ++i__) {
/*<           POINTS(I) = IDM >*/
	shotc_1.points[i__] = 29;
/*< 30      CONTINUE >*/
/* L30: */
    }
/*<         DO 35 I = 1,IDM >*/
    for (i__ = 1; i__ <= 29; ++i__) {
/*<           J = ICHAR(DUP(I:I)) >*/
	j = *(unsigned char *)&dup[i__ - 1];
/*<           POINTS(J) = I >*/
	shotc_1.points[j] = i__;
/*< 35      CONTINUE >*/
/* L35: */
    }
/*<       END >*/
    return 0;
} /* inits_ */

/*<    >*/
integer chainl_(integer *relpg, integer *lngthg, integer *lnbr, integer *rnbr,
	 integer *ngels, integer *nconts, integer *idbsiz, integer *iin)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__, j;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LEFT GEL NUMBER OR ZERO FOR ERROR */
/*<       INTEGER RELPG(IDBSIZ) >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       I = IIN >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *iin;
/*<       J = I >*/
    j = i__;
/*<       CHAINL = 0 >*/
    ret_val = 0;
/*< 10    CONTINUE >*/
L10:
/*<       IF(I.NE.0)THEN >*/
    if (i__ != 0) {
/*<         J = I >*/
	j = i__;
/*<         I = LNBR(I) >*/
	i__ = lnbr[i__];
/*<         IF(I.EQ.IIN)RETURN >*/
	if (i__ == *iin) {
	    return ret_val;
	}
/*<         GO TO 10 >*/
	goto L10;
/*<       END IF >*/
    }
/*<       CHAINL = J >*/
    ret_val = j;
/*<       END >*/
    return ret_val;
} /* chainl_ */

/*<    >*/
/* Subroutine */ int shiftc_(integer *relpg, integer *lngthg, integer *lnbr, 
	integer *rnbr, integer *ngels, integer *nconts, integer *idevr, 
	integer *idbsiz, integer *ign, integer *ncont, integer *dist)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, l;
    extern integer clen_(integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *);
    extern /* Subroutine */ int writec_(integer *, integer *, integer *, 
	    integer *, integer *), writeg_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);

/*  AUTHOR: RODGER STADEN */
/*  SHIFTS PART OF A CONTIG FORM GEL IGN TO RIGHT END */
/*  CONTIG LINE NUMBER IF NCONT */
/*<       INTEGER RELPG(IDBSIZ) >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       INTEGER DIST,CLEN >*/
/*<       EXTERNAL CLEN >*/
/*<       I = IGN >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *ign;
/*< 10    CONTINUE >*/
L10:
/*<       IF(I.NE.0)THEN >*/
    if (i__ != 0) {
/*<         RELPG(I) = RELPG(I) + DIST >*/
	relpg[i__] += *dist;
/*<         CALL WRITEG(IDEVR,I,RELPG(I),LNGTHG(I),LNBR(I),RNBR(I)) >*/
	writeg_(idevr, &i__, &relpg[i__], &lngthg[i__], &lnbr[i__], &rnbr[i__]
		);
/*<         I = RNBR(I) >*/
	i__ = rnbr[i__];
/*<         GO TO 10 >*/
	goto L10;
/*<       END IF >*/
    }
/*  UPDATE CONTIG LENGTH */
/*<    >*/
    l = clen_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, ign);
/*<       RELPG(NCONT) = L >*/
    relpg[*ncont] = l;
/*<    >*/
    i__1 = *idbsiz - *ncont;
    writec_(idevr, &i__1, &relpg[*ncont], &lnbr[*ncont], &rnbr[*ncont]);
/*<       END >*/
    return 0;
} /* shiftc_ */

/*<       SUBROUTINE FNDCON(SEQ,IDIM,CENDS,NENDS,IDCEND,MAXCON) >*/
/* Subroutine */ int fndcon_(char *seq, integer *idim, integer *cends, 
	integer *nends, integer *idcend, integer *maxcon, ftnlen seq_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static char dc[1*8];
    extern integer indexa_(char *, integer *, char *, ftnlen, ftnlen), 
	    jfromc_(char *, integer *, ftnlen);
    extern /* Subroutine */ int erromf_(char *, ftnlen);

/*   AUTHOR: RODGER STADEN */
/*   STORES THEIR POSITIONS IN CENDS AND THEIR LEFT LINE NUMBERS IN NENDS */
/*<       PARAMETER (MAXDG = 8) >*/
/*<       CHARACTER SEQ(IDIM),DC(MAXDG) >*/
/*<       INTEGER CENDS(MAXCON) >*/
/*<       INTEGER NENDS(MAXCON) >*/
/*<       EXTERNAL JFROMC,INDEXA >*/
/*<       IDCEND=0 >*/
    /* Parameter adjustments */
    --seq;
    --nends;
    --cends;

    /* Function Body */
    *idcend = 0;
/*<       DO 10 I=1,IDIM >*/
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         IF(SEQ(I).NE.'<')GO TO 10 >*/
	if (*(unsigned char *)&seq[i__] != '<') {
	    goto L10;
	}
/*<         IDCEND=IDCEND+1 >*/
	++(*idcend);
/*       PUT POSITION OF LEFT END OF CONTIG IN CENDS */
/*<         CENDS(IDCEND)=I >*/
	cends[*idcend] = i__;
/*<         K = INDEXA(SEQ(I),20,'.') >*/
	k = indexa_(seq + i__, &c__20, ".", (ftnlen)1, (ftnlen)1);
/*<         IF(K.EQ.0) GO TO 20 >*/
	if (k == 0) {
	    goto L20;
	}
/*<         K = K + I >*/
	k += i__;
/*        IF (.NOT.((SEQ(K+MAXDG).EQ.'-').OR.(SEQ(K+MAXDG).EQ.'>'))) */
/*     +     GOTO20 */
/*<         DO 5 J=1,MAXDG >*/
	for (j = 1; j <= 8; ++j) {
/*<           IF ((SEQ(K).EQ.'-').OR.(SEQ(K).EQ.'>')) THEN >*/
	    if (*(unsigned char *)&seq[k] == '-' || *(unsigned char *)&seq[k] 
		    == '>') {
/*<              GO TO 6 >*/
		goto L6;
/*<           END IF >*/
	    }
/*<           DC(J)=SEQ(K) >*/
	    *(unsigned char *)&dc[j - 1] = *(unsigned char *)&seq[k];
/*<           K=K+1 >*/
	    ++k;
/*< 5       CONTINUE >*/
/* L5: */
	}
/*< 6       CONTINUE >*/
L6:
/*<         NENDS(IDCEND)=JFROMC(DC,J-1) >*/
	i__2 = j - 1;
	nends[*idcend] = jfromc_(dc, &i__2, (ftnlen)1);
/*< 10    CONTINUE >*/
L10:
	;
    }
/*     STORE POSITION OF LAST CHAR +1 TO SIMPLIFY DISPLAY ROUTINES */
/*<       CENDS(IDCEND+1)=IDIM+1 >*/
    cends[*idcend + 1] = *idim + 1;
/*<       RETURN >*/
    return 0;
/*<  20   CONTINUE >*/
L20:
/*<       CALL ERROMF('Error in FNDCON: illegal consensus header') >*/
    erromf_("Error in FNDCON: illegal consensus header", (ftnlen)41);
/*<       IDCEND = 0 >*/
    *idcend = 0;
/*<       END >*/
    return 0;
} /* fndcon_ */

/*<    >*/
integer gclin_(integer *relpg, integer *lngthg, integer *lnbr, integer *rnbr, 
	integer *ngels, integer *nconts, integer *idbsiz, integer *iin)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer j, n;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LINE NUMBER OR ZERO FOR ERROR */
/*<       INTEGER RELPG(IDBSIZ) >*/
/*<       INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) >*/
/*<       GCLIN = 0 >*/
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    ret_val = 0;
/*<       N=IDBSIZ-NCONTS >*/
    n = *idbsiz - *nconts;
/*<       DO 10 J=N,IDBSIZ-1 >*/
    i__1 = *idbsiz - 1;
    for (j = n; j <= i__1; ++j) {
/*<         IF(LNBR(J).EQ.IIN) THEN >*/
	if (lnbr[j] == *iin) {
/*<           GCLIN = J >*/
	    ret_val = j;
/*<           RETURN >*/
	    return ret_val;
/*<         END IF >*/
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       END >*/
    return ret_val;
} /* gclin_ */

/*     SQCOM */
/*<       SUBROUTINE SQCOMM(SEQ,IDIM) >*/
/* Subroutine */ int sqcomm_(char *seq, integer *idim, ftnlen seq_len)
{
    /* Initialized data */

    static char list1[1*12] = "C" "T" "A" "G" "c" "t" "a" "g" "e" "d" "f" 
	    "i";
    static char list2[1*12] = "G" "A" "T" "C" "g" "a" "t" "c" "i" "f" "d" 
	    "e";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static char temp[1];

/*   AUTHOR: RODGER STADEN */
/*<       PARAMETER (MAXLST = 12) >*/
/*<       CHARACTER SEQ(IDIM),LIST1(MAXLST),LIST2(MAXLST),TEMP >*/
/*<       SAVE LIST1,LIST2 >*/
/*<    >*/
    /* Parameter adjustments */
    --seq;

    /* Function Body */
/*<    >*/
/*<       DO 100 I=1,IDIM >*/
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         TEMP = SEQ(I) >*/
	*(unsigned char *)temp = *(unsigned char *)&seq[i__];
/*<         DO 50 J=1,MAXLST >*/
	for (j = 1; j <= 12; ++j) {
/*<           IF(TEMP.EQ.LIST1(J))THEN >*/
	    if (*(unsigned char *)temp == *(unsigned char *)&list1[j - 1]) {
/*<             SEQ(I)=LIST2(J) >*/
		*(unsigned char *)&seq[i__] = *(unsigned char *)&list2[j - 1];
/*<             GO TO 99 >*/
		goto L99;
/*<           END IF >*/
	    }
/*< 50      CONTINUE >*/
/* L50: */
	}
/*< 99      CONTINUE >*/
L99:
/*< 100   CONTINUE >*/
/* L100: */
	;
    }
/*<       END >*/
    return 0;
} /* sqcomm_ */

/*<       SUBROUTINE ERROMF(STRING) >*/
/* Subroutine */ int erromf_(char *string, ftnlen string_len)
{
    extern /* Subroutine */ int fverr_(integer *, char *, char *, ftnlen, 
	    ftnlen);

/*<       CHARACTER STRING*(*) >*/
/*<       CALL FVERR(0, 'Error', STRING) >*/
    fverr_(&c__0, "Error", string, (ftnlen)5, string_len);
/*<       END >*/
    return 0;
} /* erromf_ */

/*<       SUBROUTINE BUSYF() >*/
/* Subroutine */ int busyf_(void)
{
/*      WRITE(*,*)'BUSY' */
/*<       END >*/
    return 0;
} /* busyf_ */

/*  ROUTINES TO CONTROL CHARACTER LOOKUP */
/*  FOR BOTH DNA AND PROTEIN SEQUENCES */
/*  THE INITIALISING ROUTINES ARE SENT THE CHARACTERSET SIZE IDM */
/*  WHICH DETERMINES WHICH CHARACTERSET IS USED */
/*<       SUBROUTINE INITLU(IDM) >*/
/* Subroutine */ int initlu_(integer *idm)
{
    /* Initialized data */

    static char dup[16] = "TCAG-RYWSMKHBVDN";
    static char pup[26] = "CSTPAGNDEQBZHRKMILVFYW-X? ";
    static char dlow[16] = "tcag-rywsmkhbvdn";
    static char plow[26] = "cstpagndeqbzhrkmilvfyw-x? ";

    static integer i__, j;

/*  AUTHOR RODGER STADEN */
/*<       INTEGER POINT1(0:255),POINT2(0:255) >*/
/*<       CHARACTER DUP*16,DLOW*16,PUP*26,PLOW*26 >*/
/*<       COMMON /IASCI1/POINT1 >*/
/*<       COMMON /IASCI2/POINT2 >*/
/*<       SAVE /IASCI1/ >*/
/*<       SAVE /IASCI2/ >*/
/*<       SAVE DUP,PUP,DLOW,PLOW >*/
/*<       DATA DUP/'TCAG-RYWSMKHBVDN'/ >*/
/*<       DATA PUP/'CSTPAGNDEQBZHRKMILVFYW-X? '/ >*/
/*<       DATA DLOW/'tcag-rywsmkhbvdn'/ >*/
/*<       DATA PLOW/'cstpagndeqbzhrkmilvfyw-x? '/ >*/
/*  ICHAR RETURNS THE COLLATING SEQUENCE NUMBER */
/*  I WANT 1-5 FOR ACGT OR 1-26 FOR AMINO ACIDS BY USING ICHAR. */
/*  THE ACTUAL VALUE RETURNED BY ICHAR IS NOT PORTABLE */
/*  SO I NEED TO INITIALIZE POINTR SO THAT THE CORRECT */
/*  ELEMENTS CONTAIN VALUES 1 - 5, OR 1 - 26 */
/*  WORKS ON UPPER AND LOWER CASE - REMOVE DLOW,PLOW AND LOOPS 41 AND 51 */
/*  IF LOWERCASE NOT ALLOWED */

/*<       IF(IDM.EQ.5)THEN >*/
    if (*idm == 5) {
/*<         DO 30 I = 0,255 >*/
	for (i__ = 0; i__ <= 255; ++i__) {
/*<           POINT1(I) = IDM >*/
	    iasci1_1.point1[i__] = *idm;
/*<           POINT2(I) = 17 >*/
	    iasci2_1.point2[i__] = 17;
/*< 30      CONTINUE >*/
/* L30: */
	}
/*<         DO 35 I = 1,5 >*/
	for (i__ = 1; i__ <= 5; ++i__) {
/*<           J = ICHAR(DUP(I:I)) >*/
	    j = *(unsigned char *)&dup[i__ - 1];
/*<           POINT1(J) = I >*/
	    iasci1_1.point1[j] = i__;
/*< 35      CONTINUE >*/
/* L35: */
	}
/*<         DO 36 I = 1,5 >*/
	for (i__ = 1; i__ <= 5; ++i__) {
/*<           J = ICHAR(DLOW(I:I)) >*/
	    j = *(unsigned char *)&dlow[i__ - 1];
/*<           POINT1(J) = I >*/
	    iasci1_1.point1[j] = i__;
/*< 36      CONTINUE >*/
/* L36: */
	}
/*<         DO 40 I = 1,16 >*/
	for (i__ = 1; i__ <= 16; ++i__) {
/*<           J = ICHAR(DUP(I:I)) >*/
	    j = *(unsigned char *)&dup[i__ - 1];
/*<           POINT2(J) = I >*/
	    iasci2_1.point2[j] = i__;
/*< 40      CONTINUE >*/
/* L40: */
	}
/*  DEAL WITH U */
/*<           J = ICHAR('U') >*/
	j = 'U';
/*<           POINT1(J) = 1   >*/
	iasci1_1.point1[j] = 1;
/*<           POINT2(J) = 1   >*/
	iasci2_1.point2[j] = 1;
/*<         DO 41 I = 1,16 >*/
	for (i__ = 1; i__ <= 16; ++i__) {
/*<           J = ICHAR(DLOW(I:I)) >*/
	    j = *(unsigned char *)&dlow[i__ - 1];
/*<           POINT2(J) = I >*/
	    iasci2_1.point2[j] = i__;
/*< 41      CONTINUE >*/
/* L41: */
	}
/*  DEAL WITH U */
/*<           J = ICHAR('u') >*/
	j = 'u';
/*<           POINT1(J) = 1   >*/
	iasci1_1.point1[j] = 1;
/*<           POINT2(J) = 1   >*/
	iasci2_1.point2[j] = 1;
/*<         ELSE IF(IDM.EQ.26)THEN >*/
    } else if (*idm == 26) {
/*<           DO 45 I = 0,255 >*/
	for (i__ = 0; i__ <= 255; ++i__) {
/*<             POINT1(I) = IDM >*/
	    iasci1_1.point1[i__] = *idm;
/*< 45        CONTINUE >*/
/* L45: */
	}

/*<         DO 50 I = 1,26 >*/
	for (i__ = 1; i__ <= 26; ++i__) {
/*<           J = ICHAR(PUP(I:I)) >*/
	    j = *(unsigned char *)&pup[i__ - 1];
/*<           POINT1(J) = I >*/
	    iasci1_1.point1[j] = i__;
/*< 50      CONTINUE >*/
/* L50: */
	}
/*<         DO 51 I = 1,26 >*/
	for (i__ = 1; i__ <= 26; ++i__) {
/*<           J = ICHAR(PLOW(I:I)) >*/
	    j = *(unsigned char *)&plow[i__ - 1];
/*<           POINT1(J) = I >*/
	    iasci1_1.point1[j] = i__;
/*< 51      CONTINUE >*/
/* L51: */
	}
/*<         DO 60 I = 0,255 >*/
	for (i__ = 0; i__ <= 255; ++i__) {
/*<           POINT2(I) = POINT1(I) >*/
	    iasci2_1.point2[i__] = iasci1_1.point1[i__];
/*< 60      CONTINUE >*/
/* L60: */
	}
/*<       END IF >*/
    }
/*<       END >*/
    return 0;
} /* initlu_ */

/*<       INTEGER FUNCTION INDEXA(STRING,ID,CHAR) >*/
integer indexa_(char *string, integer *id, char *char__, ftnlen string_len, 
	ftnlen char_len)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;

/*<       CHARACTER STRING(ID),CHAR >*/
/*  FUNCTION TO FIND FIRST OCCURRENCE OF CHAR IN STRING */
/*<       DO 10 I = 1,ID >*/
    /* Parameter adjustments */
    --string;

    /* Function Body */
    i__1 = *id;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         IF(STRING(I).EQ.CHAR)THEN >*/
	if (*(unsigned char *)&string[i__] == *(unsigned char *)char__) {
/*<           INDEXA = I >*/
	    ret_val = i__;
/*<           RETURN >*/
	    return ret_val;
/*<         END IF >*/
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       INDEXA = 0 >*/
    ret_val = 0;
/*<       END >*/
    return ret_val;
} /* indexa_ */

/*     SQCOPY */
/*   SEQUENCE COPYING PROGRAM */
/*<       SUBROUTINE SQCOPY(SEQNCE,COMSEQ,IDIM) >*/
/* Subroutine */ int sqcopy_(char *seqnce, char *comseq, integer *idim, 
	ftnlen seqnce_len, ftnlen comseq_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*   AUTHOR: RODGER STADEN */
/*<       CHARACTER SEQNCE(IDIM),COMSEQ(IDIM) >*/
/*<       DO 100 I=1,IDIM >*/
    /* Parameter adjustments */
    --comseq;
    --seqnce;

    /* Function Body */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         COMSEQ(I)=SEQNCE(I) >*/
	*(unsigned char *)&comseq[i__] = *(unsigned char *)&seqnce[i__];
/*< 100   CONTINUE >*/
/* L100: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* sqcopy_ */

/*<       INTEGER FUNCTION CTONUM(CHAR) >*/
integer ctonum_(char *char__, ftnlen char_len)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer icol;

/*  AUTHOR RODGER STADEN */
/*<       INTEGER POINT1(0:255) >*/
/*<       CHARACTER CHAR >*/
/*<       COMMON /IASCI1/POINT1 >*/
/*<       SAVE /IASCI1/ >*/

/*  GET COLLATING SEQUENCE VALUE */
/*<       ICOL = ICHAR(CHAR) >*/
    icol = *(unsigned char *)char__;
/*  THIS POINTS TO A VALUE IN POINTR */
/*<       CTONUM = POINT1(ICOL) >*/
    ret_val = iasci1_1.point1[icol];
/*<       END >*/
    return ret_val;
} /* ctonum_ */

/*<       SUBROUTINE BUB3AS(LIST,LIST2,LIST3,IDIM) >*/
/* Subroutine */ int bub3as_(integer *list, integer *list2, integer *list3, 
	integer *idim)
{
    static integer i__, j, itemp;

/*   AUTHOR: RODGER STADEN */
/*<       INTEGER LIST(IDIM),LIST2(IDIM),LIST3(IDIM) >*/
/*<       I=0 >*/
    /* Parameter adjustments */
    --list3;
    --list2;
    --list;

    /* Function Body */
    i__ = 0;
/*<       J=0 >*/
    j = 0;
/*< 10    CONTINUE >*/
L10:
/*   SET I=J IF WE HAVE JUST CORRECTLY POSITIONED AN ELEMENT */
/*<       IF(J.GT.I)I=J >*/
    if (j > i__) {
	i__ = j;
    }
/*<       I=I+1 >*/
    ++i__;
/*<       IF(I.EQ.IDIM)RETURN >*/
    if (i__ == *idim) {
	return 0;
    }
/*< 20    CONTINUE >*/
L20:
/*<       IF(LIST(I).LE.LIST(I+1))GO TO 10 >*/
    if (list[i__] <= list[i__ + 1]) {
	goto L10;
    }
/*   FIRST MOVE THIS ELEMENT? IF SO SET POINTER TO ITS INITIAL POSITION */
/*<       IF(J.LT.I)J=I >*/
    if (j < i__) {
	j = i__;
    }
/*<       ITEMP=LIST(I) >*/
    itemp = list[i__];
/*<       LIST(I)=LIST(I+1) >*/
    list[i__] = list[i__ + 1];
/*<       LIST(I+1)=ITEMP >*/
    list[i__ + 1] = itemp;
/*<       ITEMP=LIST2(I) >*/
    itemp = list2[i__];
/*<       LIST2(I)=LIST2(I+1) >*/
    list2[i__] = list2[i__ + 1];
/*<       LIST2(I+1)=ITEMP >*/
    list2[i__ + 1] = itemp;
/*<       ITEMP=LIST3(I) >*/
    itemp = list3[i__];
/*<       LIST3(I)=LIST3(I+1) >*/
    list3[i__] = list3[i__ + 1];
/*<       LIST3(I+1)=ITEMP >*/
    list3[i__ + 1] = itemp;
/*   DECREMENT BACK THRU LIST WITH THIS ELEMENT */
/*<       IF(I.GT.1)I=I-1 >*/
    if (i__ > 1) {
	--i__;
    }
/*<       GO TO 20 >*/
    goto L20;
/*<       END >*/
} /* bub3as_ */

