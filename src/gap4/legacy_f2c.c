/* legacy.f -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer points[256];
} shotc_ = {};

#define shotc_1 shotc_

struct {
    integer point1[256];
} iasci1_ = {};

#define iasci1_1 iasci1_

struct {
    integer point2[256];
} iasci2_ = {};

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
static integer c__9 = 9;

/* note to kfs: add extra option to assembly menu "Ignore previous data"? */
/* duplicate dialogue from "normal shotgun assembly" but */
/* add extra argument NOPT to argument list for dbauto. */
/* set to 0 for all existing calls to dbauto and to 1 for new option. */
/* Subroutine */ int dbauto_(relpg, lngthg, lnbr, rnbr, maxdb, idbsiz, ngels, 
	nconts, clist, maxgel, seq1, seq2, seq3, seq4, seq5, seqc2, seqg2, 
	seqg3, seqc3, rnames, maxseq, maxglm, savps, savpg, savl, maxsav, 
	cends, nends, maxcon, idev1, namarc, nampro, percd, iokent, ishow, 
	minmat, maxpg, permax, irepsc, iopt, nopt, ansjok, ansfe, iwing, nbad,
	 list, iok, minovr, seq1_len, seq2_len, seq3_len, seq4_len, seq5_len, 
	seqc2_len, seqg2_len, seqg3_len, seqc3_len, rnames_len, namarc_len, 
	nampro_len, list_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *maxdb, *idbsiz, *ngels, *nconts, *
	clist, *maxgel;
char *seq1, *seq2, *seq3, *seq4, *seq5, *seqc2, *seqg2, *seqg3, *seqc3, *
	rnames;
integer *maxseq, *maxglm, *savps, *savpg, *savl, *maxsav, *cends, *nends, *
	maxcon, *idev1;
char *namarc, *nampro;
real *percd;
integer *iokent, *ishow, *minmat, *maxpg;
real *permax;
integer *irepsc, *iopt, *nopt, *ansjok, *ansfe, *iwing, *nbad;
char *list;
integer *iok, *minovr;
ftnlen seq1_len;
ftnlen seq2_len;
ftnlen seq3_len;
ftnlen seq4_len;
ftnlen seq5_len;
ftnlen seqc2_len;
ftnlen seqg2_len;
ftnlen seqg3_len;
ftnlen seqc3_len;
ftnlen rnames_len;
ftnlen namarc_len;
ftnlen nampro_len;
ftnlen list_len;
{
    /* Format strings */
    static char fmt_1006[] = "(\002Processing \002,i8,\002 in batch\002)";
    static char fmt_1007[] = "(\002File name \002,a)";
    static char fmt_1800[] = "(\002Reading length \002,i6)";
    static char fmt_1022[] = "(a,f5.1,2i6)";
    static char fmt_10077[] = "(\002 Contig line for contig\002,i8,\002 not \
found!\002)";
    static char fmt_1014[] = "(\002New reading overlaps contig\002,i8)";
    static char fmt_1013[] = "(\002Overlap between contigs\002,i8,\002 an\
d\002,i8)";
    static char fmt_1012[] = "(\002Entering the new reading into contig\002,\
i8)";
    static char fmt_1002[] = "(\002Length of overlap between the contigs\002\
,i6)";
    static char fmt_1020[] = "(\002Could not join contigs\002,i8,\002 and\
\002,i8)";
    static char fmt_1021[] = "(\002Reading has been entered into contig\002,\
i8)";
    static char fmt_1017[] = "(\002Editing contig\002,i8)";
    static char fmt_1018[] = "(\002Completing the join between contigs\002,i\
8,\002 and\002,i8)";
    static char fmt_1030[] = "(i8,\002 sequences processed\002)";
    static char fmt_1031[] = "(i8,\002 sequences entered into database\002)";
    static char fmt_1032[] = "(i8,\002 joins made\002)";
    static char fmt_1033[] = "(i8,\002 joins failed\002)";

    /* System generated locals */
    integer seqc2_dim1, seqc2_offset, seqg2_dim1, seqg2_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    extern /* Subroutine */ int ccta_();
    static integer jobc, jgel;
    static char csen[1];
    static integer lreg, ilcr;
    extern /* Subroutine */ int info_();
    static integer mask, leno, rreg, ilct, ierr, cnum, idim1, idim2, i__, j, 
	    iladd[1], n, iradd[1], ifail[2], igelc, idim22[2], kfail, x;
    extern /* Subroutine */ int aline_();
    static real z__;
    static char infod[80];
    extern integer gclin_();
    extern /* Subroutine */ int sindb_();
    static integer jngel, imatc, igood, joinf, idsav, maxpc, itask, llino[2], 
	    iposc[2], iposg[2], joint[2], idout[2], ioptc, klass, iover, 
	    itype[2], lmost;
    extern /* Subroutine */ int sqrev_();
    static integer rmost;
    extern /* Subroutine */ int ajoin2_(), ajoin3_(), precn1_(), dbchek_(), 
	    abedin_();
    static integer pl[2], kt;
    extern integer gnread_();
    static integer pr[2];
    extern /* Subroutine */ int delcon_(), addtit_(), dbautp_();
    static integer ifcomp;
    extern /* Subroutine */ int aenter_();
    static integer lincon[2], ilefts[2], isense[2], itotpc[2];
    static real permis[2];
    static char infoud[80];
    static integer itotpg[2];
    extern integer cmpseq_();
    static integer minsli, iempty, jnjoin, maxovr;
    extern /* Subroutine */ int precon_(), erromf_();
    static integer ngelsl, ncontl;
    extern /* Subroutine */ int updout_(), arrfio_(), aerror_(), sqcopy_(), 
	    autocn_(), tolist_(), updcon_(), cmplmt_(), sqcomm_();
    static integer ilc[2], idm, ltl, ltr;

    /* Fortran I/O blocks */
    static icilist io___26 = { 0, infod, 0, fmt_1006, 80, 1 };
    static icilist io___27 = { 0, infod, 0, fmt_1007, 80, 1 };
    static icilist io___28 = { 0, infod, 0, fmt_1800, 80, 1 };
    static icilist io___48 = { 0, infoud, 0, fmt_1022, 80, 1 };
    static icilist io___55 = { 0, infod, 0, fmt_10077, 80, 1 };
    static icilist io___56 = { 0, infod, 0, fmt_1014, 80, 1 };
    static icilist io___57 = { 0, infod, 0, fmt_1013, 80, 1 };
    static icilist io___59 = { 0, infod, 0, fmt_1012, 80, 1 };
    static icilist io___60 = { 0, infod, 0, fmt_1012, 80, 1 };
    static icilist io___65 = { 0, infod, 0, fmt_1002, 80, 1 };
    static icilist io___66 = { 0, infod, 0, fmt_1012, 80, 1 };
    static icilist io___67 = { 0, infod, 0, fmt_1020, 80, 1 };
    static icilist io___68 = { 0, infod, 0, fmt_1021, 80, 1 };
    static icilist io___72 = { 0, infod, 0, fmt_1012, 80, 1 };
    static icilist io___76 = { 0, infod, 0, fmt_1017, 80, 1 };
    static icilist io___77 = { 0, infod, 0, fmt_1017, 80, 1 };
    static icilist io___80 = { 0, infod, 0, fmt_1018, 80, 1 };
    static icilist io___82 = { 0, infod, 0, fmt_1030, 80, 1 };
    static icilist io___83 = { 0, infod, 0, fmt_1031, 80, 1 };
    static icilist io___84 = { 0, infod, 0, fmt_1032, 80, 1 };
    static icilist io___85 = { 0, infod, 0, fmt_1033, 80, 1 };



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
    seqg2_offset = 1 + seqg2_dim1 * 1;
    seqg2 -= seqg2_offset;
    seqc2_dim1 = *maxglm;
    seqc2_offset = 1 + seqc2_dim1 * 1;
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
    idm = 5;
    minsli = 3;
    ifail[0] = 0;
    iempty = 0;
    maxpc = *maxpg;
/* note by KFS (22/5/95) - use of CLIST and CNUM by auto_assemble. */
/* auto_assemble is a special case where all contigs are always used */
/* so although CLIST = NULL, CNUM = NCONTS and this fact is catered for in */
/* get_contig_list */
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
    if (*ngels < 1 || *nopt == 1) {
	iempty = 1;
    }
    dbchek_(idev1, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &idm, idbsiz, 
	    ngels, nconts, &ierr);
    if (ierr > 1) {
	return 0;
    }
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


    jgel = 0;
    jngel = 0;
    jnjoin = 0;
    joinf = 0;
    mask = 0;
    imatc = 0;
    if (*iopt == 3 || *iopt == 4) {
	*minmat = 1;
    }
    if (*iopt == 2) {
	mask = 3;
    }
    if (*iopt == 5) {
	mask = 4;
    }
    if (*iopt == 1 || *iopt == 2 || *iopt == 5) {
	idim1 = 0;
	maxovr = *maxgel - max(maxpc,*maxpg) * 3;
	imatc = 0;
/* get ready for precon */

/* set task (normal consensus+title) */

	itask = 5;

/* set masking if required */

	if (mask == 3) {
	    itask += 32;
	}
	if (mask == 4) {
	    itask = 9;
	}
	if (iempty == 0) {
	    precon_(seq1 + 1, nampro, percd, idbsiz, &cnum, &clist[1], &itask,
		     idev1, &idim1, maxgel, maxseq, iwing, nbad, iladd, iradd,
		     ifail, (ftnlen)1, nampro_len);
	    if (ifail[0] != 0) {
		erromf_("Error calculating consensus", (ftnlen)27);
		goto L901;
	    }
/*        IF(IDIM1.GT.0)CALL FMTDB(SEQ1,IDIM1,1,IDIM1,60,30) */
/*        IF (0.EQ.0) RETURN */
	}
/*        WRITE(*,*)'INITIAL IDIM1',IDIM1,ITASK */

/* init hashing constants */

	ioptc = 1;
	idsav = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1], 
		maxsav, seq1 + 1, seq2 + 1, maxseq, maxgel, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
	if (idsav != 0) {
	    ioptc = 6;
	    idsav = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[
		    1], &idsav, seq1 + 1, seq2 + 1, maxseq, maxgel, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
	    return 0;
	}
    }

/* set intitial values for contig count and number of readings */
/* just to get thru the first consensus calc */
    ngelsl = *ngels + 2;
    ncontl = *nconts + 1;


/*                          MAIN LOOP */


L1:
    updout_();


    idim2 = *maxgel;
    *iok = gnread_(namarc, namarc_len);
    if (*iok == 1) {
	goto L900;
    }
    if (*iok != 0) {
	goto L1;
    }
    info_(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", (ftnlen)46);
/*      WRITE(*,*)'>>>>>>>>>>>',NAMARC */
/*      WRITE(6,*)'IDIM1',IDIM1,NGELS */
/*         CALL FMTDB1(SEQ1,IDIM1,1,IDIM1,60,6) */
/*      WRITE(*,*)'MAIN' */
/*      WRITE(*,*)(SEQ1(JJJJ),JJJJ=1,IDIM1) */
    ++jgel;
    s_wsfi(&io___26);
    do_fio(&c__1, (char *)&jgel, (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
/* L1007: */
    s_wsfi(&io___27);
    do_fio(&c__1, namarc, namarc_len);
    e_wsfi();
    info_(infod, (ftnlen)80);
/* Added by Simon 23-March-1993 */
    arrfio_(namarc, seq2 + 1, &idim2, &c__1, iok, namarc_len, (ftnlen)1);
/*      CALL OPENRS(IDEV4,NAMARC,IOK,LRECL,2) */
    if (*iok != 0) {
/*        IF(INF.EQ.1) RETURN */
	aerror_(list, namarc, &c__0, list_len, namarc_len);
	goto L1;
    }
    s_wsfi(&io___28);
    do_fio(&c__1, (char *)&idim2, (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);



    if (idim2 < *minmat) {
	aerror_(list, namarc, &c__1, list_len, namarc_len);
	goto L1;
    }
    if (*iopt == 3 || *iopt == 4) {
	dbautp_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq2 + 1, namarc, joint, itype, isense, seqc2 + seqc2_offset, 
		itotpc, &idim2, idout, llino, lincon, ifail, idbsiz, maxdb, 
		idev1, maxgel, &imatc, &iempty, rnames + rnames_len, iopt, (
		ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
	if (ifail[0] != 0) {
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
	} else {
	    ++jngel;
	}
	goto L1;
    }
    sqcopy_(seq2 + 1, seq3 + 1, &idim2, (ftnlen)1, (ftnlen)1);
    ifcomp = 0;
    imatc = 0;
    ifail[0] = 0;
    ifail[1] = 0;
    jobc = 2;
    if (ngelsl < *ngels && ncontl < *nconts) {
	jobc = 1;
    } else if (ngelsl == *ngels) {
	jobc = 0;
    }
    ngelsl = *ngels;
    ncontl = *nconts;
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
    if (*irepsc == 0) {
	if (ifcomp != 0) {
	    aerror_(list, namarc, &c__2, list_len, namarc_len);
	    goto L1;
	}
	if (imatc <= 0) {
	    permis[0] = (float)0.;
	    leno = 0;
	}
	s_wsfi(&io___48);
	do_fio(&c__1, namarc, namarc_len);
	do_fio(&c__1, (char *)&permis[0], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&idim2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&leno, (ftnlen)sizeof(integer));
	e_wsfi();
	tolist_(list, infoud, list_len, (ftnlen)80);
    }
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
    if (ifcomp != 0) {
	aerror_(list, namarc, &c__2, list_len, namarc_len);
	goto L1;
    }
    sqcopy_(seq3 + 1, seq2 + 1, &idim2, (ftnlen)1, (ftnlen)1);


/* No overlap below mismatch cutoff */

    if (imatc == 0) {

/* if masking then count as failure (code 5) unless we are entering all reads */
/* ie if we are masking and there is no match we do not enter unless "enter */
/* all reads (ANSFE=2) is set. */

	if (mask != 0) {
	    if (*ansfe == 1) {
		aerror_(list, namarc, &c__5, list_len, namarc_len);
		goto L1;
	    }
	}

/*                     NO OVERLAP NEW CONTIG */

	if (ifail[0] != 0) {
	    if (*ansfe == 1) {
		aerror_(list, namarc, &c__2, list_len, namarc_len);
		goto L1;
	    }
	    info_("New reading overlaps poorly: start a new contig", (ftnlen)
		    47);
	} else {
	    info_("New reading does not overlap: start a new contig", (ftnlen)
		    48);
	}
/*     ITYPE 0 = NO OVERLAP */
/*     ISENSE 1 = SAME SENSE AS ARCHIVE */
	itype[0] = 0;
	isense[0] = 1;
	idout[0] = *maxgel;
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq2 + 1, namarc, &x, itype, isense, seqc2 + (seqc2_dim1 + 1),
		 itotpc, &idim2, idout, llino, lincon, ifail, idbsiz, idev1, 
		maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (ftnlen)1,
		 rnames_len);
	if (ifail[0] != 0) {
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
	    goto L1;
	}
	iempty = 0;
	++idim1;

/* new start */

	lreg = 1;
	rreg = idim2;
	lincon[0] = *idbsiz - *nconts;
	precn1_(seq1 + 1, nampro, percd, idbsiz, lincon, &lreg, &rreg, &itask,
		 idev1, &idim1, maxgel, maxseq, iwing, nbad, iladd, iradd, 
		ifail, (ftnlen)1, nampro_len);
	if (ifail[0] != 0) {
	    erromf_("Error calculating consensus", (ftnlen)27);
	    goto L900;
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
	++jngel;
	goto L1;
    }


/*         OVERLAP SO TRY TO ENTER THE READING */


    i__1 = imatc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	n = *idbsiz - *nconts;
	i__2 = *idbsiz - 1;
	for (j = n; j <= i__2; ++j) {
	    if (lnbr[j] != llino[i__ - 1]) {
		goto L99;
	    }
	    lincon[i__ - 1] = j;
	    goto L100;
L99:
	    ;
	}
	s_wsfi(&io___55);
	do_fio(&c__1, (char *)&llino[i__ - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	erromf_(infod, (ftnlen)80);
	goto L900;
L100:
	;
    }

    if (imatc == 1) {


/*                     SINGLE OVERLAP */



	s_wsfi(&io___56);
	do_fio(&c__1, (char *)&llino[0], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	if (itotpg[0] > 0) {
	    ccta_(seqg2 + (seqg2_dim1 + 1), idim22, (ftnlen)1);
	}
/*      WRITE(*,*)'BEFORE entry' */
/*      WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDIM1) */
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (seqg2_dim1 + 1), namarc, joint, itype, isense, seqc2 
		+ (seqc2_dim1 + 1), itotpc, idim22, idout, llino, lincon, 
		ifail, idbsiz, idev1, maxgel, rnames + rnames_len, (ftnlen)1, 
		namarc_len, (ftnlen)1, rnames_len);
	if (ifail[0] != 0) {
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
	    goto L1;
	}
	++jngel;
/*      WRITE(*,*)'after entry' */
/*      WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDIM1) */
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, ilefts, ilc, lincon, nampro,
		 seq2 + 1, idev1, ifail, maxgel, &idm, percd, &mask, &clist[1]
		, (ftnlen)1, nampro_len, (ftnlen)1);
	if (ifail[0] != 0) {
	    erromf_("Error calculating consensus", (ftnlen)27);
	    goto L900;
	}
	if (kfail != 0) {
/*          CALL AERROR(LIST,NAMARC,4) */
	}
	goto L1;
    }


/*                     DOUBLE OVERLAP */


    s_wsfi(&io___57);
    do_fio(&c__2, (char *)&llino[0], (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    if (*ansjok != 0) {

/*  read overlaps 2 contigs but joins are forbidden */
/*  stick in in one of the contigs but do not join */

	igood = 1;
	info_("Read overlaps 2 contigs: entering it at best site", (ftnlen)49)
		;
	if (itotpg[igood - 1] > 0) {
	    ccta_(seqg2 + (igood * seqg2_dim1 + 1), &idim22[igood - 1], (
		    ftnlen)1);
	}
	s_wsfi(&io___59);
	do_fio(&c__1, (char *)&llino[igood - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (igood * seqg2_dim1 + 1), namarc, &joint[igood - 1], &
		itype[igood - 1], &isense[igood - 1], seqc2 + (igood * 
		seqc2_dim1 + 1), &itotpc[igood - 1], &idim22[igood - 1], &
		idout[igood - 1], &llino[igood - 1], &lincon[igood - 1], &
		ifail[igood - 1], idbsiz, idev1, maxgel, rnames + rnames_len, 
		(ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
	if (ifail[igood - 1] != 0) {
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
	    goto L1;
	}
	++jngel;
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, &ilefts[igood - 1], &ilc[
		igood - 1], &lincon[igood - 1], nampro, seq2 + 1, idev1, 
		ifail, maxgel, &idm, percd, &mask, &clist[1], (ftnlen)1, 
		nampro_len, (ftnlen)1);
	if (ifail[0] != 0) {
	    erromf_("Error calculating consensus", (ftnlen)27);
	    goto L900;
	}
	goto L1;
    }



    if (llino[0] == llino[1]) {

/*  read overlaps twice in one contig - stick it in the best place */

	igood = 1;
	info_("Read overlaps twice in one contig: entering it at best site", (
		ftnlen)59);
	if (itotpg[igood - 1] > 0) {
	    ccta_(seqg2 + (igood * seqg2_dim1 + 1), &idim22[igood - 1], (
		    ftnlen)1);
	}
	s_wsfi(&io___60);
	do_fio(&c__1, (char *)&llino[igood - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (igood * seqg2_dim1 + 1), namarc, &joint[igood - 1], &
		itype[igood - 1], &isense[igood - 1], seqc2 + (igood * 
		seqc2_dim1 + 1), &itotpc[igood - 1], &idim22[igood - 1], &
		idout[igood - 1], &llino[igood - 1], &lincon[igood - 1], &
		ifail[igood - 1], idbsiz, idev1, maxgel, rnames + rnames_len, 
		(ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
	if (ifail[igood - 1] != 0) {
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
	    goto L1;
	}
	++jngel;
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, &ilefts[igood - 1], &ilc[
		igood - 1], &lincon[igood - 1], nampro, seq2 + 1, idev1, 
		ifail, maxgel, &idm, percd, &mask, &clist[1], (ftnlen)1, 
		nampro_len, (ftnlen)1);
	if (ifail[0] != 0) {
	    erromf_("Error calculating consensus", (ftnlen)27);
	    goto L900;
	}
	goto L1;
    }

/* is the overlap between the contigs too large? */

    ajoin3_(&relpg[1], idbsiz, lincon, itype, isense, joint, idim22, &klass, &
	    iover, pl, pr);
    s_wsfi(&io___65);
    do_fio(&c__1, (char *)&iover, (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    if (iover > maxovr) {
	info_("Overlap too large: entry only", (ftnlen)29);

/* cannot align the two contigs, so try to enter the reading into one of them */

	ifail[1] = 1;
	igood = 0;
	if (ifail[0] == 0) {
	    igood = 1;
	}
	if (ifail[1] == 0) {
	    igood = 2;
	}
	if (igood == 0) {
	    aerror_(list, namarc, &c__2, list_len, namarc_len);
	    ++joinf;
	    goto L1;
	}
	if (itotpg[igood - 1] > 0) {
	    ccta_(seqg2 + (igood * seqg2_dim1 + 1), &idim22[igood - 1], (
		    ftnlen)1);
	}
	s_wsfi(&io___66);
	do_fio(&c__1, (char *)&llino[igood - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seqg2 + (igood * seqg2_dim1 + 1), namarc, &joint[igood - 1], &
		itype[igood - 1], &isense[igood - 1], seqc2 + (igood * 
		seqc2_dim1 + 1), &itotpc[igood - 1], &idim22[igood - 1], &
		idout[igood - 1], &llino[igood - 1], &lincon[igood - 1], &
		ifail[igood - 1], idbsiz, idev1, maxgel, rnames + rnames_len, 
		(ftnlen)1, namarc_len, (ftnlen)1, rnames_len);
	if (ifail[igood - 1] != 0) {
	    aerror_(list, namarc, &c__3, list_len, namarc_len);
	    ++joinf;
	    goto L1;
	}
	++jngel;
	updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, 
		nconts, seq1 + 1, maxseq, &idim1, &ilefts[igood - 1], &ilc[
		igood - 1], &lincon[igood - 1], nampro, seq2 + 1, idev1, 
		ifail, maxgel, &idm, percd, &mask, &clist[1], (ftnlen)1, 
		nampro_len, (ftnlen)1);
	if (ifail[0] != 0) {
	    erromf_("Error calculating consensus", (ftnlen)27);
	    goto L900;
	}
	s_wsfi(&io___67);
	do_fio(&c__2, (char *)&llino[0], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
/* L1021: */
	s_wsfi(&io___68);
	do_fio(&c__1, (char *)&llino[igood - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	++joinf;
	goto L1;
    }
/*   WHICH CONTIG IS LEFTMOST? */
    lmost = 1;
    rmost = 2;
    if (pl[0] > pl[1]) {
	lmost = 2;
	rmost = 1;
    }
/*   SAVE LENGTH OF RMOST CONTIG FOR DELETION STEP LATER */
    ilcr = ilc[rmost - 1];
    if (itotpg[lmost - 1] > 0) {
	ccta_(seqg2 + (lmost * seqg2_dim1 + 1), &idim22[lmost - 1], (ftnlen)1)
		;
    }
    s_wsfi(&io___72);
    do_fio(&c__1, (char *)&llino[lmost - 1], (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, seqg2 + 
	    (lmost * seqg2_dim1 + 1), namarc, &joint[lmost - 1], &itype[lmost 
	    - 1], &isense[lmost - 1], seqc2 + (lmost * seqc2_dim1 + 1), &
	    itotpc[lmost - 1], &idim22[lmost - 1], &idout[lmost - 1], &llino[
	    lmost - 1], &lincon[lmost - 1], &ifail[lmost - 1], idbsiz, idev1, 
	    maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (ftnlen)1, 
	    rnames_len);
    if (ifail[lmost - 1] != 0) {
	aerror_(list, namarc, &c__3, list_len, namarc_len);
	++joinf;
	goto L1;
    }
    ++jngel;
    updcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, ngels, nconts, 
	    seq1 + 1, maxseq, &idim1, &ilefts[lmost - 1], &ilc[lmost - 1], &
	    lincon[lmost - 1], nampro, seq2 + 1, idev1, ifail, maxgel, &idm, 
	    percd, &mask, &clist[1], (ftnlen)1, nampro_len, (ftnlen)1);
    if (ifail[0] != 0) {
	erromf_("Error calculating consensus", (ftnlen)27);
	goto L900;
    }
    if (itype[lmost - 1] == 1) {
	llino[lmost - 1] = *ngels;
    }
    if (ilefts[lmost - 1] < ilefts[rmost - 1]) {
	ilefts[rmost - 1] += relpg[lincon[lmost - 1]] - ilc[lmost - 1];
    }
    ilc[lmost - 1] = relpg[lincon[lmost - 1]];
    for (i__ = 1; i__ <= 2; ++i__) {
	if (isense[i__ - 1] == -1) {
	    cmplmt_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    &lincon[i__ - 1], &llino[i__ - 1], seq2 + 1, idbsiz, 
		    idev1, maxgel, (ftnlen)1);
	    sqrev_(seq1 + ilefts[i__ - 1], &ilc[i__ - 1], (ftnlen)1);
	    sqcomm_(seq1 + ilefts[i__ - 1], &ilc[i__ - 1], (ftnlen)1);
	    kt = idim1;
	    addtit_(seq1 + (ilefts[i__ - 1] - 20), nampro, &lnbr[lincon[i__ - 
		    1]], &kt, (ftnlen)1, nampro_len);
	}
/* L500: */
    }
/*   NEED TO KNOW POSITION OF OVERLAP RELATIVE TO CONTIG, TO CONSENSUS */
/*   WHICH BITS TO SEND TO ALIGNMENT ROUTINES */
/*   SET UP FOR ALINE (NOTE RMOST IS EQUIVALENT TO THE GEL READING AND */
/*   SO IS SLID ALONG THE LMOST CONTIG. THE SECTION SENT TO ALINE MUST */
/*   BE OF LENGTH < MAXGEL-2*MAX(MAXPC,MAXPG) */
/*   IT MUST START AT POSITION 1 IN THE RMOST CONTIG AND EXTEND */
    iposc[lmost - 1] = pl[rmost - 1] + relpg[*ngels] - 1;
    ilct = relpg[lincon[lmost - 1]] - relpg[*ngels] - pl[rmost - 1] + 2 + 
	    maxpc;

/* change 5-6-95 line below to line above */

/*      ILCT = RELPG(LINCON(LMOST)) - RELPG(NGELS) - PL(RMOST) + 2 */
/* Computing MIN */
    i__1 = ilct, i__2 = ilc[rmost - 1];
    ilc[rmost - 1] = min(i__1,i__2);
/*      WRITE(*,*)'ILC(LMOST)',ILC(LMOST) */
/*      WRITE(*,*)'RELPG(LINCON(LMOST))',RELPG(LINCON(LMOST)) */
/*      WRITE(*,*)'RELPG(NGELS)',RELPG(NGELS) */
/*      WRITE(*,*)'PL(RMOST)',PL(RMOST) */
/*      WRITE(*,*)'LENGTH OF OVERLAP',ILC(RMOST) */
    iposc[rmost - 1] = 1;
    idout[lmost - 1] = *maxgel;
    idout[rmost - 1] = *maxgel;
    idsav = *maxsav;
/*  ON INPUT TO ALINE ILC(RMOST) CONTAINS THE OVERLAP LENGTH */
/*  ON OUTPUT IT CONTAINS THE LENGTH OF THE ALIGNED SECTION (IE INCLUDING */
/*  PADS) */
    info_("Trying to align the two contigs", (ftnlen)31);
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
    if (ifail[0] != 0) {
	info_("Failed to align the two overlapping contigs", (ftnlen)43);
/*        CALL AERROR(LIST,NAMARC,4) */
	++joinf;
	goto L1;
    }
    if (itotpc[lmost - 1] > 0) {
	s_wsfi(&io___76);
	do_fio(&c__1, (char *)&llino[lmost - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq3 + 1, &lincon[lmost - 1], &joint[lmost - 1], seqc2 + (
		lmost * seqc2_dim1 + 1), &itotpc[lmost - 1], &idout[lmost - 1]
		, idbsiz, idev1, maxgel, (ftnlen)1, (ftnlen)1);
    }
    joint[rmost - 1] = 1;
    idout[rmost - 1] = ilc[rmost - 1];
    if (itotpc[rmost - 1] > 0) {
	s_wsfi(&io___77);
	do_fio(&c__1, (char *)&llino[rmost - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq3 + 1, &lincon[rmost - 1], &joint[rmost - 1], seqc2 + (
		rmost * seqc2_dim1 + 1), &itotpc[rmost - 1], &idout[rmost - 1]
		, idbsiz, idev1, maxgel, (ftnlen)1, (ftnlen)1);
    }
    ilc[rmost - 1] = ilcr;
    ltl = lnbr[lincon[lmost - 1]];
    ltr = lnbr[lincon[rmost - 1]];
    s_wsfi(&io___80);
    do_fio(&c__1, (char *)&lnbr[lincon[lmost - 1]], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&lnbr[lincon[rmost - 1]], (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    ajoin2_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    &joint[lmost - 1], &ltl, &ltr, &lincon[lmost - 1], &lincon[rmost 
	    - 1], idev1);
    llino[0] = ltl;
    if (ilefts[lmost - 1] > ilefts[rmost - 1]) {
	delcon_(seq1 + 1, &ilefts[lmost - 1], &ilc[lmost - 1], &idim1, (
		ftnlen)1);
	delcon_(seq1 + 1, &ilefts[rmost - 1], &ilc[rmost - 1], &idim1, (
		ftnlen)1);
    }
    if (ilefts[rmost - 1] >= ilefts[lmost - 1]) {
	delcon_(seq1 + 1, &ilefts[rmost - 1], &ilc[rmost - 1], &idim1, (
		ftnlen)1);
	delcon_(seq1 + 1, &ilefts[lmost - 1], &ilc[lmost - 1], &idim1, (
		ftnlen)1);
    }
    lreg = 1;
    rreg = joint[lmost - 1];
    igelc = llino[0];
/*      JOB = 1 */

/* assume itask set above */

/* also need to add 1 to consensus position */

    ++idim1;

/* for precon need to find the contig line number after the join */

    lincon[lmost - 1] = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], 
	    ngels, nconts, idbsiz, &igelc);
    if (lincon[lmost - 1] == 0) {
	erromf_("Cannot find contig line! Quitting", (ftnlen)33);
	goto L900;
    }
    precn1_(seq1 + 1, nampro, percd, idbsiz, &lincon[lmost - 1], &lreg, &rreg,
	     &itask, idev1, &idim1, maxgel, maxseq, iwing, nbad, iladd, iradd,
	     ifail, (ftnlen)1, nampro_len);
    if (ifail[0] != 0) {
	erromf_("Error calculating consensus", (ftnlen)27);
	goto L900;
    }
    ++jnjoin;
    if (kfail != 0) {
/*        CALL AERROR(LIST,NAMARC,4) */
/*        JOINF = JOINF + 1 */
    }
    goto L1;
L900:
    if (*iopt == 1 || *iopt == 2 || *iopt == 5) {
	ioptc = 6;
	idsav = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1], 
		&idsav, seq1 + 1, seq2 + 1, maxseq, maxgel, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
    }
L901:
    info_("Batch finished", (ftnlen)14);
    s_wsfi(&io___82);
    do_fio(&c__1, (char *)&jgel, (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    s_wsfi(&io___83);
    do_fio(&c__1, (char *)&jngel, (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    s_wsfi(&io___84);
    do_fio(&c__1, (char *)&jnjoin, (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    s_wsfi(&io___85);
    do_fio(&c__1, (char *)&joinf, (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
} /* dbauto_ */

/* Subroutine */ int dbautp_(relpg, lngthg, lnbr, rnbr, ngels, nconts, seq2, 
	namarc, joint, itype, isense, seqc2, itotpc, idim2, idout, llino, 
	lincon, ifail, idbsiz, maxdb, idev1, maxgel, imatc, iempty, rnames, 
	iopt, seq2_len, namarc_len, seqc2_len, rnames_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts;
char *seq2, *namarc;
integer *joint, *itype, *isense;
char *seqc2;
integer *itotpc, *idim2, *idout, *llino, *lincon, *ifail, *idbsiz, *maxdb, *
	idev1, *maxgel, *imatc, *iempty;
char *rnames;
integer *iopt;
ftnlen seq2_len;
ftnlen namarc_len;
ftnlen seqc2_len;
ftnlen rnames_len;
{
    extern /* Subroutine */ int aenter_();

/*  deals with entering all readings into contig 1 (IOPT=3) */
/*  or all readings into new contigs (IOPT=4) */
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
	if (*imatc == 0) {
	    *itype = 0;
	    *isense = 1;
	    *idout = *maxgel;
	    aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    seq2 + 1, namarc, joint, itype, isense, seqc2 + 1, itotpc,
		     idim2, idout, llino, lincon, ifail, idbsiz, idev1, 
		    maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (
		    ftnlen)1, rnames_len);
	    if (*ifail != 0) {
		return 0;
	    }
	    *iempty = 0;
	    *imatc = 1;
	} else {
	    *itype = -1;
	    *isense = 1;
	    *joint = 1;
	    *llino = *ngels;
	    *lincon = *idbsiz - *nconts;
	    *itotpc = 0;
	    *idout = *maxgel;
	    aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    seq2 + 1, namarc, joint, itype, isense, seqc2 + 1, itotpc,
		     idim2, idout, llino, lincon, ifail, idbsiz, idev1, 
		    maxgel, rnames + rnames_len, (ftnlen)1, namarc_len, (
		    ftnlen)1, rnames_len);
	    if (*ifail != 0) {
		return 0;
	    }
	}
    } else if (*iopt == 4) {
	*itype = 0;
	*isense = 1;
	*idout = *maxgel;
	aenter_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		seq2 + 1, namarc, joint, itype, isense, seqc2 + 1, itotpc, 
		idim2, idout, llino, lincon, ifail, idbsiz, idev1, maxgel, 
		rnames + rnames_len, (ftnlen)1, namarc_len, (ftnlen)1, 
		rnames_len);
	if (*ifail != 0) {
	    return 0;
	}
    }
} /* dbautp_ */

/*   SUBROUTINE TO ENTER NEW GEL SEQUENCES INTO DATA BASE. */
/*   IT READS IN AN ARCHIVE VERSION AND WRITES OUT A WORKING VERSION. */
/*   IT ALSO SETS UP ANY RELATIONSHIPS WITH OTHER DATA IN THE DATABASE */
/*   BOTH BY POSITION IN A CONTIG AND POINTERS TO LEFT AND RIGHT */
/*   NEIGHBOURS. */
/* Subroutine */ int aenter_(relpg, lngthg, lnbr, rnbr, ngels, nconts, gel, 
	namarc, x, itype, isense, seqc2, itotpc, idim, idc, ncontc, lincon, 
	ifail, idbsiz, idevr, maxgel, rnames, gel_len, namarc_len, seqc2_len, 
	rnames_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts;
char *gel, *namarc;
integer *x, *itype, *isense;
char *seqc2;
integer *itotpc, *idim, *idc, *ncontc, *lincon, *ifail, *idbsiz, *idevr, *
	maxgel;
char *rnames;
ftnlen gel_len;
ftnlen namarc_len;
ftnlen seqc2_len;
ftnlen rnames_len;
{
    /* Format strings */
    static char fmt_1013[] = "(\002New reading already in database with numb\
er\002,i8,\002 Entry aborted\002)";
    static char fmt_1003[] = "(\002This gel reading has been given the numbe\
r \002,i8)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    extern integer indb_();
    extern /* Subroutine */ int info_();
    static integer itmp, j, k, n, y;
    static char namid[40], infod[80];
    extern /* Subroutine */ int sindb_();
    static integer idevn;
    extern /* Subroutine */ int abedin_(), idline_(), erromf_(), writec_(), 
	    writeg_(), shiftt_(), stikit_(), writrn_();
    static integer iok;

    /* Fortran I/O blocks */
    static icilist io___89 = { 0, infod, 0, fmt_1013, 80, 1 };
    static icilist io___91 = { 0, infod, 0, fmt_1003, 80, 1 };


/*   AUTHOR: RODGER STADEN */
/*      WRITE(*,*)'IN ENTER',NGELS */
/*      WRITE(*,*)'X,ITYPE,ISENSE,IDIM,IDC' */
/*      WRITE(*,*)X,ITYPE,ISENSE,IDIM,IDC */
/*   SET FAIL FLAG */
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
    if (*idbsiz - (*ngels + *nconts) > 2) {
	goto L5;
    }
/*   FULL */
    erromf_("Database full!", (ftnlen)14);
    *ifail = 7;
    return 0;
L5:
/*   NEED TO CHECK TO SEE IF GEL ALREADY IN DB */
/*   LOOK THRU ARC FILE */
    idline_(namarc, namid, namarc_len, (ftnlen)40);
    j = indb_(ngels, rnames + 40, namid, (ftnlen)40, (ftnlen)40);
    if (j != 0) {
/*   FOUND */
	s_wsfi(&io___89);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	e_wsfi();
	erromf_(infod, (ftnlen)80);
	*ifail = 6;
	return 0;
    }
/*   INCREMENT NUMBER OF GELS */
    ++(*ngels);

/* set dummy int for idevn */

    idevn = 0;
    sindb_(&idevn, ngels, rnames + 40, namid, &c__2, (ftnlen)40, (ftnlen)40);
/*   SET LENGTH THIS GEL */
    lngthg[*ngels] = *idim * *isense;
    s_wsfi(&io___91);
    do_fio(&c__1, (char *)&(*ngels), (ftnlen)sizeof(integer));
    e_wsfi();
/*      WRITE(*,1003)NGELS */
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
    if (*itype != 0) {
	goto L100;
    }

/*   DOES NOT OVERLAP SO IT STARTS A CONTIG OF ITS OWN */

/*   SET CONTIG POINTERS AND GENERAL VALUES */
/*   INCREMENT NUMBER OF CONTIGS */
    ++(*nconts);
/*   POINTER TO THIS CONTIG */
    n = *idbsiz - *nconts;
/*   POINTER TO LEFT GEL THIS CONTIG */
    lnbr[n] = *ngels;
/*   POINTER TO RIGHT GEL THIS CONTIG */
    rnbr[n] = *ngels;
/*   LENGTH OF CONTIG */
    relpg[n] = *idim;
/*   WRITE CONTIG DESCRIPTOR */
    i__1 = *idbsiz - n;
    writec_(idevr, &i__1, &relpg[n], &lnbr[n], &rnbr[n]);
/*     Setup tags, original positions, conf values, vectors etc */
    i__1 = *idbsiz - n;
    stikit_(idevr, namarc, ngels, &lngthg[*ngels], gel + 1, maxgel, &iok, &
	    i__1, &c__1, namarc_len, (ftnlen)1);
    if (iok != 0) {
	--(*nconts);
	--(*ngels);
	*ifail = 1;
	return 0;
    }

/*   Create gel info */
/*   SET LEFT AND RIGHT POINTERS TO ZERO,RELPG TO 1 */
    lnbr[*ngels] = 0;
    rnbr[*ngels] = 0;
    relpg[*ngels] = 1;
/*   WRITE NEW GEL LINE */
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   WRITE DB DESCRIPTOR */
    writrn_(idevr, ngels, nconts);
    return 0;

L100:


/*     Shift tags if this new gel adjusts the left end */

    if (*itype == 1) {
	i__1 = *idbsiz - *lincon;
	i__2 = *x - 1;
	shiftt_(idevr, &i__1, &c__1, &i__2);
	itmp = 1;
    } else {
	itmp = *x;
    }
/*     Setup tags, original positions, conf values, vectors etc */
    i__1 = *idbsiz - *lincon;
    stikit_(idevr, namarc, ngels, &lngthg[*ngels], gel + 1, maxgel, &iok, &
	    i__1, &itmp, namarc_len, (ftnlen)1);
    if (iok != 0) {
	--(*ngels);
	*ifail = 1;
	return 0;
    }

/*   DOES OVERLAP */
/* L150: */

/*   LEFT END OR RIGHT OVERLAP? */
    if (*itype == 1) {
	goto L400;
    }
/*   RIGHT END OR INTERNAL OVERLAP */

/* L160: */
/*   NEED TO SEARCH THRU THIS CONTIG TO FIND LEFT AND RIGHT */
/*   NEIGHBOURS FOR THIS NEW GEL */
/*   LINE NUMBER OF LEFT END OF CONTIG */
    n = *ncontc;
/*   LOOK THRU UNTIL CURRENT IS >= THEN IT MUST BE THE PREVIOUS ONE */
L200:
    if (relpg[n] > *x) {
	goto L250;
    }
/*   IS THIS THE LAST GEL IN CONTIG? */
    if (rnbr[n] == 0) {
	goto L350;
    }
/*   NO SO LOOK AT NEXT */
    n = rnbr[n];
    goto L200;
L250:
/*   GEL LIES BETWEEN N AND LNBR(N) */
/*   NEED TO EDIT DB HERE */
    if (*itotpc > 0) {
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, x, seqc2 + 1, itotpc, idc, idbsiz, idevr, maxgel,
		 (ftnlen)1, (ftnlen)1);
    }


/*   SET POINTERS IN NEW GEL */
    lnbr[*ngels] = lnbr[n];
    rnbr[*ngels] = n;
    relpg[*ngels] = *x;
/*   WRITE NEW GEL LINE */
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   SET POINTERS  IN LEFT AND RIGHT NEIGHBOURS */
    k = lnbr[n];
    rnbr[k] = *ngels;
/*      RNBR(LNBR(N))=NGELS */
/*   WRITE LEFT AND RIGHT NEIGHBOURS */
    writeg_(idevr, &k, &relpg[k], &lngthg[k], &lnbr[k], &rnbr[k]);
    lnbr[n] = *ngels;
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
/*   WRITE NGELS NCONTS */
    writrn_(idevr, ngels, nconts);
/*   HAVE WE INCREASED LENGTH OF CONTIG? */
/*   ITS LINE NUMBER IS LINCON */
/*   NEED TO UPDATE IDIM IN CASE OF EDITS */
    *idim = (i__1 = lngthg[*ngels], abs(i__1));
    y = *x + *idim - 1;
    if (y <= relpg[*lincon]) {
	return 0;
    }
    relpg[*lincon] = y;
/*   WRITE NEW CONTIG LINE */
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
    return 0;
L350:
/*   MUST BE A RIGHT END OVERLAP */
/*   NEED TO EDIT DB HERE */
    if (*itotpc > 0) {
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, x, seqc2 + 1, itotpc, idc, idbsiz, idevr, maxgel,
		 (ftnlen)1, (ftnlen)1);
    }


/*   SET POINTERS FOR NEW GEL */
    lnbr[*ngels] = n;
    rnbr[*ngels] = 0;
    relpg[*ngels] = *x;
/*   WRITE NEW GEL LINE */
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   OLD RIGHT END */
    rnbr[n] = *ngels;
/*   WRITE NEW RIGHT LINE */
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
/*   RESET RIGHT NAME IN CONTIG */
/*   ITS LINE NUMBER IS LINCON */
    rnbr[*lincon] = *ngels;
/*   HAVE WE INCREASED LENGTH OF CONTIG? */
/*   NEED TO UPDATE LENGTH OF GEL IN CASE OF EDITS */
    *idim = (i__1 = lngthg[*ngels], abs(i__1));
    y = *x + *idim - 1;
/* Computing MAX */
    i__1 = relpg[*lincon];
    relpg[*lincon] = max(i__1,y);
/*   WRITE HERE */
/*   WRITE CONTIG DESCRIPTOR */
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
    writrn_(idevr, ngels, nconts);
    return 0;

L400:

/*   ADDING TO LEFT END */
/* L410: */
/*   NEED TO EDIT DB HERE */
    if (*itotpc > 0) {
	abedin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, &c__1, seqc2 + 1, itotpc, idc, idbsiz, idevr, 
		maxgel, (ftnlen)1, (ftnlen)1);
    }

/* L420: */
/*   SET POINTERS IN NEW GEL */
    relpg[*ngels] = 1;
    rnbr[*ngels] = *ncontc;
    lnbr[*ngels] = 0;
/*   WRITE NEW GEL LINE */
    writeg_(idevr, ngels, &relpg[*ngels], &lngthg[*ngels], &lnbr[*ngels], &
	    rnbr[*ngels]);
/*   SET POINTERS IN OLD LEFT END */
    lnbr[*ncontc] = *ngels;
    relpg[*ncontc] = *x;
/*   WRITE NEW LEFT END */
    writeg_(idevr, ncontc, &relpg[*ncontc], &lngthg[*ncontc], &lnbr[*ncontc], 
	    &rnbr[*ncontc]);
/*   NEW LENGTH OF CONTIG */
    relpg[*lincon] = relpg[*lincon] + *x - 1;
/*   MAY HAVE JUST ADDED A GEL LONGER THAN CONTIG */
    *idim = (i__1 = lngthg[*ngels], abs(i__1));
    y = *idim;
    if (y > relpg[*lincon]) {
	relpg[*lincon] = y;
    }
/*   NEW NAME OF LEFT END OF CONTIG */
    lnbr[*lincon] = *ngels;
/*   WRITE CONTIG DESCRIPTOR */
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
    writrn_(idevr, ngels, nconts);
/*   NOW GO THRU AND CHANGE ALL RELATIVE POSITIONS */
    n = *ncontc;
L440:
    if (rnbr[n] == 0) {
	return 0;
    }
    n = rnbr[n];
    relpg[n] = relpg[n] + *x - 1;
/*   WRITE NEW LINE */
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
    goto L440;
} /* aenter_ */

/*      ABEDIN */

/*   ROUTINE TO EDIT THE DB USING A PADDED SEQ */
/*   HAVE AN ARRAY SEQC2 LENGTH IDC OF PADDED SECTION OF CONTIG LINCON */
/*  THE LEFT END OF THE PADDED CONTIG STARTS AT X */
/*   THERE ARE ITOTPC PADS TO MAKE */

/* Subroutine */ int abedin_(relpg, lngthg, lnbr, rnbr, ngels, nconts, gel, 
	lincon, x, seqc2, itotpc, idc, idbsiz, idevr, maxgel, gel_len, 
	seqc2_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts;
char *gel;
integer *lincon, *x;
char *seqc2;
integer *itotpc, *idc, *idbsiz, *idevr, *maxgel;
ftnlen gel_len;
ftnlen seqc2_len;
{
    /* Initialized data */

    static char p[1+1] = ",";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ipad, posn, j, idone;
    extern /* Subroutine */ int padcon_(), erromf_();
    static integer iat;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seqc2;
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */

/*   POINT TO CONTIG */
    posn = *x - 1;
/*   POINT TO SEQC2 */
    iat = 0;
/*   COUNT PADS DONE */
    idone = 0;
/*   LOOP FOR ALL SEQC2 */
    i__1 = *idc;
    for (j = 1; j <= i__1; ++j) {
	++posn;
	++iat;
	ipad = 0;
/*   IS THIS A PADDING CHAR? */
	if (*(unsigned char *)&seqc2[iat] != *(unsigned char *)&p[0]) {
	    goto L100;
	}
L50:
/*   COUNT PADS */
	++ipad;
	++iat;
	if (*(unsigned char *)&seqc2[iat] == *(unsigned char *)&p[0]) {
	    goto L50;
	}
/*   END OF THIS STRETCH OF PADS,DO INSERT */
/*   HAVE IPAD INSERTS TO MAKE AT POSN */
/*      WRITE(*,*)'LINCON,POSN,IPAD',LINCON,POSN,IPAD */
	padcon_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, gel 
		+ 1, lincon, &posn, &ipad, idbsiz, idevr, maxgel, (ftnlen)1);
/*   MOVE POINTER TO CONTIG */
	posn += ipad;
/*   COUNT PADS DONE */
	idone += ipad;
/*   ANY MORE TO DO? */
	if (idone == *itotpc) {
	    goto L101;
	}
L100:
	;
    }
/*   ERROR SHOULD HAVE DONE ALL PADS */
    erromf_("Problem: some pads were not done!", (ftnlen)33);
L101:
    ;
} /* abedin_ */

/* Subroutine */ int addtit_(seq1, nampro, ngels, idim1, seq1_len, nampro_len)
char *seq1, *nampro;
integer *ngels, *idim1;
ftnlen seq1_len;
ftnlen nampro_len;
{
    extern /* Subroutine */ int cadtit_();

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seq1;

    /* Function Body */
    cadtit_(seq1 + 1, nampro, ngels, (ftnlen)1, nampro_len);
    *idim1 += 20;
} /* addtit_ */

/* Subroutine */ int adism3_(isavps, savpg, cends, nends, idcend, maxcon, 
	ilefts, ilc, iposc, iposg, isense, llino, imatc, istran, nextc, maxc, 
	jj, isavl, lmatch)
integer *isavps, *savpg, *cends, *nends, *idcend, *maxcon, *ilefts, *ilc, *
	iposc, *iposg, *isense, *llino, *imatc, *istran, *nextc, *maxc, *jj, *
	isavl, *lmatch;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, savps;
    extern /* Subroutine */ int erromf_();
    static integer lcl, lcr;

/*   AUTHOR: RODGER STADEN */
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
    *jj = 1;

/* we have a match isavps, isavpg, isavl (pos in consensus, gel, length) */
/* 1. which contig is it in? --> JJ */
/* 2. save pos of contig, pos in contig, contig length, contig number, sense */

/*      WRITE(*,*)'ENTER ADISM3, IMATC',IMATC */
    i__1 = *idcend;
    for (j = 2; j <= i__1; ++j) {
	if (savps > cends[j]) {
	    goto L5;
	}
	*jj = j - 1;
	goto L6;
L5:
	;
    }
    *jj = *idcend;
L6:
    --savps;
    lcl = savps - cends[*jj];
    lcr = cends[*jj + 1] - *isavps - 1;
    *nextc = cends[*jj + 1] + 20;
    if (*imatc <= *maxc) {
	ilefts[*imatc] = cends[*jj] + 20;
	ilc[*imatc] = lcl + lcr + 1;
	iposc[*imatc] = lcl + 1;
	iposg[*imatc] = *savpg;
	llino[*imatc] = nends[*jj];
	isense[*imatc] = 1;
	if (*istran == 2) {
	    isense[*imatc] = -1;
	}
	*lmatch = *isavl;
/*        WRITE(INFOD,1000)LLINO(IMATC),IPOSC(IMATC),ISTRAN, */
/*     +  IPOSG(IMATC) */
/* 1000   FORMAT */
/*     +  ('Contig',I8,' position',I8,' matches strand',I2, */
/*     +  ' at position',I8) */
/*        CALL INFO(INFOD) */
    } else {
	erromf_("Warning: too many overlaps", (ftnlen)26);
    }
} /* adism3_ */

/* Subroutine */ int adism4_(idim, idimg, savps, savpg, savl, idsav, cends, 
	nends, idcend, maxcon, ilefts, ilc, iposc, iposg, isense, llino, 
	imatc, istran, maxc)
integer *idim, *idimg, *savps, *savpg, *savl, *idsav, *cends, *nends, *idcend,
	 *maxcon, *ilefts, *ilc, *iposc, *iposg, *isense, *llino, *imatc, *
	istran, *maxc;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer lend, i__, lastc, nextc;
    extern /* Subroutine */ int bub3as_(), adism3_();
    static integer lmatch;

/*   AUTHOR: RODGER STADEN */
/*      CHARACTER INFOD*80 */
/*      WRITE(*,*)'ENTER ADISM4    , IMATC',IMATC */
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

    bub3as_(&savps[1], &savpg[1], &savl[1], idsav);
/*      DO 123 II = 1,IDSAV */
/*        WRITE(*,*)II,SAVPS(II),SAVPG(II),SAVL(II) */
/* 123    CONTINUE */
    ++(*imatc);
/*      WRITE(*,*)'IN ADISM4, UPDATED IMATC',IMATC */

/* get the contig info for the first match */

/*        WRITE(*,*)'sav1',SAVL(1) */
/* we have a match savps, savpg, savl (pos in consensus, gel, length) */
/* 2. save pos of contig, pos in contig, contig length, contig number, sense */
    adism3_(&savps[1], &savpg[1], &cends[1], &nends[1], idcend, maxcon, &
	    ilefts[1], &ilc[1], &iposc[1], &iposg[1], &isense[1], &llino[1], 
	    imatc, istran, &nextc, maxc, &lastc, &savl[1], &lmatch);

/* now decide when a match is with a new contig and get the relevant info. */
/* Decide its the same overlap if it is covered by the previous gel position */
/* If we want to record the longest match for each overlap we should test */
/* here to see if overlapping ones are longer than the one weve recorded */

    lend = *idimg - savpg[1] + savps[1];
    i__1 = *idsav;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*         WRITE(*,*)SAVPS(I),SAVPG(I) */
/*         WRITE(*,*)'SAVPS(I)-SAVPG(I)',SAVPS(I)-SAVPG(I) */
/*        WRITE(*,*)'savl(I),lend',SAVL(I),LEND */
	if (savps[i__] < lend && savps[i__] < nextc) {

/* test here if this match is longer */

	    if (savl[i__] > lmatch) {

/* next test added 22-11-94 because the trap in adism3 is insufficient */

		if (*imatc <= *maxc) {
		    iposc[*imatc] = savps[i__] - cends[lastc] - 19;
		    iposg[*imatc] = savpg[i__];
		    lmatch = savl[i__];
/*              WRITE(*,*)'new best g,c,l',IPOSG(IMATC),IPOSC(IMATC),LMATCH */
		}
	    }
	    goto L10;
	}
/*         WRITE(*,*)'2SAVPS(I)-SAVPG(I)',SAVPS(I)-SAVPG(I) */
/*         WRITE(*,*)IPOSC(IMATC),IPOSG(IMATC),SAVPS(I),SAVPG(I) */
	++(*imatc);
/*      WRITE(*,*)'IN ADISM4, UPDATED AGAIN IMATC',IMATC */
/* we have a match savps, savpg, savl (pos in consensus, gel, length) */
/* 2. save pos of contig, pos in contig, contig length, contig number, sense */
	adism3_(&savps[i__], &savpg[i__], &cends[1], &nends[1], idcend, 
		maxcon, &ilefts[1], &ilc[1], &iposc[1], &iposg[1], &isense[1],
		 &llino[1], imatc, istran, &nextc, maxc, &lastc, &savl[i__], &
		lmatch);
	lend = *idimg - savpg[i__] + savps[i__];
/*        RSTART = SAVPS(I) - SAVPG(I) */
L10:
	;
    }
    *imatc = min(*imatc,*maxc);
/*      WRITE(*,*)'IN ADISM4, LAST IMATC',IMATC */
} /* adism4_ */

/* C    AJOIN2 */
/* C   COMPLETES JOIN AND RETURNS LENGTH OF NEW CONTIG IN LLINOR */
/* Subroutine */ int ajoin2_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 relx, llinol, llinor, lnconl, lnconr, idevr)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *relx, *
	llinol, *llinor, *lnconl, *lnconr, *idevr;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer n;
    static real x;
    extern /* Subroutine */ int merge_(), remcnl_(), mrgtag_(), writec_(), 
	    writeg_(), mrgnot_();

/*   AUTHOR: RODGER STADEN */
/*   RELX IS THE POSITION OF THE JOINT */
/*   LLINOL IS THE LEFT GEL NUMBER OF THE LEFT CONTIG */
/*   LLINOR IS THE LEFT GEL OF THE RIGHT CONTIG */
/*   LNCONL IS THE LEFT CONTIG LINE NUMBER */
/*   LNCONR IS THE RIGHT CONTIG LINE NUMBER */

/*   ADJUST ALL RELATIVE POSITIONS IN RIGHT CONTIG */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    n = *llinor;
    relpg[n] = *relx;
L50:
    if (rnbr[n] == 0) {
	goto L60;
    }
    n = rnbr[n];
    relpg[n] = relpg[n] + *relx - 1;
    goto L50;
L60:

/*   FIX UP NEW GEL LINE FOR OLD LEFT OF RIGHT CONTIG */
    lnbr[*llinor] = rnbr[*lnconl];
/*   FIX UP RIGHT GEL OF LEFT CONTIG */
    n = rnbr[*lnconl];
    rnbr[n] = *llinor;
/*   MERGE WILL SORT OUT THE CORRECT NEIGHBOURS */

    merge_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], lnconl, idbsiz);
/*   MERGE DOES NOT WRITE TO DISK */
    n = lnbr[*lnconl];
L65:
/*      WRITE(IDEVR,REC=N)RELPG(N),LNGTHG(N),LNBR(N),RNBR(N) */
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
    n = rnbr[n];
    if (n != 0) {
	goto L65;
    }
/*   Merge annotation lists. */
    i__1 = *idbsiz - *lnconl;
    i__2 = *idbsiz - *lnconr;
    i__3 = *relx - 1;
    mrgtag_(idevr, &i__1, &i__2, &i__3);
    i__1 = *idbsiz - *lnconr;
    i__2 = *idbsiz - *lnconl;
    mrgnot_(idevr, &i__1, &i__2);
/*   CONTIG LINES */
    x = (real) (relpg[*lnconr] + *relx - 1);
/*   LENGTH MAY NOT HAVE INCREASED! */
    if (x > (real) relpg[*lnconl]) {
	relpg[*lnconl] = x;
    }
/*   SAVE LENGTH OF NEW CONTIG */
    *relx = relpg[*lnconl];
/*      WRITE(IDEVR,REC=LNCONL)RELPG(LNCONL),LNGTHG(LNCONL),LNBR(LNCONL), */
/*     1RNBR(LNCONL) */
    i__1 = *idbsiz - *lnconl;
    writec_(idevr, &i__1, &relpg[*lnconl], &lnbr[*lnconl], &rnbr[*lnconl]);
/*   Now remove the old contig. We must use the C routine for this so that */
/*   it can update the contig order, tag lists, etc. */
    remcnl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    lnconr, idevr);
    return 0;
} /* ajoin2_ */

/*     SUBROUTINE AJOIN3 */
/* Subroutine */ int ajoin3_(relpg, idbsiz, lincon, itype, isense, joint, 
	idim22, klass, iover, pl, pr)
integer *relpg, *idbsiz, *lincon, *itype, *isense, *joint, *idim22, *klass, *
	iover, *pl, *pr;
{
    static integer i__;

/*   AUTHOR: RODGER STADEN */

/*   CALC POSITIONS OF CONTIGS RELATIVE TO FIXED GEL */
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
	if (itype[i__] != -1 || isense[i__] != 1) {
	    goto L11;
	}
	pl[i__] = -joint[i__] + 2;
	pr[i__] = pl[i__] + relpg[lincon[i__]] - 1;
	goto L20;
/*   L+ */
L11:
	if (itype[i__] != 1 || isense[i__] != 1) {
	    goto L12;
	}
	pl[i__] = joint[i__];
	pr[i__] = pl[i__] + relpg[lincon[i__]] - 1;
	goto L20;
/*   R- */
L12:
	if (itype[i__] != -1 || isense[i__] != -1) {
	    goto L13;
	}
	pr[i__] = joint[i__] + idim22[i__] - 1;
	pl[i__] = pr[i__] - relpg[lincon[i__]] + 1;
	goto L20;
/*   L- */
L13:
	pr[i__] = idim22[i__] - joint[i__] + 1;
	pl[i__] = pr[i__] - relpg[lincon[i__]] + 1;
L20:
	;
    }
/*  LENGTH OF OVERLAP */
    *iover = min(pr[1],pr[2]) - max(pl[1],pl[2]) + 1;

/*  CLASS NUMBER 1-16 */
    *klass = 1;
    if (itype[1] == 1) {
	*klass += 8;
    }
    if (isense[1] == -1) {
	*klass += 4;
    }
    if (itype[2] == 1) {
	*klass += 2;
    }
    if (isense[2] == -1) {
	++(*klass);
    }
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

/* Subroutine */ int aline_(seq1, seq2, seqg2, seqc2, isav1, isav2, isav3, 
	idsav, idc, idim2, idout, ic1, ig1, minsli, joint, itotpc, itotpg, 
	ifail, itype, maxpc, maxpg, permax, seq3, maxgel, percm, leno, ishow, 
	mask, jrorc, seq1_len, seq2_len, seqg2_len, seqc2_len, seq3_len)
char *seq1, *seq2, *seqg2, *seqc2;
integer *isav1, *isav2, *isav3, *idsav, *idc, *idim2, *idout, *ic1, *ig1, *
	minsli, *joint, *itotpc, *itotpg, *ifail, *itype, *maxpc, *maxpg;
real *permax;
char *seq3;
integer *maxgel;
real *percm;
integer *leno, *ishow, *mask, *jrorc;
ftnlen seq1_len;
ftnlen seq2_len;
ftnlen seqg2_len;
ftnlen seqc2_len;
ftnlen seq3_len;
{
    extern /* Subroutine */ int maskc_();
    static integer idim2i;
    extern /* Subroutine */ int bub3as_(), dalign_(), tpchek_(), upchek_(), 
	    slides_(), lineup_(), removl_();
    static integer minslt;
    extern /* Subroutine */ int mstlkl_(), sqcopy_();
    static integer ipp;

/*   AUTHOR: RODGER STADEN */
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
    idim2i = *idim2;

/* need to unmask both(for contig joins) sequences */

/*        CALL FMTDB(SEQ2,IDIM2,1,IDIM2,60,6) */
    if (*mask != 0) {
/*        WRITE(*,*)'SEQ1 B' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
	maskc_(seq1 + 1, idc, &c__2, (ftnlen)1);
/*        WRITE(*,*)'SEQ1 A' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
/*        WRITE(*,*)'SEQ2 B' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2) */
	maskc_(seq2 + 1, idim2, &c__2, (ftnlen)1);
/*        WRITE(*,*)'SEQ2 A' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2) */
    }
/*        CALL FMTDB(SEQ2,IDIM2,1,IDIM2,60,6) */
/*   SAVE SEQ2 */
    sqcopy_(seq2 + 1, seq3 + 1, idim2, (ftnlen)1, (ftnlen)1);
    mstlkl_(seq3 + 1, idim2, (ftnlen)1);
/*        CALL FMTDB(SEQ3,IDIM2,1,IDIM2,60,6) */
    *ifail = 1;
/*   FIND MATCHES */
    ipp = *idsav;
/*      WRITE(*,*)'IC1,IG1',IC1,IG1,MAXPG,MAXPC,MINSLT */
    slides_(seq1 + 1, idc, seq3 + 1, idim2, ic1, ig1, maxpg, maxpc, &minslt, &
	    isav1[1], &isav2[1], &isav3[1], &ipp, (ftnlen)1, (ftnlen)1);
/*      WRITE(*,*)'IPP',IPP,IDSAV */
    if (ipp > *idsav) {
	goto L50;
    }
    if (ipp < 1) {
	goto L50;
    }
    removl_(&isav2[1], &isav3[1], &isav1[1], &ipp);
/*      WRITE(*,*)'IPP',IPP,IDSAV */
    bub3as_(&isav2[1], &isav3[1], &isav1[1], &ipp);
/*   DO TOPOLOGICAL CHECK */
    tpchek_(&isav2[1], &isav3[1], &isav1[1], &ipp);

/* added next routine 27-2-93 */

/*      WRITE(*,*)'IPP',IPP,IDSAV */
    upchek_(&isav2[1], &isav3[1], &isav1[1], &ipp);
/*      WRITE(*,*)'IPP',IPP,IDSAV */
    lineup_(seq2 + 1, seq1 + 1, seqg2 + 1, seqc2 + 1, idc, idim2, idout, &
	    isav3[1], &isav2[1], &isav1[1], &ipp, itotpc, itotpg, joint, 
	    itype, maxgel, ifail, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/*       WRITE(*,*)'ITOTPC,ITOTPG',ITOTPC,ITOTPG,IFAIL */
/*      IF(ITOTPC.GT.MAXPC)IFAIL=1 */
/*      IF(ITOTPG.GT.MAXPG)IFAIL=1 */
    if (*ifail != 0) {
	goto L50;
    }
/*   IDIM2 IS NOW LENGTH OF ALIGNED GEL */
    dalign_(seqc2 + 1, seqg2 + 1, seq3 + 1, maxgel, idout, idim2, joint, 
	    itype, percm, ifail, leno, permax, ishow, maxpg, maxpc, itotpg, 
	    itotpc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
L50:

/* need to remask both(for contig joins) sequences */

    if (*mask != 0) {
/*        WRITE(*,*)'SEQ1 B1' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
	maskc_(seq1 + 1, idc, &c__3, (ftnlen)1);
/*        WRITE(*,*)'SEQ1 A1' */
/*        WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */

/* only mask consensus data which is only lowercase where marked */

/*        WRITE(*,*)'SEQ2 B1' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2I) */
	if (*jrorc == 1) {
	    maskc_(seq2 + 1, &idim2i, &c__3, (ftnlen)1);
	}
/*        WRITE(*,*)'SEQ2 A1' */
/*        WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2I) */
    }
} /* aline_ */

/* Subroutine */ int autocn_(seq1, idim, gel, idimg, ilefts, ilc, iposc, 
	iposg, isense, llino, imatc, ifcomp, minmat, maxgel, maxglm, gelcop, 
	savps, savpg, savl, maxsav, cends, nends, maxcon, seqg2, seqc2, seq4, 
	idout, idim22, itotpg, itotpc, joint, ifail, itype, maxpc, maxpg, 
	permax, minsli, seqg3, seqc3, kfail, jobc, permis, leno, ishow, mask, 
	minovr, seq1_len, gel_len, gelcop_len, seqg2_len, seqc2_len, seq4_len,
	 seqg3_len, seqc3_len)
char *seq1;
integer *idim;
char *gel;
integer *idimg, *ilefts, *ilc, *iposc, *iposg, *isense, *llino, *imatc, *
	ifcomp, *minmat, *maxgel, *maxglm;
char *gelcop;
integer *savps, *savpg, *savl, *maxsav, *cends, *nends, *maxcon;
char *seqg2, *seqc2, *seq4;
integer *idout, *idim22, *itotpg, *itotpc, *joint, *ifail, *itype, *maxpc, *
	maxpg;
real *permax;
integer *minsli;
char *seqg3, *seqc3;
integer *kfail, *jobc;
real *permis;
integer *leno, *ishow, *mask, *minovr;
ftnlen seq1_len;
ftnlen gel_len;
ftnlen gelcop_len;
ftnlen seqg2_len;
ftnlen seqc2_len;
ftnlen seq4_len;
ftnlen seqg3_len;
ftnlen seqc3_len;
{
    /* Format strings */
    static char fmt_1000[] = "(\002Total matches found\002,i6)";
    static char fmt_1002[] = "(\002Contig\002,i8,\002 position\002,i8,\002 m\
atches strand \002,i2,\002 at position\002,i8)";
    static char fmt_1001[] = "(\002Trying to align with contig \002,i8)";

    /* System generated locals */
    integer seqg2_dim1, seqg2_offset, seqc2_dim1, seqc2_offset, i__1, i__2, 
	    i__3, i__4, i__5;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    static char csen[1];
    extern /* Subroutine */ int info_();
    static integer i__, jfail, jdim22;
    extern /* Subroutine */ int aline_();
    static integer jmatc;
    static char infod[80];
    static integer idsav, jposc[100], jposg[100], ioptc;
    extern /* Subroutine */ int sqcom_();
    static integer jdout;
    extern /* Subroutine */ int busyf_();
    static integer jtype, start;
    static real perms;
    extern /* Subroutine */ int copym_(), sqrev_(), adism4_();
    static integer idcend;
    extern /* Subroutine */ int fndcon_();
    static integer jlefts[100], jsense[100], jllino[100];
    extern integer cmpseq_();
    static integer istran, ksense;
    extern /* Subroutine */ int mstlkl_();
    static integer lenovr, jjoint, jtotpc, jtotpg;
    extern /* Subroutine */ int sqcopy_();
    static integer jlc[100];

    /* Fortran I/O blocks */
    static icilist io___132 = { 0, infod, 0, fmt_1000, 80, 1 };
    static icilist io___134 = { 0, infod, 0, fmt_1002, 80, 1 };
    static icilist io___139 = { 0, infod, 0, fmt_1001, 80, 1 };


/*   AUTHOR: RODGER STADEN */
/*   changed 29-11-90 to make first in list of alignments the best */

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
    seqc2_offset = 1 + seqc2_dim1 * 1;
    seqc2 -= seqc2_offset;
    seqg2_dim1 = *maxglm;
    seqg2_offset = 1 + seqg2_dim1 * 1;
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

    ifail[1] = 1;
    ifail[2] = 1;
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
    sqcopy_(gel + 1, gelcop + 1, idimg, (ftnlen)1, (ftnlen)1);
/*  COUNT NUMBER OF CONTIGS THAT MATCH */
    *imatc = 0;
    idcend = *maxcon;
    busyf_();
/*      WRITE(*,*)'IDIM',IDIM,IDCEND */
    fndcon_(seq1 + 1, idim, &cends[1], &nends[1], &idcend, maxcon, (ftnlen)1);
    if (*jobc != 0) {
	start = 1;
	if (*jobc == 1) {
	    start = cends[idcend];
	}
	ioptc = 2;
	*ifcomp = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1]
		, &idsav, seq1 + 1, gel + 1, idim, idimg, (ftnlen)1, (ftnlen)
		1, (ftnlen)1);
	if (*ifcomp != 0) {
	    return 0;
	}
    }
/* L1: */
    istran = 1;
L2:
    mstlkl_(gel + 1, idimg, (ftnlen)1);
    idsav = *maxsav;
    ioptc = 3;
    *ifcomp = cmpseq_(&ioptc, csen, minmat, &savps[1], &savpg[1], &savl[1], &
	    idsav, seq1 + 1, gel + 1, idim, idimg, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1);
    if (*ifcomp < 0) {
	return 0;
    }
    idsav = *ifcomp;
    if (idsav != 0) {
	adism4_(idim, idimg, &savps[1], &savpg[1], &savl[1], &idsav, &cends[1]
		, &nends[1], &idcend, maxcon, jlefts, jlc, jposc, jposg, 
		jsense, jllino, imatc, &istran, &c__100);
    }
    ++istran;
    if (istran == 2) {
	sqcopy_(gelcop + 1, gel + 1, idimg, (ftnlen)1, (ftnlen)1);
	sqrev_(gel + 1, idimg, (ftnlen)1);
	sqcom_(gel + 1, idimg, (ftnlen)1);
	goto L2;
    }
    sqcopy_(gelcop + 1, gel + 1, idimg, (ftnlen)1, (ftnlen)1);
    ksense = 0;
    s_wsfi(&io___132);
    do_fio(&c__1, (char *)&(*imatc), (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)80);
    if (*imatc == 0) {
	ifail[1] = 0;
	*ifcomp = 0;
	return 0;
    }
    i__1 = *imatc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfi(&io___134);
	do_fio(&c__1, (char *)&jllino[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jposc[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jsense[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jposg[i__ - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
/* L99: */
    }
    jmatc = 0;
    i__1 = *imatc;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*  3-10-95 New idea! have minimum overlap before allowing entry. */
/*  Simplest place to apply is here (though not where it should */
/*  be if we had started afresh) */
/* Computing MIN */
	i__2 = jposg[i__ - 1], i__3 = jposc[i__ - 1];
/* Computing MIN */
	i__4 = *idimg - jposg[i__ - 1], i__5 = jlc[i__ - 1] - jposc[i__ - 1];
	lenovr = min(i__2,i__3) + min(i__4,i__5);
	if (lenovr < *minovr) {
/*         WRITE(*,*)'SHORT OVERLAP', */
/*     +   LENOVR,JPOSG(I),IDIMG,JPOSC(I),JLC(I) */
	    goto L100;
	}



/*         WRITE(*,*)'*******LONG OVERLAP', */
/*     +   LENOVR,JPOSG(I),IDIMG,JPOSC(I),JLC(I) */
	if (jsense[i__ - 1] == -1) {
	    if (ksense == 0) {
		sqrev_(gel + 1, idimg, (ftnlen)1);
		sqcom_(gel + 1, idimg, (ftnlen)1);
		ksense = 1;
	    }
	}
	jdim22 = *idimg;
	jdout = *maxgel;
	idsav = *maxsav;
	s_wsfi(&io___139);
	do_fio(&c__1, (char *)&jllino[i__ - 1], (ftnlen)sizeof(integer));
	e_wsfi();
	info_(infod, (ftnlen)80);
	aline_(seq1 + jlefts[i__ - 1], gel + 1, seqg3 + 1, seqc3 + 1, &savps[
		1], &savpg[1], &savl[1], &idsav, &jlc[i__ - 1], &jdim22, &
		jdout, &jposc[i__ - 1], &jposg[i__ - 1], minsli, &jjoint, &
		jtotpc, &jtotpg, &jfail, &jtype, maxpc, maxpg, permax, seq4 + 
		1, maxgel, &perms, leno, ishow, mask, &c__0, (ftnlen)1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
	if (jfail == 0) {
	    ++jmatc;
	    if (jmatc == 1) {
/*    Save in elements 1 */
		copym_(&jlefts[i__ - 1], &ilefts[1], &jlc[i__ - 1], &ilc[1], &
			jposc[i__ - 1], &iposc[1], &jsense[i__ - 1], &isense[
			1], &jllino[i__ - 1], &llino[1], &jjoint, &joint[1], &
			jtotpc, &itotpc[1], &jtotpg, &itotpg[1], &jtype, &
			itype[1], &jdout, &idout[1], &jdim22, &idim22[1], 
			seqg3 + 1, seqg2 + (seqg2_dim1 + 1), seqc3 + 1, seqc2 
			+ (seqc2_dim1 + 1), &perms, &permis[1], (ftnlen)1, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
		ifail[1] = 0;
	    } else if (jmatc == 2) {
		if (perms < permis[1]) {
/*    Better match so save in elements 1, so copy 1 to 2 first */
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
		    ifail[2] = 0;
/*    Now save in 1 */
		    copym_(&jlefts[i__ - 1], &ilefts[1], &jlc[i__ - 1], &ilc[
			    1], &jposc[i__ - 1], &iposc[1], &jsense[i__ - 1], 
			    &isense[1], &jllino[i__ - 1], &llino[1], &jjoint, 
			    &joint[1], &jtotpc, &itotpc[1], &jtotpg, &itotpg[
			    1], &jtype, &itype[1], &jdout, &idout[1], &jdim22,
			     &idim22[1], seqg3 + 1, seqg2 + (seqg2_dim1 + 1), 
			    seqc3 + 1, seqc2 + (seqc2_dim1 + 1), &perms, &
			    permis[1], (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			    ftnlen)1);
		} else {
/*    Save in element 2 */
		    copym_(&jlefts[i__ - 1], &ilefts[2], &jlc[i__ - 1], &ilc[
			    2], &jposc[i__ - 1], &iposc[2], &jsense[i__ - 1], 
			    &isense[2], &jllino[i__ - 1], &llino[2], &jjoint, 
			    &joint[2], &jtotpc, &itotpc[2], &jtotpg, &itotpg[
			    2], &jtype, &itype[2], &jdout, &idout[2], &jdim22,
			     &idim22[2], seqg3 + 1, seqg2 + ((seqg2_dim1 << 1)
			     + 1), seqc3 + 1, seqc2 + ((seqc2_dim1 << 1) + 1),
			     &perms, &permis[2], (ftnlen)1, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
		    ifail[2] = 0;
		}
	    } else {
		if (perms < permis[1]) {
/*    Better match so save in elements 1, so copy 1 to 2 first */
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
		    ifail[2] = 0;
/*    Now save in 1 */
		    copym_(&jlefts[i__ - 1], &ilefts[1], &jlc[i__ - 1], &ilc[
			    1], &jposc[i__ - 1], &iposc[1], &jsense[i__ - 1], 
			    &isense[1], &jllino[i__ - 1], &llino[1], &jjoint, 
			    &joint[1], &jtotpc, &itotpc[1], &jtotpg, &itotpg[
			    1], &jtype, &itype[1], &jdout, &idout[1], &jdim22,
			     &idim22[1], seqg3 + 1, seqg2 + (seqg2_dim1 + 1), 
			    seqc3 + 1, seqc2 + (seqc2_dim1 + 1), &perms, &
			    permis[1], (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			    ftnlen)1);
		} else if (perms < permis[2]) {
/*    Save in element 2 */
		    copym_(&jlefts[i__ - 1], &ilefts[2], &jlc[i__ - 1], &ilc[
			    2], &jposc[i__ - 1], &iposc[2], &jsense[i__ - 1], 
			    &isense[2], &jllino[i__ - 1], &llino[2], &jjoint, 
			    &joint[2], &jtotpc, &itotpc[2], &jtotpg, &itotpg[
			    2], &jtype, &itype[2], &jdout, &idout[2], &jdim22,
			     &idim22[2], seqg3 + 1, seqg2 + ((seqg2_dim1 << 1)
			     + 1), seqc3 + 1, seqc2 + ((seqc2_dim1 << 1) + 1),
			     &perms, &permis[2], (ftnlen)1, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
		}
	    }
	} else {
	    *kfail = 1;
	}
L100:
	;
    }
    *imatc = min(2,jmatc);
    *ifcomp = 0;
} /* autocn_ */

/*     BUBBL3 */
/*   SUBROUTINE TO SORT INTEGER ARRAY (LIST) INTO ASCENDING  ORDER */

/* Subroutine */ int ccta_(seq, id, seq_len)
char *seq;
integer *id;
ftnlen seq_len;
{
    /* Initialized data */

    static char com[1+1] = ",";
    static char as[1+1] = "*";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --seq;

    /* Function Body */
    i__1 = *id;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&seq[i__] == *(unsigned char *)&com[0]) {
	    *(unsigned char *)&seq[i__] = *(unsigned char *)&as[0];
	}
/* L10: */
    }
} /* ccta_ */

/* Subroutine */ int copym_(jlefts, ilefts, jlc, ilc, jposc, iposc, jsense, 
	isense, jllino, llino, jjoint, joint, jtotpc, itotpc, jtotpg, itotpg, 
	jtype, itype, jdout, idout, jdim22, idim22, seqg3, seqg2, seqc3, 
	seqc2, perms, permis, seqg3_len, seqg2_len, seqc3_len, seqc2_len)
integer *jlefts, *ilefts, *jlc, *ilc, *jposc, *iposc, *jsense, *isense, *
	jllino, *llino, *jjoint, *joint, *jtotpc, *itotpc, *jtotpg, *itotpg, *
	jtype, *itype, *jdout, *idout, *jdim22, *idim22;
char *seqg3, *seqg2, *seqc3, *seqc2;
real *perms, *permis;
ftnlen seqg3_len;
ftnlen seqg2_len;
ftnlen seqc3_len;
ftnlen seqc2_len;
{
    extern /* Subroutine */ int sqcopy_();

    /* Parameter adjustments */
    --seqc2;
    --seqc3;
    --seqg2;
    --seqg3;

    /* Function Body */
    *ilefts = *jlefts;
    *ilc = *jlc;
    *iposc = *jposc;
    *isense = *jsense;
    *llino = *jllino;
    *joint = *jjoint;
    *itotpc = *jtotpc;
    *itotpg = *jtotpg;
    *itype = *jtype;
    *idout = *jdout;
    *idim22 = *jdim22;
    sqcopy_(seqg3 + 1, seqg2 + 1, jdim22, (ftnlen)1, (ftnlen)1);
    sqcopy_(seqc3 + 1, seqc2 + 1, jdout, (ftnlen)1, (ftnlen)1);
    *permis = *perms;
} /* copym_ */

/*     SUBROUTINE DALIGN */

/*   COUNTS MISMATCHES AND DISPLAYS OVERLAP. */
/* Subroutine */ int dalign_(seqc2, seqg2, seq3, maxgel, idout, idim2, joint, 
	itype, x, ifail, lo, permax, ishow, maxpg, maxpc, itotpg, itotpc, 
	seqc2_len, seqg2_len, seq3_len)
char *seqc2, *seqg2, *seq3;
integer *maxgel, *idout, *idim2, *joint, *itype;
real *x;
integer *ifail, *lo;
real *permax;
integer *ishow, *maxpg, *maxpc, *itotpg, *itotpc;
ftnlen seqc2_len;
ftnlen seqg2_len;
ftnlen seq3_len;
{
    /* Format strings */
    static char fmt_1052[] = "(\002Percent mismatch \002,f4.1,\002, pads in \
contig\002,i3,\002, pads in gel\002,i3)";
    static char fmt_1000[] = "(a)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi(), i_len();

    /* Local variables */
    extern /* Subroutine */ int info_();
    static char name1[15], name2[15];
    static integer i__, j, k, iendc, iendg;
    static real y;
    static char infod[80];
    extern integer forta_();
    static integer kc, lg;
    extern /* Subroutine */ int erromf_();
    extern integer ctonum_();
    extern /* Subroutine */ int mstlkl_(), sqcopy_();

    /* Fortran I/O blocks */
    static icilist io___157 = { 0, infod, 0, fmt_1052, 80, 1 };
    static icilist io___158 = { 0, infod, 0, fmt_1052, 80, 1 };
    static icilist io___159 = { 0, infod, 0, fmt_1052, 80, 1 };
    static icilist io___161 = { 0, name2, 0, fmt_1000, 15, 1 };
    static icilist io___163 = { 0, name1, 0, fmt_1000, 15, 1 };


/*   AUTHOR: RODGER STADEN */
/*      CHARACTER PAD,DASH */
/*      SAVE PAD,DASH */
/*      DATA PAD,DASH/',','-'/ */
    /* Parameter adjustments */
    --seq3;
    --seqg2;
    --seqc2;

    /* Function Body */
    iendg = 1;
    iendc = *joint;
/*   ONLY LOOK AT OVERLAP WHICH IS FROM JOINT FOR LEFT TYPE JOIN */
    if (*itype == 1) {
	iendg = *joint;
	iendc = 1;
    }
/* L100: */
/*   LENGTH OF OVERLAP? */
    lg = *idim2 - iendg + 1;
    *lo = min(*idout,lg);
/*   SAVE RAW DATA */
    sqcopy_(seqg2 + 1, seq3 + 1, idim2, (ftnlen)1, (ftnlen)1);
    mstlkl_(seq3 + 1, idim2, (ftnlen)1);
    *x = (real) (*lo);
    y = *x;
    k = iendg + *lo - 1;
/*   POINT TO CONSENSUS */
    j = 0;
/*   CHECK FOR OVERFLOW */
    if (k > *maxgel) {
	erromf_("DALIGN: matching region too long", (ftnlen)32);
	*ifail = 1;
	return 0;
    }
    i__1 = k;
    for (i__ = iendg; i__ <= i__1; ++i__) {
	++j;
	if (ctonum_(seqc2 + j, (ftnlen)1) == ctonum_(seq3 + i__, (ftnlen)1)) {
	    if (ctonum_(seqc2 + j, (ftnlen)1) < 5) {
		goto L200;
	    }
/*         SO FAR THEY ARE = AND ACGT, WHAT IS LEFT IS = AND 5 */
	    if ((*(unsigned char *)&seqc2[j] == '*' || *(unsigned char *)&
		    seqc2[j] == ',') && (*(unsigned char *)&seq3[i__] == '*' 
		    || *(unsigned char *)&seq3[i__] == ',')) {
		goto L200;
	    }
	}
	*x += (float)-1.;
L200:
	;
    }
    *x = (y - *x) * (float)100. / y;
    *ifail = 0;
    if (*x > *permax) {
	*ifail = 1;
    }
    if (*itotpc > *maxpc) {
	*ifail = 1;
    }
    if (*itotpg > *maxpg) {
	*ifail = 1;
    }

/* ISHOW 1 hide all alignments */
/*       2 show passes */
/*       3 show all alignments */
/*       4 show failures only */

/*          WRITE(*,*)X,ITOTPC,ITOTPG */
    if (*ishow == 1) {
	if (*ifail == 0) {
	    s_wsfi(&io___157);
	    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*itotpc), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itotpg), (ftnlen)sizeof(integer));
	    e_wsfi();
	    info_(infod, (ftnlen)80);
	}
	return 0;
    } else if (*ishow == 2) {
	if (*ifail != 0) {
	    return 0;
	}
    } else if (*ishow == 4) {
	if (*ifail == 0) {
	    s_wsfi(&io___158);
	    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*itotpc), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itotpg), (ftnlen)sizeof(integer));
	    e_wsfi();
	    info_(infod, (ftnlen)80);
	    return 0;
	}
    }
    s_wsfi(&io___159);
    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&(*itotpc), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*itotpg), (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___161);
    do_fio(&c__1, "     Consensus", (ftnlen)14);
    e_wsfi();
    s_wsfi(&io___163);
    do_fio(&c__1, "       Reading", (ftnlen)14);
    e_wsfi();
    i__1 = i_len(name1, (ftnlen)15);
    kc = forta_(seqc2 + 1, seqg2 + iendg, lo, name2, name1, &i__1, &iendc, &
	    iendg, infod, &c__80, (ftnlen)1, (ftnlen)1, (ftnlen)15, (ftnlen)
	    15, (ftnlen)80);
} /* dalign_ */

/*     DELCON */

/*   DELETES CONTIG FROM CONSENSUS SEQUENCE */
/* Subroutine */ int delcon_(seq1, ileft, ilc, idim1, seq1_len)
char *seq1;
integer *ileft, *ilc, *idim1;
ftnlen seq1_len;
{
    static integer i1, i2, id;
    extern /* Subroutine */ int sqcopy_();

/*   AUTHOR: RODGER STADEN */
/*   FIRST CHAR TO REPLACE */
    /* Parameter adjustments */
    --seq1;

    /* Function Body */
    i1 = *ileft - 20;
/*   FIRST CHAR TO MOVE */
    i2 = *ileft + *ilc;
/*   IS THIS RIGHTMOST CONTIG ANYWAY? */
    if (i2 > *idim1) {
	goto L10;
    }
/*   NUMBER TO MOVE */
    id = *idim1 - i2 + 1;
/*   MOVE */
    sqcopy_(seq1 + i2, seq1 + i1, &id, (ftnlen)1, (ftnlen)1);
/*   RESET LENGTH */
    *idim1 = i1 + id - 1;
    return 0;
L10:
/*   RIGHTMOST CONTIG SO DONT MOVE */
    *idim1 = i1 - 1;

    return 0;
} /* delcon_ */

/*     LINEUP */

/*   TAKES 2 SEQS SET OF MATCHES AND PRODUCES LINED UP SEQS */
/*   FINDS IF WE HAVE A LEFT OVERLAP */
/*   RETURNS POSITION OF JOINT. THIS IS RELATIVE TO THE CONTIG */
/*   FOR MOST MATCHES BUT I RELATIVE TO THE GEL FOR A LEFT OVERLAP */
/* Subroutine */ int lineup_(seqg, seqc, seqg2, seqc2, idc, idg, idout, matg, 
	matc, matl, ip, itotpc, itotpg, joint, itype, maxgel, ifail, seqg_len,
	 seqc_len, seqg2_len, seqc2_len)
char *seqg, *seqc, *seqg2, *seqc2;
integer *idc, *idg, *idout, *matg, *matc, *matl, *ip, *itotpc, *itotpg, *
	joint, *itype, *maxgel, *ifail;
ftnlen seqg_len;
ftnlen seqc_len;
ftnlen seqg2_len;
ftnlen seqc2_len;
{
    /* Initialized data */

    static char pad[1+1] = ",";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, l, m, nmtch, l5;
    extern /* Subroutine */ int padcop_(), erromf_();
    static integer ic1, ic2, lc1, ig1, ig2, lg1, lg2, lc2, is1, is2;
    extern /* Subroutine */ int sqcopy_();

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seqc;
    --seqg;
    --seqc2;
    --seqg2;
    --matl;
    --matc;
    --matg;

    /* Function Body */
    *ifail = 0;
/*   ZERO PADDING CHARS IN CONTIG (GEL DONE AT END BY DIFFERENCE */
/*   IN INPUT AND OUTPUT LENGTHS) */
    *itotpc = 0;
/*   FILL OUTPUT WITH PADDING */
    i__1 = *idout;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)&seqg2[i__] = *(unsigned char *)&pad[0];
	*(unsigned char *)&seqc2[i__] = *(unsigned char *)&pad[0];
/* L10: */
    }
    nmtch = 0;
/*   SET INITIAL POINTERS TO OUTPUT */
/*   CONSENSUS */
    is1 = 1;
/*   GEL */
    is2 = 1;
/*   FIND DISTANCE FROM LEFT MATCH IN GEL TO LEFT OF GEL */
    ig2 = matg[1] - 1;
    if (ig2 == 0) {
/*       THE LEFT END OF THE GEL MATCHES SO THIS IS NOT A LEFT OVERLAP */
/*       SET TYPE */
	*itype = -1;
/*       SET JOINT */
	*joint = matc[1];
/*       SKIP NEXT SECTION */
	goto L50;
    }
/*   FIND DISTANCE FROM LEFT MATCH IN CONTIG TO LEFT OF CONTIG */
    ic2 = matc[1] - 1;
/*   GET DISTANCE FROM FIRST MATCH IN CONTIG TO FIRST MATCH IN GEL. */
/*   IF THIS DISTANCE <0 THEN WE HAVE A LEFT OVERLAP */
    ic1 = ic2 - ig2 + 1;
    if (ic1 > 0) {
/*       THIS IS NOT A LEFT OVERLAP */
/*       SET TYPE */
	*itype = -1;
/*       SET LEFT END */
	*joint = ic1;
/*       COPY THE GEL UPTO THE FIRST MATCH, INTO THE OUTPUT ARRAY */
/*       CHECK FOR OVERFLOW */
	if (ig2 > *maxgel) {
	    goto L700;
	}
	sqcopy_(seqg + 1, seqg2 + 1, &ig2, (ftnlen)1, (ftnlen)1);
/*       COPY THE CONTIG FOR THE SAME REGION */
	if (ig2 > *maxgel) {
	    goto L700;
	}
	sqcopy_(seqc + ic1, seqc2 + 1, &ig2, (ftnlen)1, (ftnlen)1);
	is1 += ig2;
	is2 += ig2;
	goto L50;
    }
/*   MUST BE LEFT END OVERLAP */
/*   SET TYPE */
    *itype = 1;
/*   SET POSITION OF JOINT RELATIVE TO GEL */
    *joint = abs(ic1) + 2;
/*   COPY OVER THE GEL UPTO THE JOINT */
/*   CHECK FOR OVERFLOW */
    if (ig2 > *maxgel) {
	goto L700;
    }
    sqcopy_(seqg + 1, seqg2 + 1, &ig2, (ftnlen)1, (ftnlen)1);
    is2 += ig2;
/*   WE MAY ALSO HAVE MISMATCHING */
/*   DATA AT THE JOIN SO DEAL WITH THAT NOW */
/*   IF IC2 >0 THE LEFT END OF THE CONTIG MATCHES THE GEL BUT OTHERWISE */
/*   WE HAVE SOME MISMATCHED DATA TO DEAL WITH - WE NEED TO TRANSFER */
/*   THE MISMATCHED REGION OF THE CONTIG TO THE OUTPUT ARRAY */
    if (ic2 > 0) {
	if (ic2 > *maxgel) {
	    goto L700;
	}
	sqcopy_(seqc + 1, seqc2 + 1, &ic2, (ftnlen)1, (ftnlen)1);
	is1 += ic2;
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

L50:
/*   POINT TO NEXT MATCH */
    ++nmtch;
/*   COPY NEXT MATCH */
    ig1 = matg[nmtch];
    ic1 = matc[nmtch];
    l = matl[nmtch];
/*   CHECK FOR OVERFLOW */
    if (is2 + l - 1 > *maxgel) {
	goto L700;
    }
    sqcopy_(seqg + ig1, seqg2 + is2, &l, (ftnlen)1, (ftnlen)1);
/*   CHECK FOR OVERFLOW */
    if (is1 + l - 1 > *maxgel) {
	goto L700;
    }
    sqcopy_(seqc + ic1, seqc2 + is1, &l, (ftnlen)1, (ftnlen)1);
/*   POINT TO NEXT OUTPUT POSITIONS */
    is1 += l;
    is2 += l;
/*   END OF CURRENT MATCH */
    lg1 = ig1 + l;
    lc1 = ic1 + l;
/*   ANY MORE MATCHES */
    if (nmtch == *ip) {
	goto L500;
    }
    k = nmtch + 1;
    lg2 = matg[k] - lg1;
    lc2 = matc[k] - lc1;
/*   ANY DIFFERENCE IN LENGTH? IF SO WE HAVE TO PAD SO THEY BECOME THE SAME */
    l5 = (i__1 = lg2 - lc2, abs(i__1));
/*   COUNT PADDING CHARS IN CONTIG */
    if (lg2 > lc2) {
	*itotpc += l5;
    }
/*   IF DIFFERENCE INCREMENT SHORTER */
    if (lg2 > lc2) {
	is1 += l5;
    }
/*   IF GEL NEEDS PADDING TRY TO PUT PADS NEXT TO DOUBLE CODES */
    if (lc2 > lg2) {
	padcop_(seqg + 1, seqg2 + 1, &lg1, &matg[k], &l5, &is2, &lg2, maxgel, 
		ifail, seqc + 1, idc, &lc1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    }
/*   CHECK FOR OVERFLOW */
    if (*ifail == 1) {
	goto L700;
    }
/*   NOW COPY MISSMATCHED REGION */
/*   CHECK FOR OVERFLOW */
    if (is2 + lg2 - 1 > *maxgel) {
	goto L700;
    }
    if (lg2 > 0) {
	sqcopy_(seqg + lg1, seqg2 + is2, &lg2, (ftnlen)1, (ftnlen)1);
    }
/*   CHECK FOR OVERFLOW */
    if (is1 + lc2 - 1 > *maxgel) {
	goto L700;
    }
    if (lc2 > 0) {
	sqcopy_(seqc + lc1, seqc2 + is1, &lc2, (ftnlen)1, (ftnlen)1);
    }
/*   POINT TO NEXT OUTPUT POSITIONS */
    is1 += lc2;
    is2 += lg2;
/*   GET NEXT MATCH */
    goto L50;
L500:

/*   FINISH RIGHT ENDS */
/*   ONLY COPY TO END OF GEL IN GEL AND TO THE SAME RELATIVE POSITION */
/*   IN THE CONTIG FOR DISPLAY PURPOSES AND FOR COUNTING MISMATCH */
/*   CURRENT ENDS AT LG1,LC1 */
/*   HOW FAR TO END OF GEL? */
/*   SET M */
    m = 0;
    l = *idg - lg1 + 1;
    if (l < 1) {
	goto L600;
    }
/*   CHECK FOR OVERFLOW */
    if (is2 + l - 1 > *maxgel) {
	goto L700;
    }
    sqcopy_(seqg + lg1, seqg2 + is2, &l, (ftnlen)1, (ftnlen)1);
/*   NEED TO COPY TO END OF GEL IN CONTIG FOR DISPLAY */
/*   POINT TO POSN IN CONTIG LEVEL WITH END OF GEL */
    m = lc1 + l - 1;
/*   IS THIS OVER END OF CONTIG? */
    if (m > *idc) {
	m = *idc;
    }
/*   NUMBER TO COPY */
    m = m - lc1 + 1;
/*   CHECK FOR OVERFLOW */
    if (is1 + m - 1 > *maxgel) {
	goto L700;
    }
    if (m > 0) {
	sqcopy_(seqc + lc1, seqc2 + is1, &m, (ftnlen)1, (ftnlen)1);
    }
L600:
/*   COUNT PADDING IN GEL */
    *itotpg = is2 + l - 1 - *idg;
/*   SET NEW LENGTHS FOR RETURN TO CALLING ROUTINE */
    *idout = is1 + m - 1;
    *idg = is2 + l - 1;
    *ifail = 0;
    return 0;
L700:
    erromf_("Matching region too long in lineup: alignment aborted", (ftnlen)
	    53);
    *ifail = 1;
    return 0;
} /* lineup_ */

/*     MERGE */

/*   ROUTINE SENT CONTIG WHOSE GELS MAY BE OUT OF ORDER */
/*   REORDERS GELS ON POSITION OF LEFT ENDS AND SETS LEFT */
/*   GEL NUMBER FOR THE REORDERED CONTIG */

/* Subroutine */ int merge_(relpg, lngthg, lnbr, rnbr, lincon, idbsiz)
integer *relpg, *lngthg, *lnbr, *rnbr, *lincon, *idbsiz;
{
    static integer m, n, i1, i2, nr;

/*   AUTHOR: RODGER STADEN */

/*   START AT LEFT END */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    n = lnbr[*lincon];
    goto L22;
L21:
/*   SET POINTER TO NEXT GEL TO RIGHT IN LIST */
    n = nr;
    if (i1 > 0) {
	n = i2;
    }
L22:
/*   SET POINTER TO NEXT GEL TO RIGHT */
    nr = rnbr[n];
    if (nr == 0) {
	goto L30;
    }
/*   HAVENT REACHED END YET */
    i1 = 0;
L23:
/*   ARE THESE 2 IN CORRECT ORDER IE N<=NR ? */
    if (relpg[n] <= relpg[nr]) {
	goto L21;
    }
/*   NOT IN ORDER SO CHAIN LEFT UNTIL CORRECTLY POSITIONED */
/*   THEN COME BACK TO THIS POINT AND CONTINUE */
/*   IF FIRST MOVE SAVE POSITION */
    if (i1 == 0) {
	i2 = n;
    }
    i1 = 1;
/*   EXCHANGE NEIGHBOURS */
    m = rnbr[nr];
    if (m != 0) {
	lnbr[m] = n;
    }
    m = lnbr[n];
    if (m != 0) {
	rnbr[m] = nr;
    }
    rnbr[n] = rnbr[nr];
    rnbr[nr] = n;
    lnbr[nr] = lnbr[n];
    lnbr[n] = nr;
/*   CHAIN BACK THRU LIST */
    n = lnbr[nr];
    if (n == 0) {
	goto L21;
    }
/*   END NOT REACHED */
    goto L23;
L30:
/*  ALL DONE POINTER AT RIGHT GEL */
    rnbr[*lincon] = n;
    return 0;
} /* merge_ */

/* Subroutine */ int removl_(matc, matg, matl, ip)
integer *matc, *matg, *matl, *ip;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, idelt, nmtch, k1, k2, k3, k4, k5, k6;
    extern /* Subroutine */ int bubbl3_();
    static integer ipp;

/*   AUTHOR: RODGER STADEN */

/*   SET POINTER TO FIRST MATCH */
    /* Parameter adjustments */
    --matl;
    --matg;
    --matc;

    /* Function Body */
    nmtch = 0;
L10:
/*   POINT TO NEXT MATCH */
    ++nmtch;
/*   SORT MATCHES ON LENGTH */
    ipp = *ip - nmtch + 1;
    bubbl3_(&matl[nmtch], &matg[nmtch], &matc[nmtch], &ipp);
/*   LOOK FOR END OF POSITIVES */
    i__1 = *ip;
    for (i__ = nmtch; i__ <= i__1; ++i__) {
	j = i__;
/* L20: */
	if (matl[i__] < 1) {
	    goto L30;
	}
    }
    ++j;
L30:
    *ip = j - 1;
/*   END OF POSITIVES AT IP */
    if (nmtch >= *ip) {
	return 0;
    }
    k1 = matc[nmtch];
    k2 = k1 + matl[nmtch] - 1;
    k3 = matg[nmtch];
    k4 = k3 + matl[nmtch] - 1;
/*   POINT TO FIRST MATCH TO TEST */
    k6 = nmtch + 1;
    i__1 = *ip;
    for (i__ = k6; i__ <= i__1; ++i__) {
/*   DO CONSENSUS FIRST */
/*   OVERLAP? */
	if (matc[i__] > k2) {
	    goto L100;
	}
	k5 = matc[i__] + matl[i__] - 1;
	if (k5 < k1) {
	    goto L100;
	}
/*   DOES OVERLAP */
/*   WHICH END */
	if (k5 <= k2) {
	    goto L80;
	}
/*   LENGTH TO REDUCE MATCH BY IS IDELT */
	idelt = k2 - matc[i__] + 1;
/*   NEW LENGTH */
	matl[i__] -= idelt;
/*  MOVE LEFT ENDS */
	matc[i__] += idelt;
	matg[i__] += idelt;
	goto L100;
L80:
/*   LENGTH */
	matl[i__] = k1 - matc[i__];
L100:
/*   NOW LOOK FOR OVERLAPS WITH GEL */
/*   OVERLAP? */
	if (matg[i__] > k4) {
	    goto L200;
	}
	k5 = matg[i__] + matl[i__] - 1;
	if (k5 < k3) {
	    goto L200;
	}
/*   DOES OVERLAP */
/*   WHICH END? */
	if (k5 <= k4) {
	    goto L180;
	}
/*   LENGTH TO REDUCE MATCH BY IS IDELT */
	idelt = k4 - matg[i__] + 1;
/*   NEW LENGTH */
	matl[i__] -= idelt;
/*   MOVE LEFT ENDS */
	matc[i__] += idelt;
	matg[i__] += idelt;
	goto L200;
L180:
/*   LENGTH */
	matl[i__] = k3 - matg[i__];
L200:
	;
    }
    goto L10;
} /* removl_ */

/* Subroutine */ int tpchek_(pc, pg, l, n)
integer *pc, *pg, *l, *n;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j1, k1;
    extern /* Subroutine */ int ml_();

/*     AUTHOR RODGER STADEN */
/*     IF OVERLAPPING BLOCKS ARE FOUND REMOVE THE SHORTER ONE */
/*     THEN REMOVE LARGE GAPS AT ENDS (THOSE AS LARGE AS THE END BLOCK) */
    /* Parameter adjustments */
    --l;
    --pg;
    --pc;

    /* Function Body */
    k1 = 2;
L1:
    i__1 = *n;
    for (i__ = k1; i__ <= i__1; ++i__) {
	j1 = i__;
	if (pc[i__] <= pc[i__ - 1]) {
	    goto L20;
	}
	if (pg[i__] <= pg[i__ - 1]) {
	    goto L20;
	}
/* L10: */
    }
/*     REMOVE LARGE GAPS FROM ENDS */
/*     THIS RULE OF THUMB COULD BE CHANGED TO USE A DIFFERENCE */
/*     BETWEEN THE NUMBERS OF MISMATCHING CHARACTERS */
    if (*n > 1) {
	k1 = pc[2] - pc[1] - l[1];
	j1 = pg[2] - pg[1] - l[1];
	if (max(k1,j1) > l[1]) {
	    ml_(&pc[1], &pg[1], &l[1], n, &c__1);
	    --(*n);
	}
	if (*n > 1) {
	    k1 = pc[*n] - pc[*n - 1] - l[*n - 1];
	    j1 = pg[*n] - pg[*n - 1] - l[*n - 1];
	    if (max(k1,j1) > l[*n]) {
		ml_(&pc[1], &pg[1], &l[1], n, n);
		--(*n);
	    }
	}
    }
    return 0;
L20:
    if (l[j1 - 1] > l[j1]) {
	ml_(&pc[1], &pg[1], &l[1], n, &j1);
    } else {
	i__1 = j1 - 1;
	ml_(&pc[1], &pg[1], &l[1], n, &i__1);
    }
/*  Until 25-11-90 next line was k1=j1 but this does not deal with all */
/*  cases: when a line is deleted we must compare it with the previous */
/*  one before dealing with the rest, because it could be left of that */
/*   one as well! */
/* Computing MAX */
    i__1 = 2, i__2 = j1 - 1;
    k1 = max(i__1,i__2);
    --(*n);
    goto L1;
} /* tpchek_ */

/* Subroutine */ int updcon_(relpg, lngthg, lnbr, rnbr, idbsiz, ngels, nconts,
	 seq, maxseq, idim1, cstart, cleno, lincon, nampro, seq2, idevr, 
	ifail, maxgel, idm, percd, mask, clist, seq_len, nampro_len, seq2_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *idbsiz, *ngels, *nconts;
char *seq;
integer *maxseq, *idim1, *cstart, *cleno, *lincon;
char *nampro, *seq2;
integer *idevr, *ifail, *maxgel, *idm;
real *percd;
integer *mask, *clist;
ftnlen seq_len;
ftnlen nampro_len;
ftnlen seq2_len;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer nbad, lreg, rreg, l, iladd[1], iradd[1], igelc;
    extern integer chnrp_();
    static integer itask, iwing, b1, s1;
    extern /* Subroutine */ int precn1_();
    static integer ld;
    extern /* Subroutine */ int makhca_(), addtit_(), erromf_();

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
    if (relpg[*lincon] - *cleno > 0) {

/* it is longer so we probably need to shift */

	if (s1 == 0) {

/* no readings start to the right of the new data */

	    if (*cstart + *cleno - 1 < *idim1) {

/* there are other contigs to the right */

/*           WRITE(*,*)'CSTART,CLENO',CSTART,CLENO */
		b1 = *cstart + *cleno;
/*            WRITE(*,*)'B1',B1 */
		i__1 = relpg[*lincon] - *cleno;
		makhca_(seq + 1, maxseq, &b1, &i__1, idim1, ifail, (ftnlen)1);
		if (*ifail != 0) {
		    erromf_("Error: consensus too long", (ftnlen)25);
		    return 0;
		}
	    } else {

/* there are no contigs to the right and no readings start to the right of */
/* the new one so nothing to shift */

	    }
	} else {

/* there are readings starting to the right of the new one */

/* shift from start of next reading to right */

	    l = *cstart + *cleno - 1;
/*           WRITE(*,*)'CSTART,CLENO,L',CSTART,CLENO,L */
	    ld = relpg[*lincon] - relpg[s1] + 1;
/*           WRITE(*,*)'LD',LD */
	    b1 = l - ld + 1;
/*            WRITE(*,*)'B1',B1 */
	    i__1 = relpg[*lincon] - *cleno;
	    makhca_(seq + 1, maxseq, &b1, &i__1, idim1, ifail, (ftnlen)1);
	    if (*ifail != 0) {
		erromf_("Error: consensus too long", (ftnlen)25);
		return 0;
	    }
	}
    }

/* now make new consensus (where do we put it,  do we need */
/* to give it a header, and what region do we make it for ? */
/* in the simplest case make it for relpg(ngels) to relpg(s1) -1 */
/* if s1=0 make it for relpg(ngels) to end of contig (relpg(lincon)) */
/* we give it a header if it is at the left end of the contig ie lnbr(ngels)=0 */

/* we always start at the left end of the new reading */

    lreg = relpg[*ngels];

/* we end at the next reading to the right or the end of the contig */

    if (s1 != 0) {
	rreg = relpg[s1] - 1;
    } else {
	rreg = relpg[*lincon];
    }

/* where do we put the new consensus ? */

    b1 = *cstart + relpg[*ngels] - 1;
/*      WRITE(*,*)'LREG,RREG',LREG,RREG */
/*            WRITE(*,*)'B1',B1 */

/* do we need to add a title */

    if (lnbr[*ngels] == 0) {
	b1 = *cstart - 20;
/*        WRITE(*,*)'ADD NEW TIT AT',B1 */
	addtit_(seq + b1, nampro, ngels, &b1, (ftnlen)1, nampro_len);
    }
    igelc = lnbr[*lincon];

/* note aconsn will chain along until it find the first useful reading!! */

/*      JOB = 2 */

/* set dummy values for precon (and iladd,iradd above) */

    nbad = 0;
    iwing = 0;

/* set task (normal consensus) */

    itask = 4;

/* add masking if required */

    if (*mask == 3) {
	itask += 32;
    }
    if (*mask == 4) {
	itask = 8;
    }
/*      WRITE(*,*)'BEFORE' */
/*      WRITE(*,*)(SEQ(JJJ),JJJ=1,IDIM1+RELPG(LINCON)-CLENO) */
/*         CALL FMTDB1(SEQ,IDIM1,1,IDIM1,60,6) */
/*      WRITE(*,*)'NOCONT,LREG,RREG,ITASK,B1' */
/*      WRITE(*,*)NOCONT,LREG,RREG,ITASK,B1 */
    precn1_(seq + 1, nampro, percd, idbsiz, lincon, &lreg, &rreg, &itask, 
	    idevr, &b1, maxgel, maxseq, &iwing, &nbad, iladd, iradd, ifail, (
	    ftnlen)1, nampro_len);
    if (*ifail != 0) {
	erromf_("Error calculating consensus", (ftnlen)27);
	return 0;
    }

/* before we leave we must make the overall consensus length correct */
/*  so add on the extra length (if any) which is the new length - old length */

/*      WRITE(*,*)'OLD IDIM1',IDIM1 */
    *idim1 = *idim1 + relpg[*lincon] - *cleno;
/*      WRITE(*,*)'after NEW IDIM1/2',IDIM1 */
/*      WRITE(*,*)(SEQ(JJJ),JJJ=1,IDIM1) */
/*         CALL FMTDB1(SEQ,IDIM1,1,IDIM1,60,6) */
} /* updcon_ */

/* Subroutine */ int makhca_(string, maxar, from, hsize, asize, ifail, 
	string_len)
char *string;
integer *maxar, *from, *hsize, *asize, *ifail;
ftnlen string_len;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;


/* make a hole of size hsize in character array size asize */

    /* Parameter adjustments */
    --string;

    /* Function Body */
    j = *asize + *hsize;
    if (j > *maxar) {
	*ifail = 1;
	return 0;
    }
    i__1 = *from;
    for (i__ = *asize; i__ >= i__1; --i__) {
	*(unsigned char *)&string[j] = *(unsigned char *)&string[i__];
	--j;
/* L10: */
    }
    *ifail = 0;
} /* makhca_ */

integer chnrp_(relpg, lngthg, lnbr, rnbr, idbsiz, lgel, ncont, lreg)
integer *relpg, *lngthg, *lnbr, *rnbr, *idbsiz, *lgel, *ncont, *lreg;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;


/* find first reading starting past lreg (0=none found) */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *lgel;
    ret_val = 0;
L10:
    if (i__ != 0) {
	if (relpg[i__] <= *lreg) {
	    i__ = rnbr[i__];
	    goto L10;
	}
	ret_val = i__;
	return ret_val;
    }
    return ret_val;
} /* chnrp_ */

integer chnrp1_(relpg, lngthg, lnbr, rnbr, idbsiz, lgel, lreg)
integer *relpg, *lngthg, *lnbr, *rnbr, *idbsiz, *lgel, *lreg;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;


/* find first reading with data covering or past lreg (0=none found) */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *lgel;
    ret_val = 0;
L10:
    if (i__ != 0) {
	if (relpg[i__] + (i__1 = lngthg[i__], abs(i__1)) - 1 < *lreg) {
	    i__ = rnbr[i__];
	    goto L10;
	}
	ret_val = i__;
	return ret_val;
    }
    return ret_val;
} /* chnrp1_ */

/* Subroutine */ int aerror_(list, name__, ierr, list_len, name_len)
char *list, *name__;
integer *ierr;
ftnlen list_len;
ftnlen name_len;
{
    /* Format strings */
    static char fmt_1000[] = "(a,i2)";
    static char fmt_1010[] = "(a,a,a)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(), s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    extern /* Subroutine */ int info_();
    static integer j, l;
    static char infod[60];
    extern /* Subroutine */ int erromf_();
    static char errmsg[333];
    extern /* Subroutine */ int tolist_();

    /* Fortran I/O blocks */
    static icilist io___223 = { 0, infod, 0, fmt_1000, 60, 1 };
    static icilist io___225 = { 0, errmsg, 0, fmt_1010, 333, 1 };



/* handle errors for assembly */

/* errors are: */
/* 0 file not found or file is of invalid format */
/* 1 read too short */
/* 2 failed to align and not entered */
/* 3 failed on entry */
/* 4 failed to align but entered */
/* 5 no match found during masked assembly */
    l = 1;
    i__1 = i_len(name__, name_len);
    for (j = 1; j <= i__1; ++j) {
	l = j;
	if (*(unsigned char *)&name__[j - 1] == ' ') {
	    goto L6;
	}
/* L5: */
    }
L6:
    s_wsfi(&io___223);
    do_fio(&c__1, name__, l);
    do_fio(&c__1, (char *)&(*ierr), (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___225);
    do_fio(&c__1, "Failed file ", (ftnlen)12);
    do_fio(&c__1, name__, l);
    do_fio(&c__1, "written to error file", (ftnlen)21);
    e_wsfi();
    erromf_(errmsg, (ftnlen)333);
    tolist_(list, infod, list_len, (ftnlen)60);
    info_(errmsg, (ftnlen)333);
} /* aerror_ */

/* Subroutine */ int sindb_(idevn, ngels, rnames, name__, job, rnames_len, 
	name_len)
integer *idevn, *ngels;
char *rnames, *name__;
integer *job;
ftnlen rnames_len;
ftnlen name_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int readn_();

    /* Parameter adjustments */
    rnames -= rnames_len;

    /* Function Body */
    if (*job == 1) {
	i__1 = *ngels;
	for (j = 1; j <= i__1; ++j) {
	    readn_(idevn, &j, rnames + j * rnames_len, rnames_len);
/*          WRITE(*,*)'INITIALISING RNAMES ',RNAMES(J) */
/* L10: */
	}
	return 0;
    } else if (*job == 2) {
	s_copy(rnames + *ngels * rnames_len, name__, rnames_len, name_len);
/*          WRITE(*,*)' ADDING TO RNAMES ',RNAMES(NGELS) */
    }
} /* sindb_ */

integer indb_(ngels, rnames, name__, rnames_len, name_len)
integer *ngels;
char *rnames, *name__;
ftnlen rnames_len;
ftnlen name_len;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer s_cmp();

    /* Local variables */
    static integer j;

    /* Parameter adjustments */
    rnames -= 40;

    /* Function Body */
    i__1 = *ngels;
    for (j = 1; j <= i__1; ++j) {
/*        WRITE(*,*)'CHECKING RNAMES ',NAME,' ',RNAMES(J) */
	if (s_cmp(name__, rnames + j * 40, name_len, (ftnlen)40) == 0) {
	    ret_val = j;
	    return ret_val;
	}
/* L10: */
    }
    ret_val = 0;
    return ret_val;
} /* indb_ */

/* Subroutine */ int slides_(seq1, idc, seq2, idim2, ms1, ms2, maxpg, maxpc, 
	minsli, matl, matc, matg, ip, seq1_len, seq2_len)
char *seq1;
integer *idc;
char *seq2;
integer *idim2, *ms1, *ms2, *maxpg, *maxpc, *minsli, *matl, *matc, *matg, *ip;
ftnlen seq1_len;
ftnlen seq2_len;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, n;
    extern /* Subroutine */ int savit_();
    static integer p1, p2;
    extern integer ctonum_();
    static integer ip1, p1s;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seq1;
    --seq2;
    --matg;
    --matc;
    --matl;

    /* Function Body */
    ip1 = *ip;
    *ip = 0;
/*   LEFT END S2 RELATIVE S1 - MAX PADS -2 READY FOR LOOP */
    p1s = *ms1 - *ms2 - *maxpc - 1;
/*   TRY NSLIDE START POSNS FOR SEQ2 */
/*      WRITE(*,*)'IDC,IDIM2',IDC,IDIM2 */
/*      WRITE(*,*)'SEQ1' */
/*      WRITE(*,*)(SEQ1(JJJ),JJJ=1,IDC) */
/*      WRITE(*,*)'SEQ2' */
/*      WRITE(*,*)(SEQ2(JJJ),JJJ=1,IDIM2) */
    i__1 = *maxpg + *maxpc + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       POINT TO SEQ1 START */
	++p1s;
/*       POINT TO CURRENT SEQ1 POSN */
	p1 = p1s;
	n = 0;
/*       COMPARE WHOLE LENGTH OF SEQ2 (IF P1 WITHIN RANGE) */
	i__2 = *idim2;
	for (j = 1; j <= i__2; ++j) {
	    p2 = j;
	    ++p1;
	    if (p1 < 1) {
		goto L50;
	    }
/*         OFF RIGHT END? IF SO MAY HAVE BEEN A MATCH */
	    if (p1 > *idc) {
		goto L40;
	    }
	    if (ctonum_(seq1 + p1, (ftnlen)1) == ctonum_(seq2 + p2, (ftnlen)1)
		    ) {
		goto L45;
	    }
L40:
	    if (n >= *minsli) {
		savit_(&n, &p1, &p2, ip, &matl[1], &matc[1], &matg[1], &ip1);
	    }
	    n = 0;
	    goto L50;
L45:
	    ++n;
L50:
	    ;
	}
/*       GOOD SCORE AT END? NEED TO INCREMENT POINTERS FOR SAVIT */
	++p1;
	++p2;
	if (n >= *minsli) {
	    savit_(&n, &p1, &p2, ip, &matl[1], &matc[1], &matg[1], &ip1);
	}
/* L100: */
    }
} /* slides_ */

/* Subroutine */ int padcon_(relpg, lngthg, lnbr, rnbr, ngels, nconts, gel, 
	lincon, posn, nc, idbsiz, idevr, maxgel, gel_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts;
char *gel;
integer *lincon, *posn, *nc, *idbsiz, *idevr, *maxgel;
ftnlen gel_len;
{
    /* Initialized data */

    static char pad[1+1] = "*";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_sign();

    /* Local variables */
    static integer i__, k, x;
    extern /* Subroutine */ int readw_();
    static integer llino;
    extern /* Subroutine */ int insbas_(), writec_(), writeg_(), shiftt_();

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */
/*   NOW FIND FIRST CHAR THAT OVERLAPS REGION */
    llino = lnbr[(0 + (0 + (*lincon << 2))) / 4];
L30:
    x = relpg[llino] + (i__1 = lngthg[llino], abs(i__1)) - 1;
    if (x >= *posn) {
	goto L40;
    }
/*   NOT IN REGION */
    llino = rnbr[llino];
    goto L30;
L40:
/*   NOW GET THIS GEL FROM DISK */
    readw_(idevr, &llino, gel + 1, maxgel, (ftnlen)1);
/*   CALC POSN IN THIS GEL TO EDIT */
    x = *posn - relpg[llino] + 1;
    k = x;
    lngthg[llino] += i_sign(nc, &lngthg[llino]);
    if ((i__1 = lngthg[llino], abs(i__1)) > *maxgel) {
	lngthg[llino] = i_sign(maxgel, &lngthg[llino]);
    }
/*   INSBAS TAKES CARE OF OVERFLOW OF READ IN ITS OWN WAY */
    i__1 = *nc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	insbas_(idevr, &llino, &k, pad, (ftnlen)1);
/* L61: */
    }
/*   WRITE NEW LINE */
    writeg_(idevr, &llino, &relpg[llino], &lngthg[llino], &lnbr[llino], &rnbr[
	    llino]);
L65:
/*   NOW GET NEXT GEL */
    llino = rnbr[llino];
/*   LAST GEL? */
    if (llino == 0) {
	goto L70;
    }
/*   DOES IT HAVE DATA IN REGION? */
/*   IE DO RELPG  AND RELPG+LNGTHG-1 LIE EITHER SIDE OF POSN? */
    if (relpg[llino] > *posn) {
	goto L70;
    }
    x = relpg[llino] + (i__1 = lngthg[llino], abs(i__1)) - 1;
    if (x < *posn) {
	goto L65;
    }
/*  WITHIN */
    goto L40;
L70:
/*   INSERTS FINISHED SO NEED TO INCREMENT ALL THOSE GELS TO RIGHT */
    llino = lnbr[*lincon];
L75:
    if (relpg[llino] > *posn) {
	goto L80;
    }
L76:
    llino = rnbr[llino];
    if (llino == 0) {
	goto L90;
    }
    goto L75;
L80:
    relpg[llino] += *nc;
/*   WRITE NEW LINE */
    writeg_(idevr, &llino, &relpg[llino], &lngthg[llino], &lnbr[llino], &rnbr[
	    llino]);
    goto L76;
L90:
/*   NEED TO INCREMENT CONTIG LINE */
    relpg[*lincon] += *nc;
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
/*   Now move tags along on the consensus */
    i__1 = *idbsiz - *lincon;
    shiftt_(idevr, &i__1, posn, nc);
} /* padcon_ */

/* Subroutine */ int upchek_(pc, pg, l, n)
integer *pc, *pg, *l, *n;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j1, k1, dc, dg;
    extern /* Subroutine */ int ml_();

/*     AUTHOR RODGER STADEN */

/* only allow gaps that are shorter than the next block of identity */

    /* Parameter adjustments */
    --l;
    --pg;
    --pc;

    /* Function Body */
    k1 = 2;
L1:
    i__1 = *n;
    for (i__ = k1; i__ <= i__1; ++i__) {
	j1 = i__;
	dc = pc[i__] - pc[i__ - 1] - l[i__ - 1];
	dg = pg[i__] - pg[i__ - 1] - l[i__ - 1];
	if ((i__2 = dc - dg, abs(i__2)) >= l[i__]) {
	    goto L20;
	}
/* L10: */
    }
    return 0;
L20:
/*      WRITE(*,*)'REMOVING!!' */
    ml_(&pc[1], &pg[1], &l[1], n, &j1);
/*      IF(L(J1-1).GT.L(J1)) THEN */
/*        CALL ML(PC,PG,L,N,J1) */
/*      ELSE */
/*        CALL ML(PC,PG,L,N,J1-1) */
/*      END IF */
/* Computing MAX */
    i__1 = 2, i__2 = j1 - 1;
    k1 = max(i__1,i__2);
    --(*n);
    goto L1;
} /* upchek_ */

/* Subroutine */ int padcop_(seqg, seqg2, lg1, mg, l5, is2, lg2, maxgel, 
	ifail, seqc, idc, ic1, seqg_len, seqg2_len, seqc_len)
char *seqg, *seqg2;
integer *lg1, *mg, *l5, *is2, *lg2, *maxgel, *ifail;
char *seqc;
integer *idc, *ic1;
ftnlen seqg_len;
ftnlen seqg2_len;
ftnlen seqc_len;
{
    /* Initialized data */

    static char dubbl[1*4+1] = "DBVH";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int info_();
    static integer j, m, idone, jc1, maxreq, mgm1;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seqg2;
    --seqg;
    --seqc;

    /* Function Body */
    jc1 = *ic1;
/* Make seqg2 from seqg placing L5 padding chars before position MG */
/* which is the start of the next block of identity. Try to put the */
/* padding either in line with consensus pads, or next to double */
/* codes. The positions in seqg are LG1 to MG-1. seqg2 needs to be long */
/* enough to be extended from IS2 to IS2 + L5 -1 + MGM1-LG1 +1 */
/* ie we add L5 pads, plus the chars between and including  LG1 and MGM1 */
    idone = 0;
/*   POINT TO END OF MISMATCH */
    mgm1 = *mg - 1;
/*   MAY BE NO CHARS TO COPY */
    if (mgm1 < *lg1) {
	goto L111;
    }
/*  Next check added 26-2-91 */
    maxreq = *is2 + *l5 - 1 + mgm1 - *lg1 + 1;
    if (mgm1 > *maxgel || maxreq > *maxgel) {
	info_("Matching region too large in padcop: alignment aborted", (
		ftnlen)54);
	*ifail = 1;
	return 0;
    }
    i__1 = mgm1;
    for (j = *lg1; j <= i__1; ++j) {
	if (idone < *l5) {
	    if (jc1 > 0 && jc1 < *idc) {
		if (*(unsigned char *)&seqc[jc1] == '*') {
		    ++(*is2);
		    ++jc1;
		    ++idone;
		    goto L109;
		}
	    }
	    for (m = 1; m <= 4; ++m) {
		if (*(unsigned char *)&seqg[j] == *(unsigned char *)&dubbl[m 
			- 1]) {
		    ++(*is2);
		    ++jc1;
		    ++idone;
		    goto L109;
		}
/* L108: */
	    }
L109:
	    ;
	}
	*(unsigned char *)&seqg2[*is2] = *(unsigned char *)&seqg[j];
	++(*is2);
	++jc1;
/* L110: */
    }
L111:
/*   ALL CHARS COPIED. ENOUGH PADDING? */
    if (idone < *l5) {
	*is2 = *is2 + *l5 - idone;
    }
/*   IS2 SHOULD NOW BE POINTING AT NEXT CHAR */
/*   ZERO LG2 TO SHOW CALLING ROUTINE COPYING DONE */
    *lg2 = 0;
    *ifail = 0;
} /* padcop_ */

/* Subroutine */ int ml_(pc, pg, l, n, j)
integer *pc, *pg, *l, *n, *j;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --l;
    --pg;
    --pc;

    /* Function Body */
    i__1 = *n - 1;
    for (i__ = *j; i__ <= i__1; ++i__) {
	pc[i__] = pc[i__ + 1];
	pg[i__] = pg[i__ + 1];
	l[i__] = l[i__ + 1];
/* L10: */
    }
} /* ml_ */

/* Subroutine */ int mstlkl_(seq, idim, seq_len)
char *seq;
integer *idim;
ftnlen seq_len;
{
    /* System generated locals */
    integer i__1;
    char ch__1[1];

    /* Local variables */
    static integer i__, j, k;
    extern /* Character */ VOID charsu_();
    extern integer indexs_();

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seq;

    /* Function Body */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = indexs_(seq + i__, &k, (ftnlen)1);
	charsu_(ch__1, (ftnlen)1, &j);
	*(unsigned char *)&seq[i__] = *(unsigned char *)&ch__1[0];
/* L100: */
    }
} /* mstlkl_ */

integer indexs_(c__, s, c_len)
char *c__;
integer *s;
ftnlen c_len;
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

/*      DATA DUP/'CTAG1234DVBHKLMNRY5678ctag*,-'/ */
/*  changed 28-7-91 to give 10 to old zeroes and 100 to lowercase */
    i__ = *(unsigned char *)c__;
    i__ = shotc_1.points[i__];
    *s = scores[i__ - 1];
    ret_val = ind[i__ - 1];
    return ret_val;
} /* indexs_ */

/* Character */ VOID charsu_(ret_val, ret_val_len, i__)
char *ret_val;
ftnlen ret_val_len;
integer *i__;
{
    /* Initialized data */

    static char c__[6+1] = "CTAG*-";

    /* Builtin functions */
    /* Subroutine */ int s_copy();

    s_copy(ret_val, c__ + (0 + (0 + (*i__ - 1))), (ftnlen)1, *i__ - (*i__ - 1)
	    );
} /* charsu_ */

/*     SAVIT */

/* Subroutine */ int savit_(n, j, k, ip, s1, s2, s3, ip1)
integer *n, *j, *k, *ip, *s1, *s2, *s3, *ip1;
{
/*   AUTHOR: RODGER STADEN */

    /* Parameter adjustments */
    --s3;
    --s2;
    --s1;

    /* Function Body */
    ++(*ip);
/*   TEST FOR OVERFLOW */
    if (*ip > *ip1) {
	return 0;
    }
    s1[*ip] = *n;
    s2[*ip] = *j - *n;
    s3[*ip] = *k - *n;

    return 0;
} /* savit_ */

integer clen_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz, iin)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *iin;
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, len;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LEFT GEL NUMBER OR ZERO FOR ERROR */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *iin;
    ret_val = 0;
    len = 0;
L10:
    if (i__ != 0) {
/* Computing MAX */
	i__2 = len, i__3 = relpg[i__] + (i__1 = lngthg[i__], abs(i__1)) - 1;
	len = max(i__2,i__3);
	i__ = rnbr[i__];
	if (i__ == *iin) {
	    return ret_val;
	}
	goto L10;
    }
    ret_val = len;
    return ret_val;
} /* clen_ */

/* Subroutine */ int gllino_(relpg, lngthg, lnbr, rnbr, idbsiz, nconts, llino,
	 lincon)
integer *relpg, *lngthg, *lnbr, *rnbr, *idbsiz, *nconts, *llino, *lincon;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n, mxt;


/* routine to get the left gel number and contig line of the longest contig */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    *llino = 0;
    *lincon = 0;
    n = *idbsiz - *nconts;
    mxt = 0;
    i__1 = *idbsiz - 1;
    for (i__ = n; i__ <= i__1; ++i__) {
	if (relpg[i__] > mxt) {
	    mxt = relpg[i__];
	    *llino = lnbr[i__];
	    *lincon = i__;
	}
/* L4: */
    }
} /* gllino_ */

integer clinno_(lnbr, idbsiz, nconts, iin)
integer *lnbr, *idbsiz, *nconts, *iin;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer j, n;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LINE NUMBER OR ZERO FOR ERROR */
    /* Parameter adjustments */
    --lnbr;

    /* Function Body */
    ret_val = 0;
    n = *idbsiz - *nconts;
    i__1 = *idbsiz - 1;
    for (j = n; j <= i__1; ++j) {
	if (lnbr[j] == *iin) {
	    ret_val = j;
	    return ret_val;
	}
/* L10: */
    }
    return ret_val;
} /* clinno_ */

/* Subroutine */ int shftla_(string, maxar, froms, to, frome, string_len)
char *string;
integer *maxar, *froms, *to, *frome;
ftnlen string_len;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;


/* shift an array left from froms to to */

    /* Parameter adjustments */
    --string;

    /* Function Body */
    j = *to;
    i__1 = *frome;
    for (i__ = *froms; i__ <= i__1; ++i__) {
	*(unsigned char *)&string[j] = *(unsigned char *)&string[i__];
	++j;
/* L10: */
    }
} /* shftla_ */

integer randc_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz, iin, lincon, 
	llino)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *iin, *
	lincon, *llino;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;
    extern integer gclin_(), chainl_();


/* return reading and contig number for contig containing IIN */
/* -1 = ERROR, 0 = OK */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = chainl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, iin);
    if (i__ == 0) {
	ret_val = -1;
	return ret_val;
    }
    *llino = i__;
    i__ = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, llino);
    if (i__ == 0) {
	ret_val = -2;
	return ret_val;
    }
    *lincon = i__;
    ret_val = 0;
    return ret_val;
} /* randc_ */

/*      CMPLMT */

/*   SUBROUTINE TO REVERSE AND COMPLEMENT GELS AND DATA BASE */
/*   THE POSITIONS OF THE RIGHT ENDS OF GELS ARE FIRST STORED */
/*   IN RELPG THEN WE DO A BUBBLE SORT ON THESE POSITIONS */
/*   UPDATING RELATIONSHIPS AS WE GO */
/*   ALSO SEQUENCES ARE COMPLEMENTED, SIGNS OF LENGTH ARE */
/*   MULTIPLIED BY -1 AND THE CONTIG LINE IS ALTERED */
/* Subroutine */ int cmplmt_(relpg, lngthg, lnbr, rnbr, ngels, nconts, lincon,
	 llino, gel, idbsiz, idevr, maxgel, gel_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *lincon, *llino;
char *gel;
integer *idbsiz, *idevr, *maxgel;
ftnlen gel_len;
{
    /* Format strings */
    static char fmt_1000[] = "(\002Complementing contig\002,i8)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    extern /* Subroutine */ int info_();
    static integer m, n, x;
    static char infod[30];
    static integer i1, i2, nl;
    extern /* Subroutine */ int comtag_(), cplseq_(), writec_(), writeg_();

    /* Fortran I/O blocks */
    static icilist io___271 = { 0, infod, 0, fmt_1000, 30, 1 };


/*   AUTHOR: RODGER STADEN */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */
    s_wsfi(&io___271);
    do_fio(&c__1, (char *)&(*llino), (ftnlen)sizeof(integer));
    e_wsfi();
    info_(infod, (ftnlen)30);
/*   CHAIN THRU AND PUT RIGHT ENDS IN RELPG */
    n = *llino;
L10:
    relpg[n] = relpg[n] + (i__1 = lngthg[n], abs(i__1)) - 1;
    if (rnbr[n] == 0) {
	goto L20;
    }
    n = rnbr[n];
    goto L10;
L20:

/*   NOW EFFECTIVELY BUBBLE SORT ON RELPG */
    n = rnbr[*lincon];
    goto L22;
L21:
    n = nl;
    if (i1 > 0) {
	n = i2;
    }
L22:
    nl = lnbr[n];
    if (nl == 0) {
	goto L30;
    }
    i1 = 0;
L23:
    if (relpg[n] >= relpg[nl]) {
	goto L21;
    }
/*   NOT IN CORRECT ORDER SO CHAIN ALONG UNTIL CORRECT,THEN COME */
/*   BACK TO THIS POINT AND CONTINUE */
/*   IF FIRST MOVE THIS LINE SET POINTER TO CURRENT POSITION */
    if (i1 == 0) {
	i2 = n;
    }
    i1 = 1;

/*   EXCHANGE NEIGHBOURS. CURRENTLY LOOKING AT N AND ITS LEFT */
/*   NBR, AND THE LEFT NBR IS FURTHER RIGHT THAN N */
/*   FIX UP POINTERS TO LEFT AND RIGHT OF THESE TWO */
    m = lnbr[nl];
    if (m != 0) {
	rnbr[m] = n;
    }
    m = rnbr[n];
    if (m != 0) {
	lnbr[m] = nl;
    }
    lnbr[n] = lnbr[nl];
    lnbr[nl] = n;
    rnbr[nl] = rnbr[n];
    rnbr[n] = nl;
/*   CHAIN BACK THRU LIST WITH THIS LINE */
    n = rnbr[nl];
    if (n == 0) {
	goto L21;
    }
/*   IE END MET */
    goto L23;
L30:
/*   FINISH WITH LEFT END IN N */
L40:
/*   NOW REVERSE NBRS SO CHAIN BACK RIGHT */
    nl = rnbr[n];
    if (nl == 0) {
	goto L50;
    }
    rnbr[n] = lnbr[n];
    lnbr[n] = nl;
    n = nl;
    goto L40;
L50:
/*   NEED TO FIX UP NEW LEFT END */
    rnbr[n] = lnbr[n];
    lnbr[n] = 0;
/*   ALL POINTERS FIXED NOW DO RELATIVE POSITION */
/*   FINISH WITH LEFT END IN N */
/*   SO CHAIN BACK RIGHT */
/*   SAVE RIGHT LINE NUMBER */
    nl = n;
    x = relpg[n];
L60:
    relpg[n] = -(relpg[n] - x) + 1;
    if (rnbr[n] == 0) {
	goto L70;
    }
    n = rnbr[n];
    goto L60;
L70:
/*   NOW FIX CONTIG LINE */
    lnbr[*lincon] = nl;
    rnbr[*lincon] = n;
/*   WRITE NEW CONTIG LINE */
    i__1 = *idbsiz - *lincon;
    writec_(idevr, &i__1, &relpg[*lincon], &lnbr[*lincon], &rnbr[*lincon]);
/*   NOW REVERSE AND COMPLEMENT GELS */
    n = nl;
L80:
/* Added by Simon 17-March-1993 */
    cplseq_(idevr, &n, maxgel);
/*      READ(IDEVW,REC=N)GEL */
/*      CALL READW(IDEVW,N,GEL,MAXGEL) */
/*      M=ABS(LNGTHG(N)) */
/*      CALL SQREV(GEL,M) */
/*      CALL SQCOM(GEL,M) */
/*      CALL WRITEW(IDEVW,N,GEL,MAXGEL) */
/*   CHANGE SIGNS */
    lngthg[n] = -lngthg[n];
/*   WRITE NEW GEL LINE */
    writeg_(idevr, &n, &relpg[n], &lngthg[n], &lnbr[n], &rnbr[n]);
/*   ANY MORE? */
    n = rnbr[n];
    if (n != 0) {
	goto L80;
    }
/*   NO MORE */
/*   Now update consensus tag list */
    i__1 = *idbsiz - *lincon;
    comtag_(idevr, &i__1, &relpg[*lincon]);
    return 0;
} /* cmplmt_ */

/* Subroutine */ int inits_()
{
    /* Initialized data */

    static char dup[29+1] = "CTAG1234DVBHKLMNRY5678ctag*,-";

    static integer i__, j;

/*  AUTHOR RODGER STADEN */
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

    for (i__ = 0; i__ <= 255; ++i__) {
	shotc_1.points[i__] = 29;
/* L30: */
    }
    for (i__ = 1; i__ <= 29; ++i__) {
	j = *(unsigned char *)&dup[i__ - 1];
	shotc_1.points[j] = i__;
/* L35: */
    }
} /* inits_ */

integer chainl_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz, iin)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *iin;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__, j;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LEFT GEL NUMBER OR ZERO FOR ERROR */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *iin;
    j = i__;
    ret_val = 0;
L10:
    if (i__ != 0) {
	j = i__;
	i__ = lnbr[i__];
	if (i__ == *iin) {
	    return ret_val;
	}
	goto L10;
    }
    ret_val = j;
    return ret_val;
} /* chainl_ */

/* Subroutine */ int shiftc_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idevr, 
	idbsiz, ign, ncont, dist)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idevr, *idbsiz, *ign,
	 *ncont, *dist;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern integer clen_();
    static integer i__, l;
    extern /* Subroutine */ int writec_(), writeg_();

/*  AUTHOR: RODGER STADEN */
/*  SHIFTS PART OF A CONTIG FORM GEL IGN TO RIGHT END */
/*  CONTIG LINE NUMBER IF NCONT */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    i__ = *ign;
L10:
    if (i__ != 0) {
	relpg[i__] += *dist;
	writeg_(idevr, &i__, &relpg[i__], &lngthg[i__], &lnbr[i__], &rnbr[i__]
		);
	i__ = rnbr[i__];
	goto L10;
    }
/*  UPDATE CONTIG LENGTH */
    l = clen_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, ign);
    relpg[*ncont] = l;
    i__1 = *idbsiz - *ncont;
    writec_(idevr, &i__1, &relpg[*ncont], &lnbr[*ncont], &rnbr[*ncont]);
} /* shiftc_ */

/* Subroutine */ int fndcon_(seq, idim, cends, nends, idcend, maxcon, seq_len)
char *seq;
integer *idim, *cends, *nends, *idcend, *maxcon;
ftnlen seq_len;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static char dc[1*8];
    extern integer indexa_(), jfromc_();
    extern /* Subroutine */ int erromf_();

/*   AUTHOR: RODGER STADEN */
/*   STORES THEIR POSITIONS IN CENDS AND THEIR LEFT LINE NUMBERS IN NENDS */
    /* Parameter adjustments */
    --seq;
    --nends;
    --cends;

    /* Function Body */
    *idcend = 0;
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&seq[i__] != '<') {
	    goto L10;
	}
	++(*idcend);
/*       PUT POSITION OF LEFT END OF CONTIG IN CENDS */
	cends[*idcend] = i__;
	k = indexa_(seq + i__, &c__20, ".", (ftnlen)1, (ftnlen)1);
	if (k == 0) {
	    goto L20;
	}
	k += i__;
/*        IF (.NOT.((SEQ(K+MAXDG).EQ.'-').OR.(SEQ(K+MAXDG).EQ.'>'))) */
/*     +     GOTO20 */
	for (j = 1; j <= 8; ++j) {
	    if (*(unsigned char *)&seq[k] == '-' || *(unsigned char *)&seq[k] 
		    == '>') {
		goto L6;
	    }
	    *(unsigned char *)&dc[j - 1] = *(unsigned char *)&seq[k];
	    ++k;
/* L5: */
	}
L6:
	i__2 = j - 1;
	nends[*idcend] = jfromc_(dc, &i__2, (ftnlen)1);
L10:
	;
    }
/*     STORE POSITION OF LAST CHAR +1 TO SIMPLIFY DISPLAY ROUTINES */
    cends[*idcend + 1] = *idim + 1;
    return 0;
L20:
    erromf_("Error in FNDCON: illegal consensus header", (ftnlen)41);
    *idcend = 0;
} /* fndcon_ */

integer gclin_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz, iin)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *iin;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer j, n;

/*  AUTHOR: RODGER STADEN */
/*  RETURNS CONTIG LINE NUMBER OR ZERO FOR ERROR */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    ret_val = 0;
    n = *idbsiz - *nconts;
    i__1 = *idbsiz - 1;
    for (j = n; j <= i__1; ++j) {
	if (lnbr[j] == *iin) {
	    ret_val = j;
	    return ret_val;
	}
/* L10: */
    }
    return ret_val;
} /* gclin_ */

/*     SQCOM */
/* Subroutine */ int sqcomm_(seq, idim, seq_len)
char *seq;
integer *idim;
ftnlen seq_len;
{
    /* Initialized data */

    static char list1[1*12+1] = "CTAGctagedfi";
    static char list2[1*12+1] = "GATCgatcifde";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static char temp[1];
    static integer i__, j;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seq;

    /* Function Body */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)temp = *(unsigned char *)&seq[i__];
	for (j = 1; j <= 12; ++j) {
	    if (*(unsigned char *)temp == *(unsigned char *)&list1[j - 1]) {
		*(unsigned char *)&seq[i__] = *(unsigned char *)&list2[j - 1];
		goto L99;
	    }
/* L50: */
	}
L99:
/* L100: */
	;
    }
} /* sqcomm_ */

/* Subroutine */ int arrfim_(idev, seqnce, j, seqnce_len)
integer *idev;
char *seqnce;
integer *j;
ftnlen seqnce_len;
{
    /* Initialized data */

    static char endchr[1+1] = "@";
    static char space[1+1] = " ";
    static char titchr[1+1] = ";";

    /* Format strings */
    static char fmt_1001[] = "(80a1)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsfe(), do_fio(), e_rsfe();

    /* Local variables */
    static integer idmx;
    static char temp[1*80];
    static integer i__;
    extern /* Subroutine */ int erromf_();

    /* Fortran I/O blocks */
    static cilist io___300 = { 1, 0, 1, fmt_1001, 0 };


    /* Parameter adjustments */
    --seqnce;

    /* Function Body */
    idmx = *j;
    *j = 0;
L1:
    io___300.ciunit = *idev;
    i__1 = s_rsfe(&io___300);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_fio(&c__80, temp, (ftnlen)1);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = e_rsfe();
L100001:
    if (i__1 < 0) {
	goto L30;
    }
    if (i__1 > 0) {
	goto L40;
    }
    if (*(unsigned char *)&temp[0] == *(unsigned char *)&titchr[0]) {
	goto L1;
    }
/* L10: */
    for (i__ = 1; i__ <= 80; ++i__) {
	if (*(unsigned char *)&temp[i__ - 1] != *(unsigned char *)&space[0]) {
	    if (*(unsigned char *)&temp[i__ - 1] == *(unsigned char *)&endchr[
		    0]) {
		return 0;
	    }
	    if (*j == idmx) {
		erromf_("Too much data, so input truncated", (ftnlen)33);
		return 0;
	    }
	    ++(*j);
	    *(unsigned char *)&seqnce[*j] = *(unsigned char *)&temp[i__ - 1];
	}
/* L20: */
    }
    goto L1;
L30:
    return 0;
L40:
    erromf_("Error reading file", (ftnlen)18);
    *j = 0;
} /* arrfim_ */

integer jfromc_(chars, length, chars_len)
char *chars;
integer *length;
ftnlen chars_len;
{
    /* Format strings */
    static char fmt_1002[] = "(i10)";

    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();
    integer s_rsfi(), do_fio(), e_rsfi();

    /* Local variables */
    static integer lens, list;
    static char number[10];
    extern /* Subroutine */ int erromf_(), rjstfy_();

    /* Fortran I/O blocks */
    static icilist io___305 = { 1, number, 0, fmt_1002, 10, 1 };


/*   AUTHOR: RODGER STADEN */
/*   INTEGER FUNCTION TO CONVERT CHARACTER ARRAYS OF */
/*   NUMERALS TO BINARY FORM */
/*   LENGTH OF STRING NUMBER */
    /* Parameter adjustments */
    --chars;

    /* Function Body */
    lens = 10;
    s_copy(number, " ", (ftnlen)10, (ftnlen)1);
    rjstfy_(chars + 1, number, &lens, length, (ftnlen)1, (ftnlen)10);
    i__1 = s_rsfi(&io___305);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_fio(&c__1, (char *)&list, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = e_rsfi();
    if (i__1 != 0) {
	goto L100;
    }
    ret_val = list;
    return ret_val;
L100:
    ret_val = 0;
    erromf_("Error in internal read, value set to zero", (ftnlen)41);
    return ret_val;
} /* jfromc_ */

/* Subroutine */ int erromf_(string, string_len)
char *string;
ftnlen string_len;
{
    extern /* Subroutine */ int fverr_();

    fverr_(&c__0, "Error", string, (ftnlen)5, string_len);
} /* erromf_ */

/* Subroutine */ int busyf_()
{
/*      WRITE(*,*)'BUSY' */
} /* busyf_ */

/* Subroutine */ int remgbc_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 gel, maxgel, idevr, iok, array, iall, iopt, gel_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz;
char *gel;
integer *maxgel, *idevr, *iok, *array, *iall, *iopt;
ftnlen gel_len;
{
    /* Format strings */
    static char fmt_1000[] = "(\002Reading name not found: \002,a)";
    static char fmt_1001[] = "(\002Renumbering reading\002,i8,\002 to\002,i8)"
	    ;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    extern /* Subroutine */ int info_();
    static integer idum1, idum2, i__, iredc;
    extern integer gclin_();
    static char infod[100];
    static integer remme, icont, ifrom;
    extern /* Subroutine */ int dbchek_(), delgel_();
    extern integer chainl_(), gnread_();
    static char namarc[40];
    extern logical crucal_();
    extern integer nameno_();
    extern /* Subroutine */ int remcnl_(), rmgtag_(), movenc_(), movgel_(), 
	    erromf_(), unlnkr_(), writrn_();

    /* Fortran I/O blocks */
    static icilist io___310 = { 0, infod, 0, fmt_1000, 100, 1 };
    static icilist io___315 = { 0, infod, 0, fmt_1000, 100, 1 };
    static icilist io___318 = { 0, infod, 0, fmt_1001, 100, 1 };


/*      SUBROUTINE REMGBC(RELPG,LNGTHG,LNBR,RNBR,NGELS,NCONTS,IDBSIZ, */
/*     +KBIN,KBOUT,GEL,MAXGEL,IDEVR,IDEV2,FILNAM, */
/*     +IHELPS,IHELPE,HELPF,IDEVH,IOK,ARRAY) */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --array;
    --gel;

    /* Function Body */
    dbchek_(idevr, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &c__5, idbsiz, 
	    ngels, nconts, iok);
    if (*iok > 1) {
	return 0;
    }

/* IALL = 1 all reads, = 2 non-crucial only */
/* IOPT = 1 remove, = 2 move */

/* assumes db is logically consistent */

/* here we start a new contig with the selected readings */

    if (*iopt == 2) {
L40:

/* use list of names */

	*iok = gnread_(namarc, (ftnlen)40);
	if (*iok == 1) {
	    goto L200;
	}
	if (*iok != 0) {
	    goto L40;
	}
	remme = nameno_(namarc, ngels, idevr, (ftnlen)40);
	if (remme == 0) {
	    s_wsfi(&io___310);
	    do_fio(&c__1, namarc, (ftnlen)40);
	    e_wsfi();
	    erromf_(infod, (ftnlen)100);
	    goto L40;
	}
	i__ = chainl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, &remme);
	icont = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, &i__);
	if (icont == 0) {
	    erromf_("No contig line for this reading", (ftnlen)31);
	    *iok = 1;
	    goto L200;
	}

/* special case: read is already a single contig (do nothing) */

	if (lnbr[remme] == 0 && rnbr[remme] == 0) {
	    goto L40;
	}
	if (*iall == 2) {
	    if (crucal_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, &
		    remme, &lnbr[icont], &array[1], maxgel, &c__3)) {
		goto L40;
	    }
	}
	unlnkr_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		idbsiz, &remme, &icont, idevr, &idum1, &idum2, &c__1, iok);
	if (*iok == 0) {
	    movenc_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    idbsiz, idevr, &remme, iok);
	    goto L40;
	}
	erromf_("Escape from REMGBC", (ftnlen)18);
L200:
	dbchek_(idevr, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &c__5, 
		idbsiz, ngels, nconts, iok);
	if (*iok > 1) {
	    return 0;
	}
	return 0;

/* here we remove reads from the database */

    } else if (*iopt == 1) {
L30:

/* use list of names */

	*iok = gnread_(namarc, (ftnlen)40);
	if (*iok == 1) {
	    goto L100;
	}
	if (*iok != 0) {
	    goto L30;
	}
	remme = nameno_(namarc, ngels, idevr, (ftnlen)40);
	if (remme == 0) {
	    s_wsfi(&io___315);
	    do_fio(&c__1, namarc, (ftnlen)40);
	    e_wsfi();
	    erromf_(infod, (ftnlen)100);
	    goto L30;
	}
	i__ = chainl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, &remme);
	icont = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, &i__);
	if (icont == 0) {
	    erromf_("No contig line for this reading", (ftnlen)31);
	    *iok = 1;
	    goto L100;
	}

/* allow test for being crucial */

	if (*iall == 2) {
	    if (crucal_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, &
		    remme, &lnbr[icont], &array[1], maxgel, &c__3)) {
		goto L30;
	    }
	}

/* special case: read is already a single contig so dont unlink */
/* but must reduce number of contigs which is normally handled */
/* by remcnl in unlnkr */

	iredc = 0;
	if (lnbr[remme] != 0 || rnbr[remme] != 0) {
	    iredc = 1;
	    unlnkr_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    idbsiz, &remme, &icont, idevr, &idum1, &idum2, &c__1, iok)
		    ;
	    if (*iok != 0) {
		erromf_("Escape from REMGBC", (ftnlen)18);
		dbchek_(idevr, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &
			c__5, idbsiz, ngels, nconts, iok);
		if (*iok > 1) {
		    return 0;
		}
		return 0;
	    }
	}
	ifrom = *ngels;
	rmgtag_(idevr, &remme, &c__0, &c__0);
	delgel_(idevr, &remme);
	if (remme != ifrom) {
	    s_wsfi(&io___318);
	    do_fio(&c__1, (char *)&ifrom, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&remme, (ftnlen)sizeof(integer));
	    e_wsfi();
	    info_(infod, (ftnlen)100);
	    movgel_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    idbsiz, gel + 1, &ifrom, &remme, idevr, maxgel, (ftnlen)1)
		    ;
	}
	if (iredc == 0) {

/* remove contig line icont */

/*          WRITE(*,*)'REMOVING CONTIG LINE',ICONT */
	    remcnl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		    idbsiz, &icont, idevr);
	}
	--(*ngels);
	writrn_(idevr, ngels, nconts);
	goto L30;
L100:
	dbchek_(idevr, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &c__5, 
		idbsiz, ngels, nconts, iok);
	if (*iok > 1) {
	    return 0;
	}
    }
} /* remgbc_ */

/* Subroutine */ int breakc_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 idevr, ir, iok)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *idevr, *ir, 
	*iok;
{
    extern integer gclin_();
    extern /* Subroutine */ int dbchek_();
    static integer il;
    extern /* Subroutine */ int cbreak_();
    extern integer chainl_();
    extern /* Subroutine */ int erromf_();
    static integer lconto, nconto, ncontr, idm, ilo;

/*   AUTHOR: RODGER STADEN */

/* note IR is the number of the read that will become a left end */

/* ROUTINE TO BREAK A CONTIG INTO 2 */
/* LEFT GEL OF NEW RIGHT CONTIG IS IR */
/* RIGHT GEL OF NEW LEFT CONTIG IS IL */
/* LEFT GEL OF OLD LEFT CONTIG IS ILO */
/* CONTIG LINE OF OLD CONTIG IS NCONTO */
/* CONTIG LINE OF NEW RIGHT CONTIG IS NCONTR */
/* CONTIG LINE OF NEW LEFT CONTIG IS NCONTO */
/* LENGTH OF OLD CONTIG IS LCONTO */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    il = lnbr[*ir];
    if (il == 0) {
	erromf_("Reading is already a left end", (ftnlen)29);
	*iok = 1;
	return 0;
    }
    idm = 5;
    dbchek_(idevr, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &idm, idbsiz, 
	    ngels, nconts, iok);
    if (*iok > 1) {
	return 0;
    }
    *iok = 1;
    ncontr = *idbsiz - *nconts - 1;
    if (ncontr <= *ngels) {
	erromf_("Insufficient space for new contig line", (ftnlen)38);
	return 0;
    }
    ilo = chainl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, ir);
    nconto = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, &ilo);
    lconto = relpg[nconto];
    cbreak_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    idevr, ir, &il, &ilo, &nconto, &ncontr, iok);
    dbchek_(idevr, &relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &idm, idbsiz, 
	    ngels, nconts, iok);
} /* breakc_ */

/* Subroutine */ int cbreak_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 idevr, ir, il, ilo, nconto, ncontr, iok)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *idevr, *ir, 
	*il, *ilo, *nconto, *ncontr, *iok;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern integer clen_();
    extern /* Subroutine */ int info_();
    static integer i__, l;
    extern /* Subroutine */ int movec_(), busyf_();
    static integer iremme;
    extern /* Subroutine */ int shiftc_(), erromf_(), spltag_(), writec_(), 
	    writeg_(), dupnot_(), unlnkr_(), writrn_();
    static integer irr;

/*   AUTHOR: RODGER STADEN */
/* ROUTINE TO BREAK A CONTIG INTO 2 */
/* LEFT GEL OF NEW RIGHT CONTIG IS IR */
/* RIGHT GEL OF NEW LEFT CONTIG IS IL */
/* LEFT GEL OF OLD LEFT CONTIG IS ILO */
/* CONTIG LINE OF OLD CONTIG IS NCONTO */
/* CONTIG LINE OF NEW RIGHT CONTIG IS NCONTR */
/* CONTIG LINE OF NEW LEFT CONTIG IS NCONTO */
/* LENGTH OF OLD CONTIG IS LCONTO */

/* new strategy: do the right and left contigs then use unlnkr to handle */
/* the consequences for the right contig - it may make a single */
/* new contig or it might have to be broken into several depending */
/* on whether the right end reads of the left contig held together */
/* some of the left end reads of the righthand contig. We send unlnkr */
/* a job=2 for this case (job=1 is for unlinking and dealing with the */
/* consequences) */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    *iok = 1;
    busyf_();
    ++(*nconts);
/*  MAKE NEW CONTIG A COPY OF OLD */
    relpg[*ncontr] = relpg[*nconto];
    lngthg[*ncontr] = lngthg[*nconto];
    lnbr[*ncontr] = *ir;
    rnbr[*ncontr] = rnbr[*nconto];
/*      CALL INFO('Writing new right contig line') */
    i__1 = *idbsiz - *ncontr;
    writec_(idevr, &i__1, &relpg[*ncontr], &lnbr[*ncontr], &rnbr[*ncontr]);
/*  WRITE LAST LINE OF DB */
/*      CALL INFO('Increasing number of contigs by 1') */
    writrn_(idevr, ngels, nconts);
/* Move new contig next to old contig */
    i__1 = *idbsiz - *ncontr;
    i__2 = *idbsiz - *nconto;
    movec_(idevr, &i__1, &i__2);

/*  NEED LENGTH FOR OLD LEFT CONTIG */
    rnbr[*il] = 0;
/*      L = CLEN(RELPG,LNGTHG,LNBR,RNBR,NGELS,NCONTS, */
/*     +IDBSIZ,IL) */
/*  Change 24/6/93 jkb */
    l = clen_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
	    idbsiz, ilo);
    if (l < 1) {
	erromf_("New left contig has zero length. Break not made", (ftnlen)47)
		;
	return 0;
    }
/*  Split tags */
    i__1 = *idbsiz - *nconto;
    i__2 = *idbsiz - *ncontr;
    spltag_(idevr, &i__1, &i__2, &relpg[*ir], &l);
    i__1 = *idbsiz - *nconto;
    i__2 = *idbsiz - *ncontr;
    dupnot_(idevr, &i__1, &i__2);
    relpg[*nconto] = l;
    rnbr[*nconto] = *il;
/*  DO CONTIG LINE FOR NEW LEFT CONTIG */
    i__1 = *idbsiz - *nconto;
    writec_(idevr, &i__1, &relpg[*nconto], &lnbr[*nconto], &rnbr[*nconto]);
/*  DO GEL LINE FOR RIGHT GEL OF NEW LEFT CONTIG */
    writeg_(idevr, il, &relpg[*il], &lngthg[*il], &lnbr[*il], &rnbr[*il]);
/*  DO GEL LINE FOR NEW RIGHT CONTIG */
    lnbr[*ir] = 0;
    writeg_(idevr, ir, &relpg[*ir], &lngthg[*ir], &lnbr[*ir], &rnbr[*ir]);
/*  NOW SHIFT */
    i__ = 1 - relpg[*ir];
    info_("Shifting readings in right contig", (ftnlen)33);
    shiftc_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idevr, 
	    idbsiz, ir, ncontr, &i__);
    irr = rnbr[*ncontr];
    unlnkr_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    &iremme, ncontr, idevr, ir, &irr, &c__2, iok);
    info_("Break completed", (ftnlen)15);
    *iok = 0;
} /* cbreak_ */

/* Subroutine */ int unlnkr_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 remme, icont, idevr, leftr, rightr, job, iok)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *remme, *
	icont, *idevr, *leftr, *rightr, *job, *iok;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern integer clen_();
    static integer l, iread, redge;
    extern /* Subroutine */ int movec_();
    static integer jcont, ocont;
    extern /* Subroutine */ int unlnk_(), busyf_();
    static integer nextr, lo;
    extern /* Subroutine */ int rmctag_(), shiftc_(), remcnl_();
    static integer ishift;
    extern /* Subroutine */ int erromf_(), spltag_();
    extern integer nextcl_();
    extern /* Subroutine */ int writec_(), writeg_();
    static integer ispltl;
    extern /* Subroutine */ int dupnot_();
    static integer ispltr;
    extern /* Subroutine */ int writrn_();


/* This routine has 2 functions: originally it was part of removing a */
/* read from a contig, now it also tidies up after breaking a contig. */

/* Original function: job = 1 */
/* we want to unlink a read from a contig so we do that then we have */
/* to tidy up the consequences which can be complicated. So always start */
/* at the left end of the contig and check contiguity for each read, starting */
/* new contigs if required. When we leave here the read is unlinked and */
/* unless an error occurs everything is tidy. We always start a new contig */
/* line for the contig we begin with. Before we leave we delete the original */
/* contig line. The number of reads is unchanged but we expect a different */
/* number of contigs if the unlinked read was crucial. */
/* unlink a read REMME from a contig */
/* unlink REMME making links to its nbrs (if they exist) */
/* ICONT contig line for remme */

/* New function: job = 2 */
/* We have broken a contig and hence created two contigs. Weve done all */
/* the work except the right hand contig may not be contiguous because */
/* reads at the right end of the left contig may have provided the contiguity. */
/* So we come here to sort out the right hand contig. */
/* ICONT is contig line for this right hand contig */

/* starting from the left end of the contig chain thru checking contiguity */
/* if contiguity ok continue chaining */
/* if contiguity fails, end the current contig and start a new one */
/* use shiftc to write out the relationships (gel lines and contig line) as */
/* well as shift if necessary. */

/* LEFTR left read current cont */
/* RIGHTR right read current contig */
/* IREAD is current read */
/* NEXTR is next read */

/* unlink remme */

/*      WRITE(*,*)'UNLINK',REMME,ICONT */
/*      WRITE(*,*)'LEFT,RIGHT',LNBR(REMME),RNBR(REMME) */

/* if we come in to remove a reading then */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    busyf_();
    if (*job == 1) {
	unlnk_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, remme, 
		icont, leftr, rightr, iok);
	if (*iok == 0) {
	    return 0;
	}
    }

/* the rest of the routine takes a contig and turns it into */
/* contiguous contigs */

    ocont = *icont;
    ispltr = relpg[*leftr];
    ispltl = ispltr;
L10:


/* create a new contig if theres space */

    jcont = nextcl_(idbsiz, ngels, nconts);
/*      WRITE(*,*)'NEW CONTIG LINE NUMBER',JCONT */
    if (jcont == 0) {
	erromf_("Not enough space to create new contig", (ftnlen)37);
	erromf_("We are in trouble now - revert to copy!", (ftnlen)39);
	*iok = 1;
	return 0;
    } else {
	++(*nconts);
	lnbr[jcont] = *leftr;
	rnbr[jcont] = *rightr;
	lngthg[jcont] = 0;

/* write contig record to ensure it's allocation */

	i__1 = *idbsiz - jcont;
	writec_(idevr, &i__1, &relpg[jcont], &lnbr[jcont], &rnbr[jcont]);

/* reorder contig order */

/*        WRITE(*,*)'current reading ',REMME */
	i__1 = *idbsiz - *icont;
	movec_(idevr, nconts, &i__1);

/* reset the left neighbour for the left read */

	lnbr[*leftr] = 0;

/* Duplicate the tags using spltag_() */

/*        WRITE(*,*)'Splitting at ',ISPLTL,ISPLTR */
	i__1 = *idbsiz - *icont;
	i__2 = *idbsiz - jcont;
	spltag_(idevr, &i__1, &i__2, &ispltr, &ispltr);
	i__1 = *idbsiz - *icont;
	i__2 = *idbsiz - jcont;
	dupnot_(idevr, &i__1, &i__2);
	if (ispltl != ispltr) {
	    i__1 = *idbsiz - *icont;
	    i__2 = ispltl + 1;
	    i__3 = ispltr + 1;
	    rmctag_(idevr, &i__1, &i__2, &i__3);
	}

/* always use shiftc to sort out relpgs and contig length and line and */
/* write it all to disk */

	ishift = 1 - relpg[*leftr];
/* Remember old length */
	lo = relpg[*icont];
	shiftc_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		idevr, idbsiz, leftr, &jcont, &ishift);
	l = clen_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, 
		idbsiz, leftr);
/*  If it's shrunk then remove consensus tags */
	if (l < lo) {
	    i__1 = *idbsiz - jcont;
	    i__2 = l + 1;
	    i__3 = lo + 1;
	    rmctag_(idevr, &i__1, &i__2, &i__3);
	}
    }

    iread = *leftr;

/* set right edge for contiguity check (its as far right as weve reached) */

    redge = (i__1 = lngthg[iread], abs(i__1));
L20:
/*      WRITE(*,*)'IREAD',IREAD */

    if (rnbr[iread] != 0) {
	nextr = rnbr[iread];
/*        WRITE(*,*)'NEXTR',NEXTR */

/* contiguous ? */

	if (relpg[nextr] > redge) {

/* not contiguous so start new contig (first correct last one! having */
/* saved the right read number) and also set rnbr(iread)=0 */

/*          WRITE(*,*)'NOT CONTIGUOUS' */
	    *rightr = rnbr[jcont];
	    rnbr[jcont] = iread;
	    rnbr[iread] = 0;
	    relpg[jcont] = clen_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], 
		    ngels, nconts, idbsiz, &lnbr[jcont]);
/*          WRITE(*,*)'NEW CONTIG LENGTH',RELPG(JCONT) */
	    writeg_(idevr, &iread, &relpg[iread], &lngthg[iread], &lnbr[iread]
		    , &rnbr[iread]);
	    i__1 = *idbsiz - jcont;
	    writec_(idevr, &i__1, &relpg[jcont], &lnbr[jcont], &rnbr[jcont]);
	    ispltr = relpg[nextr];
	    ispltl = relpg[iread] + (i__1 = lngthg[iread], abs(i__1)) - 1;
	    *leftr = nextr;
	    *icont = jcont;
	    goto L10;
	} else {
/*          WRITE(*,*)'CONTIGUOUS' */

/* reset redge */

/* Computing MAX */
	    i__2 = redge, i__3 = relpg[nextr] + (i__1 = lngthg[nextr], abs(
		    i__1)) - 1;
	    redge = max(i__2,i__3);
	    iread = nextr;
	    goto L20;
	}
    }

/* reached end of linked list so sort out the number of contigs and gels */
/* counts on disk */
/*       WRITE(*,*)'WRITING NGELS,NCONTS',NGELS,NCONTS */
    writrn_(idevr, ngels, nconts);

/* remove contig line ocont */

    remcnl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    &ocont, idevr);
    *iok = 0;
} /* unlnkr_ */

/* Subroutine */ int unlnk_(relpg, lngthg, lnbr, rnbr, idbsiz, remme, icont, 
	leftr, rightr, iok)
integer *relpg, *lngthg, *lnbr, *rnbr, *idbsiz, *remme, *icont, *leftr, *
	rightr, *iok;
{
    static integer j, k;

/*      WRITE(*,*)'UNLINK',REMME,ICONT */
/*      WRITE(*,*)'LEFT,RIGHT',LNBR(REMME),RNBR(REMME) */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    *iok = 1;
    j = lnbr[*remme];
    k = rnbr[*remme];
    if (j != 0) {
	rnbr[j] = k;
    }
    if (k != 0) {
	lnbr[k] = j;
    }

/* initialise leftr and rightr (making sure they are not remme) */

    if (lnbr[*icont] != *remme) {
	*leftr = lnbr[*icont];
    } else {
	if (rnbr[*remme] != 0) {
	    *leftr = rnbr[*remme];
	} else {

/* single read contig surely theres nothing to do? CHECK THIS with an example */

	    *iok = 0;
	    return 0;
	}
    }
    if (rnbr[*icont] != *remme) {
	*rightr = rnbr[*icont];
    } else {
	if (lnbr[*remme] != 0) {
	    *rightr = lnbr[*remme];
	} else {

/* single read contig weve done above */

	    *iok = 0;
	    return 0;
	}
    }
} /* unlnk_ */

/* Subroutine */ int movgel_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 gel, from, to, idevr, maxgel, gel_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz;
char *gel;
integer *from, *to, *idevr, *maxgel;
ftnlen gel_len;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern integer gclin_();
    static logical lefte;
    extern /* Subroutine */ int movnm_();
    extern integer chainl_();
    static logical righte;
    extern /* Subroutine */ int erromf_(), writec_();
    static integer nconto;
    extern /* Subroutine */ int writeg_();

/*   Subroutine to move a gel from line from to line to */
    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */
    lefte = FALSE_;
    righte = FALSE_;

/* left end ? */

    if (lnbr[*from] == 0) {
	lefte = TRUE_;
    }

/* right end ? */

    if (rnbr[*from] == 0) {
	righte = TRUE_;
    }

/* if both true remove the contig line, then overwrite the gel */

    if (lefte && righte) {
	nconto = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, from);
	if (nconto == 0) {
	    erromf_("No contig line for this reading   ", (ftnlen)34);
	} else {
	    lnbr[nconto] = *to;
	    rnbr[nconto] = *to;
	    i__1 = *idbsiz - nconto;
	    writec_(idevr, &i__1, &relpg[nconto], &lnbr[nconto], &rnbr[nconto]
		    );
	}
    } else if (lefte) {
	nconto = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, from);
	if (nconto == 0) {
	    erromf_("No contig line for this reading   ", (ftnlen)34);
	} else {
	    lnbr[nconto] = *to;
	    i__1 = *idbsiz - nconto;
	    writec_(idevr, &i__1, &relpg[nconto], &lnbr[nconto], &rnbr[nconto]
		    );
	}
    } else if (righte) {
	i__ = chainl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, from);
	nconto = gclin_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, 
		nconts, idbsiz, &i__);
	if (nconto == 0) {
	    erromf_("No contig line for this reading   ", (ftnlen)34);
	} else {
	    if (rnbr[nconto] != *from) {
		erromf_("No contig line for this reading   ", (ftnlen)34);
	    } else {
		rnbr[nconto] = *to;
		i__1 = *idbsiz - nconto;
		writec_(idevr, &i__1, &relpg[nconto], &lnbr[nconto], &rnbr[
			nconto]);
	    }
	}
    }
    relpg[*to] = relpg[*from];
    lngthg[*to] = lngthg[*from];
    lnbr[*to] = lnbr[*from];
    rnbr[*to] = rnbr[*from];
/* Added by Simon 17-March-1993 */
    movnm_(idevr, from, to);
/*      CALL READW(IDEVW,FROM,GEL,MAXGEL) */
/*      CALL WRITEW(IDEVW,TO,GEL,MAXGEL) */
/*      CALL READN(IDEVN,FROM,NAMGEL) */
/*      CALL WRITEN(IDEVN,TO,NAMGEL) */
    writeg_(idevr, to, &relpg[*to], &lngthg[*to], &lnbr[*to], &rnbr[*to]);
/*   Do neighbours */
    if (lnbr[*from] != 0) {
	i__ = lnbr[*from];
	rnbr[i__] = *to;
	writeg_(idevr, &i__, &relpg[i__], &lngthg[i__], &lnbr[i__], &rnbr[i__]
		);
    }
    if (rnbr[*from] != 0) {
	i__ = rnbr[*from];
	lnbr[i__] = *to;
	writeg_(idevr, &i__, &relpg[i__], &lngthg[i__], &lnbr[i__], &rnbr[i__]
		);
    }
/*      CALL MOVTAG(FROM,TO) */
} /* movgel_ */

logical crucal_(relpg, lngthg, lnbr, rnbr, idbsiz, remme, llino, array, 
	maxgel, job)
integer *relpg, *lngthg, *lnbr, *rnbr, *idbsiz, *remme, *llino, *array, *
	maxgel, *job;
{
    /* System generated locals */
    integer i__1, i__2;
    logical ret_val;

    /* Local variables */
    static integer rreg, last, irno, i__, j;
    extern logical fandl_();
    extern /* Subroutine */ int filli_();
    static integer first;
    extern /* Subroutine */ int multi_();
    extern integer chnrp1_();
    static integer lo;
    extern /* Subroutine */ int erromf_();


/* Routine to see if a reading is crucial. It is crucial if */
/* a. it is not completely covered by data on both strands job=3 */
/* b. it is not completely covered by data on the plus strand job=2 */
/* c. it is not completely covered by data on the minus strand job=3 */

/* we test by filling an array of length=readlength with 1's then */
/* multiplying each element covered by data on the plus strand by 2 */
/* and each lement on the minus strand by 3. Then if any element is */
/* not zero given the test mod(array(i),6) we know the read is crucial. */
/* note only for job=3 is contiguity assured (others could be touching) */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --array;

    /* Function Body */
    ret_val = TRUE_;
    i__2 = (i__1 = lngthg[*remme], abs(i__1));
    filli_(&array[1], &i__2, &c__1);
    irno = chnrp1_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], idbsiz, llino, &
	    relpg[*remme]);
    if (irno == 0) {
	erromf_("Scream 1: error in crucal", (ftnlen)25);
	return ret_val;
    }
    rreg = relpg[*remme] + (i__1 = lngthg[*remme], abs(i__1)) - 1;
L10:

/* get the first and last useful data positions for this reading */

/*      WRITE(*,*)'IRNO',IRNO */
    i__2 = relpg[*remme] + (i__1 = lngthg[*remme], abs(i__1)) - 1;
    if (fandl_(&relpg[irno], &lngthg[irno], &relpg[*remme], &i__2, &first, &
	    last)) {
/*        WRITE(*,*)'FIRST,LAST',FIRST,LAST */
	if (irno != *remme) {
	    lo = last - first + 1;
	    if (lngthg[irno] > 0) {
		multi_(&array[relpg[irno] + first - relpg[*remme]], &lo, &
			c__2);
	    } else {
		multi_(&array[relpg[irno] + first - relpg[*remme]], &lo, &
			c__3);
	    }
	}
    }
    irno = rnbr[irno];
    if (irno > 0) {
	if (relpg[irno] <= rreg) {
	    goto L10;
	}
    }
/*      WRITE(*,*)'END OF FIRST LOOP',IRNO */

/* all relevant readings processed. Is the reading covered */

    if (*job == 1) {

/* want plus strand covered */

	j = 2;
	i__2 = (i__1 = lngthg[*remme], abs(i__1));
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (array[i__] % j != 0) {
		return ret_val;
	    }
/* L20: */
	}
	ret_val = FALSE_;
    } else if (*job == 2) {

/* want minus strand covered */

	j = 3;
	i__1 = (i__2 = lngthg[*remme], abs(i__2));
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (array[i__] % j != 0) {
		return ret_val;
	    }
/* L30: */
	}
	ret_val = FALSE_;
    } else if (*job == 3) {

/* want plus and minus strand covered */

	j = 6;
/*        WRITE(*,*)(ARRAY(K),K=1,ABS(LNGTHG(REMME))) */
	i__2 = (i__1 = lngthg[*remme], abs(i__1));
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (array[i__] % j != 0) {
		return ret_val;
	    }
/* L40: */
	}
	ret_val = FALSE_;
/*        WRITE(*,*)CRUCAL */
    } else {
	erromf_("Scream 2: in crucal", (ftnlen)19);
    }
    return ret_val;
} /* crucal_ */

logical fandl_(relpg, lngthg, lreg, rreg, first, last)
integer *relpg, *lngthg, *lreg, *rreg, *first, *last;
{
    /* System generated locals */
    logical ret_val;


/* Function to decide which parts of a reading are needed */
/* for processing. The read starts at relpg and has length */
/* lngthg. The active region of the contig is lreg to rreg. */
/* This cannot be used as a check if there is more relevent data */
/* because it will return false if a read does not reach lreg and */
/* the next read might reach. */
/* First check if theres any useful data */

    if (*relpg > *rreg || *relpg + abs(*lngthg) <= *lreg) {
	ret_val = FALSE_;
	return ret_val;
    }
    if (*relpg >= *lreg) {
	*first = 1;
    } else {
	*first = *lreg - *relpg + 1;
    }
    if (*relpg + abs(*lngthg) - 1 <= *rreg) {
	*last = abs(*lngthg);
    } else {
	*last = *rreg - *relpg + 1;
    }
    ret_val = TRUE_;
    return ret_val;
} /* fandl_ */

/* Subroutine */ int multi_(array, size, factor)
integer *array, *size, *factor;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* multiply each element of ARRAY by value FACTOR */

    /* Parameter adjustments */
    --array;

    /* Function Body */
    i__1 = *size;
    for (i__ = 1; i__ <= i__1; ++i__) {
	array[i__] *= *factor;
/* L10: */
    }
} /* multi_ */

/* Subroutine */ int movenc_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 idevr, moveme, iok)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *idevr, *
	moveme, *iok;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer jcont;
    extern /* Subroutine */ int erromf_();
    extern integer nextcl_();
    extern /* Subroutine */ int writec_(), writeg_(), writrn_();


/* routine to move a read MOVEME to start a new contig */
/* we can fail if there is not enough space */

    /* Parameter adjustments */
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;

    /* Function Body */
    lnbr[*moveme] = 0;
    rnbr[*moveme] = 0;
    relpg[*moveme] = 1;

/* leave orientation the same */

    writeg_(idevr, moveme, &relpg[*moveme], &lngthg[*moveme], &lnbr[*moveme], 
	    &rnbr[*moveme]);

/* start a new contig */

    jcont = nextcl_(idbsiz, ngels, nconts);
    if (jcont == 0) {
	erromf_("Not enough space to create new contig", (ftnlen)37);
	erromf_("We are in trouble now - revert to copy!", (ftnlen)39);
	*iok = 1;
	return 0;
    }
    lnbr[jcont] = *moveme;
    rnbr[jcont] = *moveme;
    lngthg[jcont] = 0;
    relpg[jcont] = (i__1 = lngthg[*moveme], abs(i__1));
    i__1 = *idbsiz - jcont;
    writec_(idevr, &i__1, &relpg[jcont], &lnbr[jcont], &rnbr[jcont]);
    ++(*nconts);
    writrn_(idevr, ngels, nconts);
    *iok = 0;
} /* movenc_ */

/*      SUBROUTINE REMCNL(RELPG,LNGTHG,LNBR,RNBR,NGELS,NCONTS,IDBSIZ, */
/*     +REMME,IDEVR) */
/*      INTEGER RELPG(IDBSIZ),LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ) */
/*      INTEGER REMME */
/* Routine to remove a contig line from a db */
/* Loop deals with case of remove top contig */
/* We also remove all annotations associated with this contig. */
/*      CALL RMCTAG(IDEVR, IDBSIZ-REMME, 0, 0) */
/* Move down all lines from above */
/*      DO 10 I = REMME,IDBSIZ-NCONTS+1,-1 */
/*        RELPG(I) = RELPG(I-1) */
/*        LNGTHG(I) = LNGTHG(I-1) */
/*        LNBR(I) = LNBR(I-1) */
/*        RNBR(I) = RNBR(I-1) */
/*        CALL GETCTG(IDEVR, IDBSIZ-(I-1), IANNO) */
/*        CALL PUTCTG(IDEVR, IDBSIZ-I, IANNO) */
/*        CALL WRITEC(IDEVR,IDBSIZ-I,RELPG(I),LNBR(I),RNBR(I)) */
/* 10    CONTINUE */
/*      NCONTS = NCONTS - 1 */
/*      CALL WRITRN(IDEVR,NGELS,NCONTS) */
/*      END */
integer nextcl_(idbsiz, ngels, nconts)
integer *idbsiz, *ngels, *nconts;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;


/* return next free contig line or zero if full */

    ret_val = 0;
    i__ = *idbsiz - *nconts - 1;
    if (i__ > *ngels) {
	ret_val = i__;
    }
    return ret_val;
} /* nextcl_ */

/* Subroutine */ int remcon_(relpg, lngthg, lnbr, rnbr, ngels, nconts, idbsiz,
	 icont, gel, llino, idevr, maxgel, temp, gel_len)
integer *relpg, *lngthg, *lnbr, *rnbr, *ngels, *nconts, *idbsiz, *icont;
char *gel;
integer *llino, *idevr, *maxgel, *temp;
ftnlen gel_len;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ndel, i__, j;
    extern /* Subroutine */ int delgel_(), remcnl_(), rmgtag_(), movgel_();

/*   AUTHOR: RODGER STADEN */

/* problem is to remove a contig. Strategy is to find all its read numbers */
/* and store them in temp. Then go thru and replace them by reads from the */
/* end of the list of reads. An earlier version tried to chain thru but */
/* ran into complications about replacing reads by reads they pointed to! */

    /* Parameter adjustments */
    --temp;
    --rnbr;
    --lnbr;
    --lngthg;
    --relpg;
    --gel;

    /* Function Body */
    i__1 = *idbsiz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp[i__] = 0;
/* L10: */
    }

    ndel = 0;
    i__ = *llino;
L20:
    temp[i__] = 1;
    ++ndel;
    if (rnbr[i__] != 0) {
	i__ = rnbr[i__];
	goto L20;
    }

/* now in temp all the reads in the contig are 1, the rest of temp is zero */
/* let i be the read to move and j the last read left in the list of reads */
/* so if temp(i) is 1 and temp(j) is zero move read j to i and set j = j - 1 */
/* then deal with the next i. */
/* if temp(i) is 1 and temp(j) is also 1 simply set j = j - 1 and try to move */
/* that one */
/* we stop when weve gone so far along the list that the next read to */
/* delete is equal to the number of reads left in the list */

/* what are the difficult cases? */
/* 1. only one contig */
/* 2. all reads to right of llino are are near the end of the list of reads */


/* Firstly we remove all tags and sequence (etc) */

    i__1 = *ngels;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (temp[i__] == 1) {
	    rmgtag_(idevr, &i__, &c__0, &c__0);
	    delgel_(idevr, &i__);
	}
/* L40: */
    }

/* And secondly shift the readings down. */

    i__ = 1;
    j = *ngels;
L30:
    if (temp[i__] == 1) {
	if (temp[j] == 0) {
/*            WRITE(KBOUT,*)'MOVE ',J,' TO ',I */
	    movgel_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], &j, nconts, 
		    idbsiz, gel + 1, &j, &i__, idevr, maxgel, (ftnlen)1);
	    --j;
	} else {
	    --j;
	    if (i__ < j) {
		goto L30;
	    }
	}
    }
    if (i__ < j) {
	++i__;
	goto L30;
    }

/* fix up number of reads and the contig record */

    *ngels -= ndel;
    remcnl_(&relpg[1], &lngthg[1], &lnbr[1], &rnbr[1], ngels, nconts, idbsiz, 
	    icont, idevr);
} /* remcon_ */

/* Subroutine */ int fmtdb_(seq1, idim, isw, ise, linlen, idev, seq1_len)
char *seq1;
integer *idim, *isw, *ise, *linlen, *idev;
ftnlen seq1_len;
{
    /* Format strings */
    static char fmt_1003[] = "()";
    static char fmt_1001[] = "(\002 \002,12(5x,i6))";
    static char fmt_1002[] = "(\002  \002,12(10a1,1x))";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfi(), e_wsfi(), do_fio();

    /* Local variables */
    extern /* Subroutine */ int info_();
    static integer isww, j, k;
    static char infod[100];
    static integer ie, kl[12], is, kkk;

    /* Fortran I/O blocks */
    static icilist io___361 = { 0, infod, 0, fmt_1003, 100, 1 };
    static icilist io___365 = { 0, infod, 0, fmt_1001, 100, 1 };
    static icilist io___367 = { 0, infod, 0, fmt_1002, 100, 1 };


/*   NOTE SAME AS FMTSEP! */
/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seq1;

    /* Function Body */
    isww = *isw - 1;
    ie = *isw - 1;
L1:
    s_wsfi(&io___361);
    e_wsfi();
    info_(infod, (ftnlen)100);
/*   SET UP DECIMAL COUNTERS */
    i__1 = *linlen / 10;
    for (j = 1; j <= i__1; ++j) {
	isww += 10;
	kl[j - 1] = isww;
/* L50: */
    }
    is = ie + 1;
    ie += *linlen;
    if (ie > *ise) {
	ie = *ise;
    }
    s_wsfi(&io___365);
/* Computing MIN */
    i__2 = ie - is + 1;
    i__1 = min(i__2,*linlen) / 10;
    for (kkk = 1; kkk <= i__1; ++kkk) {
	do_fio(&c__1, (char *)&kl[kkk - 1], (ftnlen)sizeof(integer));
    }
    e_wsfi();
    info_(infod, (ftnlen)100);
    s_wsfi(&io___367);
    i__2 = ie;
    for (k = is; k <= i__2; ++k) {
	do_fio(&c__1, seq1 + k, (ftnlen)1);
    }
    e_wsfi();
    info_(infod, (ftnlen)100);
    if (ie == *ise) {
	return 0;
    }
    goto L1;
} /* fmtdb_ */

/*     SQCOM */
/* Subroutine */ int sqcom_(seq, idim, seq_len)
char *seq;
integer *idim;
ftnlen seq_len;
{
    /* Initialized data */

    static char list1[1*31+1] = "ACGTN-acgtnBDHKMRSVWYbdhkmrsvwy";
    static char list2[1*31+1] = "TGCAN-tgcanVHDMKYSBWRvhdmkysbwr";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static char temp[1];
    static integer i__, j;

/*   AUTHOR: RODGER STADEN */
/*   CHANGE TO IUB CODES 21-1-98 */
    /* Parameter adjustments */
    --seq;

    /* Function Body */
/*      DATA LIST1/ */
/*     +'C','T','A','G', */
/*     +'c','t','a','g', */
/*     +'D','V','B','H', */
/*     +'d','v','b','h', */
/*     +'K','L','M','N', */
/*     +'k','l','m','n', */
/*     +'R','Y','U', */
/*     +'r','y','u', */
/*     +'1','2','3','4', */
/*     +'5','6','7','8'/ */
/*      DATA LIST2/ */
/*     +'G','A','T','C', */
/*     +'g','a','t','c', */
/*     +'H','B','V','D', */
/*     +'h','b','v','d', */
/*     +'N','M','L','K', */
/*     +'n','m','l','k', */
/*     +'Y','R','A', */
/*     +'y','r','a', */
/*     +'4','3','2','1', */
/*     +'6','5','7','8'/ */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)temp = *(unsigned char *)&seq[i__];
	for (j = 1; j <= 31; ++j) {
	    if (*(unsigned char *)temp == *(unsigned char *)&list1[j - 1]) {
		*(unsigned char *)&seq[i__] = *(unsigned char *)&list2[j - 1];
		goto L99;
	    }
/* L50: */
	}
L99:
/* L100: */
	;
    }
} /* sqcom_ */

/*   SQREV */
/* Subroutine */ int sqrev_(seqnce, idim, seqnce_len)
char *seqnce;
integer *idim;
ftnlen seqnce_len;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iend;
    static char temp[1];
    static integer i__;

/*   AUTHOR: RODGER STADEN */
/*   REVERSE THE SEQUENCE */
    /* Parameter adjustments */
    --seqnce;

    /* Function Body */
    iend = *idim / 2;
    i__1 = iend;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)temp = *(unsigned char *)&seqnce[i__];
	*(unsigned char *)&seqnce[i__] = *(unsigned char *)&seqnce[*idim + 1 
		- i__];
	*(unsigned char *)&seqnce[*idim + 1 - i__] = *(unsigned char *)temp;
/* L100: */
    }
    return 0;
} /* sqrev_ */

/*  ROUTINES TO CONTROL CHARACTER LOOKUP */
/*  FOR BOTH DNA AND PROTEIN SEQUENCES */
/*  THE INITIALISING ROUTINES ARE SENT THE CHARACTERSET SIZE IDM */
/*  WHICH DETERMINES WHICH CHARACTERSET IS USED */
/* Subroutine */ int initlu_(idm)
integer *idm;
{
    /* Initialized data */

    static char dup[16+1] = "TCAG-RYWSMKHBVDN";
    static char pup[26+1] = "CSTPAGNDEQBZHRKMILVFYW-X? ";
    static char dlow[16+1] = "tcag-rywsmkhbvdn";
    static char plow[26+1] = "cstpagndeqbzhrkmilvfyw-x? ";

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer i__, j;

    /* Fortran I/O blocks */
    static cilist io___383 = { 0, 6, 0, 0, 0 };


/*  AUTHOR RODGER STADEN */
/*  ICHAR RETURNS THE COLLATING SEQUENCE NUMBER */
/*  I WANT 1-5 FOR ACGT OR 1-26 FOR AMINO ACIDS BY USING ICHAR. */
/*  THE ACTUAL VALUE RETURNED BY ICHAR IS NOT PORTABLE */
/*  SO I NEED TO INITIALIZE POINTR SO THAT THE CORRECT */
/*  ELEMENTS CONTAIN VALUES 1 - 5, OR 1 - 26 */
/*  WORKS ON UPPER AND LOWER CASE - REMOVE DLOW,PLOW AND LOOPS 41 AND 51 */
/*  IF LOWERCASE NOT ALLOWED */

    if (*idm == 5) {
	for (i__ = 0; i__ <= 255; ++i__) {
	    iasci1_1.point1[i__] = *idm;
	    iasci2_1.point2[i__] = 17;
/* L30: */
	}
	for (i__ = 1; i__ <= 5; ++i__) {
	    j = *(unsigned char *)&dup[i__ - 1];
	    iasci1_1.point1[j] = i__;
/* L35: */
	}
	for (i__ = 1; i__ <= 5; ++i__) {
	    j = *(unsigned char *)&dlow[i__ - 1];
	    iasci1_1.point1[j] = i__;
/* L36: */
	}
	for (i__ = 1; i__ <= 16; ++i__) {
	    j = *(unsigned char *)&dup[i__ - 1];
	    iasci2_1.point2[j] = i__;
/* L40: */
	}
/*  DEAL WITH U */
	j = 'U';
	iasci1_1.point1[j] = 1;
	iasci2_1.point2[j] = 1;
	for (i__ = 1; i__ <= 16; ++i__) {
	    j = *(unsigned char *)&dlow[i__ - 1];
	    iasci2_1.point2[j] = i__;
/* L41: */
	}
/*  DEAL WITH U */
	j = 'u';
	iasci1_1.point1[j] = 1;
	iasci2_1.point2[j] = 1;
    } else if (*idm == 26) {
	for (i__ = 0; i__ <= 255; ++i__) {
	    iasci1_1.point1[i__] = *idm;
/* L45: */
	}

	for (i__ = 1; i__ <= 26; ++i__) {
	    j = *(unsigned char *)&pup[i__ - 1];
	    iasci1_1.point1[j] = i__;
/* L50: */
	}
	for (i__ = 1; i__ <= 26; ++i__) {
	    j = *(unsigned char *)&plow[i__ - 1];
	    iasci1_1.point1[j] = i__;
/* L51: */
	}
	for (i__ = 0; i__ <= 255; ++i__) {
	    iasci2_1.point2[i__] = iasci1_1.point1[i__];
/* L60: */
	}
    } else {
	s_wsle(&io___383);
	do_lio(&c__9, &c__1, "ERROR INITIALISING CHARACTER LOOKUP POINTERS", (
		ftnlen)44);
	e_wsle();
    }
} /* initlu_ */

integer indexa_(string, id, char__, string_len, char_len)
char *string;
integer *id;
char *char__;
ftnlen string_len;
ftnlen char_len;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;

/*  FUNCTION TO FIND FIRST OCCURRENCE OF CHAR IN STRING */
    /* Parameter adjustments */
    --string;

    /* Function Body */
    i__1 = *id;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&string[i__] == *(unsigned char *)char__) {
	    ret_val = i__;
	    return ret_val;
	}
/* L10: */
    }
    ret_val = 0;
    return ret_val;
} /* indexa_ */

/* Subroutine */ int fillc_(seq, idim, ch, seq_len, ch_len)
char *seq;
integer *idim;
char *ch;
ftnlen seq_len;
ftnlen ch_len;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seq;

    /* Function Body */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)&seq[i__] = *(unsigned char *)ch;
/* L10: */
    }
    return 0;
} /* fillc_ */

/*   FILLI */
/* Subroutine */ int filli_(seq, idim, ch)
integer *seq, *idim, *ch;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --seq;

    /* Function Body */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	seq[i__] = *ch;
/* L10: */
    }
    return 0;
} /* filli_ */

/*     SQCOPY */
/*   SEQUENCE COPYING PROGRAM */
/* Subroutine */ int sqcopy_(seqnce, comseq, idim, seqnce_len, comseq_len)
char *seqnce, *comseq;
integer *idim;
ftnlen seqnce_len;
ftnlen comseq_len;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --comseq;
    --seqnce;

    /* Function Body */
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)&comseq[i__] = *(unsigned char *)&seqnce[i__];
/* L100: */
    }
    return 0;
} /* sqcopy_ */

integer ctonum_(char__, char_len)
char *char__;
ftnlen char_len;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer icol;

/*  AUTHOR RODGER STADEN */

/*  GET COLLATING SEQUENCE VALUE */
    icol = *(unsigned char *)char__;
/*  THIS POINTS TO A VALUE IN POINTR */
    ret_val = iasci1_1.point1[icol];
    return ret_val;
} /* ctonum_ */

/* Subroutine */ int rjstfy_(array, string, lens, length, array_len, 
	string_len)
char *array, *string;
integer *lens, *length;
ftnlen array_len;
ftnlen string_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer i__, k, k1, k3;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --array;

    /* Function Body */
    s_copy(string, " ", string_len, (ftnlen)1);
/*   LOOK FOR FIRST NON SPACE CHAR */
    k = *length + 1;
    i__1 = *length;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--k;
/* L1: */
	if (*(unsigned char *)&array[k] != ' ') {
	    goto L2;
	}
    }
/*   ALL SPACES! */
    return 0;
L2:
    k1 = k;
/*  POINT TO RIGHT END OF STRING */
    k3 = *lens + 1;
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--k3;
	*(unsigned char *)&string[k3 - 1] = *(unsigned char *)&array[k];
/* L3: */
	--k;
    }
    return 0;
} /* rjstfy_ */

/* Subroutine */ int bub3as_(list, list2, list3, idim)
integer *list, *list2, *list3, *idim;
{
    static integer i__, j, itemp;

/*   AUTHOR: RODGER STADEN */
    /* Parameter adjustments */
    --list3;
    --list2;
    --list;

    /* Function Body */
    i__ = 0;
    j = 0;
L10:
/*   SET I=J IF WE HAVE JUST CORRECTLY POSITIONED AN ELEMENT */
    if (j > i__) {
	i__ = j;
    }
    ++i__;
    if (i__ == *idim) {
	return 0;
    }
L20:
    if (list[i__] <= list[i__ + 1]) {
	goto L10;
    }
/*   FIRST MOVE THIS ELEMENT? IF SO SET POINTER TO ITS INITIAL POSITION */
    if (j < i__) {
	j = i__;
    }
    itemp = list[i__];
    list[i__] = list[i__ + 1];
    list[i__ + 1] = itemp;
    itemp = list2[i__];
    list2[i__] = list2[i__ + 1];
    list2[i__ + 1] = itemp;
    itemp = list3[i__];
    list3[i__] = list3[i__ + 1];
    list3[i__ + 1] = itemp;
/*   DECREMENT BACK THRU LIST WITH THIS ELEMENT */
    if (i__ > 1) {
	--i__;
    }
    goto L20;
} /* bub3as_ */

integer notrl_(text, itext, word, text_len, word_len)
char *text;
integer *itext;
char *word;
ftnlen text_len;
ftnlen word_len;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;

/*   AUTHOR: RODGER STADEN */
/*   LOOKS RIGHT TO LEFT THRU TEXT FOR FIRST ELEMENT THAT IS NOT WORD */
/*   RETURNS ELEMENT NUMBER OR ZERO IF ALL ELEMENTS ARE WORD */
    for (i__ = *itext; i__ >= 1; --i__) {
	if (*(unsigned char *)&text[i__ - 1] != *(unsigned char *)word) {
	    ret_val = i__;
	    return ret_val;
	}
/* L1: */
    }
    ret_val = 0;
    return ret_val;
} /* notrl_ */

