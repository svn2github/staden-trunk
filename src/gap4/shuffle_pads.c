#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "qual.h"
#include "IO.h"
#include "align.h"
#include "dna_utils.h"
#include "align_lib.h"

#define DEBUG
#define CHECKING

/*
 * Converts a basecall into 0-4 inclusive for A,C,G,T,*,-.
 * Anything else is treated as '-'.
 */
static int base2ind(char base) {
    switch(base) {
    case 'a':
    case 'A':
	return 0;

    case 'c':
    case 'C':
	return 1;

    case 'g':
    case 'G':
	return 2;

    case 't':
    case 'T':
	return 3;

    case '*':
	return 4;
    }

    return 5;
}


/*
 * Adds a single sequence to the vector array
 */
static void add_to_vector(FastInt (*vec)[6], char *seq, int pos, int len) {
    int i;
    for (i = 0; i < len; i++, pos++) {
	vec[pos][base2ind(seq[i])]++;
    }
}

/*
 * Removes a single sequence from the vector array
 */
static void remove_from_vector(FastInt (*vec)[6], char *seq, int pos, int len)
{
    int i;
    for (i = 0; i < len; i++, pos++) {
	if (--vec[pos][base2ind(seq[i])] < 0) {
	    printf("ERROR: negative item in vector at pos %d\n", pos);
	}
    }
}

void dump_vector(FastInt (*vec)[6], int len) {
    int i, j;
    char line[6][80];

    for (i = 0; i < len; i++) {
	if ((i % 60) == 0) {
	    puts("");
	    for (j = 0; j < 6; j++) {
		if (i)
		    printf("%.60s\n", line[j]);
		memset(line[j], 0, 80);
	    }
	}
	for (j = 0; j < 6; j++)
	    line[j][i%60] = vec[i][j]+'0';
    }

    puts("");
    for (j = 0; j < 6; j++)
	printf("%.*s\n", len%60, line[j]);
}


/*
 * Builds an array of 6-wide vectors of A,C,G,T,*,- for each consensus
 * position in the specified contig.
 *
 * Returns vectors on success.
 *         NULL on failure.
 */
FastInt (*build_vector(int contig,
		       int (*info_func)(int         job,
					void       *mydata,
					info_arg_t *theirdata),
		       void *info_data))[6] {
    info_arg_t info, iseq;
    FastInt (*c6)[6];
    int c6_len;

    /* Get contig info */
    info.contig_info.contig = contig;
    info_func(GET_CONTIG_INFO, info_data, &info);

    c6_len = info.contig_info.length;
    c6 = (FastInt (*)[6])xcalloc(c6_len, 6*sizeof(FastInt));
#ifdef DEBUG
    printf("Contig %d, len %d\n",
	   info.contig_info.contig,
	   info.contig_info.length);
#endif

    /* Loop through all sequences in the contig */
    info.gel_info.gel = info.contig_info.leftgel;
    do {
	/* Get sequence and add it to the consensus vector */
	info_func(GET_GEL_INFO, info_data, &info);

	iseq.gel_seq.gel = info.gel_info.gel;
	info_func(GET_SEQ, info_data, &iseq);
	
	add_to_vector(c6, &iseq.gel_seq.gel_seq[iseq.gel_seq.gel_start],
		      info.gel_info.position-1, info.gel_info.length);
	
	info_func(DEL_SEQ, info_data, &iseq);
	
    } while (info.gel_info.gel = info.gel_info.next_right);

    /* dump_vector(c6, c6_len); */

    return c6;
}

int edit_mseqs_old(MALIGN *malign, CONTIGL *cl, MOVERLAP *o, int cons_pos) {
    int i, j, poso, posn, len;
    char *oseq;
    char *cp;
    int leno, lenn;

    printf("EDIT_MSEQS, s1_len=%d s2_len=%d\n",
	   o->s1_len, o->s2_len);

    /* Cons vector */
    for (poso = i = 0; i < o->s1_len; i++) {
	if (o->S1[i] < 0) {
	    printf("S1:Ins %d pads as pos %d\n",
		   -o->S1[i], poso);
	} else {
	    poso += o->S1[i];
	}
    }

    /* sequence */
    xfree(cl->mseg->seq);
    cl->mseg->seq = strdup(o->seq2_out);

#if 0
    /* NOT NEEDED. we just replace verbatim with the newly aligned seq. */

    /*
     * Now compare sequence to old sequence working out whether the pads added
     * are new or not.
     */
    for (i = 0; i < o->seq_out_len; i++)
	if (o->seq2_out[i] == '.')
	    o->seq2_out[i] = '*';

    /*    for (i = j = 0; i < poso && j < posn; i++, j++) {*/
    leno = cl->mseg->length;
    lenn = o->seq_out_len;
    cl->mseg->seq = (char *)xrealloc(cl->mseg->seq, MAX(lenn, leno)+2);
    for (i = leno-1, j = lenn-1; i >= 0 && j >= 0; i--, j--) {
	while (o->seq2_out[j] != cl->mseg->seq[i] && i >= 0 && j >= 0) {
	    if (cl->mseg->seq[i] == '*') {
		printf("EDIT: deletion at %d\n", j);
		memmove(&cl->mseg->seq[j], &cl->mseg->seq[j+1], leno-(j+1));
		leno--;
		i--;
		printf("mseg now %s\n", cl->mseg->seq);
	    } else if (o->seq2_out[j] == '*') {
		printf("EDIT: insertion at %d\n", j);
		memmove(&cl->mseg->seq[j+1], &cl->mseg->seq[j], leno-j);
		cl->mseg->seq[j] = '+';
		leno++;
		j--;
		printf("mseg now %s\n", cl->mseg->seq);
	    } else {
		printf("EDIT: unknown diff at %d, %c/%c\n", j,
		       o->seq2_out[j], cl->mseg->seq[i]);
		break;
	    }
	}
    }

    printf("Nseq = %s %d %d\n", cl->mseg->seq, leno, lenn);
#endif

    return 0;
}

/*
 * Insert 'size' pads into a contig at position 'pos'.
 */
void malign_padcon(MALIGN *malign, int pos, int size) {
    CONTIGL *cl = malign->contigl;

    for (; cl; cl = cl->next) {
	/* We do one of three things: nothing, insert, or shift */
	/* Nothing: */
	if (cl->mseg->offset+cl->mseg->length-1 < pos)
	    continue;

	/* Shift right: */
	if (cl->mseg->offset >= pos) {
	    cl->mseg->offset += size;
	    continue;
	}

	/* Insert */
	cl->mseg->length += size;
	cl->mseg->seq = (char *)realloc(cl->mseg->seq, cl->mseg->length+1);
	memmove(&cl->mseg->seq[pos - cl->mseg->offset + size],
		&cl->mseg->seq[pos - cl->mseg->offset],
		cl->mseg->length-size - (pos - cl->mseg->offset));
	memset(&cl->mseg->seq[pos - cl->mseg->offset], '*', size);
	cl->mseg->seq[cl->mseg->length] = 0;
    }

    /* Update malign... */
    /* TODO */
    /* malign->length += size; */
}

/*
 * Returns the number of consensus pads added or -1 for error.
 */
int edit_mseqs(MALIGN *malign, CONTIGL *cl, MOVERLAP *o, int cons_pos) {
    int i, j, poso, posn, len;
    char *oseq, *cp;
    int leno, lenn, npads;

    printf("EDIT_MSEQS, s1_len=%d s2_len=%d\n",
	   o->s1_len, o->s2_len);
    /* Cons vector */
    npads = 0;
    for (poso = i = 0; i < o->s1_len; i++) {
	if (o->S1[i] < 0) {
	    printf("S1:Ins %d pads at pos %d+%d=%d\n",
		   -o->S1[i], poso, cons_pos, poso+cons_pos);
	    malign_padcon(malign, poso+cons_pos+npads, -o->S1[i]);
	    npads += -o->S1[i];
	} else {
	    poso += o->S1[i];
	}
    }

    /* sequence */
    xfree(cl->mseg->seq);
    cl->mseg->seq = strdup(o->seq2_out);
    for (cp = cl->mseg->seq; *cp; cp++) {
	if (*cp == '.')
	    *cp = '*';
    }

    return npads;
}

/*
 * Given an alignment of the depadded sequence to the consensus vector, this
 * compares the edits back ot the unpadded sequence and makes any edits to the
 * sequence and/or consensus as appropriate.
 *
 * seq1 and seq2, with lengths len1 and len2, were the inputs to the
 * alignment algorithm, which returned S as the edit buffer.
 * gel_info and gel_seq contain information about the original sequence
 * before unpadded (eg as it is currently held in the database).
 *
 * Returns the number of additional bases in the consensus.
 */
int edit_seqs(int contig,
	      int cons_pos,
	      FastInt (**seq2_ptr)[6],
	      int *seq2_len,
	      char *seq1,
	      FastInt (*seq2)[6],
	      int len1,
	      int len2,
	      align_int *S,
	      info_arg_t *gel_info,
	      info_arg_t *gel_seq,
	      int (*info_func)(int         job,
			       void       *mydata,
			       info_arg_t *theirdata),
	      void *info_data) {
    int s1pos, s2pos;
    int o1pos;
    int op = 0;
    char *oseq;
    int offset = 0;
    info_arg_t info;
    int cons_pads = 0;

    s1pos = s2pos = 0;
    o1pos = gel_info->gel_info.start;
    oseq = gel_seq->gel_seq.gel_seq;

    while (s1pos < len1 && s2pos < len2) {
	if (op == 0)
	    op = *S++;

	if (op == 0) {
	    /*
	     * No insertion between seq and consensus, but if there was a pad
	     * in the original sequence then we need to delete it.
	     */
	    while (oseq[o1pos] == '*') {
#ifdef DEBUG
		printf("Del at pos %d\n", o1pos-(gel_info->gel_info.start-1));
#endif
		info.seq_del.gel = gel_info->gel_info.gel;
		info.seq_del.length = 1;
		info.seq_del.position = o1pos+offset-gel_info->gel_info.start+1;
		info_func(SEQ_DEL, info_data, &info);
		offset--;
		o1pos++;
	    }
	    /*
#ifdef DEBUG
	    printf("%d %d %c(%c) %d%d%d%d%d%d\n",
		   s1pos, s2pos, seq1[s1pos], oseq[o1pos],
		   seq2[s2pos][0], seq2[s2pos][1], seq2[s2pos][2], 
		   seq2[s2pos][3], seq2[s2pos][4], seq2[s2pos][5]);
#endif
	    */
	    o1pos++;
	    s1pos++;
	    s2pos++;
	} else if (op > 0) {
	    /*
	     * Insertion to sequence. If there's not already a pad at this
	     * point in the original sequence then we add one.
	     */

	    /*
#ifdef DEBUG
	    printf("%d %d +(%c) %d%d%d%d%d%d\n",
		   s1pos, s2pos, oseq[o1pos],
		   seq2[s2pos][0], seq2[s2pos][1], seq2[s2pos][2], 
		   seq2[s2pos][3], seq2[s2pos][4], seq2[s2pos][5]);
#endif
	    */
	    s2pos++;
	    if (oseq[o1pos] != '*') {
#ifdef DEBUG
		printf("Ins at pos %d\n", o1pos-(gel_info->gel_info.start-1));
#endif
		info.seq_ins.gel = gel_info->gel_info.gel;
		info.seq_ins.bases = "*";
		info.seq_ins.length = 1;
		info.seq_ins.position = o1pos+offset-gel_info->gel_info.start+1;
		info_func(SEQ_INS, info_data, &info);
		offset++;
	    } else {
		o1pos++;
	    }
	    op--;
	} else {
	    /*
	     * Insertion to consensus.
	     * This is tricky to handle - we need to pad the consensus (moving
	     * sequences as required, but that is handled for us), but delete
	     * the pad from the sequence being aligned.
	     * We also need to keep the 6-wide consensus vector up to date.
	     */
	    int depth;
	    int cpos;
	    int b;
	    off_t dist;

#ifdef DEBUG
	    printf("%d %d %c(%c) ++++++\n",
		   s1pos, s2pos, seq1[s1pos], oseq[o1pos]);
#endif

	    info.cons_ins.contig = contig;
	    info.cons_ins.bases = "*";
	    info.cons_ins.length = 1;
	    info.cons_ins.position = (cpos = s2pos + cons_pos+1) ;
	    info_func(CONS_INS, info_data, &info);
	    cons_pads++;

	    /* Also remove pad from sequence - we were padding all others */
	    info.seq_del.gel = gel_info->gel_info.gel;
	    info.seq_del.length = 1;
	    info.seq_del.position = cpos-gel_info->gel_info.position+1;
	    info_func(SEQ_DEL, info_data, &info);

	    /* Also update 'c6' with insertion */
	    *seq2_len = *seq2_len+1;
	    dist = seq2 - *seq2_ptr;
	    *seq2_ptr = (FastInt (*)[6])xrealloc(*seq2_ptr,
						 *seq2_len * 6 *
						 sizeof(FastInt));
	    cpos--; /* our arrays count from zero */
	    seq2 = *seq2_ptr + dist;
	    for (b = depth = 0; b < 6; b++)
		depth += (*seq2_ptr)[cpos][b];
	    memmove(&((*seq2_ptr)[cpos+1]),
		    &((*seq2_ptr)[cpos]),
		    (*seq2_len-1 - cpos) * 6 * sizeof(FastInt));
	    for (b = 0; b < 6; b++)
		(*seq2_ptr)[cpos][b] = 0;
	    (*seq2_ptr)[cpos][4] = depth;
	    

	    o1pos++;
	    s1pos++;
	    op++;
	}
    }

#ifdef DEBUG
    printf("Added %d pads\n", cons_pads);
#endif
    return cons_pads;
}

/*
 * Iterates through all sequences in a contig realigning them against the
 * consensus vector.
 *
 * It then adds the newly aligned sequence back into the consensus, editing the
 * sequence and tag positions/lengths too.
 * To do this we may need to shuffle the start position of sequences
 * downstream, and hence also move consensus tags.
 */
void realign_seqs_old(int contig,
		  FastInt (*c6)[6], 
		  int (*info_func)(int         job,
				   void       *mydata,
				   info_arg_t *theirdata),
		  void *info_data,
		  MALIGN *malign) {
    info_arg_t info, iseq, cinfo;
    int c6_len;
    int len;
    int cons_pos, cons_len;
    int pads;
    static int Wmat[128][128];
    static int Wmat_init = 0;
    static int **nt_matrix = NULL;
    static char *nt_order = "ACGTURYMWSKDHVB-*";
    CONTIGL *cl = malign->contigl;

    if (!Wmat_init) {
	char buf[1024];
	char *env = getenv("STADTABL");

	sprintf(buf, "%s/shuffle_pads_matrix", env);
	nt_matrix = create_matrix(buf, nt_order);
	if (nt_matrix)
	    init_align_mat(nt_matrix, nt_order, 0, Wmat);
	else
	    verror(ERR_FATAL, "init_globals",
		   "%s: file not found", buf);
	
	Wmat_init = 1;
    }

    /* Get contig info */
    cinfo.contig_info.contig = contig;
    info_func(GET_CONTIG_INFO, info_data, &cinfo);

    c6_len = cinfo.contig_info.length;

    /* Loop through all sequences in the contig */
    info.gel_info.gel = cinfo.contig_info.leftgel;
    do {
	align_int *alignment;
	char *depadded_seq;
	int *depad_to_pad;
	
	/* Get and depad the sequence */
	info_func(GET_GEL_INFO, info_data, &info);
	if (info.gel_info.next_right == 0 &&
	    cinfo.contig_info.leftgel == info.gel_info.gel) {
#ifdef DEBUG
	    printf("Gel %d is a single read contig; skipping.\n",
		   info.gel_info.gel);
#endif
	    break;
	}

#ifdef DEBUG
	printf("Gel %d, pos %d, len %d, seq ",
	       info.gel_info.gel,
	       info.gel_info.position,
	       info.gel_info.length);
#endif

	iseq.gel_seq.gel = info.gel_info.gel;
	info_func(GET_SEQ, info_data, &iseq);
	{
	    int i;
	    for (i = 0; i < info.gel_info.unclipped_len; i++) {
		iseq.gel_seq.gel_seq[i] = toupper(iseq.gel_seq.gel_seq[i]);
	    }
	}
#ifdef DEBUG
	printf("%.*s\n",
	       MIN(20, iseq.gel_info.length),
	       &iseq.gel_seq.gel_seq[iseq.gel_seq.gel_start]);
#endif
	
	len = info.gel_info.length;
#if 0
	cons_pos = MAX(info.gel_info.position-3,0);
	cons_len = info.gel_info.position + len-1 + 3;
	if (cons_len >= c6_len)
	    cons_len = c6_len-1;
	cons_len = cons_len - cons_pos;
#endif
	cons_pos = info.gel_info.position-1;
	cons_len = len;

	depad_to_pad = (int *)xmalloc((len+1)*sizeof(int));
	depadded_seq = (char *)xmalloc(len+1);
	strncpy(depadded_seq,
		&iseq.gel_seq.gel_seq[iseq.gel_seq.gel_start], len);
	depadded_seq[len] = 0;
	depad_seq(depadded_seq, &len, depad_to_pad);

	/* Take seq out of consensus vector and realign it */
	remove_from_vector(c6, &iseq.gel_seq.gel_seq[iseq.gel_seq.gel_start],
			   info.gel_info.position-1, info.gel_info.length);

	alignment = (align_int *)xmalloc((cons_len*2+1)*sizeof(align_int));
	calignm(depadded_seq, &c6[cons_pos], len, cons_len,
		NULL, NULL, NULL, NULL,
		-10, +10, 1, 8,
		ALIGN_J_SV | ALIGN_GAP_S1 | ALIGN_GAP_E1 | ALIGN_GAP_S2 |
		ALIGN_GAP_E2,
		0, alignment, Wmat);

	/*
	 * If the first value in the alignment is padding, then we need to
	 * shift the seq left or right appropriately instead of adding a
	 * bunch of pads.
	 * 
	 * FIXME: still need to implement this.
	 */

#ifdef DEBUG
	vmessage("Seq %d at pos %d\n", info.gel_info.gel, cons_pos);
	cdisplay(depadded_seq, &c6[cons_pos], len, cons_len, ALIGN_J_SV,
		 alignment, 0, 0);
#endif

	pads = edit_seqs(contig, cons_pos, &c6, &c6_len,
			 depadded_seq, &c6[cons_pos], len, cons_len,
			 alignment, &info, &iseq, info_func, info_data);

	xfree(alignment);
	xfree(depadded_seq);
	xfree(depad_to_pad);
	info_func(DEL_SEQ, info_data, &iseq);

	/* Contig length may have changed */
	cinfo.contig_info.contig = contig;
	info_func(GET_CONTIG_INFO, info_data, &cinfo);
	if (c6_len < cinfo.contig_info.length+1) {
	    c6_len = cinfo.contig_info.length+1;
	    c6 = (FastInt (*)[6])xrealloc(c6, c6_len * 6 * sizeof(FastInt));
	}


	/* edit_seqs may have changed it, so read new before adding back */
	info_func(GET_GEL_INFO, info_data, &info);
	iseq.gel_seq.gel = info.gel_info.gel;
	info_func(GET_SEQ, info_data, &iseq);
	add_to_vector(c6, &iseq.gel_seq.gel_seq[iseq.gel_seq.gel_start],
		      info.gel_info.position-1, info.gel_info.length);
	info_func(DEL_SEQ, info_data, &iseq);

	/* FIXME - hideously slow, but a test for where the bug is */
	if (1) {
	    info_arg_t c;
	    FastInt (*c6v2)[6];

	    c.contig_info.contig = contig;
	    info_func(GET_CONTIG_INFO, info_data, &c);

	    c6_len = c.contig_info.length;
	    c6v2 = build_vector(contig, info_func, info_data);

#ifdef CHECKING	    
	    {
		int i, j;
		for (i = 0; i < c6_len; i++) {
		    for (j = 0; j < 6; j++) {
			if (c6[i][j] != c6v2[i][j]) {
			    printf("Diff at pos %d: ", i);
			    for (j = 0; j < 6; j++) {
				printf("%d/%d ",
				       c6[i][j],
				       c6v2[i][j]);
			    }
			    putchar('\n');
			    break;
			}
		    }
		}
	    }
#endif
	    xfree(c6);
	    c6 = c6v2;
	}

	info_func(IF_FLUSH, info_data, NULL);

	cl = cl->next;
	
    } while (info.gel_info.gel = info.gel_info.next_right);
}

/*
 * Iterates through all sequences in a contig realigning them against the
 * consensus vector.
 *
 * It then adds the newly aligned sequence back into the consensus, editing the
 * sequence and tag positions/lengths too.
 * To do this we may need to shuffle the start position of sequences
 * downstream, and hence also move consensus tags.
 */
void realign_seqs(int contig,
		  MALIGN *malign) {
    CONTIGL *lastl = NULL, *contigl;
    int nsegs;
    int seg_num = 0;

    for (contigl = malign->contigl, nsegs = 0; contigl; nsegs++)
	contigl = contigl->next;

    /* Loop through all sequences in the contig */
    contigl = malign->contigl;
    while (contigl) {
	int *depad_to_pad;
	char *depadded_seq;
	int len;
	MOVERLAP *o;
	ALIGN_PARAMS *p;
	int cons_pos;
	int npads, band;

	printf("Seq %d/%d, contigl = %p\n", seg_num++, nsegs, contigl);

	/* Obtain a depadded copy of this mseg */
	len = contigl->mseg->length;
	depad_to_pad = (int *)xmalloc((len+1)*sizeof(int));
	depadded_seq = (char *)xmalloc(len+1);
	strncpy(depadded_seq, contigl->mseg->seq, len);
	depadded_seq[len] = 0;
	depad_seq(depadded_seq, &len, depad_to_pad);
	npads = contigl->mseg->length - len;
	band = MAX(20, npads);
	band = MIN(1000, band);

#if 0
	/* Remove sequence from malign */
	if (lastl) {
	    lastl->next = contigl->next;
	} else {
	    malign->contigl = contigl->next;
	}
#endif

	/* Recalc scores (create_malign_counts) and rescale
	 * (scale_malign_scores) over the region we have changed.
	 * FIXME: TODO
	 */

	/* Align sequence to malign */
	p = create_align_params();
	set_align_params (p,
			  band, /*band*/
			  1, /*gap_open*/
			  8, /*gap_extend*/
			  EDGE_GAPS_COUNT,
			  /* EDGE_GAPS_ZERO | BEST_EDGE_TRACE, */
			  RETURN_EDIT_BUFFERS | RETURN_SEQ |
			  RETURN_NEW_PADS,
			  0,  /*seq1_start*/
			  0,  /*seq2_start*/
			  0,  /*old pad sym*/
			  0,  /*new pad sym*/
			  0   /*set_job*/);

	o = create_moverlap();
	init_moverlap(o, malign, depadded_seq, malign->length, len);

	cons_pos = contigl->mseg->offset;
	o->malign_len = malign->length - cons_pos;
	if (o->malign_len > contigl->mseg->length)
	    o->malign_len = contigl->mseg->length;
	malign->consensus += cons_pos;
	malign->scores += cons_pos;
	    
	printf("RESULT=%d\n", affine_malign(o, p));
	printf("score= %f\n", o->score);

	
	malign->consensus -= cons_pos;
	malign->scores -= cons_pos;

	edit_mseqs(malign, contigl, o, cons_pos);

	print_moverlap(malign, o, cons_pos);
	print_malign(malign);

#if 0
	/* Put sequence back */
	if (lastl) {
	    lastl->next = contigl;
	} else {
	    malign->contigl = contigl;
	}
#endif

	/* Memory leak, slow, etc etc etc, but a test */
	{
	    if ( malign->msegs ) xfree ( malign->msegs );
	    if ( malign->consensus ) xfree ( malign->consensus );
	    destroy_malign_counts(malign->matrix,
				  malign->charset_size,
				  malign->charset_size);
	    destroy_malign_counts(malign->scores,
				  malign->length,
				  malign->charset_size);
	    xfree ( malign );
	}
	malign = contigl_to_malign(malign->contigl);

	if (0) {
	    /* new consensus - is this used? */
	    malign->length = contigl_length(malign->contigl);
	    get_malign_consensus(malign);
	}

	destroy_moverlap(o);
	destroy_alignment_params(p); 

	xfree(depadded_seq);
	xfree(depad_to_pad);

	lastl = contigl;
	contigl = contigl->next;
    }

    print_malign(malign);
}

/**
 * Builds and returns MALIGN from a Gap4 IO handle for the contig 'cnum'.
 */
MALIGN *build_malign(GapIO *io, int cnum) {
    CONTIGL *contig, *first_contig, *last_contig = NULL;
    GContigs c;
    GReadings r;
    int rnum;

    /* Generate contigl linked list */
    contig_read(io, cnum, c);
    for (rnum = c.left; rnum; rnum = r.right) {
	char *seq;
	gel_read(io, rnum, r);
	contig = create_contig_link();
	contig->mseg = create_mseg();
	seq = TextAllocRead(io, r.sequence);
	seq[r.start + r.sequence_length] = 0;
	init_mseg(contig->mseg, strdup(seq+r.start),
		  r.sequence_length, r.position-1);
	xfree(seq);
	if (last_contig) {
	    last_contig->next = contig;
	} else {
	    first_contig = contig;
	}
	last_contig = contig;
    }

    return contigl_to_malign(first_contig);
}

#define LLEN 60
struct clist {
    char *seq;
    int len;
    char line[LLEN];
};

void print_malign(MALIGN *malign) {
    int i, j;
    struct clist *depth = NULL;
    int ndepth = 0;
    CONTIGL *cl = malign->contigl;

    for (i = 0; i < malign->length; i++) {
	/* Maintain a list of CONTIGLs covering this point */

	/* ... adding new items to the list */
	while (cl && cl->mseg->offset <= i) {
	    ndepth++;
	    /* runaway loops completely kills deskpros */
	    if (ndepth > 1000)
		abort();
	    depth = (struct clist *)realloc(depth, ndepth * sizeof(*depth));
	    depth[ndepth-1].seq = cl->mseg->seq;
	    depth[ndepth-1].len = cl->mseg->length;
	    memset(depth[ndepth-1].line, ' ', LLEN);
	    cl = cl->next;
	}

	for (j = 0; j < ndepth; j++) {
	    depth[j].line[i%LLEN] = (depth[j].seq) ? *depth[j].seq++ : ' ';
	    if (depth[j].len > 0 && --depth[j].len == 0) {
		depth[j].seq = NULL;
	    }
	}

	/* Print line, and remove items from depth as and when needed */
	if (i%LLEN == LLEN-1) {
	    for (j = LLEN * (int)(i/LLEN); j < i; j+=10)
		printf("%10d", j+10);
	    printf("\n");
	    for (j = 0; j < ndepth; j++) {
		printf("%.*s\n", LLEN, depth[j].line);
		if (!depth[j].seq) {
		    memmove(&depth[j], &depth[j+1],
			    (ndepth-(j+1)) * sizeof(depth[j]));
		    ndepth--;
		    j--;
		}
	    }
	    printf("\n");
	}
    }

    /* Print remainder of lines */
    if ((i-1)%LLEN != LLEN-1) {
	for (j = LLEN * (int)(i/LLEN); j < i; j+=10)
	    printf("%10d", j+10);
	printf("\n");
	for (j = 0; j < ndepth; j++) {
	    printf("%.*s\n", i - LLEN * (int)(i/LLEN), depth[j].line);
	}
	printf("\n");
    }

    free(depth);
}

void print_moverlap(MALIGN *malign, MOVERLAP *o, int offset) {
    int i, j;
    struct clist *depth = NULL;
    int ndepth = 0;
    CONTIGL *cl = malign->contigl;
    int s1op = 0, s2op = 0;
    int *S1 = o->S1;
    int *S2 = o->S2;
    char *seq = o->seq2;
    int cins = 0;

    for (i = offset; i < malign->length+offset; i++) {
	/* Maintain a list of CONTIGLs covering this point */

	/* ... adding new items to the list */
	for (; cl && cl->mseg->offset+cins <= i; cl = cl->next) {
	    if (cl->mseg->offset+cins + cl->mseg->length-1 < i)
		continue;
	    ndepth++;
	    /* runaway loops completely kills deskpros */
	    if (ndepth > 1000)
		abort();
	    depth = (struct clist *)realloc(depth, ndepth * sizeof(*depth));
	    depth[ndepth-1].seq = cl->mseg->seq + i-(cl->mseg->offset+cins);
	    depth[ndepth-1].len = cl->mseg->length - (i-(cl->mseg->offset+cins));
	    memset(depth[ndepth-1].line, ' ', LLEN);
	}

	if (!s1op) {
	    s1op = *S1++;
	    if (S1-o->S1 > o->s1_len)
		break;
	}
	if (!s2op) {
	    s2op = *S2++;
	    if (S2-o->S2 > o->s2_len)
		break;
	}

	printf("%4d: ", i);

	if (s1op < 0) {
	    /* Ins to consensus */
	    s1op++;
	    printf("%c\n", *seq++);
	    cins++;
	    continue;
	} else if (s2op > 0) {
	    /* Match/mismatch */
	    printf("%c ", *seq++);
	    s2op--;
	} else if (s2op < 0) {
	    /* Ins to sequence */
	    printf("  ");
	    s2op++;
	}

	s1op--;
	for (j = 0; j < ndepth; j++) {
	    printf("%c", *depth[j].seq++);
	    if (--depth[j].len == 0) {
		depth[j].seq = NULL;
		memmove(&depth[j], &depth[j+1],
			(ndepth-(j+1)) * sizeof(depth[j]));
		ndepth--;
		j--;
	    }
	}
	printf("\n");
    }

    free(depth);
}

int shuffle_contigs_io(GapIO *io) {
    int cnum;
    FastInt (*c6)[6];
    MALIGN *malign;
    
    set_malign_lookup(5);
    set_alignment_matrix("nuc_matrix", "ACGTURYMWSKDHVB-*");

#if 0
    for (cnum = 1; cnum <= NumContigs(io); cnum++) {
	malign = build_malign(io, cnum);
	print_malign_seqs(malign);
	/* print_malign(malign); */
	c6 = build_vector(cnum, database_info, (void *)io);
	realign_seqs(cnum, c6,  database_info, (void *)io, malign);

	/* xfree(c6);*/
    }
#endif

    for (cnum = 1; cnum <= NumContigs(io); cnum++) {
	malign = build_malign(io, cnum);
	/* print_malign_seqs(malign); */
	puts("==PASS 1==");
	realign_seqs(cnum, malign);
    }

    return 0;
}
