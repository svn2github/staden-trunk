#include <ctype.h>
#include <stdlib.h>

#include "parse_ft.h"
#include "edUtils.h"
#include "edStructs.h"
#include "array.h"
#include "misc.h"
#include "dna_utils.h"
#include "genetic_code.h"
#include "dstring.h"
#include "tman_interface.h"

/*
 * HTML templates for generation of mutation reports.
 */
static char *html_template0 =
"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"
"<html>\n"
"<head>\n"
"  <title>Mutation Report</title>\n"
"  <meta http-equiv=\"content-type\" content=\"text/html; charset=ISO-8859-1\">\n"
"</head>\n"
"<body bgcolor=\"#ffffff\">\n"
"Date: %s, Database %s, Contig %s\n"
"<h1>Summary Mutation Report</h1>\n";

static char *html_template1 =
"<h2>Sorted by %s</h2>\n"
"<table cellspacing=0>\n"
"  <thead>\n"
"    <tr>\n"
"      <th align=left>Name</th>\n"
"      <th width=20></th>\n"
"      <th align=left>Mutation</th>\n"
"      <th width=20></th>\n"
"      <th align=left>Effect</th>\n"
"      <th width=20></th>\n"
"      <th align=left>Fwd</th>\n"
"      <th width=5></th>\n"
"      <th align=left>Rev</th>\n"
"      <th width=25></th>\n"
"      <th align=left>Fwd range</th>\n"
"      <th width=25></th>\n"
"      <th align=left>Rev range</th>\n"
"    </tr>\n"
"    <tr bgcolor=\"green\"><td colspan=13></td></tr>\n"
"  </thead>\n"
"  <tbody>\n";

static char *html_template2 =
"  </tbody>\n"
"</table>\n";


static Array generate_exons(EdStruct *xx, int pos, int width);

/*
 * Feature table information. Only used locally for computing the translation
 * status lines.
 */
typedef struct {
    ft_entry *entry;
    Array ranges; /* An array of ft_range pointers (not linked up) */
    char name[DB_NAMELEN+1];
    int num; /* Feature number */
    int sense;
} ft_entry_ptr;

/*
 * Recreates the ft_entry structures from sequence annotations.
 *
 * The cons_locations flag, when true, indicates that the feature positions are
 * to be computed from the observed position on the consensus rather than
 * the actual value in the feature. This is useful when editing has taken
 * place (eg padding).
 */
static Array get_ft_entries(EdStruct *xx, int seq, int cons_locations) {
    tagStruct *t;
    Array fta;
    int i;


    fta = ArrayCreate(sizeof(ft_entry_ptr), 0);

    /* Exons are marked by FCDS tags */
    DBgetSeq(DBI(xx), seq);
    t = DBgetTags(DBI(xx), seq);
    while (t && (t = t->next)) {
	int tpos;
	int ft_num;
	int ele_num;
	int ind;
	ft_range *range;

	if (!(t->tagrec.type.c[0] == 'F' &&
	      t->tagrec.type.c[1] == 'C' &&
	      t->tagrec.type.c[2] == 'D' &&
	      t->tagrec.type.c[3] == 'S'))
	    continue;

	/* Found one */
	tpos = normalisePos2(xx, seq, t->tagrec.position, t->tagrec.length);

	/* Load up the comment, if none then it's an illegal FCDS */
	force_comment(DBI_io(xx), t);
	if (!t->newcomment) {
	    verror(ERR_WARN, "find_exons", "No comment for FCDS tag");
	    continue;
	}

	/* Parse the first line to get the FT number and range element */
	if (2 != sscanf(t->newcomment, "#FEATURE %d ELEMENT %d\n%n",
			&ft_num, &ele_num, &ind)) {
	    verror(ERR_WARN, "find_exons",
		   "Illegal comment format for FCDS tag: '%s'", t->newcomment);
	    continue;
	}


	/* Is this a duplicate of one we've seen already? */
	for (i = 0; i < ArrayMax(fta); i++) {
	    if (arr(ft_entry_ptr, fta, i).num == ft_num) {
		break;
	    }
	}

	if (i == ArrayMax(fta)) {
	    int name_len;
	    char *name;
	    ft_entry *e;
	    ft_value_element *ele;

	    /* Not a duplicate, so add it to the fta array */
	    i = ArrayMax(fta);
	    ArrayRef(fta, i);
	    arr(ft_entry_ptr, fta, i).num = ft_num;
	    arr(ft_entry_ptr, fta, i).ranges =
		ArrayCreate(sizeof(ft_range *), 0);
	    arr(ft_entry_ptr, fta, i).sense = t->tagrec.sense;
	    
	    /* Parse the qualifiers, we'll need /gene, /transl_table, etc */
	    e = arr(ft_entry_ptr, fta, i).entry = new_ft_entry();
	    init_ft_qual_hash(e, t->newcomment);

	    /* Extract the /gene qualifier */
	    ele = search_ft_qual_hash(e, "gene");
	    name = ele ? ele->value : NULL;

	    if (name) {
		name_len = MIN(strlen(name), DB_NAMELEN);
		sprintf(arr(ft_entry_ptr, fta, i).name, "%.*s",
			name_len, name);
	    } else {
		sprintf(arr(ft_entry_ptr, fta, i).name, "CDS %d", ft_num);
	    }
	}

	/* Do we want to update the locations to be consensus positions? */
	if (!cons_locations) {
	    continue;
	}

	/* Allocate a new range and link to linked list */
	{
	    /* Extend the array, filling it out with NULLs */
	    Array a = arr(ft_entry_ptr, fta, i).ranges;
	    int max = ArrayMax(a);
	    int ind;

	    ArrayRef(a, ele_num);
	    for (ind = max; ind < ArrayMax(a); ind++) {
		arr(ft_range *, a, ind) = NULL;
	    }
	}
	range = arr(ft_range *, arr(ft_entry_ptr, fta, i).ranges, ele_num) = 
	    new_ft_range();

	/* Add locations for left and right */
	range->left = new_ft_location();
	range->left->min = tpos;
	range->right = new_ft_location();
	range->right->min = tpos + t->tagrec.length-1;
	range->complemented = t->tagrec.sense;
    }

    return fta;
}

/*
 * Stitches together exons to produce an mRNA sequence.
 * For CDSs on the complementary strand the sequence return is automatically
 * reversed and complemented.
 * 'include_pads' indicates whether the mRNA will contain padding characters
 * (which clearly need to be stripped out during translation). It's primary
 * use is to detect when a pad in the reference is within an exon and when it
 * is not, in order to determine whether insertions in other sequences fall
 * within an exon or not.
 *
 * Returns a malloced mRNA string and mapping, to be freed by the caller.
 *        or NULL for failure.
 *
 *        'len' will contain (if not NULL) the length of the returned string.
 *
 *        'mapping' will contain a translation of base position in the mRNA
 *	  (counting from zero), to positions in the contig (counting from 1).
 */
static char *create_mRNA(EdStruct *xx, int seq, Array ranges,
			 int codon_start, int include_pads,
			 int **mapping_p, int *len) {
    char *bases;
    char *mRNA = NULL;
    int *mapping = NULL;
    int mRNA_len = 0, exon_len;
    int pos = 0;
    int seq_offset;
    int ind;

    if (seq) {
	bases = DBgetSeq(DBI(xx), seq);
    } else {
	if (!(bases = (char *)xmalloc(DB_Length(xx, 0) + 1)))
	    return NULL;
	/* Ask for all consensus. This is cached (usually) */
	DBcalcConsensus(xx, 1, DB_Length(xx, 0), bases,
			NULL, BOTH_STRANDS);
    }

    if (codon_start < 1 || codon_start > 3)
	codon_start = 1;
    codon_start--;

    seq_offset = DB_RelPos(xx, seq);

    for (ind = 0; ind < ArrayMax(ranges); ind++) {
	int i, j;
	int reverse;
	ft_range *r = arr(ft_range *, ranges, ind);

	if (!r) {
	    verror(ERR_WARN, "create_mRNA", "missing exon tag");
	    continue;
	}

	exon_len = r->right->min - r->left->min + 1;
	mRNA_len += exon_len;

	if (NULL == (mRNA = (char *)xrealloc(mRNA, mRNA_len+1)))
	    goto error;
	if (NULL == (mapping = (int *)xrealloc(mapping,
					       mRNA_len * sizeof(int))))
	    goto error;

	if ((r->complemented && DB_Comp(xx, seq) != COMPLEMENTED) ||
	    (!r->complemented && DB_Comp(xx, seq) == COMPLEMENTED)) {
	    reverse = 1;
	} else {
	    reverse = 0;
	}

	/* Copy sequence, in whatever direction is appropriate */
	exon_len -= codon_start;
	for (j = i = 0; i < exon_len; i++) {
	    int seq_pos = reverse
		? r->left->min-1 + exon_len-1-i - DB_Start(xx, seq)
		: r->left->min-1+i + codon_start - DB_Start(xx, seq);
	    char base;

	    if (seq == 0) {
		if (seq_pos >= 0 && seq_pos < DB_Length(xx, 0))
		    base = bases[seq_pos];
		else
		    base = '-';
	    } else {
		base =  bases[seq_pos];
	    }

	    if (include_pads || (!include_pads && base != '*')) {
		if (reverse)
		    mRNA[pos+j] = complement_base(base);
		else
		    mRNA[pos+j] = base;
		mapping[pos+j++] = seq_pos + seq_offset;
	    }
	}
	mRNA_len -= i-j; /* difference is no. of pads */
	pos += j; /* unpadded length */
    }
    mRNA[pos] = 0;

    if (len)
	*len = mRNA_len;
    if (mapping_p)
	*mapping_p = mapping;
    else
	xfree(mapping);

    if (!seq && bases)
	xfree(bases);

    return mRNA;

 error:
    if (!seq && bases)
	xfree(bases);

    return NULL;
}

static mutation_t *new_mutation_t(void) {
    mutation_t *m = (mutation_t *)xcalloc(sizeof(mutation_t), 1);
    m->type = non_coding;
    return m;
}

/*
 * Translates an mRNA sequence.
 * Any pads in the mRNA are automatically stripped.
 *
 * Memory is allocated by this routine using malloc. The calling function
 * has the responsibility to free it.
 *
 * Returns: buffer on success plus fills out the mapping array (if not NULL)
 *	    NULL on failure
 */
static char *translate_mRNA(char *mRNA, int mRNA_len, int *mapping) {
    char *translation = (char *)xmalloc((mRNA_len+2)/3+2);
    int i, j, k;

    if (!translation)
	return NULL;

    for (j = 0, i = 0; i < mRNA_len;) {
	char codon[3];
	for (k = 0; k < 3; k++) {
	    while (i < mRNA_len && mRNA[i] == '*') {
		if (mapping)
		    mapping[i] = j;
		i++;
	    }
	    codon[k] = (i < mRNA_len) ? mRNA[i] : '-';
	    if (mapping)
		mapping[i] = j;
	    i++;
	}
	translation[j++] = codon_to_acid1(codon);
    }
    translation[j++] = 0;

    return translation;
}

/*
 * As per translate_mRNA
 */
static aa_pair_t *translate_mRNA2(char *mRNA, int mRNA_len, int *mapping) {
    aa_pair_t *translation = (aa_pair_t *)xmalloc(((mRNA_len+2)/3+2) *
						  sizeof(aa_pair_t));
    int i, j, k;

    if (!translation)
	return NULL;

    for (j = 0, i = 0; i < mRNA_len;) {
	char codon[3];
	int ambig[4];
	int total_ambigs;
	int ambig_pos;

	for (k = 0; k < 3; k++) {
	    while (i < mRNA_len && mRNA[i] == '*') {
		if (mapping)
		    mapping[i] = j;
		i++;
	    }
	    codon[k] = (i < mRNA_len) ? mRNA[i] : '-';
	    if (mapping)
		mapping[i] = j;
	    i++;
	}

	/*
	 * Count ambiguities.
	 * 0 ambiguities is easy - just print it
	 * 1 ambiguity, through looping, can be displays as [??], but we
	 * only have a data structure to store 2 amino acid characters.
	 * anything else we report as "-"
	 */
	total_ambigs = 0;
	ambig_pos = 0;
	for (k = 0; k < 3; k++) {
	    if (strchr("ACGTacgt", codon[k]) == 0) {
		total_ambigs++;
		ambig_pos = k;
	    }
	}

	if (total_ambigs == 1) { 
	    ambiguity2bases(toupper(codon[ambig_pos]),
			    &ambig[0], &ambig[1], &ambig[2], &ambig[3]);

	    if (ambig[0] + ambig[1] + ambig[2] + ambig[3] != 2)
		total_ambigs = 100; /* can't handle this code */
	}

	switch (total_ambigs) {
	case 0:
	    /* Simple case - it's a normal base */
	    translation[j].AA1 = codon_to_acid1(codon);
	    translation[j].AA2 = codon_to_acid1(codon);
	    break;

	case 1:
	    /* Ambiguity for 2 bases (eg R) */
	    {
		int a;

		for (a = 0; a < 4; a++) {
		    if (ambig[a]) {
			codon[ambig_pos] = "ACGT"[a];
			break;
		    }
		}
		translation[j].AA1 = codon_to_acid1(codon);

		for (a++; a < 4; a++) {
		    if (ambig[a]) {
			codon[ambig_pos] = "ACGT"[a];
			break;
		    }
		}
		translation[j].AA2 = codon_to_acid1(codon);
	    }
	    break;

	default:
	    /* Anything else - cannot decide to call it an X */
	    translation[j].AA1 = 'X';
	    translation[j].AA2 = 0;
	}

	j++;
    }
    translation[j].AA1 = 0;
    translation[j].AA2 = 0;
    j++;

    return translation;
}

/*
 * With a known exon, this function compares the individual sequences against
 * the reference sequence to find differences. From these differences we then
 * check whether the difference causes a silent mutation.
 */
static void check_mutations(EdStruct *xx, int ref_len, 
			    char *mRNA, int mRNA_len, int *mapping,
			    mutation_t ***muts) {
    int seq;
    int i;
    char *bases = NULL;
    char *mRNA_copy = NULL;
    char *translation = NULL; 
    aa_pair_t *translation_copy = NULL;
    int *translation_map = NULL;
    int *translation_copy_map = NULL;
    int gene_plus;
    int offset = DB_RelPos(xx, DBI(xx)->reference_seq);
    int reflen = DB_Length(xx, DBI(xx)->reference_seq);
   
    gene_plus = mapping[0] <= mapping[mRNA_len-1] ? 1 : 0;
   
    if (NULL == (mRNA_copy = (char *)xmalloc(mRNA_len+1)))
	goto error;
    
    if (NULL == (translation_map = (int *)xcalloc(mRNA_len+1, sizeof(int))))
	goto error;

    if (NULL == (translation_copy_map = (int *)xcalloc(mRNA_len+1,
						       sizeof(int))))
	goto error;

    /* Translate the original mRNA */
    translation = translate_mRNA(mRNA, mRNA_len, translation_map);

    /* Iterate around sequences */
    for (seq = 1; seq <= DBI_gelCount(xx); seq++) {
	int edited;

	/* Does sequence overlap mRNA somewhere? Crude search */
	if (mapping[0] < mapping[mRNA_len-1]) {
	    if (mapping[mRNA_len-1] < DB_RelPos(xx, seq) ||
		mapping[0] > DB_RelPos(xx, seq) + DB_Length(xx, seq))
		continue;
	} else {
	    if (mapping[0] < DB_RelPos(xx, seq)	||
		mapping[mRNA_len-1] > DB_RelPos(xx, seq) + DB_Length(xx, seq))
		continue;
	}

	/* Get sequence and produce a new mRNA */
	strcpy(mRNA_copy, mRNA);
	edited = 0;
	bases = DBgetSeq(DBI(xx), seq);
	for (i = 0; i < mRNA_len; i++) {
	    int pos = mapping[i];
	    if (pos >= DB_RelPos(xx, seq) &&
		pos <= DB_RelPos(xx, seq) + DB_Length(xx, seq)-1) {
		if (mapping[i]-offset < 0 ||
		    mapping[i]-offset >= reflen ||
		    !muts[seq][mapping[i]-offset])
		    continue;
		/* mRNA_copy[i] = bases[pos - DB_RelPos(xx, seq)]; */
		mRNA_copy[i] = muts[seq][mapping[i]-offset]->nucleotide_to;
		if (!gene_plus) {
		    mRNA_copy[i] = complement_base(mRNA_copy[i]);
		}
		edited = 1;
	    }
	}
	if (!edited)
	    continue;

	/* Translate the new sequence */
	translation_copy = translate_mRNA2(mRNA_copy, mRNA_len,
					   translation_copy_map);

	/* Loop through mapping array identifying location of mutations */
	for (i = 0; i < mRNA_len; i++) {
	    int pos = mapping[i] - offset;
	    int cpos;

	    if (pos < 0 || pos >= ref_len || !muts[seq][pos])
		continue;

	    muts[seq][pos]->AA_from = translation[translation_map[i]];
	    cpos = translation_copy_map[i];
	    muts[seq][pos]->AA_to = translation_copy[cpos];
	    muts[seq][pos]->type =
		(muts[seq][pos]->AA_to.AA1 == muts[seq][pos]->AA_to.AA2 &&
		 muts[seq][pos]->AA_from == muts[seq][pos]->AA_to.AA1)
		? silent : expressed;
	}

	xfree(translation_copy);
    }

 error:
    if (mRNA_copy)
	xfree(mRNA_copy);
    if (translation)
	xfree(translation);
    if (translation_map)
	xfree(translation_map);
    if (translation_copy_map)
	xfree(translation_copy_map);

    return;
}

static void dump_single_mutation(EdStruct *xx, mutation_t ***muts,
				 int refcmp, int seq, int pos, int *transpos) {
    char nuc_from, nuc_to;

    if (refcmp) {
	nuc_from = complement_base(muts[seq][pos]->nucleotide_from);
	nuc_to   = complement_base(muts[seq][pos]->nucleotide_to);
    } else {
	nuc_from = muts[seq][pos]->nucleotide_from;
	nuc_to   = muts[seq][pos]->nucleotide_to;
    }

    if (muts[seq][pos]->type == no_mutation) {
	/* insertion */
	vmessage("%s (No mutations found)\n",
		 io_rname(DBI_io(xx), DB_Number(xx, seq)));
	return;
    }

    if (nuc_from == '*' && nuc_to != '*') {
	/* insertion */
	vmessage("%s %5dins%c",
		 io_rname(DBI_io(xx), DB_Number(xx, seq)),
		 transpos[pos], nuc_to);
    } else if (nuc_from != '*' && nuc_to == '*') {
	/* deletion */
	vmessage("%s %5ddel%c",
		 io_rname(DBI_io(xx), DB_Number(xx, seq)),
		 transpos[pos], nuc_from);
    } else {
	/* substitution */
	vmessage("%s %5d%c>%c",
		 io_rname(DBI_io(xx), DB_Number(xx, seq)),
		 transpos[pos], nuc_from, nuc_to);
    }

    switch(muts[seq][pos]->type) {
    case no_mutation:
	vmessage(" (No mutations found)");
	break;

    case non_coding:
	vmessage(" (noncoding)");
	break;

    case silent:
	vmessage(" (silent %c)",
		 muts[seq][pos]->AA_from);
	break;

    case expressed:
	if (muts[seq][pos]->AA_to.AA1 == muts[seq][pos]->AA_to.AA2) {
	    vmessage(" (expressed %c>%c)",
		     muts[seq][pos]->AA_from,
		     muts[seq][pos]->AA_to.AA1);
	} else {
	    vmessage(" (expressed %c>[%c%c])",
		     muts[seq][pos]->AA_from,
		     muts[seq][pos]->AA_to.AA1,
		     muts[seq][pos]->AA_to.AA2);
	}

	break;
    }

    if (muts[seq][pos]->strands == 3) {
	vmessage(" (double stranded)");
    } else {
	vmessage(" (strand %c only)",
		 muts[seq][pos]->strands == 1 ? '+' : '-');
    }

    if (muts[seq][pos]->conflict) {
	vmessage(" (strand conflict)");
    }

    vmessage("\n");
}

static void html_mutation_summary(dstring_t *html,
				  EdStruct *xx, mutation_t ***muts,
				  mutation_cov_t *mcov,
				  int refcmp, int seq, int pos,
				  int *transpos,
				  int bg_toggle,
				  int print_ranges) {
    char nuc_from, nuc_to;
    char fwd_str, rev_str;
    int covered;
    int fseq;
    int rseq;

    if (mcov) {
	fseq = mcov[seq].fwd ? seq : mcov[seq].sibling;
	rseq = mcov[seq].fwd ? mcov[seq].sibling : seq;
    } else {
	fseq = seq;
	rseq = 0;
	print_ranges = 0;
    }

    dstring_appendf(html,
		    "    <tr bgcolor=\"%s\">\n"
		    "      <td><a href=\"#Sample_%s_%d\">%s</a></td>\n"
		    "      <td></td>\n",
		    bg_toggle ? "#ffffff" : "#e0e0e0",
		    io_rname(DBI_io(xx), DB_Number(xx, seq)),
		    transpos[pos],
		    io_rname(DBI_io(xx), DB_Number(xx, seq)));

    /* Check for dummy 'no_mutation' mutations, and display if appropriate */
    if (muts[seq][pos]->type == no_mutation) {
	dstring_appendf(html,
			"      <td>(None)</td>\n"
			"      <td></td>\n"
			"      <td>-</td>\n"
			"      <td></td>\n"
			"      <td align=center>-</td>\n"
			"      <td></td>\n"
			"      <td align=center>-</td>\n"
			"      <td></td>\n"
			"      <td>%d - %d</td>\n"
			"      <td></td>\n"
			"      <td>%d - %d</td>\n"
			"    </tr>\n",
			mcov ? mcov[fseq].ref_start : 0,
			mcov ? mcov[fseq].ref_end : 0,
			mcov ? mcov[rseq].ref_start : 0,
			mcov ? mcov[rseq].ref_end : 0);
	return;
    }

    if (refcmp) {
	nuc_from = complement_base(muts[seq][pos]->nucleotide_from);
	nuc_to   = complement_base(muts[seq][pos]->nucleotide_to);
    } else {
	nuc_from = muts[seq][pos]->nucleotide_from;
	nuc_to   = muts[seq][pos]->nucleotide_to;
    }

    if (nuc_from == '*' && nuc_to != '*') {
	/* insertion */
	dstring_appendf(html, "      <td>%dins%c</td>\n",
			transpos[pos], nuc_to);
    } else if (nuc_from != '*' && nuc_to == '*') {
	/* deletion */
	dstring_appendf(html, "      <td>%ddel%c</td>\n",
			transpos[pos], nuc_from);
    } else {
	/* substitution */
	dstring_appendf(html, "      <td>%d%c&gt;%c</td>\n",
			transpos[pos], nuc_from, nuc_to);
    }
    dstring_append(html, "      <td></td>\n");

    switch(muts[seq][pos]->type) {
    case no_mutation:
	dstring_append(html, "      <td>-</td>\n");
	break;

    case non_coding:
	dstring_append(html, "      <td>noncoding</td>\n");
	break;

    case silent:
	dstring_appendf(html, "      <td>silent %c</td>\n",
			muts[seq][pos]->AA_from);
		 
	break;

    case expressed:
	if (muts[seq][pos]->AA_to.AA1 == muts[seq][pos]->AA_to.AA2) {
	    dstring_appendf(html, "      <td>expressed %c&gt;%c</td>\n",
			    muts[seq][pos]->AA_from,
			    muts[seq][pos]->AA_to.AA1);
	} else {
	    dstring_appendf(html, "      <td>expressed %c>[%c%c]</td>\n",
			    muts[seq][pos]->AA_from,
			    muts[seq][pos]->AA_to.AA1,
			    muts[seq][pos]->AA_to.AA2);
	}

	break;
    }
    dstring_append(html, "      <td></td>\n");

			
    covered = mcov
	? (pos >= mcov[fseq].ref_start && pos <= mcov[fseq].ref_end)
	: 0;
    fwd_str = (muts[seq][pos]->strands & 1) ? "yY"[covered] : "-N"[covered];
	
    covered = mcov
	? (pos >= mcov[rseq].ref_start && pos <= mcov[rseq].ref_end)
	: 0;
    rev_str = (muts[seq][pos]->strands & 2) ? "yY"[covered] : "-N"[covered];

    if (muts[seq][pos]->conflict) {
	fwd_str = rev_str = 'X';
    }
    dstring_appendf(html,
		    "      <td align=center>%c</td>\n"
		    "      <td></td>\n"
		    "      <td align=center>%c</td>\n",
		    fwd_str, rev_str);
    
    if (print_ranges) {
	dstring_appendf(html,
			"      <td></td>\n"
			"      <td>%d - %d</td>\n"
			"      <td></td>\n"
			"      <td>%d - %d</td>\n",
			mcov[fseq].ref_start,
			mcov[fseq].ref_end,
			mcov[rseq].ref_start,
			mcov[rseq].ref_end);
    } else {
	dstring_append(html,
		       "      <td></td>\n"
		       "      <td></td>\n"
		       "      <td></td>\n"
		       "      <td></td>\n");
    }

    dstring_append(html, "    </tr>\n");
}

static void html_mutation_detailed(dstring_t *html,
				   EdStruct *xx, mutation_t ***muts,
				   int refcmp, int seq, int pos,
				   int *transpos, int page_break,
				   char *dir) {
    int refpos = DB_RelPos(xx, DBI(xx)->reference_seq);
    int spos;

    /*
     * Pos is a position in the ref sequence.
     * Convert it to an original position within the trace file.
     */
    spos = pos + refpos; /* consensus pos */
    save_trace_images(html, xx, seq, spos, muts[seq][pos], transpos[pos],
		      page_break, dir);
}

/*
 * qsort callback function; used by get_sorted_seqids().
 */
static EdStruct *seqids_sort_xx = NULL;
static int seqids_sort_func(const void *v1, const void *v2) {
    const int *id1 = (const int *)v1;
    const int *id2 = (const int *)v2;
    char *name1 = io_rname(DBI_io(seqids_sort_xx),
			   DB_Number(seqids_sort_xx, *id1));
    char *name2 = io_rname(DBI_io(seqids_sort_xx),
			   DB_Number(seqids_sort_xx, *id2));

    return strcmp(name1, name2);
}

/*
 * Returns an alphabetical sorted version of the sequences in this database.
 * The memory returned should be freed using xfree().
 *
 * Returns a malloced buffer of seqid integers on success
 *         NULL on failure.
 */
static int *get_sorted_seqids(EdStruct *xx) {
    int *ids;
    int i;

    ids = (int *)xmalloc(DBI_gelCount(xx) * sizeof(int));
    if (!ids)
	return NULL;

    /* Seed with 1 to N */
    for (i = 0; i < DBI_gelCount(xx); i++) {
	ids[i] = i+1;
    }

    /* Sort - alas we can't pass over client data, so we use a global. */
    seqids_sort_xx = xx;
    qsort(ids, DBI_gelCount(xx), sizeof(int), seqids_sort_func);

    return ids;
}

static void dump_mutations(dstring_t *html, EdStruct *xx,
			   mutation_t ***muts, mutation_cov_t *mcov,
			   int sort_by_position, char *dir, int detail) {
    int i, j, reflen, refoff, refcmp;
    char *bases;
    int *transpos;
    int toggle;
    int *sorted_seqs;

    reflen = DB_Length(xx, DBI(xx)->reference_seq);
    refoff = DBI(xx)->reference_offset;
    refcmp = (DB_Comp(xx, DBI(xx)->reference_seq) == COMPLEMENTED);

    /*
     * Allocate a translation buffer for positions in the reference sequence
     * to compensate for pads, complementing and cyclic sequences.
     */
    bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
    transpos = (int *)xcalloc(reflen+1, sizeof(int));
    if (!transpos)
	return;

    if (DB_Comp(xx, DBI(xx)->reference_seq) == UNCOMPLEMENTED) {
	for (i = 0, j = 1; i < reflen; i++) {
	    if (DBI(xx)->reference_len)
		transpos[i] = ((j-1+refoff-1) % DBI(xx)->reference_len) + 1;
	    else
		transpos[i] = ((j-1+refoff-1)) + 1;
	    if (bases[i] != '*')
		j++;
	}
    } else {
	for (i = reflen-1, j = 1; i >= 0; i--) {
	    if (DBI(xx)->reference_len)
		transpos[i] = ((j-1+refoff-1) % DBI(xx)->reference_len) + 1;
	    else
		transpos[i] = ((j-1+refoff-1)) + 1;
	    if (bases[i] != '*')
		j++;
	}
    }

    /* Textual report */
    vfuncheader("Report Mutations");
    if (sort_by_position) {
	int found_one;
	for (j = 0; j < reflen; j++) {
	    found_one = 0;
	    for (i = 1; i <= DBI_gelCount(xx); i++) {
		if (muts[i][j]) {
		    dump_single_mutation(xx, muts, refcmp, i, j, transpos);
		    found_one = 1;
		}
	    }
	    if (found_one)
		vmessage("\n");
	}
    } else {
	int found_one;
	for (i = 1; i <= DBI_gelCount(xx); i++) {
	    found_one = 0;
	    for (j = 0; j < reflen; j++) {
		if (muts[i][j]) {
		    dump_single_mutation(xx, muts, refcmp, i, j, transpos);
		    found_one = 1;
		}
	    }
	    if (found_one)
		vmessage("\n");
	}
    }

    /* HTML report */
    if  (html && detail >= 1) {
	char date[1024];
	time_t t;

	strftime(date, sizeof(date)-1, "%c %Z", localtime((time(&t),&t)));
	dstring_appendf(html, html_template0,
			date, io_name(DBI_io(xx)),
			io_rname(DBI_io(xx), DBI_contigNum(xx)));

	/* Summary HTML */
	dstring_appendf(html, html_template1, "Position");
	toggle = 0;
	for (j = 0; j < reflen; j++) {
	    int to_toggle = 0;
	    for (i = 1; i <= DBI_gelCount(xx); i++) {
		if (muts[i][j] && muts[i][j]->type != no_mutation) {
		    html_mutation_summary(html, xx, muts, mcov, refcmp, i, j,
					  transpos, toggle, 1);
		    to_toggle = 1;
		}
	    }
	    if (to_toggle)
		toggle ^= 1;
	}
	dstring_append(html, html_template2);
	
	dstring_appendf(html, html_template1, "Name");
	toggle = 0;
	sorted_seqs = get_sorted_seqids(xx);
	for (i = 0; i < DBI_gelCount(xx); i++) {
	    int seq = sorted_seqs[i];
	    int to_toggle = 0;
	    int ranges = 1;
	    for (j = 0; j < reflen; j++) {
		if (muts[seq][j]) {
		    html_mutation_summary(html, xx, muts, mcov, refcmp, seq, j,
					  transpos, toggle, ranges);
		    to_toggle = 1;
		    ranges = 0;
		}
	    }
	    if (to_toggle)
		toggle ^= 1;
	}
	xfree(get_sorted_seqids(xx));
	dstring_append(html, html_template2);

	if (detail == 2) {
	    /* Detailed HTML */
	    dstring_append(html,
			   "<p><hr><p><h1 style=\"page-break-before: always\">"
			   "Detailed Mutation Report</h1>\n");
	    if (sort_by_position) {
		int first = 1;
		for (j = 0; j < reflen; j++) {
		    for (i = 1; i <= DBI_gelCount(xx); i++) {
			if (muts[i][j] && muts[i][j]->type != no_mutation) {
			    html_mutation_detailed(html, xx, muts, refcmp,
						   i, j, transpos,
						   first ? 0 : 1,
						   dir);
			    first = 0;
			}
		    }
		}
	    } else {
		int first = 1;
		for (i = 1; i <= DBI_gelCount(xx); i++) {
		    for (j = 0; j < reflen; j++) {
			if (muts[i][j] && muts[i][j]->type != no_mutation) {
			    html_mutation_detailed(html, xx, muts, refcmp,
						   i, j, transpos,
						   first ? 0 : 1,
						   dir);
			    first = 0;
			}
		    }
		}
	    }
	}

	dstring_append(html, "</body></html>\n");
    }

    xfree(transpos);
}

/*
 * Sort's and combines mutations by fwd/rev pairs.
 */
static void sort_mutations(EdStruct *xx, mutation_t ***muts,
			   mutation_cov_t *mcov) {
    int i, j, pos, reflen;
    int *processed;

    reflen = DB_Length(xx, DBI(xx)->reference_seq);

    processed = (int *)xcalloc(DBI_gelCount(xx)+1, sizeof(int));
    if (!processed)
	return;

    for (i = 1; i <= DBI_gelCount(xx); i++) {
	int i_has_mut = 0, j_has_mut = 0;
	int fwd_seq = i;
	int temp = DBI_DB(xx)[i].template;

	if (processed[i])
	    continue;

	/* Is there a mutation on this sequence anywhere? */
	for (pos = 0; pos < reflen; pos++) {
	    if (muts[i][pos]) {
		i_has_mut = 1;
		break;
	    }
	}

	/* Find 'sibling' sequences - from other direction */
	for (j = i+1; j <= DBI_gelCount(xx); j++) {
	    int fwd, rev, primeri, primerj;
	    GReadings r;

	    if (DBI_DB(xx)[j].template != temp)
		continue;

	    /*
	     * Set fwd/rev to be i/j or j/i depending on which is which
	     * It's possible that we have 2 forwards or 2 reverses, so this
	     * is just a best guess. Either way both i and j are contained
	     * in fwd and rev regardless of whether it's a real pair.
	     */
	    gel_read(DBI_io(xx), DB_Number(xx, i), r);
	    primeri = PRIMER_TYPE(r);
	    gel_read(DBI_io(xx), DB_Number(xx, j), r);
	    primerj = PRIMER_TYPE(r);

	    if (primerj == 2 || primerj == 4) {
		fwd_seq = fwd = i;
		rev = j;
	    } else {
		fwd_seq = fwd = j;
		rev = i;
	    }

	    processed[j] = 1;

	    if (mcov) {
		mcov[i].sibling = j;
		mcov[j].sibling = i;
		mcov[i].fwd = (i == fwd);
		mcov[j].fwd = (j == fwd);
	    }

	    /*
	     * Found a pair - combine results.
	     * In combining we move or merge the 'rev' sequence into the 'fwd'
	     * sequence and delete the 'rev' records.
	     */
	    for (pos = 0; pos < reflen; pos++) {
		/* Neither have a mutation at this point - skip to next base */
		if (!muts[fwd][pos] && !muts[rev][pos])
		    continue;

		if (muts[j][pos])
		    j_has_mut = 1;

		/* Both have mutations - merge */
		if (muts[fwd][pos] && muts[rev][pos]) {
		    /* Are they the same? */
		    if ((muts[fwd][pos]->nucleotide_from ==
			 muts[rev][pos]->nucleotide_from) &&
			(muts[fwd][pos]->nucleotide_to ==
			 muts[rev][pos]->nucleotide_to)) {
			muts[fwd][pos]->conflict = 0;
		    } else {
			/*
			 * Perhaps they're compatible MUTA + HETE. This often
			 * happens due to missed HETE tags.
			 */
			int A[2], C[2], G[2], T[2];
			int het[2];
			ambiguity2bases(muts[fwd][pos]->nucleotide_to,
					&A[0], &C[0], &G[0], &T[0]);
			ambiguity2bases(muts[rev][pos]->nucleotide_to,
					&A[1], &C[1], &G[1], &T[1]);

			/* One 'To' needs to be a MUTA A,C,G or T */
			het[0] = (A[0] + C[0] + G[0] + T[0] != 1);
			het[1] = (A[1] + C[1] + G[1] + T[1] != 1);
			if ( (!het[0] || !het[1]) &&
			     ((A[0] & A[1]) |
			      (C[0] & C[1]) |
			      (G[0] & G[1]) |
			      (T[0] & T[1]))) {
			    if (het[1]) /* Reverse is het, copy to fwd */
				muts[fwd][pos]->nucleotide_to =
				    muts[rev][pos]->nucleotide_to;
			    muts[fwd][pos]->conflict = 0;
			} else {
			    muts[fwd][pos]->conflict = 1;
			}
		    }

		    /* Set strand coverage */
		    muts[fwd][pos]->strands = 0;
		    muts[fwd][pos]->strands |=
			DB_Comp(xx, fwd) == UNCOMPLEMENTED ? 1 : 2;
		    muts[fwd][pos]->strands |=
			DB_Comp(xx, rev) == UNCOMPLEMENTED ? 1 : 2;

		    memcpy(muts[fwd][pos]->tag_type_bot,
			   muts[rev][pos]->tag_type_top,
			   4);
		    muts[fwd][pos]->tag_text_bot =
			muts[rev][pos]->tag_text_top;

		    xfree(muts[rev][pos]);
		    muts[rev][pos] = NULL;

		} else if (muts[rev][pos]) {
		    /* 2nd sequence only has a mutation */
		    muts[fwd][pos] = muts[rev][pos];
		    memcpy(muts[fwd][pos]->tag_type_bot,
			   muts[fwd][pos]->tag_type_top,
			   4);
		    muts[fwd][pos]->tag_type_top[0] = 0;
		    muts[fwd][pos]->tag_text_bot =
			muts[fwd][pos]->tag_text_top;
		    muts[fwd][pos]->tag_text_top = NULL;
		    muts[rev][pos] = NULL;
		    muts[fwd][pos]->strands |=
			DB_Comp(xx, rev) == UNCOMPLEMENTED ? 1 : 2;
		} else {
		    /* 1st sequence only */
		    muts[fwd][pos]->strands |=
			DB_Comp(xx, fwd) == UNCOMPLEMENTED ? 1 : 2;
		}

		muts[fwd][pos]->seq_bot = rev;
		muts[fwd][pos]->seq_top = fwd;
	    }
	}

	/*
	 * If we did not find a mutation on either the forward or reverse
	 * sequence (or if we only have one sequence and it doesn't contain
	 * a mutation) then create a dummy mutation of type 'none' in order
	 * for subsequent algorithms to be able to iterate through all
	 * sequences by iterating through all mutations.
	 */
	if (i_has_mut == 0 && j_has_mut == 0) {
	    muts[fwd_seq][0] = new_mutation_t();
	    muts[fwd_seq][0]->type = no_mutation;
	}
    }

    xfree(processed);
}

static void free_mutations(EdStruct *xx, mutation_t ***muts) {
    int i, j, reflen;

    reflen = DB_Length(xx, DBI(xx)->reference_seq);

    for (i = 1; i <= DBI_gelCount(xx); i++) {
	for (j = 0; j < reflen; j++) {
	    if (muts[i][j]) {
		if (muts[i][j]->tag_text_top)
		    xfree(muts[i][j]->tag_text_top);
		if (muts[i][j]->tag_text_bot)
		    xfree(muts[i][j]->tag_text_bot);
		xfree(muts[i][j]);
	    }
	}
	xfree(muts[i]);
    }
    xfree(muts);
}

/*
 * Having found an exon, this function stores the snippet of translation
 * that is visible within the contig editor window into a status_line
 * structure.
 */
static void store_translation(EdStruct *xx, int pos, int width,
			      char *mRNA, int mRNA_len, int *mapping,
			      int dir, char *name) {
    int i, blank = 1;
    char line[MAX_DISPLAY_WIDTH+1];
    int p, l;

    /* Compute a translation */
    line[MAX_DISPLAY_WIDTH]=0;
    memset(line, ' ', MAX_DISPLAY_WIDTH);

    for (i = 1; i < mRNA_len; i+=3) {
	if (mapping[i] >= pos && mapping[i] < pos+width) {
	    line[mapping[i]-pos] = codon_to_acid1(&mRNA[i-1]);
	    blank = 0;
	}
    }

    /* Did we translate anything within pos to pos+width? */
    if (blank)
	return;

    /* Allocate a status_line structure */
    l = xx->status_depth++;
    xx->status_lines =
	(EdStatus *)xrealloc(xx->status_lines,
			     xx->status_depth * sizeof(EdStatus));
    if (!xx->status_lines) {
	xx->status_depth = 0;
	return;
    }
   
    /* Set the line contents */
    memcpy(xx->status_lines[l].line, line, width);

    /* Set the line colour */
    for (p = 0; p < width; p++) {
	xx->status_lines[l].colours[p].sh = sh_default;
    } 

    /* Set the line name */
    sprintf(xx->status_lines[l].name, " %*c %-*s",
	    DB_GELNOLEN, dir == 0 ? '+' : '-',
	    DB_NAMELEN, name);
}

/*
 * Looks for exons overlapping the editor between pos and pos+width.
 * Exons are defined by annotations on the reference sequence.
 *
 * TODO: If no reference sequence exists then we fabricate 6 overlapping exons
 * spanning the entire contig, to simulate the old style display of
 * translation display.
 */
void find_exons(EdStruct *xx, int pos, int width, int generate) {
    Array fta;
    int i;
    int seq;

    if (generate) {
	/* Create artificial CDS structures */
	fta = generate_exons(xx, pos, width);
	seq = 0;
    } else {
	seq = DBI(xx)->reference_seq;
	if (!seq) {
	    return;
	}

	/* Regenerate the feature table structures from tag comments */
	fta = get_ft_entries(xx, seq, 1);
    }

    if (!fta) {
	verror(ERR_WARN, "find_exons", "Couldn't find any CDS lines");
	return;
    }

    /*
     * Scan through consensus ranges convered by CDS exons looking for items
     * which overlap pos to pos+width. We assume here that the ranges
     * are sorted in real exon order (ie reverse order if complemented).
     */
    for (i = 0; i < ArrayMax(fta); i++) {
	ft_range *r;
	char *mRNA;
	int *mapping;
	int mRNA_len;
	int overlap = 0;
	Array ranges;
	int ind;
	ft_entry *e;
	int codon_start = 1;
	int transl_table;

	ranges = arr(ft_entry_ptr, fta, i).ranges;
	/* Does CDS contain an exon overlapping pos to pos+width? */
	for (ind = 0; ind < ArrayMax(ranges); ind++) {
	    int left, right;
	    if (NULL == (r = arr(ft_range *, ranges, ind)))
		continue;

	    left = r->left->min - DB_Start(xx, seq) + (DB_RelPos(xx, seq)-1);
	    right = r->right->min - DB_Start(xx, seq) + (DB_RelPos(xx, seq)-1);

	    if (left <= pos+width && right >= pos)
		overlap = 1;
	}

	if (!overlap)
	    continue;

	/* Look for a /codon_start qualifier */
	if ((e = arr(ft_entry_ptr, fta, i).entry)) {
	    ft_value_element *ele;
	    ele = search_ft_qual_hash(e, "codon_start");
	    if (ele && ele->value)
		codon_start = atoi(ele->value);
	}

	/* Look for a genetic code qualifier /transl_table */
	transl_table = 1; /* Standard */
	if ((e = arr(ft_entry_ptr, fta, i).entry)) {
	    ft_value_element *ele;
	    ele = search_ft_qual_hash(e, "transl_table");
	    if (ele && ele->value)
		transl_table = atoi(ele->value);
	}
	if (load_genetic_code_number(transl_table) == -1) {
	    verror(ERR_WARN, "load_genetic_code_number",
		   "Failed to load code %d; using standard code",
		   transl_table);
	    load_genetic_code_number(1);
	}

	/* Look for /codon qualifiers and edit table accordingly */
	if ((e = arr(ft_entry_ptr, fta, i).entry)) {
	    ft_value_element *ele;
	    for (ele = search_ft_qual_hash(e, "codon");
		 ele;
		 ele = ele->next) {
		if (!ele->value)
		    continue;
		
		if (edit_genetic_code(ele->value) == -1) {
		    verror(ERR_WARN, "edit_genetic_code",
			   "Invalid /codon '%s'\n", ele->value);
		}
	    }
	}

	/* It does, so produce an mRNA */
	mRNA = create_mRNA(xx, seq, ranges, codon_start, 0,
			   &mapping, &mRNA_len); 
	if (!mRNA)
	    continue;

	/* Add the result */
	store_translation(xx, pos, width, mRNA, mRNA_len, mapping,
			  arr(ft_entry_ptr, fta, i).sense,
			  arr(ft_entry_ptr, fta, i).name);

	xfree(mRNA);
	xfree(mapping);
    }

    /* Tidy up */
    for (i = 0; i < ArrayMax(fta); i++) {
	Array r;

	r = arr(ft_entry_ptr, fta, i).ranges;
	if (r) {
	    int j;
	    for (j = 0; j < ArrayMax(r); j++)
		if (arr(ft_range *, r, j))
		    del_ft_range(arr(ft_range *, r, j));
	    ArrayDestroy(r);
	}

	if (arr(ft_entry_ptr, fta, i).entry)
	    del_ft_entry(arr(ft_entry_ptr, fta, i).entry);
    }
    ArrayDestroy(fta);
}

/*
 * Allocate an array for all sequences overlapping the reference sequence.
 * Each array is indexed by ref sequence position and contains NULL if
 * there is no mutation or a pointer to a mutation_t structure otherwise.
 *
 * Mutations here are identified by their difference to the reference sequence.
 * MUTA and HETE tags are ignored.
 */
static mutation_t ***allocate_mutations_diff(EdStruct *xx) {
    int seq, refseq, reflen;
    mutation_t ***muts;
    char *bases, *refbases;

    refseq = DBI(xx)->reference_seq;
    reflen = DB_Length(xx, refseq);

    /* Create an empty mutation array to start with */
    muts = (mutation_t ***)xcalloc(DBI_gelCount(xx)+1, sizeof(muts));
    for (seq = 1; seq <= DBI_gelCount(xx); seq++) {
	int j;
	muts[seq] = (mutation_t **)xcalloc(reflen, sizeof(muts[seq]));
	for (j = 0; j < reflen; j++)
	    muts[seq][j] = NULL;
    }

    /*
     * Now loop through all sequences spotting where they differ to the
     * reference sequence. We create mutations here that default to
     * intronic. When processing the CDS records later we change the
     * appropriate records.
     */
    refbases = DBgetSeq(DBI(xx), refseq);
    for (seq = 1; seq <= DBI_gelCount(xx); seq++) {
	int start_ref, start_seq, len; /* positions and length */

	bases = DBgetSeq(DBI(xx), seq);
	if (DB_RelPos(xx, seq) >= DB_RelPos(xx, refseq)) {
	    start_seq = 0;
	    start_ref = DB_RelPos(xx, seq) - DB_RelPos(xx, refseq);
	} else {
	    start_seq = DB_RelPos(xx, refseq) - DB_RelPos(xx, seq);
	    start_ref = 0;
	}

	if (DB_RelPos(xx, seq) + DB_Length(xx, seq) >
	    DB_RelPos(xx, refseq) + DB_Length(xx, refseq)) {
	    len = DB_Length(xx, refseq) - start_ref;
	} else {
	    len = DB_Length(xx, seq) - start_seq;
	}

	for (; len; start_ref++, start_seq++, len--) {
	    if (toupper(refbases[start_ref]) !=
		toupper(bases[start_seq])) {
		muts[seq][start_ref] = new_mutation_t();
		muts[seq][start_ref]->nucleotide_from =
		    toupper(refbases[start_ref]);
		muts[seq][start_ref]->nucleotide_to =
		    toupper(bases[start_seq]);
	    }
	}
    }

    return muts;
}

/*
 * Allocate an array for all sequences overlapping the reference sequence.
 * Each array is indexed by ref sequence position and contains NULL if
 * there is no mutation or a pointer to a mutation_t structure otherwise.
 *
 * Mutations here are identified by the presence of MUTA or HETE tags.
 * Differences to the reference sequence are ignored.
 */
static mutation_t ***allocate_mutations_tagged(EdStruct *xx,
					       mutation_cov_t **cov) {
    int seq, refseq, reflen, refcomp;
    mutation_t ***muts;
    char *bases, *refbases;

    refseq = DBI(xx)->reference_seq;
    reflen = DB_Length(xx, refseq);
    refcomp = (DB_Comp(xx, DBI(xx)->reference_seq) == COMPLEMENTED);

    /* Create an empty mutation array to start with */
    muts = (mutation_t ***)xcalloc(DBI_gelCount(xx)+1, sizeof(muts));
    *cov = (mutation_cov_t *)xcalloc(DBI_gelCount(xx)+1,
				     sizeof(mutation_cov_t));
    for (seq = 1; seq <= DBI_gelCount(xx); seq++) {
	int j;
	muts[seq] = (mutation_t **)xcalloc(reflen, sizeof(muts[seq]));
	for (j = 0; j < reflen; j++)
	    muts[seq][j] = NULL;
    }

    /*
     * Now loop through all sequences looking for HETE and MUTA tags.
     * reference sequence. We create mutations here that default to
     * intronic. When processing the CDS records later we change the
     * appropriate records.
     */
    refbases = DBgetSeq(DBI(xx), refseq);
    for (seq = 1; seq <= DBI_gelCount(xx); seq++) {
	tagStruct *t;

	(void)DBgetSeq(DBI(xx), seq);
	t = DBgetTags(DBI(xx), seq);

	while (t && (t = t->next)) {
	    int hete = 0;
	    int muta = 0;
	    int mcov = 0;
	    int refpos = 0, seqpos = 0, abspos = 0;
	    char base_to;

	    /* Look only for HETE and MUTA tags - assume tag is 1 base long */
	    if (t->tagrec.type.c[0] == 'H' &&
		t->tagrec.type.c[1] == 'E' &&
		t->tagrec.type.c[2] == 'T' &&
		t->tagrec.type.c[3] == 'E')
		hete = 1;
	    else if (t->tagrec.type.c[0] == 'M' &&
		     t->tagrec.type.c[1] == 'U' &&
		     t->tagrec.type.c[2] == 'T' &&
		     t->tagrec.type.c[3] == 'A')
		muta = 1;
	    else if (t->tagrec.type.c[0] == 'M' &&
		     t->tagrec.type.c[1] == 'C' &&
		     t->tagrec.type.c[2] == 'O' &&
		     t->tagrec.type.c[3] == 'V')
		mcov = 1;

	    if (!hete && !muta && !mcov)
		continue;

	    /* Find the position on the reference sequence. */
	    seqpos = normalisePos2(xx,
				   seq,
				   t->tagrec.position,
				   t->tagrec.length) - DB_Start(xx, seq);

	    abspos = seqpos + DB_RelPos(xx, seq);
	    refpos = abspos - DB_RelPos(xx, refseq);

	    if (mcov) {
		(*cov)[seq].ref_start = refpos;
		(*cov)[seq].ref_end   = refpos + t->tagrec.length-1;
		continue;
	    }

	    if (seqpos <= 0 || seqpos > DB_Length(xx, seq)) {
		continue;
	    }

	    if (abspos < DB_RelPos(xx, refseq) ||
		abspos > DB_RelPos(xx, refseq) + DB_Length(xx, refseq)) {
		continue;
	    }

	    /* Create a mutation structure for the tag */
	    bases = DBgetSeq(DBI(xx), seq);
	    refpos--; /* index from 0 */
	    seqpos--;

	    /* Heterozygous - read the tag comment */
	    force_comment(DBI_io(xx), t);
	    if (hete && t->newcomment && strlen(t->newcomment) >= 4) {
		base_to = ambiguity_code(t->newcomment[0], t->newcomment[3]);
		if (refcomp)
		    base_to = complement_base(base_to);
	    } else {
		base_to = toupper(bases[seqpos]);
	    }

	    muts[seq][refpos] = new_mutation_t();
	    muts[seq][refpos]->nucleotide_from = toupper(refbases[refpos]);
	    muts[seq][refpos]->nucleotide_to = base_to;
	    memcpy(muts[seq][refpos]->tag_type_top, t->tagrec.type.c, 4);
	    if (t->newcomment) {
		muts[seq][refpos]->tag_text_top = strdup(t->newcomment);
	    }
	    muts[seq][refpos]->seq_top = seq;
	}
    }

    return muts;
}

/*
 * Reports mutations occuring inside and outside of exons. When in an exon
 * we also report if it is a silent mutation.
 *
 * diffs_tagged is true when we only want to report mutations where HETE and
 * MUTA tags exist, otherwise base call differences are used instead.
 *
 * sort_by_position, when true, outputs mutations column-by-column instead of
 * sequence by sequence.
 *
 * Returns: HTML report of the mutations as a dstring pointer. This should
 *		be freed by the caller.
 *	    NULL if failed, and sets err_msg pointer to an error string.
 */
dstring_t *report_mutations(EdStruct *xx, int diffs_tagged,
			    int sort_by_position, 
			    char *dir, int detail, char **err_msg) {
    Array fta;
    int i;
    int seq;
    mutation_t ***muts;
    dstring_t *html;
    mutation_cov_t *mut_coverage = NULL;

    seq = DBI(xx)->reference_seq;
    if (!seq) {
	bell();
	if (err_msg)
	    *err_msg = "No reference sequence has been specified.";
	return NULL;
    }

    if (diffs_tagged)
	muts = allocate_mutations_tagged(xx, &mut_coverage);
    else
	muts = allocate_mutations_diff(xx);

    /* Regenerate the feature table structures from tag comments */
    fta = get_ft_entries(xx, seq, 1);

    /*
     * Scan through consensus ranges convered by CDS exons.
     */
    for (i = 0; i < ArrayMax(fta); i++) {
	char *mRNA;
	int *mapping;
	int mRNA_len;
	Array ranges;
	ft_entry *e;
	int codon_start = 1;
	int transl_table;

	ranges = arr(ft_entry_ptr, fta, i).ranges;

	/* Look for a /codon_start qualifier */
	if ((e = arr(ft_entry_ptr, fta, i).entry)) {
	    ft_value_element *ele;
	    ele = search_ft_qual_hash(e, "codon_start");
	    if (ele && ele->value)
		codon_start = atoi(ele->value);
	}

	/* Look for a genetic code qualifier /transl_table */
	transl_table = 1; /* Standard */
	if ((e = arr(ft_entry_ptr, fta, i).entry)) {
	    ft_value_element *ele;
	    ele = search_ft_qual_hash(e, "transl_table");
	    if (ele && ele->value)
		transl_table = atoi(ele->value);
	}
	if (load_genetic_code_number(transl_table) == -1) {
	    verror(ERR_WARN, "load_genetic_code_number",
		   "Failed to load code %d; using standard code",
		   transl_table);
	    load_genetic_code_number(1);
	}

	/* Look for /codon qualifiers and edit table accordingly */
	if ((e = arr(ft_entry_ptr, fta, i).entry)) {
	    ft_value_element *ele;
	    for (ele = search_ft_qual_hash(e, "codon");
		 ele;
		 ele = ele->next) {
		if (!ele->value)
		    continue;
		
		if (edit_genetic_code(ele->value) == -1) {
		    verror(ERR_WARN, "edit_genetic_code",
			   "Invalid /codon '%s'\n", ele->value);
		}
	    }
	}

	/* It does, so produce an mRNA */
	mRNA = create_mRNA(xx, seq, ranges, codon_start, 1,
			   &mapping, &mRNA_len); 
	if (!mRNA)
	    continue;

	/* Check the effect of the mutations? */
	check_mutations(xx, DB_Length(xx, seq), mRNA, mRNA_len, mapping, muts);

	xfree(mRNA);
	xfree(mapping);
    }

    html = dstring_create(NULL);
    sort_mutations(xx, muts, mut_coverage);
    dump_mutations(html, xx, muts, mut_coverage, sort_by_position,
		   dir, detail);

    free_mutations(xx, muts);
    xfree(mut_coverage);

    /* Tidy up */
    for (i = 0; i < ArrayMax(fta); i++) {
	Array r;

	r = arr(ft_entry_ptr, fta, i).ranges;
	if (r) {
	    int j;
	    for (j = 0; j < ArrayMax(r); j++)
		if (arr(ft_range *, r, j))
		    del_ft_range(arr(ft_range *, r, j));
	    ArrayDestroy(r);
	}

	if (arr(ft_entry_ptr, fta, i).entry)
	    del_ft_entry(arr(ft_entry_ptr, fta, i).entry);
    }
    ArrayDestroy(fta);

    return html;
}

/*
 * When we haven't got a reference sequence containing FCDS tags we have the
 * option of translation all or specific frames in their entirety. This is
 * achieved by create fake exons.
 */
static Array generate_exons(EdStruct *xx, int pos, int width) {
    Array fta, ranges;
    int i, j;
    ft_range *r;
    int fstart[3], fend[3];

    /*
     * Compute the unpadded base number for the left most part of the screen.
     * To do this we need to get the consensus (hopefully cached) and work
     * out the position that is divisible by 3 in unpadded form and is
     * <= pos as a padded value.
     */
    {
	/* Get consensus */
	char *con;
	int npads, i, frame;
	int count;

	if (!(con = xmalloc(DB_Length(xx, 0)+1)))
	    return NULL;
	DBcalcConsensus(xx, 1, DB_Length(xx, 0), con, NULL,
			BOTH_STRANDS);

	/* Count the number of pads from here to pos */
	for (i = npads = 0; i < pos; i++) {
	    if (i >= 0 && i < DB_Length(xx, 0) && con[i] == '*')
		npads++;
	}

	/* Back up to the start of a codon plus another 3 */
	count = 2;
	do {
	    while (--i >= 0 && i < DB_Length(xx, 0) && con[i] == '*')
		npads--;
	} while ((i - npads) % 3 || --count);

	/* Produce the 3 frame positions from here. */
	for (frame = 0; frame < 3; frame++) {
	    int count;

	    for (; ((i - npads) % 3) != frame; i++) {
		if (i >= 0 && i < DB_Length(xx, 0) && con[i] == '*')
		    npads++;
	    }
	    fstart[frame] = i+1;
	    
	    /*
	     * The frame ends we also wish to end on a codon so we can use
	     * these as the start points for the reverse strand translations.
	     * We simplify this by simply counting along until we are one
	     * codon clear beyond the end of pos+width.
	     */
	    count = 2;
	    j = i;
	    while (j < pos + width || --count) {
		int base;
		for (base = 0; base < 3; base++) {
		    while (j >= 0 && j < DB_Length(xx, 0) && con[j] == '*')
			j++;
		    j++;
		}
	    }
	    fend[frame] = j;
	}

	xfree(con);
    }
    
    

    fta = ArrayCreate(sizeof(ft_entry_ptr), 0);

    for (j = i = 0; j < 6; j++) {
	char name[10];

	if (!xx->status[EDITOR_SL_FRAME1p + j])
	    continue;

	/* Create ft_entry_ptr record */
	ArrayRef(fta, i);
	arr(ft_entry_ptr, fta, i).num = 1;
	arr(ft_entry_ptr, fta, i).ranges = ArrayCreate(sizeof(ft_range *), 0);
	arr(ft_entry_ptr, fta, i).sense = (j >= 3);
	arr(ft_entry_ptr, fta, i).entry = NULL;
	sprintf(name, "Frame %d", 1+j%3);
	sprintf(arr(ft_entry_ptr, fta, i).name, "%.*s", DB_NAMELEN, name);

	/* Create range */
	ranges = arr(ft_entry_ptr, fta, i).ranges;
	ArrayRef(ranges, 0);
	r = arr(ft_range *, ranges, 0) = new_ft_range();
	r->left = new_ft_location();
	r->left->min = fstart[j%3];
	r->right = new_ft_location();
	r->right->min = fend[j%3];
	r->complemented = j>=3;

	i++;
    }

    return fta;
}

