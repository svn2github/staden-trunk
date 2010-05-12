#include <staden_config.h>

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "os.h"
#include "getfile.h"
#include "array.h"
#include "misc.h" /* MIN MAX */
#include "sequence_formats.h"


char feat_quas[number_quas][20]={
    "/db_xref",		"/transl_table",	"/gene",
    "/allele",		"/anticodon",		"/bound_moiety",
    "/cell_line",	"/cell_type",		"/chloroplast",
    "/chromosome",	"/citation",		"/clone",
    "/clone_lib",	"/codon",		"/codon_start",
    "/cons_splice",	"/country",		"/cultivar",
    "/dev_stage",	"/direction",		"/EC_number",
    "/evidence",	"/exception",		"/focus",
    "/frequency",	"/function",		"/germline",
    "/haplotype",	"/insertion_seq",	"/isolate",
    "/label",		"/lab_host",		"/map",
    "/macronuclear",	"/mod_base",		"/note",
    "/number",		"/organelle",		"/organism",
    "/partial",		"/PCR_condition",	"/pop_variant",
    "/phenotype",	"/plasmid",		"/product",
    "/protein_id",	"/proviral",		"/pseudo",
    "/rearranged",	"/replace",		"/rpt_type",
    "/rpt_unit",	"/rpt_family",		"/sequenced_mol",
    "/serotype",	"/sex",			"/specific_host",
    "/specimen_voucher","/standard_name",	"/strain",
    "/sub_clone",	"/sub_species",		"/sub_strain",
    "/tissue_lib",	"/tissue_type",		"/transl_except",
    "/transposon",	"/usedin",		"/variety",
    "/virion"
};
	   

char feat_key[number_keys][16]={
    "CDS",		"tRNA",			"mRNA",
    "exon",		"intron",		"repeat_unit",
    "misc_feature",	"attenuator",		"C_region",
    "CAAT_signal",	"conflict",		"D-loop",
    "D_segment",	"enhancer",		"GC_signal",
    "gene",		"iDNA",			"J_segment",
    "LTR",		"mat_peptide",		"misc_binding",
    "misc_difference",	"misc_recomb",		"misc_RNA",
    "misc_signal",	"misc_structure",	"modified_base",
    "mutation",		"N_region",		"old_sequence",
    "polyA_signal",	"polyA_site",		"precursor_RNA",
    "prim_transcript",	"primer_bind",		"promoter",
    "protein_bind",	"RBS",			"repeat_region",
    "rep_origin",	"rRNA", 		"source",
    "S_region",		"satellite",	    	"scRNA",
    "sig_peptide",	"snRNA",		"stem_loop",
    "STS",		"TATA_signal",		"terminator",
    "transit_peptide",	"unsure",		"V_region",
    "V_segment",	"variation",		"3'clip",
    "3'UTR",		"5'clip",		"5'UTR",
    "-",		"-10_signal",		"-35_signal"
};

char genetic_code_ft[16][10]={
    " ",	"code_1",	"code_2",	"code_3",
    "code_4",	"code_5",	"code_6",	" ",
    " ",	"code_9",	"code_10",	"code_11",
    "code_12",	"code_13",	"code_14",	"code_15"
};

/* routines to read the following sequence file formats:

   staden
   embl
   genbank
   pir
   fasta
   gcg

   For embl, genbank and fasta an entryname can also be supplied to enable
   entries to be extracted from concatenated files. The default is to extract
   the first entry.

   seq_file_format works out the file format.
   get_seq_type returns 1 for dna, 2 for protein, 0 for anything else
   At present this facility is not used.
*/

#define STADEN 1
#define EMBL 2
#define GENBANK 3
#define CODATA 4
#define FASTA 5
#define GCG 6
#define MAX_SEQ_LINE 1024

#if 0
/* lists an array of chars, 50 per line */

void print_char_array ( FILE *file, char *array, int array_len) {
    
#define LINELENGTH 60
    int lines, i, i1, i2, j;
    lines = 1 + array_len/LINELENGTH;
    if ( (array_len % LINELENGTH) == 0 ) lines -= 1;
    for ( j=0; j <= lines; j++) {
	i1 =  j * LINELENGTH;
	i2 = MIN (i1+LINELENGTH-1,array_len-1);
	for ( i=i1; i <= i2; i++)
	    putc ( array[i], file );
	putc ( '\n', file );
    }
}
#endif

/*
 * realloc the sequence buffer when necessary
 */
int realloc_sequence(char **array,
		     int *array_size,
		     int increment)
{
    (*array_size) += increment;

    if (NULL == ((*array) = (char *)xrealloc((*array), 
					     *array_size * sizeof(char)))) {
	return -1;
    }

    return 0;
}

int dotty_gcg_format ( FILE *fp ) {

  /* if we get here we have recognised a possible format such as embl
     but it could also be a gcg mangled embl file. So here we allow
     max_line records to contain a .. which might signify it is gcg.
     Another waste of time.
  */
    char line[MAX_SEQ_LINE];
    int i, max_lines = 2;
    for ( i=0; i<max_lines; i++ ) {
      if ( fgets( line,sizeof(line),fp ) != NULL ) {

	if ( (strlen(line) > 3) && (strstr(line," .."))) return 1;
      }
    }
    return 0;
  }

int seq_file_format ( FILE *fp )

/* try to find the sequence file format */

/* note we cannot know if it is "staden"! so this is the default */

{


    char line[MAX_SEQ_LINE];

    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	if ( 0 == (strncmp("ID   ",line,5)) ) {
            if ( dotty_gcg_format ( fp )) return GCG;
	    return EMBL;
	  }
	else if ( 0 == (strncmp("LOCUS",line,5)) )
/*            if ( dotty_gcg_format ( fp ) return GCG; */
	    return GENBANK;
	else if ( 0 == (strncmp("SEQUENCE",line,8)) )
	    return CODATA;
	else if ( (strlen(line) > 3) && (strstr(line," ..")))
		return GCG;
	else if ( '>' == line[0] )
	    return FASTA;
	/* the only reason we bother with the following tests
	   is to avoid reading the whole file unnecessarily */
	else if ( ';' == line[0] )
	    return STADEN;
	else if ( '<' == line[0] )
	    return STADEN;
    }
    return STADEN;
}



int get_seq_type ( char *seq, int seq_len )

/* check sequence content to see if it looks ok: 85% a,c,g,t is DNA
   98% protein codes is protein, else is crap */

{
    char protein_chars[] = {"ARNDBCQEZGHILKMFPSTWXYV"};
    char dna_chars[] = {"ACGTUN"};
    char padding_chars[] = {"-*."};
    int i, dna = 0, protein = 0, padding = 0;
    float perc;

    if (seq_len < 1)
	return 0;

    for (i=0;i<seq_len;i++) {

	if (strchr(dna_chars,toupper(seq[i]))) dna++;
	if (strchr(protein_chars,toupper(seq[i]))) protein++;
	if (strchr(padding_chars,toupper(seq[i]))) padding++;
    }
    perc = (float) dna / (float) (seq_len - padding);
    if ( 0.85 < perc ) {
	return DNA;
    }
    perc = (float) protein / (float) (seq_len - padding);
    if ( 0.98 < perc ) {
	return PROTEIN;
    }
    return 0;
}

void
write_sequence(char *line,
	       char **seq,
	       int *seq_len,
	       int *buf_size)
{
    int j;
    int increment = 50000;

    for (j = 0; j < MAX_SEQ_LINE && line[j]; j++) {
	if ( isalpha ( (int) line[j]) || (int) line[j] == '-') {
	    if ( *seq_len+1 >= (*buf_size)) {
		realloc_sequence(seq, buf_size, increment);
	    }
	    (*seq)[*seq_len] = line[j];
	    *seq_len += 1;
	}
    }
    (*seq)[*seq_len] = 0; /* nul terminate */
}


void get_staden_format_seq ( char **seq, int max_len, int *seq_len, FILE *fp )

/* read in a staden format (yuk) sequence file */

/* Deal with 2 special line types: comments that have ";" in column 0
   and contig consensus sequence headers that have "<----abc.00001---->"
   embedded in them */
{


    char line[MAX_SEQ_LINE];
    int j;
    int buf_size = 0;
    int increment = 50000;

    *seq_len = 0;
    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	/* Check for special lines of type ";"*/

	if ( ';' != line[0] ) {

	    for (j = 0;j < MAX_SEQ_LINE && line[j]; j++) {

		if ( '<' == line[j] ) j += 20;
		if (isalpha ( (int) line[j]) || (int) line[j] == '-') {
		    if ( *seq_len >= buf_size) {
			realloc_sequence(seq, &buf_size, increment);
		    }
		    (*seq)[*seq_len] = line[j];
		    *seq_len += 1;
		}
	    }
	}
    }
}

int get_fasta_format_seq ( char **seq, int max_len, int *seq_len, FILE *fp,
			   char *entry_name, char **identifier)

/* read in a fasta format sequence file */

/* Assume entry starts with > in line[0] */

{


    char line[MAX_SEQ_LINE];
    int looking_for_sequence, looking_for_entry, expecting_sequence;
    int buf_size = 0;
    char *local_id;

    *seq_len = 0;

    if (!identifier)
	identifier = &local_id;
    if(NULL == (*identifier = (char *)xmalloc((MAX_SEQ_LINE) * sizeof(char))))
       return -1;

    if ( *entry_name ) {
	looking_for_entry = 1;
	looking_for_sequence = 0;
    }
    else {
	looking_for_entry = 0;
	looking_for_sequence = 1;
    }
    expecting_sequence = 0;

    while ( fgets( line,sizeof(line),fp ) != NULL ) {
	if ( looking_for_entry ) {
	    if ( 0 == (strncmp(">",line,1)) ) {
		char *end = line+1;
		while (!isspace(*end))
		    end++;
		*end = 0;

		if (0 == strcmp(entry_name, line+1)) {
		    looking_for_entry = 0;
		    expecting_sequence = 1;
		    strcpy((*identifier), entry_name);
		
		}
	    }
	}
	else if ( looking_for_sequence ) {
	    /* Check for header line of type ">"*/
	    if ( '>' == line[0] ) {
		if (1 != sscanf(line, ">%s\n", *identifier)) {
		    strcpy(*identifier, "MISSING_ID");
		}
		expecting_sequence = 1;
		looking_for_sequence = 0;
	    }
	}
	else if ( expecting_sequence ) {

	    /* Check for header line of type ">"*/
	    
	    if ( '>' == line[0] )
		return 0; 
	    write_sequence(line, seq, seq_len, &buf_size);
	}
    }

    if (identifier == &local_id)
	xfree(*identifier);

    return 0;  
}

/*
 * Frees memory used by a key_index structure
 */
void free_key_index(Featcds **key_index) {
    int i, j, k;
    BasePos *lthis, *lnext;
    
    if (!key_index)
	return;

    for (i = 0; i < number_keys; i++) {
	if (!key_index[i])
	    continue;

	for (j = 1; j <= key_index[i]->id; j++) {
	    /* Expression */
	    if (key_index[i][j].cdsexpr)
		xfree(key_index[i][j].cdsexpr);

	    /* Qualifiers */
	    if (key_index[i][j].qualifier) {
		for (k = 0; k < number_quas; k++) {
		    if (key_index[i][j].qualifier[k])
			xfree(key_index[i][j].qualifier[k]);
		}
	    }

	    /* BasePos list */
	    for (lthis = key_index[i][j].loca; lthis; lthis = lnext) {
		lnext = lthis->next;
		xfree(lthis);
	    }
	}

	xfree(key_index[i]);
    }

    xfree(key_index);
}

/*int purify_qual (char *qual) {
  
    char *t, *c;
    int in_q = 0;

    if (NULL == (t = (char *)malloc((strlen(qual) + 1)*sizeof(char)))) return -1;

    strcpy (t, qual);
    for (c = t; *c; c++) {
	if(*c == '"') {
	    if(!in_q) {
		in_q = 1; // in 
		c++;
	    } else {
		in_q = 0; // out 
	    }
	}
	if (in_q) {
	    *qual++ = *c;
	} else {
	    if(isspace(*c)) {
		*qual = 0;
		break;
	    }
	}
    }
    return 0;
    }*/

int purify_qual (char *qual) {
  
    int len;

    len = strlen (qual);
    while (isspace (qual[len - 1])) {	
	    qual [len - 1] = 0;
	    len --;
    }

    return 0;
}

/*
 * Parses loc_expr and remaining part of feature in file pointed to by fp.
 * Adds this to key_index[i].
 *
 * Returns 0 for success (last qualifier), 1 for next location,
 *        -1 for error
 */
static int add_feat (char *loc_expr, Featcds **key_index, int i, FILE *fp,
		     char *line) {
    int k;
    char *qua_expr = NULL;
    int current_qlen;

    if (NULL == (qua_expr = (char *)xmalloc((loc_len) * sizeof(char))))      
	return -1;

    if (NULL == (key_index[i][key_index[i]->id].cdsexpr =
		 (char *)xmalloc(loc_len * sizeof(char)))){     
	return -1;
    }
    strcpy(key_index[i][key_index[i]->id].cdsexpr, loc_expr);
    for(k = 0; k < number_quas; k++){
	if (NULL == (key_index[i][key_index[i]->id].qualifier[k] =
		     (char *)xmalloc(sizeof(char))))
	    return -1;
	*key_index[i][key_index[i]->id].qualifier[k] = 0;
    }	      
    for(;;){
      if (strncmp(line, "FT", 2) != 0 || line[21] != '/' ) 
	{
	xfree(qua_expr);
	return 1;
	} else {
	    strcpy(qua_expr, &line[21]);
	    qua_expr[strlen(qua_expr) - 1] = 0;
	    purify_qual (qua_expr);
	    while (strchr(qua_expr,'"') != NULL &&
		   qua_expr[strlen(qua_expr) - 1] != '"' &&
		   strncmp(&qua_expr[0], "/translation", 12)) {  
		fgets(line, loc_len, fp);	
		qua_expr[strlen(qua_expr) - 1] = 0;
		strcat(qua_expr, " ");
		strcat(qua_expr, &line[21]);
		qua_expr[strlen(qua_expr) - 1] = 0;
		purify_qual (qua_expr);
	      }
	    for (k = 0; k < number_quas; k++){
	      current_qlen=0;
	      while(strchr(qua_expr, '=' ) != NULL &&
		    qua_expr[current_qlen] != '=') current_qlen++;
	      if( (int)strlen(feat_quas[k]) > current_qlen)
		current_qlen = strlen(feat_quas[k]);
	      if(current_qlen != 0 &&
		 !strncmp( qua_expr, feat_quas[k], current_qlen)){		
		    char **qual = key_index[i][key_index[i]->id].qualifier;
		    qual[k] = (char *)xrealloc(qual[k],
					       strlen(qual[k]) + 1 +
					       strlen(qua_expr) +
					       1 /* nul */);
		    if (NULL == qual[k])
		      return -1;
		    if(strlen(qual[k]) > 1)
		     strcat(qual[k],"?");
		    strcat(qual[k], qua_expr);
		   
	      }
	    }
      }/* end else */
      fgets(line, loc_len, fp);
    }/*end for (;;) */
}




int get_embl_format_seq_no_ft ( char **seq, int max_len, int *seq_len, FILE *fp, 
			  char *entry_name)

/* read in an embl format sequence file */

/* Assume line before sequence starts with SQ and sequence ends with a line 
   commencing // 
   A typical line of sequence is:

     aagttgatgc agatcaatta atacgatacc tgcgtcataa ttgattattt gacgtggttt     48420

   so we have to throw away the base numbers.

   First if entry_name is filled in we must locate the appropriate entry.
   
   Returns 0 for success
          -1 for failure

*/

{


    char line[MAX_SEQ_LINE];
    int looking_for_sequence, looking_for_entry;

    int buf_size = 0;
    *seq_len = 0;


    if ( *entry_name ) {
	looking_for_entry = 1;
	looking_for_sequence = 0;
    } else {
	looking_for_entry = 0;
	looking_for_sequence = 1;
    }

    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	if ( looking_for_entry ) {

	    if ( 0 == (strncmp("ID",line,2)) ) {
		char *end = line+5;

		/*
		 * must deal with case of exp file ID line like
		 * "ID name\n" OR "ID name      \n"
		 */
		while (!isspace(*end))
		    end++;
		*end = 0;

		if (0 == strcmp(line+5, entry_name)) {
		    looking_for_entry = 0;
		    looking_for_sequence = 1;
		}
	    }
	} else {

	    if ( looking_for_sequence ) {
		
		/* look for the SQ line that precedes the sequence data */

		if ( 0 == (strncmp("SQ",line,2)) )
		    looking_for_sequence = 0;
	    } else {

		/* look for the // that follows the sequence */
		
		if ( 0 == (strncmp("//",line,2)) ) {
		    return 0;
		}
		write_sequence(line, seq, seq_len, &buf_size);
	    }
	}
    }
    return -1;
}

int purify_range (char *range) {

    char *t, *c;
 
    if (NULL == (t = (char *)xmalloc((strlen(range) + 1)*sizeof(char)))) return -1;
    strcpy (t, range);
    for (c = t; *c; c++) {
	if(!isspace(*c)) {
	    *range++ = *c;
	}
    }
    *range = 0;
    xfree (t);
    return 0;
}

/* modified for reading feature tables */
int get_embl_format_seq ( Featcds **key_index, char **seq, int max_len,
			  int *seq_len, FILE *fp, char *entry_name, 
                          char **identifier, int *err)

/* read in an embl format sequence file */

/* Assume line before sequence starts with SQ and sequence ends with a line 
   commencing // 
   A typical line of sequence is:

     aagttgatgc agatcaatta atacgatacc tgcgtcataa ttgattattt gacgtggttt     48420

   so we have to throw away the base numbers.

   First if entry_name is filled in we must locate the appropriate entry.

   Returns 0 for success,
          -1 for failure
*/

{


    /* added for reading feature tables */
    char *tmp = NULL;
    char *loc_expr = NULL;
    int i, k; 
    int inkey = 0;
   
    /*static int current_number_keys[number_keys]= */

    int current_number_keys[number_keys] =
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
     ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
     ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
     ,0,0,0}; 
    /* added for reading feature tables */
  
    char line[MAX_SEQ_LINE];
    int looking_for_sequence, looking_for_entry;
    int buf_size = 0;
    *seq_len = 0;

    /* added for reading feature tables */
    if (NULL == (loc_expr = (char *)xmalloc(sizeof(char))))    
	goto error;
    if (NULL == (tmp = (char *)xmalloc((MAX_SEQ_LINE) * sizeof(char))))
	goto error;
    if (NULL == (*identifier = (char *)xmalloc((MAX_SEQ_LINE+1)*sizeof(char))))
	goto error;

    /* need to initialise loc_expr */
    *loc_expr = 0;

    /* added for reading feature tables */
    if ( *entry_name ) {
	looking_for_entry = 1;
	looking_for_sequence = 0;
    } else {
	looking_for_entry = 0;
	looking_for_sequence = 1;
    }

    while ( fgets( line,sizeof(line),fp ) != NULL ) {
	if ( 0 == (strncmp("ID", line, 2)) ) {
	    if (1 != sscanf(line, "ID %20s\n", *identifier)) {
		strcpy(*identifier, "MISSING_ID");
	    }
       	}

	if ( looking_for_entry ) {
	    if ( 0 == (strncmp("ID", line, 2)) ) {
		char *end = line+5;

		/*
		 * must deal with case of exp file ID line like
		 * "ID name\n" OR "ID name      \n"
		 */
		while (!isspace(*end))
		    end++;
		*end = 0;

		if (0 == strcmp(line+5, entry_name)) {
		    looking_for_entry = 0;
		    looking_for_sequence = 1;
		    *identifier = strdup(entry_name);
		}
	    }

	} else {
	    if ( looking_for_sequence ) {
	     
		/* added for reading feature tables */
	    ft_cds:

		
		if(key_index && !strncmp(&line[0],"FT",2)){
		    inkey = 0;	     
		    for(i = 0; i < number_keys; i++){
			for(k = 0; k < 20; k++){
			    if(!isgraph(line[5+k]))
				break;
			}		    
			if (k != 0 && !strncmp(&line[5], feat_key[i],
					       max(k, (int)strlen(feat_key[i])))) {
			    inkey++;
			    loc_expr = (char *)xrealloc(loc_expr,
							strlen(loc_expr) + 1
							+ strlen(&line[21]) + 1);
			    strcpy(loc_expr, &line[21]);
			    loc_expr[strlen(loc_expr) - 1] = 0;
			    purify_range ( loc_expr);

			    /* to get next line to judge whether feature location 
			       reading has finished or not */
			    fgets(line, loc_len, fp);			    
			    while(line[5] == ' ' && strncmp(&line[21], "/", 1) && strncmp(&line[0],"XX",2)) {
				strcpy(tmp, &line[21]);
				tmp[strlen(tmp)-1] = 0;
				purify_range (tmp);
				loc_expr = (char *)xrealloc(loc_expr,
							    strlen(loc_expr) + 1
							    + strlen(tmp) + 1);
				strcat(loc_expr, tmp);
				fgets(line, loc_len, fp);
			    }
			    if (current_number_keys[i] <= key_index[i]->id){
				current_number_keys[i]++;  
				key_index[i] =
				    (Featcds *)xrealloc(key_index[i],
							(current_number_keys[i]+1)*
							sizeof(Featcds));
				if (NULL == key_index[i]) {
				    verror(ERR_WARN, "Load sequence", 
					   "Error reading FT line '%s'\n", line);
				    goto error;
				}
			    }
			    
			    if( parse_feat(loc_expr, key_index, i)){
				int r = add_feat(loc_expr, key_index, i, fp, line);
				if (r == 1) {
				    goto ft_cds;
				}
				if (r == -1) {
				/* error */
				    verror(ERR_WARN, "Load sequence",
					   "Error reading FT line '%s'\n", line);
				    goto error;
				}
			    }/* end if(parse */
			    else {
				goto ft_cds;
			    }
			}/* end if(k  */
		    }/* end for(i */	
		    if(inkey == 0 && line[5] !=' '){
			*err = 1;
			/*      verror(ERR_WARN, "Load Sequence",
				"Error in feature table\n");*/
		    }
		}/* if(... "FT"  */
		/* added for reading feature tables */
	    
		/* look for the SQ line that precedes the sequence data */

		if ( 0 == (strncmp("SQ",line,2)) ){
		    looking_for_sequence = 0;
		}
	    } else {		
		if ( 0 == (strncmp("//",line,2)) ) {
		    goto success;
		}
		write_sequence(line, seq, seq_len, &buf_size);	
	    }
	}    
    }

 success:
    xfree(loc_expr);
    xfree(tmp);
    return 0;

 error:
    if (loc_expr)
	xfree(loc_expr);
    if (tmp)
	xfree(tmp);
    return -1;
}
   
void get_genbank_format_seq ( char **seq, int max_len, int *seq_len, FILE *fp,
			     char *entry_name)

/* read in an genbank format sequence file */

/* Assume line before sequence starts with ORIGIN and sequence ends with a line 
   commencing // 
   A typical line of sequence is:

     aagttgatgc agatcaatta atacgatacc tgcgtcataa ttgattattt gacgtggttt     48420

   so we have to throw away the base numbers.
*/

{


    char line[MAX_SEQ_LINE];
    int looking_for_sequence, looking_for_entry;
    int buf_size = 0;

    *seq_len = 0;

    if ( *entry_name ) {
	looking_for_entry = 1;
	looking_for_sequence = 0;
    }
    else {
	looking_for_entry = 0;
	looking_for_sequence = 1;
    }

    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	if ( looking_for_entry ) {

	    if ( 0 == (strncmp("LOCUS",line,5)) ) {
		char *end = line+12;
		while (!isspace(*end))
		    end++;
		*end = 0;

		if (0 == strcmp(entry_name, line+12)) {
		    looking_for_entry = 0;
		    looking_for_sequence = 1;
		}
	    }
	}
	else {

	    if ( looking_for_sequence ) {

		/* look for the ORIGIN line that precedes the sequence data */

		if ( 0 == (strncmp("ORIGIN",line,6)) )
		    looking_for_sequence = 0;
	    }
	    else {

		/* look for the // that follows the sequence */

		if ( 0 == (strncmp("//",line,2)) )
		    return;

		write_sequence(line, seq, seq_len, &buf_size);
	    }
	}
    }
}

void get_pir_format_seq ( char **seq, int max_len, int *seq_len, FILE *fp )

/* read in an pir codata format sequence file */

/* Assume line before sequence starts with SEQUENCE and sequence ends with a line 
   commencing ///
   A typical line of sequence is:

     31 A L M G L G T L Y F L V K G M G V S D P D A K K F Y A I T T             

   so we have to throw away the base numbers.
*/

{


    char line[MAX_SEQ_LINE];
    int looking_for_sequence;
    int buf_size = 0;

    *seq_len = 0;

    looking_for_sequence = 1;

    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	if ( looking_for_sequence ) {

	    /* look for the SEQUENCE line that precedes the sequence data */

	    if ( 0 == (strncmp("SEQUENCE",line,8)) )
		looking_for_sequence = 0;
	    }

	else {

	    /* look for the /// that follows the sequence */

	    if ( 0 == (strncmp("///",line,3)) )
		return;

	    write_sequence(line, seq, seq_len, &buf_size);
	}
    }
}

void get_gcg_format_seq ( char **seq, int max_len, int *seq_len, FILE *fp )

/* read in a gcg format sequence file */

/* Assume line before sequence contains " .." and sequence ends at the end of file
   A typical few lines of sequence is:

EGFR_HUMAN  Length: 1210  November 26, 1997 16:11  Type: P  Check: 1521  ..
 
       1  MRPSGTAGAA LLALLAALCP ASRALEEKKV CQGTSNKLTQ LGTFEDHFLS 
 
      51  LQRMFNNCEV VLGNLEITYV QRNYDLSFLK TIQEVAGYVL IALNTVERIP 

*/

{


    char line[MAX_SEQ_LINE];
    int looking_for_sequence;
    int buf_size = 0;

    *seq_len = 0;

    looking_for_sequence = 1;

    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	if ( looking_for_sequence ) {

	    /* look for the ..\n line that precedes the sequence data */

	    if ( (strlen(line) > 3) && (strstr(line," .."))) {
		looking_for_sequence = 0;
	      }
	  }
	else {

	    write_sequence(line, seq, seq_len, &buf_size);
	}
    }
}

/* modified for reading feature tables */
int get_seq ( char **seq, int max_len, int *seq_len, char *file_name, char *entry_name_in)

{

/* completion codes: 
   0 OK (but still can have 0 length sequence, meaning that the required data was not found)
   1 file open failed
   2 not the expected type of sequence - i.e. not dna or protein
   3 unrecognised format
   4 seek to start of file failed!
   */

    char entry_name[256];
    FILE *file_ptr;
    int fmt;

    entry_name[0] = '\0';
    if ( entry_name_in && entry_name_in[0] ) {
	strcpy(entry_name,entry_name_in);
    }

    if ( file_ptr = my_fopen ( file_name, "r" ) ) {
 

	/* determine the file format */

	if ( fmt = seq_file_format ( file_ptr ) ) {

	    if (fseeko ( file_ptr, 0, SEEK_SET ) ) return 4;

	    if ( STADEN == fmt ) {
		(void) get_staden_format_seq (seq, max_len, seq_len, file_ptr );
		/* for this catchall, default format we had better check what weve
		   got looks like a dna or protein sequence! */
		if ( seq_len ) {
		    if ( !get_seq_type ( *seq, *seq_len )) {
			*seq_len = 0;
			return 2;
		    }
		}
	    }
	    else if ( EMBL == fmt )
	      {
		  if (get_embl_format_seq_no_ft ( seq, max_len, seq_len,
					    file_ptr, entry_name))
		      return 3;
	      }
	    else if ( CODATA == fmt )
		(void) get_pir_format_seq ( seq, max_len, seq_len, file_ptr );
	    else if ( GCG == fmt )
		(void) get_gcg_format_seq ( seq, max_len, seq_len, file_ptr );
	    else if ( GENBANK == fmt )
		(void) get_genbank_format_seq ( seq, max_len, seq_len, file_ptr, entry_name);
	    else if ( FASTA == fmt )
		(void) get_fasta_format_seq ( seq, max_len, seq_len, file_ptr, entry_name, NULL);
	}
	else {
	    /* we never get here because "staden format" is the default */
	    return 3;
	}
    }
    else {
	return 1;
    }
    fclose ( file_ptr );
    return 0;
}

int get_seq_ft (Featcds **key_index,  char **seq, int max_len, int *seq_len, char *file_name, char *entry_name_in, char **identifier, int *err)

{

/* completion codes: 
   0 OK (but still can have 0 length sequence, meaning that the required data was not found)
   1 file open failed
   2 not the expected type of sequence - i.e. not dna or protein
   3 unrecognised format
   4 seek to start of file failed!
   5 out of memory error
   */

    char entry_name[256];
    FILE *file_ptr;
    int fmt;
    int retval = 0;


    entry_name[0] = '\0';
    if ( entry_name_in && entry_name_in[0] ) {
	strcpy(entry_name,entry_name_in);
    }



    /* Open file */
    file_ptr = my_fopen ( file_name, "r" );
    if( !file_ptr )
	return 1;



    /* Determine the file format */
    fmt = seq_file_format( file_ptr );
    if( !fmt )
    {
	fclose( file_ptr );
        return 3;
    }



    /* Go to start of file */
    if( fseeko(file_ptr, 0, SEEK_SET) )
    {
        fclose( file_ptr );
	return 4;
    }



    /* Read file */
    switch( fmt )
    {
	case EMBL:
	    if(get_embl_format_seq (key_index, seq, max_len, seq_len,
				    file_ptr, entry_name, &*identifier,&*err))
		retval = 3;
	    break;


	case CODATA:
	    get_pir_format_seq ( seq, max_len, seq_len, file_ptr );
	    break;


	case GCG:
	    get_gcg_format_seq ( seq, max_len, seq_len, file_ptr );
	    break;


	case GENBANK:
	    get_genbank_format_seq ( seq, max_len, seq_len, file_ptr, entry_name );
	    break;


	case FASTA:
	    get_fasta_format_seq( seq, max_len, seq_len, file_ptr, entry_name, &*identifier );
	    break;


	default:
	    get_staden_format_seq (seq, max_len, seq_len, file_ptr );
	    /* Ensure it's a dna or protein sequence */
	    if ( seq_len ) {
		if ( !get_seq_type ( *seq, *seq_len )) {
		   *seq_len = 0;
		   retval   = 2;
		}
	    }
	    break;
    }

    /* Cleanup & exit */
    fclose( file_ptr );
    return retval;
}

f_proc_ret getseq_ ( char *seq, 
		    f_int *f_max_len, 
		    f_int *f_seq_len, 
		    char *file_name_in,
		    char *entry_name_in, 
		    f_implicit file_1, 
		    f_implicit name_1)

{

    int i, max_len, seq_len;
    char entry_name[256], file_name[51];
    entry_name[0] = '\0';

    for (i=0;i<50;i++) {
	file_name[i] = file_name_in[i];
    }
    i = 14;

    file_name[i] = '\0';

    max_len = *f_max_len;

    i = get_seq ( &seq, max_len, &seq_len, file_name, entry_name);
    *f_seq_len = seq_len;
    
    f_proc_return();
}

/*
 * 
 */
int realloc_char_array(char ***array,
		       int *array_size,
		       int word_len)
{
    int increment = 100;
    int prev_array_size;
    int i;

    prev_array_size = *array_size;
    (*array_size) += increment;

    if (NULL == ( (*array) = (char **)xrealloc( (*array), 
					   *array_size * sizeof(char *)))) {
	return -1;
    }

    for (i = prev_array_size; i < *array_size; i++) {
    
	if (NULL == ((*array)[i] = (char *)xmalloc((word_len+1) * sizeof(char)))) {
	    
	    return -1;
	}
    }
    return 0;
}

int get_identifiers (char *file_name, 
		     char ***list,
		     int *num_identifiers)
{
#define NAME_LEN 50

/* completion codes: 
   0 OK (but still can have 0 length sequence, meaning that the required data was not found)
   1 file open failed
   2 not the expected type of sequence - i.e. not dna or protein
   3 unrecognised format
   4 seek to start of file failed!
   */

    FILE *file_ptr;
    int fmt;
    int cnt = 0;
    char **identifier = NULL;
    int array_size = 0;
    char line[MAX_SEQ_LINE];

    if (NULL == (file_ptr = my_fopen ( file_name, "r" )))
	return 1;
    
    /* determine the file format */

    
    if ( fmt = seq_file_format ( file_ptr ) ) {
	
	if (fseeko ( file_ptr, 0, SEEK_SET ) ) return 4;

	if ( EMBL == fmt ) {
	    while ( fgets( line, sizeof(line), file_ptr ) != NULL ) {
		if (cnt >= array_size) 
		    realloc_char_array(&identifier, &array_size, NAME_LEN);
		if (sscanf(line, "ID %20s\n", identifier[cnt]) == 1) {
		    cnt++;
		}
	    }
	    
	} else if ( GENBANK == fmt ) {
	    while ( fgets( line, sizeof(line), file_ptr ) != NULL ) {
		if (cnt >= array_size) 
		    realloc_char_array(&identifier, &array_size, NAME_LEN);
		if (sscanf(line, "LOCUS       %14s\n", identifier[cnt]) == 1) {
		    cnt++;
		}
	    }
	} else if ( FASTA == fmt ) {
	    while ( fgets( line, sizeof(line), file_ptr ) != NULL ) {
		if (cnt >= array_size) 
		    realloc_char_array(&identifier, &array_size, NAME_LEN);
		if (sscanf(line, ">%50s\n", identifier[cnt]) == 1) {
		    cnt++;
		}
	    }

	} else if ( STADEN == fmt ) {
	    while ( fgets( line,sizeof(line),file_ptr) != NULL ) {

		if (cnt >= array_size) 
		    realloc_char_array(&identifier, &array_size, NAME_LEN);
		if (sscanf(line, "<%18s>", identifier[cnt]) == 1) {
		    cnt++;
		}
	    }
	} else {
	    return 3;
	}
    }

    fclose ( file_ptr );
    *list = identifier;
    *num_identifiers = cnt;
    return 0;
}


/* added for reading feature tables */
int parse_feat(char *locexpr,  Featcds **key_index, int i)
{  
    char *locexpr1 = NULL;  
    char *locexpr2 = NULL;  
    char *tmp = NULL;
    char type_range[2]=" ";
    int  start_pos, end_pos;
    BasePos *head=NULL;
    BasePos *tail=NULL;
    BasePos *current;
    int retcode = -1;

    if (NULL == (locexpr1=(char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	goto ret;
    if (NULL == (locexpr2=(char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	goto ret;
    if (NULL == (tmp=(char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	goto ret;

    if(!strncmp(locexpr,"complement(", 11)){
	sscanf(locexpr,"%11s%s", tmp,locexpr1);
	if(!strncmp(locexpr1,"join(", 5)){
	    /*complement(join*/	  
	    if(!read_cds_pos_join(&head,locexpr1)) {
		retcode = 0;
		goto ret;
	    }
	    else {   
		key_index[i]->id++; 
		key_index[i][key_index[i]->id].id=key_index[i]->id;
		sprintf(key_index[i][key_index[i]->id].type_loca,"%s","cj");
		key_index[i][key_index[i]->id].loca=head;   
	    }
       	}
	else {
	    /*complement*/
	    if(!read_cds_pos(locexpr1, &start_pos, &end_pos)) {
		retcode = 0;
		goto ret;
	    }
	    else {
		key_index[i]->id++;
		/* printf("%d  %d\n",start_pos,end_pos);*/
		key_index[i][key_index[i]->id].id=key_index[i]->id;
		sprintf(key_index[i][key_index[i]->id].type_loca,"%s","c");
		strcpy(type_range, "n");
		tail=add_list_item(&head,tail, start_pos, end_pos, type_range);
		key_index[i][key_index[i]->id].loca=head;
	    }
	}
    }
    /*join*/
    else if(!strncmp(locexpr,"join(", 5)){
	if(!read_cds_pos_join(&head,locexpr)) {
	    retcode = 0;
	    goto ret;
	} else {
	    key_index[i]->id++; 
       	    key_index[i][key_index[i]->id].id=key_index[i]->id;
	    sprintf(key_index[i][key_index[i]->id].type_loca,"%s","j");
	    key_index[i][key_index[i]->id].loca=head;
	    for(current=key_index[i][key_index[i]->id].loca; 
		current!=NULL; current=current->next);
	}
    }/*locexpr==join*/
    else{
	/*x..y*/
	if(!read_cds_pos(locexpr, &start_pos, &end_pos)) {
	    retcode = 0;
	    goto ret;
	} else {
	    key_index[i]->id++;
	    key_index[i][key_index[i]->id].id=key_index[i]->id;
	    sprintf(key_index[i][key_index[i]->id].type_loca,"%s","n");
	    strcpy(type_range,"n");
	    tail=add_list_item(&head,tail, start_pos, end_pos, type_range);
	    key_index[i][key_index[i]->id].loca=head;
	}
    }

    retcode = 1;
 ret:
    if (locexpr1)
	free(locexpr1);
    if (locexpr2)
	free(locexpr2);
    if (tmp)
	free(tmp);

    return retcode;
}/*end of parse_feat*/

int read_cds_pos(char *locexpr, int *start_pos, int *end_pos){

    int i=0,n=0; 
    int locleng;     
    char *a,*b;
    int j=0;

    locleng=strlen(locexpr);
    if (NULL == (a = (char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	return -1;
    if (NULL == (b = (char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	return -1;

    if (locexpr[0] == '<' || !isdigit(locexpr[i])){
	free(a);
	free(b);
	return (0);
    }
    
    while( locexpr[i] != '.' ) {
	a[n] = locexpr[i];
	i++, n++;
	if(i == locleng-1) {
	    free(a);
	    free(b);
	    return (0);
	}
    }
    a[n] = 0;
    while( locexpr[i]=='.' )
    i++, n = 0;
    for(j = 0; j < locleng; j++){
	if (locexpr[j] == '>') {
	    free(a);
	    free(b);
	    return (0);
	}
    }
    while( i < locleng )
	b[n++] = locexpr[i++];
    b[n] = 0; 
    n=0,i=0;
    *start_pos = atoi(a);
    *end_pos = atoi(b);
    free(a);
    free(b);
    return(1);
}

int read_cds_pos_join(BasePos **head,char *locexpr){

    char *locexpr1 = NULL;  
    char *locexpr2 = NULL;  
    char *tmp = NULL;
    char *tmp1 = NULL;
    int retcode = -1;
    char type_range[2]=" ";
    int start_pos, end_pos;
    BasePos *tail=NULL; 

    if (NULL == (locexpr1 = (char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	goto ret;
    if (NULL == (locexpr2 = (char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	goto ret;
    if (NULL == (tmp1 =( char*)xmalloc((strlen(locexpr)+1)*sizeof(char))))
	goto ret;
  
    sscanf(locexpr,"%5s%s", tmp1, locexpr1);
    if(!strncmp(locexpr1, "complement(", 11)){
      strcpy(type_range, "c");
      sscanf(locexpr1, "%11s%s", tmp1, locexpr2);
      strcpy(locexpr1, locexpr2);   
    } else
      {
	if(!strncmp(locexpr1, "join(", 5))
	    goto ret;
	strcpy(type_range,"n"); 
      }
    if(!read_cds_pos(locexpr1, &start_pos, &end_pos)) {
	retcode = 0;
	goto ret;
    }
    tail=add_list_item(&*head,tail, start_pos, end_pos, type_range); 
    tmp=strchr(locexpr1,',');
    if(!strncmp(tmp,",complement(", 12)){
	sscanf(tmp,"%12s%s", tmp1,locexpr2);
	strcpy(type_range,"c");
    }else
	{
	  sscanf(tmp,",%s",locexpr2);
	  strcpy(type_range,"n");
	}
    while(locexpr2!=NULL){
      if(!read_cds_pos(locexpr2, &start_pos, &end_pos)) {
	  retcode = 0;
	  goto ret;
      }
      tail=add_list_item(&*head, tail, start_pos, end_pos, type_range);
      tmp=strchr(locexpr2, ',');
      if(tmp!=NULL){
	sscanf(tmp,",%s", locexpr2);
	if(!strncmp(locexpr2, "complement(", 11)){
	  strcpy(type_range, "c");
	  sscanf(locexpr2, "%11s%s", tmp1,locexpr1); 
	  strcpy(locexpr2, locexpr1);   
	}else
	  {
	    strcpy(type_range,"n"); 
	  }
      } else {
	  if (locexpr2)
	      free(locexpr2);
	  locexpr2 = NULL;
      }
    }

    retcode = 1;

 ret:
    if (locexpr1)
	free(locexpr1);
    if (locexpr2)
	free(locexpr2);
    if (tmp1)
	free(tmp1);

    return retcode;
}

BasePos *add_list_item(BasePos **head, BasePos *entry, int start_pos, 
                      int end_pos, char *type_range){

    BasePos *new_list_item;
    new_list_item=xmalloc(sizeof(BasePos));

    if(entry==NULL) *head=new_list_item;
    else  entry->next=new_list_item;

    new_list_item->start_pos=start_pos;
    new_list_item->end_pos=end_pos;
    strcpy(new_list_item->type_range,type_range);
    new_list_item->next=NULL;
    return new_list_item;
}



int display_info(FILE *pw, Featcds **key_index){
  int i,j,k;
  BasePos *current;

 for(i=0; i<number_keys; i++){
    fprintf(pw,"The feature information for %s...\n", feat_key[i]);
    fprintf(pw,"--------------------------------------------------\n");

  for(k=1; k <= key_index[i]->id; k++){
    fprintf(pw,"%d    %s   ",k,key_index[i][k].type_loca);
   
    for(current=key_index[i][k].loca; current!=NULL; current=current->next)
      fprintf(pw," %s %d..%d    ", current->type_range,
                                 current->start_pos,current->end_pos );
   
    fprintf(pw,"\n\n");
    for(j=0; j < number_quas; j++){
      if(strlen(key_index[i][k].qualifier[j])>1) fprintf(pw,"%s %s", 
                            feat_quas[j], key_index[i][k].qualifier[j]);
    }/* end j */
      fprintf(pw,"\n");
  }/* end k */
 }/* end i */
  return 1;
}
/* added for reading feature tables */



int vmsg_info(Featcds **key_index){
  int i,j,k;
  BasePos *current;
  int num_ft=0;

  for(i=0; i<number_keys; i++){
   num_ft+=key_index[i]->id;
  }
  if(num_ft){

 for(i=0; i<number_keys; i++){
    vmessage("The feature information for %s...\n", feat_key[i]);
    vmessage("--------------------------------------------------\n");

   for(k=1; k <= key_index[i]->id; k++){
    vmessage("%d    %s   ",k,key_index[i][k].type_loca);
   
    for(current=key_index[i][k].loca; current!=NULL; current=current->next)
      vmessage(" %s %d..%d    ", current->type_range,
                                 current->start_pos,current->end_pos );
   
    vmessage("\n\n");
    for(j=0; j < number_quas; j++){
      if(strlen(key_index[i][k].qualifier[j])>1) vmessage("%s", 
                             key_index[i][k].qualifier[j]);
    }/* end j */
      vmessage("\n");
   }/* end k */
 }/* end i */
  return 1;
  }/* end if(num_ft ) */
  return -1;
}
