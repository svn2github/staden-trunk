#include "base_comp.h"
#include "text_output.h"
#include "edge.h"
#include "dna_utils.h"
#include "seq_results.h"

void sequence_info(char *seq_name,
		   char *sequence,
		   int start,
		   int end,
		   int seq_structure,
		   int seq_type)
{
    double base_comp[5];
    double aa_comp[25];
    char aa[] = {"ABCDEFGHIKLMNPQRSTVWYZX*-"};
    int i;
    double aa_mass[25];
    int seq_length;

    seq_length = end - start + 1;
    vmessage("Sequence %s: %d to %d\n", seq_name, start, end);
    
    if (seq_type == 1) {
	if (seq_structure == 0) {
	    vmessage("linear ");
	} else {
	    vmessage("circular ");
	}

	vmessage("DNA\n");
	set_char_set(DNA);
	get_base_comp(&sequence[start-1], seq_length, base_comp);
	
	vmessage("Sequence composition\n");
	vmessage("\tA %d (%.2f%%) C %d (%.2f%%) G %d (%.2f%%) T %d (%.2f%%) - %d (%.2f%%)\n", 
		 (int)base_comp[0], base_comp[0]/seq_length*100.0, 
		 (int)base_comp[1], base_comp[1]/seq_length*100.0, 
		 (int)base_comp[2], base_comp[2]/seq_length*100.0,
		 (int)base_comp[3], base_comp[3]/seq_length*100.0,
		 (int)base_comp[4], base_comp[4]/seq_length*100.0); 
	vmessage("Mass %f\n", get_base_comp_mass((int)base_comp[0], 
						 (int)base_comp[1],
						 (int)base_comp[2], 
						 (int)base_comp[3]));
    } else {
	vmessage("Protein\n");
	set_char_set(PROTEIN);
	get_aa_comp(&sequence[start-1], seq_length, aa_comp);
	get_aa_comp_mass(aa_comp, aa_mass);

	/* amino acid name */
	vmessage("AA ");
	for (i = 0; i < 13; i++) {
	    vmessage(" %-5c", aa[i]);
	}
	vmessage("\n");

	/* number of each aa */
	vmessage("N  ");
	for (i = 0; i < 13; i++) {
	    vmessage("%-6g", aa_comp[i]);
	}
	vmessage("\n");

	/* % of each aa */
	vmessage("%%  ");
	for (i = 0; i < 13; i++) {
	    vmessage("%-6.1f", aa_comp[i]/seq_length*100.0);
	}
	vmessage("\n");

	/* mass of each aa */
	vmessage("M  ");
	for (i = 0; i < 13; i++) {
	    vmessage("%-6.0f", aa_mass[i]);
	}
	vmessage("\n\n");

	/* amino acid name */
	vmessage("AA ");
	for (i = 13; i < 25; i++) {
	    vmessage(" %-5c", aa[i]);
	}
	vmessage("\n");

	/* number of each aa */
	vmessage("N  ");
	for (i = 13; i < 25; i++) {
	    vmessage("%-6g", aa_comp[i]);
	}
	vmessage("\n");

	/* % of each aa */
	vmessage("%%  ");
	for (i = 13; i < 25; i++) {
	    vmessage("%-6.1f", aa_comp[i]/seq_length*100.0);
	}
	vmessage("\n");

	/* mass of each aa */
	vmessage("M  ");
	for (i = 13; i < 25; i++) {
	    vmessage("%-6.0f", aa_mass[i]);
	}
	vmessage("\n");
    }
}

double get_seq_mass (int seq_num)
{
    double mass = 0.0;

    if (GetSeqType(seq_num) == 1) {
	double base_comp[5];

	get_base_comp(GetSeqSequence(seq_num), GetSeqLength(seq_num), 
		      base_comp);
	
	mass = get_base_comp_mass((int)base_comp[0], 
				  (int)base_comp[1],
				  (int)base_comp[2], 
				  (int)base_comp[3]);
    } else {
	double aa_mass[25];
	double aa_comp[25];
	int i;

	get_aa_comp(GetSeqSequence(seq_num), GetSeqLength(seq_num), aa_comp);
	get_aa_comp_mass(aa_comp, aa_mass);
	for (i = 0; i < 25; i++) {
	    mass += aa_mass[i];
	}
    }
    return mass;
}


