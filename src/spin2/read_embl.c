#include <stdio.h>
#include <string.h>
#include "sequence_formats.h"

#define MAX_LEN 500000

int get_embl(char *file) {
    int ret;
    int i;
    
    /* Allocated and zeroed by us and input to get_seq_ft */
    Featcds **key_index = NULL;

    /* Items returned from get_seq_ft */
    char *seq = NULL;
    char *id = NULL;
    int seq_len = 0;
    int err = 0;

    key_index = (Featcds **)xmalloc(number_keys * sizeof(Featcds*));
    if (NULL == key_index)
	goto fail;

    for (i = 0; i < number_keys; i++)
	key_index[i] = NULL;

    for (i=0; i < number_keys; i++){
	if (NULL == (key_index[i] = (Featcds *)xmalloc(sizeof(Featcds))))  
	    goto fail;
	key_index[i]->id=0;
    } 


    ret = get_seq_ft(key_index, &seq, MAX_LEN, &seq_len, file, NULL, &id,
		     &err);

    printf("ret=%d, seq_len=%d, err=%d\n", ret, seq_len, err);
    printf("id='%s'\n", id);

    /* Iterate around feature key _types_ (CDS, tRNA, etc) */
    for (i = 0; i < number_keys; i++) {
	int j;

	printf("Key %d (%s)\n", i, feat_key[i]);
	printf("    id=%d\n", key_index[i]->id);

	/* Iterate around each instance (feature) of this key type */
	for (j = 1; j <= key_index[i]->id; j++) {
	    BasePos *loca;

	    printf("    %s instance %d\n", feat_key[i], j);
	    loca = key_index[i][j].loca;

	    /* Iterate around the location for this feature */
	    while (loca) {
		printf("        loca=%d..%d type='%s'\n",
		       loca->start_pos,
		       loca->end_pos,
		       loca->type_range);
		loca = loca->next;
	    }
	    printf("        type_loca='%s'\n", key_index[i][j].type_loca);
	    if (key_index[i][j].cdsexpr) {
		printf("        cdsexpr='%s'\n", key_index[i][j].cdsexpr);
	    }

	    /* Iterate around the qualifier types */
	    if (key_index[i][j].qualifier) {
		int q;
		for (q = 0; q < number_quas; q++) {
		    char *c = key_index[i][j].qualifier[q];

		    /* Iterate around each instance of this qualifier type */
		    while (c && *c) {
			char *cp;

			if (cp = strstr(c, "?/")) {
			    printf("        qualifier[%d(%s)]='%.*s'\n",
				   q, feat_quas[q], cp - c, c);
			    c = cp + 1;
			} else {
			    printf("        qualifier[%d(%s)]='%s'\n",
				   q, feat_quas[q], c);
			    c = NULL;
			}
		    }
		}
	    }
	}
    }
    /*printf("seq='%s'\n", seq);*/

    return ret;
    
 fail:
    if (key_index) {
	for (i = 0; i < number_keys; i++) {
	    if (key_index[i])
		xfree(key_index[i]);
	}
	xfree(key_index);
    }

    return -1;
}

int main(int argc, char **argv) {
    return get_embl(argv[1]);
}
