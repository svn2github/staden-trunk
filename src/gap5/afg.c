/*
 * afg.h - afg to gap5 conversion for tg_index
 *
 * Andrew Whitwham, January 2011
 *
 * AFG contig format:
 
    {CTG
    iid:1336
    eid:1435-0
    seq:
    GCTGAGGTTGGCACAAATGTAAACATTGGTTGTGGCTCTATCACTGTCAACTACGACGGA
    AAAAATAAATTCTTAACAAAGATCGAAGACAATGCATTTATTGGATGTAATTCAAATTTA
    GTGGCTCCAGTTAC
    .
    qlt:
    GCTGAGGTTGGCACAAATGTAAACATTGGTTGTGGCTCTATCACTGTCAACTACGACGGA
    AAAAATAAATTCTTAACAAAGATCGAAGACAATGCATTTATTGGATGTAATTCAAATTTA
    GTGGCTCCAGTTAC
    .
    {TLE
    src:1624978
    off:13
    clr:0,54
    }
    {TLE
    src:6278089
    off:59
    clr:54,0
    }
    }
    
    
    AFG read format:
    
    {RED
    iid:1
    eid:1
    seq:
    GGCATTGATTGAAGAACACTCAGAAGACGAAGGCGGTTCGACTATCTATCGTCA
    .
    qlt:
    GGCATTGATTGAAGAACACTCAGAAGACGAAGGCGGTTCGACTATCTATCGTCA
    .
    }
    
 *
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/time.h> 

#include "tg_gio.h"
#include "tg_index_common.h"
#include "hache_table.h"
#include "afg.h"


// version of the perl favourite
static long chomp(char *str) {
    long length = 0;

    if (str && (length = strlen(str))) {
    	if (str[length - 1] == '\n') {
    	    length--;
    	    str[length] = '\0';
    	}
    }

    return length;
}


/*
   if string matches key, return the rest of the string as
   the value.  Removes white space at the value ends.
   Modifies the field string.
*/
static char *get_value(char *key, char *field) {
    size_t key_length, field_length;
    char *value = NULL;
    
    key_length   = strlen(key);
    field_length = strlen(field);
    
    if (strncmp(field, key, key_length) == 0) {
    	size_t i;
	
	for (i = key_length; i < field_length; i++) {
	    if (!isspace(*(field + i))) break;
	}
	
	if (i != field_length) {
	    value = field + i;
	
	    for (i = field_length - 1; i > key_length; i--) {
	    	if (!isspace(*(field + i))) break;
	    }
	    
    	    if (i >= key_length) {
	    	field[i + 1] = '\0';
    	    } else {
	    	value = NULL;
	    }
	}
    }
    
    return value;
}


static long *index_reads(FILE *fp) {
    long *ri   = NULL;
    long rsize  = 16000;
    long count = 0;
    
    char *line = NULL;
    long size = 0;
    size_t line_size;
    int err = 0;
    int line_num = 0;
    long pos = 0;

    if (NULL == (ri = (long *)malloc(rsize * sizeof(long)))) {
    	fprintf(stderr, "index_reads initialisation out of memory\n");
	return NULL;
    }
    
    while (!err && tg_get_line(&line, &size, fp) > 0) {
    	char *value;
    	line_num++;
	
	if ((line_size = strlen(line))) {
	
	    if (strncmp(line, "{CTG", 4) == 0) break;
	    
	    if ((value = get_value("iid:", line))) {
	    	long val = atol(value);
		
		if (val != (count + 1)) {
		    err = 1;
		    continue;
		}
		
		if (count == rsize) {
		    
		    long *tmp;
		
		    rsize += rsize;
		    
		    if (NULL == (tmp = (long *)realloc(ri, rsize * sizeof(long)))) {
		    	fprintf(stderr, "index_reads out of memory at line %d\n", line_num);
			err = 1;
			continue;
		    }
		    
		    ri = tmp;
		}
		
		ri[count++] = pos;
	    }
	}
	    
    	pos = ftell(fp);
    }
    
    if (err || count == 0) {
    	fprintf(stderr, "Error at line %d: %s\n", line_num, line);
	
	if (line != NULL) free(line);
	free(ri);
	
	return NULL;
    }
    
    fprintf(stderr, "Number of reads %ld\n", count);
    
    return ri;
}


typedef struct {
    long src;
    long off;
    long start;
    long end;
    char *gaps;
} tle_t;


/* 
    should be inside the TLE entry
*/ 
static void read_tle(FILE *fp, tle_t *tle) {
    char line[255];
    
    while (fgets(line, 255, fp)) {
    	char *value;
    
    	if ((value = get_value("src:", line))) {
    	    tle->src = atol(value);
	} else if ((value = get_value("off:", line))) {
    	    tle->off = atol(value);
	} else if ((value = get_value("clr:", line))) {
	    tle->start = atol(value);
	    value = strchr(value, ',') + 1;  // FIX ME, error check
	    tle->end = atol(value);
	} else if (line[0] == '}') {
	    break;
	}
    }
}


static long read_in_tles(FILE *fp, tle_t **in_tle) {
    tle_t *tle;
    long tcount = 0;
    long tle_size = 2000;
    
    char *line = NULL;
    long size = 0;

    if (NULL == (tle = (tle_t*)malloc(tle_size * sizeof(tle_t)))) {
    	fprintf(stderr, "Out of memory in initial TLE acquisition.\n");
	return 0;
    }
    
            
    // find the first TLE line
    
    do {
    	tg_get_line(&line, &size, fp);
    } while (strncmp(line, "{TLE", 4) != 0 && line[0] != '}');
    
    while (strncmp(line, "{TLE", 4) == 0) {
    	tcount++;
	
	if (tcount == tle_size) {
	    tle_t *tmp;
	    
	    tle_size += tle_size;
	    
	    if (NULL == (tmp = (tle_t *)realloc(tle, tle_size * sizeof(tle_t)))) {
	    	fprintf(stderr, "Out of memory in TLE acquisition.\n");
		free(tle);
		return 0;
	    }
	    
	    tle = tmp;
	}
	
	read_tle(fp, &tle[tcount]);
	
	// should be on a '}' 
	tg_get_line(&line, &size, fp);	
    }

    if (line) free(line);
    
    *in_tle = tle;
    
    return tcount;
}



/* doesn't retrieve any data yet
   just skips over the SCF section */

static void read_in_scf(FILE *fp) {
    char *line = NULL;
    long size = 0;
    tle_t tle;

     // find the first TLE line
    
    do {
    	tg_get_line(&line, &size, fp);
    } while (strncmp(line, "{TLE", 4) != 0 && line[0] != '}');
    
    while (strncmp(line, "{TLE", 4) == 0) {
	
	read_tle(fp, &tle);
	
	// should be on a '}' 
	tg_get_line(&line, &size, fp);	
    }

    if (line) free(line);
    
    return;
}



/* 
    add a string of arbitrary length to another string.  
    Changes alloc size to 0 on a failure.
*/
static char *add_line(char *read, char *line, long line_size, long *read_size, long *alloc_size) {
	
    if (*alloc_size == 0) {
    	*alloc_size  = 16384;
	read = malloc(*alloc_size * sizeof(char));
	
	if (!read) {
	    *alloc_size = 0;
	    return read;
	}
    } 

    while ((line_size + *read_size) >= *alloc_size) {
    	char *tmp = NULL;
	*alloc_size += *alloc_size;
	tmp = realloc(read, *alloc_size * sizeof(char));
	
	if (tmp) {
	    read = tmp;
	} else {
	    *alloc_size = 0;
	    return read;
	}
    }

    strncpy(read + *read_size, line, line_size);
    *read_size += line_size;
    read[*read_size] = '\0';

    return read;
}


static long get_read_data(FILE *fp, tle_t *tle, long *ri, char **name, 
    	    	    	    char **read, char **qual, char **gaps) {

    char *line = NULL;
    long size = 0;
    
    long alloc_len;
    long offset; 
    long length; 
    
    long seq_len = 0;

    fseek(fp, ri[tle->src - 1], SEEK_SET);
    
    while (tg_get_line(&line, &size, fp) && strncmp(line, "}", 1) != 0) {
    	char *value;
	
	if ((value = get_value("eid:", line))) {
	    *name = strdup(value);
	} else if (strncmp(line, "seq:", 4) == 0) {
	    alloc_len = 0;
	    offset = 0;
	    length = 0;
	    
	    while (tg_get_line(&line, &size, fp) && strncmp(line, ".", 1) != 0) {
	    	length = chomp(line);
	    	*read = add_line(*read, line, length, &offset, &alloc_len);
		seq_len += length;
	    }
	} else if (strncmp(line, "qlt:", 4) == 0) {
	    alloc_len = 0;
	    offset = 0;
	    length = 0;

	    while (tg_get_line(&line, &size, fp) && strncmp(line, ".", 1) != 0) {
	    	length = chomp(line);
	    	*qual = add_line(*qual, line, length, &offset, &alloc_len);
	    }
	}
    }
    
    if (line) free(line);
    
    return seq_len;
}
		

static int store_contig(FILE *fp, char *fn, GapIO *io, tg_args *a, 
    	    	    	contig_t *c, tg_pair_t *pair, tle_t *tle, long *ri, long tc) {
    long store_pos;
    long i = 0;
    
    store_pos = ftell(fp); // we need to go back to this point at the end
    
    for (i = 1; i < tc; i++) {
    	seq_t seq;
	int dir;
	char *template_name = NULL;
	tg_rec recno;
	int flags, is_pair = 0;
	char *read_name = NULL;
	char *read = NULL;
	char *qual = NULL;
	char *gaps = NULL;
	
	memset(&seq, 0, sizeof(seq_t));
	
	dir = tle[i].start < tle[i].end ? 1 : -1;
	
    	// not sure that I have the alignment info
	// correct, run it and see what I have to change
	seq.pos   = tle[i].off;
	
	if (dir == 1) {
	    seq.left  = tle[i].start;
	    seq.right = tle[i].end;
	} else {
	    seq.left  = tle[i].end;
	    seq.right = tle[i].start;
	}
	
	seq.flags = dir < 0 ? SEQ_COMPLEMENTED : 0;
	seq.mapping_qual = 50; // doesn't appear to be set in AFG files
	
	// now we need the actual read
	seq.len = get_read_data(fp, &tle[i], ri, &read_name, &read, &qual, &gaps);
	
    	// FIX ME some pos manipulation may be required here
	
	seq.format = SEQ_FORMAT_CNF1;
	
	if (a->data_type & DATA_SEQ) {
	    int i;
	    
	    for (i = 0; i < seq.len; i++) {
	    	if (read[i] == '-') {
		    read[i] = '*';
		} else if (read[i] == 'n' || read[i] == 'N') {
		    read[i] = '-';
		}
	    }
	} else {
	    memset(read, 'N', seq.len);
	}
	
	if (!(a->data_type & DATA_QUAL)) {
	    memset(qual, 0, seq.len);
	} else {
	    // convert to internal form
	    int i;
	    
	    for (i = 0; i < seq.len; i++) {
	    	qual[i] = qual[i] - 32;
	    }
	}
	
	seq.name_len       = strlen(read_name);
	seq.trace_name_len = 0;
	seq.alignment_len  = 0;
	
	template_name = read_name; // not strictly necessary but keeping for clarity
	
	seq.name = (char *)calloc(seq.name_len + 5 + 2 * seq.len, sizeof(char));

	strcpy(seq.name, read_name);
	seq.trace_name = seq.name + seq.name_len + 1;
	seq.alignment = seq.trace_name + seq.trace_name_len + 1;
	seq.seq = seq.alignment + seq.alignment_len + 1;
	seq.alignment = 0;

	memcpy(seq.seq, read, seq.len);

	seq.conf = (int8_t *) seq.seq + seq.len;
	memcpy(seq.conf, qual, seq.len);

	// seq_t struct filled in, now to save it
	
	seq.len = seq.len * dir;
	
	flags = GRANGE_FLAG_TYPE_SINGLE;

	if (seq.flags & SEQ_END_REV)
	    flags |= GRANGE_FLAG_END_REV;
	else
	    flags |= GRANGE_FLAG_END_FWD;
	if (seq.len < 0)
	    flags |= GRANGE_FLAG_COMP1;

	if (pair) is_pair = 1;
	
	recno = save_range_sequence(io, &seq, seq.mapping_qual, pair,
				    is_pair, template_name, c, a, flags, NULL,
				    NULL);

	
	if (read) free(read);
	
	if (qual) free(qual);
	
	if (read_name) free(read_name);
	
	if (gaps) free(gaps);
    
    
    
    }
    
    fseek(fp, store_pos, SEEK_SET);
    
    return 0;
}
    

static int convert_afg(FILE *fp, char *fn, GapIO *io, tg_args *a, 
    	    	    	tg_pair_t *pair, long *ri) {
    contig_t *c = NULL;
    tle_t  *tle = NULL;
    
    char *line = NULL;
    long size = 0;
    size_t line_size;
    int err = 0;
    int line_num = 0;
    
    while (!err && tg_get_line(&line, &size, fp) > 0) {
    	line_num++;
	long tc;
	
	if ((line_size = strlen(line))) {
	    // we are at the beginning of a new contig
	    char *value;
	    
	    // miss out the next line (iid)
	    if (tg_get_line(&line, &size, fp) == -1) {
	    	err = 1;
		continue;
	    }
	    
	    line_num++;
	    
	    if (NULL == (value = get_value("eid:", line))) {
	    	fprintf(stderr, "Expecting eid: in line %d:'%s'\n", line_num, line);
		err = 1;
		continue;
	    }
	    
	    // now we have the name we can create a new contig
	    create_new_contig(io, &c, value, a->merge_contigs);
	    fprintf(stderr, "Storing contig %s\n", value);
	    
	    
	    // now to load in the TLEs (alignment info)
	    tc = read_in_tles(fp, &tle);
	    
	    /* ok, now we should be ready to get the sequences
	       and finish off the contig */
	       
	    err = store_contig(fp, fn, io, a, c, pair, tle, ri, tc);
	    
	    tg_get_line(&line, &size, fp);
	    
	    if (strncmp(line, "{SCF", 4) == 0) {
	    	read_in_scf(fp); // doesn't actually return any data
		tg_get_line(&line, &size, fp);
	    }

	    if (tle) free(tle);
	}
    }
    
    return 0;
}
	    
	    
    
   
   

int parse_afg(GapIO *io, char *fn, tg_args *a) {
    FILE *fp;
    long *read_index;
    tg_pair_t *pair = NULL;
    
    
    if (NULL == (fp = fopen(fn, "r"))) {
    	fprintf(stderr, "Can't open %s\n", fn);
	return -1;
    }
    
    if (a->pair_reads) {
    	pair = create_pair(a->pair_queue);
    }
    
    if (NULL == (read_index = index_reads(fp))) {
    	fprintf(stderr, "Unable to index reads\n");
	return -1;
    }
    
    convert_afg(fp, fn, io, a, pair, read_index);        
    
    if (pair && !a->fast_mode) {    
	finish_pairs(io, pair, a->link_pairs);
    }
    
    cache_flush(io);
    fclose(fp);
 
    if (pair) delete_pair(pair);
   
    
    return 0;
    
}
