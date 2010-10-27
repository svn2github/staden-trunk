/*
 * caf.c - caf to gap5 conversion for tg_index
 *
 * Andrew Whitwham, August 2010
 * Wellcome Trust Sanger Institute
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
#include "zfio.h"
#include "string_alloc.h"
#include "caf.h"


/*
    This caf conversion works by first indexing a caf on Sequence, DNA and BaseQuality
    then loading the data into gap5 one contig at a time.
*/ 



// index entry for file positions
typedef struct {
    char *name;
    int length;
    long pos;
    char *data;
} caf_entry;

// index list
typedef struct {
    caf_entry *entry;
    long size;
    long alloc;
} caf_index;

// store for dna and quality data
typedef struct {
    char *seq;  // dna
    long s_len; // dna length;
    char *qual; // qual as bytes
} seq_qual;


// store for annotations
typedef struct {
    int  type;
    int  start;
    int  end;
    char *text;
} anno_type;

// types for node storage
typedef struct {
    long seq;
    long dna;
    long qua;
} index_data;

typedef struct caf_node_s {
    char *prefix;
    struct caf_node_s *children;
    struct caf_node_s *sibling;
    index_data *data;
} caf_node;

typedef struct {
    pool_alloc_t   *data_pool; 
    pool_alloc_t   *node_pool;
    string_alloc_t *str_pool;
} pools;

typedef enum {IS_CONTIG, IS_SEQUENCE, IS_DNA, IS_QUALITY} entry_type;
    

// Next few functions are for string and file reading

static long replace(char *str, char rep) {
    long length = 0;
	
    if (str && (length = strlen(str))) {
    	if (str[length - 1] == '\n') {
	    str[length - 1] = rep;
	}
    }

    return length;
}


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
	
    	    if (i != key_length) {
	    	field[i + 1] = '\0';
    	    } else {
	    	value = NULL;
	    }
	}
    }
    
    return value;
}


/* 
    add a string of arbitrary length to another string.  
    Changes alloc size to 0 on a failier.
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

    if ((line_size + *read_size) >= *alloc_size) {
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


/*
    a more specialised version of get_value, could
    be replaced.
*/
static int find_name(char *in_line, char **name) {
    char *pos = NULL;
    int i = 0;

    if (in_line && (pos = strchr(in_line, ':')) != NULL) {
    
    	do {
	    pos++;
	} while (*pos && isspace(*pos));
	
	while (*(pos + i) && !isspace(*(pos + i))) {
	    i++;
	}
	
	*(pos + i) = '\0';
    }
    
    *name = pos;
    
    return i;
}


/* an alternative to the Gnu specific getline */
static long get_line(char **line, long *length, zfp *fp) {
    char *in_line;
    long len;
    long offset = 0;
    long in_chars;
    
    if (line == NULL || fp == NULL || length == NULL) {
    	return -1;
    }
    
    if (*line == NULL || *length <= 0) {
    	if ((*line = malloc(256 * sizeof(char))) == NULL) {
	    return -1;
	}
	
	*length = 256;
    }
    
    in_line = *line;
    len     = *length;
    
    while ((zfgets((in_line + offset), (len - offset), fp))) {
    	char *tmp = NULL;
	
    	in_chars = strlen(in_line);
	
    	offset = in_chars;
	
	// see if we have our full line
	if (*(in_line + offset - 1) == '\n') {
	    break;
	}
    	
    	len += len;
	tmp = realloc(in_line, len * sizeof(char));
	    
    	if (tmp) {
	    in_line = tmp;
	} else {
	    fprintf(stderr, "Memory error in get_line\n");
	    return -1;
	}
    }
    
    *line = in_line;
    *length = len;
    
    return offset;
}


// following fuctions are part of the Patricia Tree/Trie implementation

/*
    test function to see what is in the entire tree
*/
void tree_walk(caf_node *caf, int depth) {
    
    fprintf(stderr, "<%d> ", depth);
    
    if (caf->prefix) {
    	fprintf(stderr, "%s ", caf->prefix);
    }
    
    if (caf->children) {
    	caf_node *sib = caf->children;
	depth++;
	
	while (sib) {
	    tree_walk(sib, depth);
	    sib = sib->sibling;
	}
    }
    
    fprintf(stderr, "|\n");
    
    return;
}


/*
    match the incoming name to the node prefix
*/
static int match_prefix(char *name, int name_pos, char *prefix) {
    int pos = 0;

    if (prefix) {
    	while ((name[name_pos + pos] != '\0') && (prefix[pos] != '\0') && (name[name_pos + pos] == prefix[pos])) {
	    pos++;
	}
    }
    
    return pos;
}


/*
    add a new node to the tree and return a pointer to the new node
*/
static caf_node *add_child_node(pools *pool, caf_node *caf, char *name) {
    caf_node *node = NULL;
    caf_node *sib;
    
    // first create new node
        
    node = pool_alloc(pool->node_pool);
    
    if (NULL == node) return NULL;
        
    node->prefix   = string_dup(pool->str_pool, name);
    node->children = NULL;
    node->data     = NULL;
    node->sibling  = NULL;
   
    if (caf->children) {
    	sib = caf->children;
	
	while (sib->sibling) {
	    sib = sib->sibling;
	}
	
	sib->sibling = node;

    } else {
    	caf->children = node;
    }
    
    return node;
}


/*
    add a caf node to the right place in the tree.  The most
    complex function in the patricia tree.
*/
static int add_caf_node(pools *pool, caf_node **caf, char *name, int name_pos) {
    caf_node *node = *caf;
    
    if (!(*caf)->prefix && !(*caf)->children) {
    	// this is a special case, we're at the base node
	node = add_child_node(pool, *caf, name);
    } else {
    	// see how much of the prefix matches
    	int pos = match_prefix(name, name_pos, (*caf)->prefix);
	
	if ((*caf)->prefix && ((*caf)->prefix[pos] != '\0')) {
	    if (name[name_pos + pos] == '\0') {
	    	/* name is shorter than prefix, so create a new node
		   for rest of prefix and add name here */
		   
		caf_node *children = (*caf)->children;   
		(*caf)->children = NULL; 
		  
		node = add_child_node(pool, *caf, &((*caf)->prefix[pos]));
		
		(*caf)->prefix = string_dup(pool->str_pool, name + name_pos);		
		node->data = (*caf)->data;
		(*caf)->data = NULL;
		node->children = children;
		
		// make node the same as caf for return value
		node = *caf;		
	    } else {
		caf_node *children = (*caf)->children;
		index_data *data = (*caf)->data;

		(*caf)->children = NULL;
		(*caf)->data = NULL;

		/* node needs to be split
		   pos should give us the differing character, so create a left
		   node holding data from the current node and a new right node */
		add_child_node(pool, *caf, (*caf)->prefix + pos);
		node = add_child_node(pool, *caf, name + name_pos + pos);

		// shorten the prefix to the right length, maybe deallocate memory later
		(*caf)->prefix[pos] = '\0';

		// copy the data to the new node left most node
		(*caf)->children->data = data;
		(*caf)->children->children = children;
	    }
	} else {
	    // just need to add a new child
	    node = add_child_node(pool, *caf, name + name_pos + pos);
	}
    }
    
    *caf = node;
    
    return 1; // need to put in error checks
}


/*
    find the the correct node, optionally create one if it
    does not exist.
*/
static int find_caf_node(pools *pool, caf_node **caf, char *name, int add) {
    int found = 0;
    int searching = 1;
    int name_len = strlen(name);
    int name_pos = 0;
    caf_node *this_node = *caf;
    
    while (searching) {
    	// if no children we're at a leaf node
	if (!this_node->children) {
	    // see if the rest of the name matches
	    if (this_node->prefix) {
	    	if (strcmp(name + name_pos, this_node->prefix) == 0) {
	    	    found = 1;
		}
	    }
	    
	    searching = 0;
	} else {
	    int check_children = 1;
	
	    // check prefix first
	    if (this_node->prefix) {
	    	if (strncmp(name + name_pos, this_node->prefix, strlen(this_node->prefix)) != 0) {
		    // doesn't match as far as it goes
		    searching = 0;
		    check_children = 0;
		} else if (strcmp(name + name_pos, this_node->prefix) == 0) {
		    // an exact match, we've found our node
		    searching = 0;
		    check_children = 0;
		    found = 1;
		}
	    }
	    
	    if (check_children) {
		int pos;

		pos = match_prefix(name, name_pos, this_node->prefix);
		searching = 0;

		if ((name_pos) < name_len) {
	    	    // search the children
		    caf_node *child;

    	    	    child = this_node->children;
		    
		    while (child) {
			if (name[name_pos + pos] == child->prefix[0]) {
		    	    this_node = child;
		    	    name_pos += pos; // found so move name along
			    searching = 1;
			    break;
			}
			
			child = child->sibling;
		    }
		}
	    }
	}
    }
    
    if (!found && add) {
    	if (add_caf_node(pool, &this_node, name, name_pos)) {
	    found = 1;
	}
    }
    
    *caf = this_node;
    
    return found;
}


/*
    add position data to the appropriate node, create a node if
    it does not exist.
*/
static int add_pos_data(pools *pool, caf_node *caf, char *name, long pos, entry_type type) {
    int ok = 1;
    caf_node *node = caf;

    if (!find_caf_node(pool, &node, name, 1)) {
    	ok = 0;
    }

    if (ok) {
    	if (!node->data) {
	    if ((node->data = pool_alloc(pool->data_pool))) {
	    	node->data->seq = -1;
		node->data->dna = -1;
		node->data->qua = -1;
	    } else {
	    	fprintf(stderr, "Unable to alloc data memory in add_pos_data for %s\n", name);
		ok = 0;
	    }
	}
	
	if (ok) {
	    switch (type) {
	    	case IS_SEQUENCE:
		    node->data->seq = pos;
		    break;
		    
		case IS_DNA:
		   node->data->dna = pos;
		   break;
		   
		case IS_QUALITY:
		   node->data->qua = pos;
		   break;
		   
		default:
		    fprintf(stderr, "add_pos_data wrong type %d for %s\n", type, name);
	    }
	}
    }
    
    return ok ? 0 : 1; 
}


// following functions are for index manipulation

/*
    find an entry in the index and return its position
*/
static long find_entry(caf_index *index, char *name, long start) {
    long found = -1;

    if (start > index->size) start = 0;
    
    if (strcmp(name, index->entry[start].name) == 0) {
    	found = start;
    } else {
    	long mid = index->size / 2;
	long begin = 0;
	long end = (index->size - 1);
	
	// a basic binary search, seems good enough
	while (found == -1 && begin <= end) {
	    int match = strcmp(name, index->entry[mid].name);
	    
	    if (match > 0) {
	    	begin = mid + 1;
		mid = begin + (end - begin) / 2;
	    } else if (match < 0) {
	    	end = mid - 1;
		mid = begin + (end - begin) / 2;
	    } else {
	    	found = mid;
	    }
	}
    }
    
    return found;
}


/*
    add an entry to the index, name must be unique
*/
static int add_entry(caf_index *index, char *name, int length, long pos) {

    // check for existing entry 
    if (index->size) {
    	long found = find_entry(index, name, index->size - 1);
	
	if (found != -1) return 0;
    }    

    if (index->size == index->alloc) {
    	caf_entry *tmp;
		
    	index->alloc += index->alloc;
	tmp = realloc(index->entry, index->alloc * sizeof(caf_entry));
	
	if (tmp) {
	    index->entry = tmp;
	} else {
	    // unable to allocate memory
	    return 1;
	}
    }
    
    if (!(index->entry[index->size].name = calloc(length + 1, sizeof(char)))) {
    	return 1;
    }
    
    strncpy(index->entry[index->size].name, name, length);
    index->entry[index->size].pos = pos;
    index->entry[index->size].length = length;
    index->entry[index->size].data = NULL;
    index->size++;
    
    return 0;
}


/*
    add extra data to the index entry, entry must exist
*/
static int add_entry_data(caf_index *index, char *data, int size, long pos) {

    if (pos >= index->size) return 1;
    
    if (data && size) {
    	index->entry[pos].data = calloc(size + 1, sizeof(char));
	
	if (index->entry[pos].data) {
	    strncpy(index->entry[pos].data, data, size);
	} else {
	    return 1;
	}
    }
    
    return 0;
}
    

static int initialise_index(caf_index *index) {
    index->size = 0;
    index->alloc = 1024;
    
    if ((index->entry = malloc(index->alloc * sizeof(caf_entry)))) {
    	return 0;
    }
    
    return 1;
}


static void clear_index(caf_index *index) {
    int i;

    for (i = 0; i < index->size; i++) {
    	free(index->entry[i].name);
	
	if (index->entry[i].data) {
	    free(index->entry[i].data);
	}
    }

    free(index->entry);
    index->size = 0;
    index->alloc = 0;
}


/*
    index a caf file on Sequence Reads, Sequence Contigs, DNA and BaseQuality
*/
static int index_caf(zfp *fp, pools *pool, caf_node *caf, caf_index *contig_entry) {
		
    char *line = NULL;
    long size = 0;
    size_t line_size;
    int err = 0;
    int line_num = 0;
    long pos = 0;
    
    int read_no   = 0;
    int contig_no = 0;
    int dna_no    = 0;
    int qual_no   = 0;
    
    while (!err && get_line(&line, &size, fp) > 0) {
	line_num++;
	char *name = NULL;
	
	if ((line_size = strlen(line))) {
	    if ((name = get_value("Sequence :", line))) {
		if (name) {
		    if (get_line(&line, &size, fp) == -1) {
		    	err = 1;
			continue;
		    }
		    
		    if (strncmp(line, "Is_read", 7) == 0) {
		    	err = add_pos_data(pool, caf, name, pos, IS_SEQUENCE);
			read_no++;
		    } else if (strncmp(line, "Is_contig", 9) == 0) {
		    	err = add_pos_data(pool, caf, name, pos, IS_SEQUENCE);
		    	err = add_entry(contig_entry, name, strlen(name), pos);
			contig_no++;
		    }
		}
	    } else if ((name = get_value("DNA :", line))) {
    	        err = add_pos_data(pool, caf, name, pos, IS_DNA);
		dna_no++;
	    } else if ((name = get_value("BaseQuality :", line))) {
    	        err = add_pos_data(pool, caf, name, pos, IS_QUALITY);
		qual_no++;
	    } else if (strncmp(line, "Unpadded", 8) == 0) {
	    	fprintf(stderr, "Error, padded data only\n");
		err = 1;
	    }
	    	    
	}

    	pos = zftello(fp);
    }
    
    if (err) {
    	fprintf(stderr, "Error at line %d: %s", line_num, line);
    }
    
    fprintf(stderr, "Input summary\nReads %d\nQuality %d\nContigs %d\nBases %d\nR + C %d\n", read_no, qual_no, contig_no, dna_no, (read_no + contig_no));
    
    if (line != NULL) {
    	free(line);
    }
    
    return err;
}


static int cmp_string(const void *p1, const void *p2) {
    caf_entry *c1 = (caf_entry *)p1;
    caf_entry *c2 = (caf_entry *)p2;
    
    return strcmp(c1->name, c2->name);
} 


static void sort_index(caf_index *index) {
    qsort(index->entry, index->size, sizeof(caf_entry), cmp_string);
}


/* The remaining functions load the indexed data into gap5 */


static int parse_annotation(anno_type **annotation, int *anno_count, int *anno_size, char *value) {
    char *anno_entry = NULL;
    size_t txt_len;

    if (*anno_size == 0) {
	*annotation = malloc(8 * sizeof(anno_type));
	*anno_size = 8;
    } else if (*anno_count == *anno_size) {
	*anno_size += *anno_size;
	*annotation = realloc(*annotation, *anno_size * sizeof(anno_type));
    }

    if (!*annotation) {
	return 1;
    }
    
    (*annotation)[*anno_count].text = NULL;

    // get the anno type
    anno_entry = strtok(value, " ");
    (*annotation)[*anno_count].type = str2type(anno_entry);

    // get the start
    anno_entry = strtok(NULL, " ");
    (*annotation)[*anno_count].start = atoi(anno_entry);

    // get the end
    anno_entry = strtok(NULL, " ");
    (*annotation)[*anno_count].end = atoi(anno_entry);

    // get the text (it's quoted with ")
    anno_entry = strtok(NULL, "\"");
    
    if (anno_entry) {
       txt_len = strlen(anno_entry);

       if (((*annotation)[*anno_count].text = calloc(txt_len + 1, sizeof(char)))) {
	   strncpy((*annotation)[*anno_count].text, anno_entry, txt_len);
       } else {
	   return 1;
       }
    }

    (*anno_count)++;

    return 0;
}


static int read_contig_section(zfp *fp, GapIO *io, contig_t **contig, tg_args *a, long pos, caf_index *reads) {
    char *line = NULL;
    long size = 0;
    anno_type *annotation = NULL;
    int anno_count  = 0;
    int anno_size = 0;
    int i;
    char *value;
    
    zfseeko(fp, pos, SEEK_SET);
 
    while (get_line(&line, &size, fp) > 0) {
    	if (isspace(line[0])) break; // blank line at end of section
	
    	// grab the read name and data from the assembly part of the section
	if (strncmp(line, "Assembled_from", 14) == 0) {
	    char *start = strchr(line, ' ');
	    char *end;
	    start++;
	    end = strchr(start, ' ');
	    *end = '\0';
	    add_entry(reads, start, (end - start), 0);
	    end++;
	    if (add_entry_data(reads, end, strlen(end), (reads->size - 1))) {
	    	fprintf(stderr, "Unable to add contig data, out of memory\n");
		return 1;
	    }
	    
	} else if ((value = get_value("Tag", line))) {
	    if (a->data_type & DATA_ANNO) {
		if (parse_annotation(&annotation, &anno_count, &anno_size, value)) {
		    fprintf(stderr, "Out of memory, unable to process tags\n");
		    return 1;
		}
	    }
	}
    }
    
    
    // handle annotations (tags), if we have them 
    
    for (i = 0; i < anno_count; i++) {
	range_t r;
	anno_ele_t *e;
	int an_pos;
	int an_len;
	bin_index_t *bin;

	if (annotation[i].end >= annotation[i].start) {
	    an_len = annotation[i].end - annotation[i].start + 1;
	} else {
	    an_len = annotation[i].start - annotation[i].end + 1;
	    annotation[i].start = annotation[i].end;
	}

	an_pos = ((*contig)->start + 1) + (annotation[i].start - 1);

	if (annotation[i].text) {
	    unescape_line(annotation[i].text);
	}

	r.start    = an_pos;
	r.end      = an_pos + (an_len ? an_len -1 : 0);
	r.mqual    = annotation[i].type;
	r.pair_rec = (*contig)->rec;
	r.flags    = GRANGE_FLAG_ISANNO;

	r.rec = anno_ele_new(io, 0, GT_Contig, (*contig)->rec, 0, annotation[i].type,
	    	    annotation[i].text);

	e = (anno_ele_t *)cache_search(io, GT_AnnoEle, r.rec);
	e = cache_rw(io, e);

	bin = bin_add_range(io, contig, &r, NULL, NULL, 0);
	e->bin = bin->rec;

	if (annotation[i].text) free(annotation[i].text);
    }

    if (annotation) free(annotation);
    
    if (line != NULL) {
    	free(line);
    }
    
    return 0;
}


/*
    DNA and BaseQuality data can be split over several lines, join them together
    in one long string.
*/
static long read_section_as_line(zfp *fp, long pos, char **line, int is_seq) {
    char *line_in = NULL;
    long size  = 0;
    long length = 0;
    char *seq_line = NULL;
    long alloc_len = 0;
    long offset = 0;


    zfseeko(fp, pos, SEEK_SET);
    
    get_line(&line_in, &size, fp); // skip over the header

    while ((length = get_line(&line_in, &size, fp)) > 0) {
    	if (isspace(line_in[0])) break; // blank line at end of section
	
	if (is_seq) {
	    length = chomp(line_in);
	} else {
	    length = replace(line_in, ' ');
	}
	
	seq_line = add_line(seq_line, line_in, length, &offset, &alloc_len);
	
	if (alloc_len == 0) {
	    fprintf(stderr, "Out of memory while reading data\n");
	    return -1;
	}
    }
    
    *line = seq_line; // this doesn't get freed till later
    
    if (line_in) {
    	free(line_in);
    }
    
    return offset;
}


/*
    read in sequence bases as a single line
*/
static long sequence_data(zfp *fp, seq_qual *seq, long pos) {
    long length = 0;
    char *seq_line = NULL;
    
    if (pos != -1) {
	if ((length = read_section_as_line(fp, pos, &seq_line, 1)) > 0) {
	    seq->seq   = seq_line;
	    seq->s_len = length;
	} else {
	    pos = -1;
	}
    }
    
    return pos;
}


/*
    read in quality as a single line and convert to be held as a single char per value
*/
static int quality_data(zfp *fp, seq_qual *qual, long pos) {
    char *qual_line = NULL;
    long length = 0;
    int err = 1;

    if (pos != -1) {
    	if ((length = read_section_as_line(fp, pos, &qual_line, 0)) > 0) {
	    char *value = NULL;
	    int i = 0;
	    
	    if (!(qual->qual = calloc(qual->s_len + 1, sizeof(char)))) {
	    	return 1;
	    }
	    
	    value = strtok(qual_line, " ");
	    
	    while (value) {
	    	qual->qual[i++] = (char) atoi(value);
		value = strtok(NULL, " ");
	    }
	    
	    err = 0;
	    
	    if (qual_line) {
	    	free(qual_line);
	    }
	}
    }
    
    return err;
}


/*
    read in sequence data, bases and quality plus any annotations and enter them into
    the gap5 db
*/
static int read_data(zfp *fp, GapIO *io, tg_args *a, contig_t **c, HacheTable *pair,
    	    	    	pools *pool, caf_node *caf, char *name, char *align) {
			
    long pos;
    seq_qual sq;
    caf_node *this_node = caf;
    
    if (!find_caf_node(pool, &this_node, name, 0)) {
    	fprintf(stderr, "Can't find entry on %s\n", name);
	return 1;
    }
			
    // get the sequence and quality data
    if ((pos = sequence_data(fp, &sq, this_node->data->dna)) == -1 ||
    	       quality_data(fp, &sq, this_node->data->qua)) {
	       
    	fprintf(stderr, "Failed to get sequnce and quality data\n");
	return 1;
    }
    
    pos = this_node->data->seq;
    
    if (pos != -1) {
	char *line_in = NULL;
	long size  = 0;
	seq_t seq;
	int cstart, cend, rstart, rend;
	int dir;
	char *trace_name = NULL;
	int tr_len = 0;
	char *template_name = NULL;
	int tm_len = 0;
	tg_rec recno;
	int flags, is_pair = 0;
	anno_type *annotation = NULL;
	int anno_count  = 0;
	int anno_size = 0;
	int i;
	int need_free = 0;
	
	memset(&seq, 0, sizeof(seq_t));
	
    	// we already have some of the data
	if (sscanf(align, "%d %d %d %d", &cstart, &cend, &rstart, &rend) != 4) {
	    fprintf(stderr, "Malformed alignment %s in %s\n", align, name);
	    return 1;
	}
	
	dir = cstart < cend ? 1 : -1;
	
	if (dir == -1) {
	    int tmp = cstart;
	    cstart = cend;
	    cend = tmp;
	}
	
	seq.pos   = cstart;
	seq.left  = rstart;
	seq.right = rend;
	
	seq.flags = dir < 0 ? SEQ_COMPLEMENTED : 0;
	
	seq.mapping_qual = 50; // not set in CAF, default to this
	
	seq.len = sq.s_len * dir;
	
	if (dir == 1) {
	    seq.pos -= seq.left - 1;
	} else {
	    seq.pos -= -seq.len - seq.right;
	}
	
	seq.format = SEQ_FORMAT_CNF1; // pack bytes
	
	if (a->data_type & DATA_SEQ) {
	    int i;
	    
	    for (i = 0; i < sq.s_len; i++) {
	    	if (sq.seq[i] == '-') {
		    sq.seq[i] = '*';
		} else if (sq.seq[i] == 'n' || sq.seq[i] == 'N') {
		    sq.seq[i] = '-';
		}
	    }
	} else {
	    memset(sq.seq, 'N', sq.s_len);
	}
	
	if (!(a->data_type & DATA_QUAL)) {
	    memset(sq.qual, 0, sq.s_len);
	}

	zfseeko(fp, pos, SEEK_SET);
	
    	get_line(&line_in, &size, fp); // skip over the header
	
	while (get_line(&line_in, &size, fp) > 0) {
	    char *value;
	    
	    if (isspace(line_in[0])) break; // blank line at end of section
	    
	    if ((value = get_value("SCF_File", line_in))) {
	    	tr_len = strlen(value);
		
		if(!(trace_name = calloc(tr_len + 1, sizeof(char)))) {
		    return 1;
		}
		
		strncpy(trace_name, value, tr_len);
		
	    } else if ((value = get_value("Template", line_in))) {
	    	tm_len = strlen(value);
		
		if(!(template_name = calloc(tm_len + 1, sizeof(char)))) {
		    return 1;
		}
		
		need_free = 1;
		
		strncpy(template_name, value, tm_len);

	    } else if ((value = get_value("Strand", line_in))) {
	    	if (strncmp(value, "Reverse", 7) == 0) {
		    seq.flags |= SEQ_END_REV;
		}
	    } else if ((value = get_value("Tag", line_in))) {
	    	if (a->data_type & DATA_ANNO) {
		    if (parse_annotation(&annotation, &anno_count, &anno_size, value)) {
		    	fprintf(stderr, "Out of memory, unable to process tags\n");
		     	return 1;
		    }
		}
	    }
	}
	
	if (line_in) {
	    free(line_in);
	}
	
	// we now should have all we need
	
	seq.name_len       = strlen(name);
	seq.trace_name_len = tr_len;
	seq.alignment_len  = 0;
	
	if (template_name == NULL) {
	    template_name = name;
	}
	
	// seq.data freed by save_range_sequence
	seq.name = seq.data = (char *) calloc(seq.name_len + 1 +
	    	    	    	    	      seq.trace_name_len + 1 +
					      1 +
					      2 * sq.s_len, sizeof(char));
					      
	strcpy(seq.name, name);

	seq.trace_name = seq.name + seq.name_len + 1;
	
	if (trace_name) {
	    strcpy(seq.trace_name, trace_name);
	}

	seq.alignment = seq.trace_name + seq.trace_name_len + 1;
	seq.seq = seq.alignment + seq.alignment_len + 1;
	seq.alignment = 0;

	memcpy(seq.seq, sq.seq, sq.s_len);

	seq.conf = seq.seq + sq.s_len;
	memcpy(seq.conf, sq.qual, sq.s_len);
	
	// seq_t struct filled in, now to save it
	
	flags = GRANGE_FLAG_TYPE_SINGLE;

	if (seq.flags & SEQ_END_REV)
	    flags |= GRANGE_FLAG_END_REV;
	else
	    flags |= GRANGE_FLAG_END_FWD;
	if (seq.len < 0)
	    flags |= GRANGE_FLAG_COMP1;

	if (pair) is_pair = 1;
	
	recno = save_range_sequence(io, &seq, seq.mapping_qual, pair,
					is_pair, template_name, *c, a, flags, NULL);
					
	if (trace_name) {
	    free(trace_name);
	}
	
	if (need_free) {    // if template_name is not using the same memory as name
	    free(template_name);
	}
					
	// now we have a recno, do the annotations (tags)
	
	for (i = 0; i < anno_count; i++) {
	    range_t r;
	    anno_ele_t *e;
	    int an_pos;
	    int an_len;
	    bin_index_t *bin;
	    
	    if (annotation[i].end >= annotation[i].start) {
	    	an_len = annotation[i].end - annotation[i].start + 1;
	    } else {
	    	an_len = annotation[i].start - annotation[i].end + 1;
		annotation[i].start = annotation[i].end;
	    }
	    
	    if (seq.len >= 0) {
	    	an_pos = seq.pos + annotation[i].start - 1;
	    } else {
	    	int tmp_pos = seq.pos - seq.len - 1;
		
		an_pos = tmp_pos - (annotation[i].start - 1) -
		    	    (an_len ? an_len -1 : 0);
	    }
	    
	    if (annotation[i].text) {
	    	unescape_line(annotation[i].text);
	    }
	    
	    r.start    = an_pos;
	    r.end      = an_pos + (an_len ? an_len -1 : 0);
	    r.mqual    = annotation[i].type;
	    r.pair_rec = recno;
	    r.flags    = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	    
	    r.rec = anno_ele_new(io, 0, GT_Seq, recno, 0, annotation[i].type,
	    	    	annotation[i].text);
			
	    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, r.rec);
	    e = cache_rw(io, e);
	    
	    bin = bin_add_range(io, c, &r, NULL, NULL, 0);
	    e->bin = bin->rec;
	    
	    if (annotation[i].text) free(annotation[i].text);
	}
	
	if (annotation) free(annotation);
	
    } else {
    	return 1;
    }
    
    if (sq.seq) {
    	free(sq.seq);
    }
    
    if (sq.qual) {
    	free(sq.qual);
    }
	
    return 0;
}
	

static int import_caf(zfp *fp, GapIO *io, tg_args *a, HacheTable *pair, pools *pool, caf_node *caf, caf_index *contig_entry) {
    contig_t *c = NULL;
    long read_no       = 0;
    long contig_no     = 0;
    long i;
    caf_index contig_reads;
    
    for (i = 0; i < contig_entry->size; i++) {
    	long j;
	int err;
	
	initialise_index(&contig_reads);

	create_new_contig(io, &c, contig_entry->entry[i].name, a->merge_contigs);
	
	if (read_contig_section(fp, io, &c, a, contig_entry->entry[i].pos, &contig_reads)) {
	    return 1;
	}
	
	contig_no++;
	
	for (j = 0; j < contig_reads.size; j++) {
	    err = read_data(fp, io, a, &c, pair, pool, caf, contig_reads.entry[j].name, contig_reads.entry[j].data);

	    if (err) {
	    	fprintf(stderr, "Unable to import data for read %s on contig %s",
		    contig_reads.entry[j].name, contig_entry->entry[i].name);
	    } else {
	    	read_no++;
	    }
	}
	
	clear_index(&contig_reads);
	
	if (((read_no + contig_no) & 0x3fff) == 0) {
	    cache_flush(io);
	} 
    }
    
    fprintf(stderr, "Done %ld reads in %ld contigs.\n", read_no, contig_no);
    
    return 0;
}
     


int parse_caf(GapIO *io, char *fn, tg_args *a) {
    zfp *fp;
    struct stat sb;
    HacheTable *pair = NULL;
    int err = 0;
    caf_index contig_entry;
    caf_node caf;
    pools pool;
    
    // some variables for stats
    struct timeval tv1, tv2;
    double dt;

    printf("Loading %s...\n", fn);

    if (-1 == stat(fn, &sb) ||
	NULL == (fp = zfopen(fn, "r"))) {
	perror(fn);
	return -1;
    }
    
    /* for pair data */
    open_tmp_file();

    if (a->pair_reads) {
	pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
	pair->name = "pair";
    }
    
    // some initialisations 
    pool.data_pool = pool_create(sizeof(index_data));
    pool.node_pool = pool_create(sizeof(caf_node));
    pool.str_pool  = string_pool_create(0);
    
    caf.prefix   = NULL;
    caf.children = NULL;
    caf.data     = NULL;
    caf.sibling  = NULL;

    err += initialise_index(&contig_entry);

    if (err) {
    	fprintf(stderr, "Unable to initialise index, out of memory\n");
	return -1;
    }

    gettimeofday(&tv1, NULL);
    
    if (index_caf(fp, &pool, &caf, &contig_entry)) {
    	fprintf(stderr, "Unable to populate index, out of memory\n");
	return -1;
    }
    
    gettimeofday(&tv2, NULL);
    dt = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;
    
    fprintf(stderr, "Indexed in %f seconds\n", dt);
    
    gettimeofday(&tv1, NULL);
    
    sort_index(&contig_entry);
    
    gettimeofday(&tv2, NULL);
    dt = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;
    
    fprintf(stderr, "Sorted in %f seconds\n", dt);

    gettimeofday(&tv1, NULL);
    
    fprintf(stderr, "Loading ...\n");
    
    if (import_caf(fp, io, a, pair, &pool, &caf, &contig_entry)) {
    	fprintf(stderr, "Unable to populate gap5 db");
	return -1;
    }
    

    if (pair && !a->fast_mode) {    
	sort_pair_file();
	
	complete_pairs(io);
	
	close_tmp_file();
    }
    
    cache_flush(io);
    zfclose(fp);

    if (pair) HacheTableDestroy(pair, 1);

    clear_index(&contig_entry);

    pool_destroy(pool.data_pool);
    pool_destroy(pool.node_pool);
    string_pool_destroy(pool.str_pool);

    gettimeofday(&tv2, NULL);
    dt = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;
    
    fprintf(stderr, "Loaded in %f seconds\n", dt);

    return 0;
}
   
