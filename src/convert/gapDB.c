#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "IO.h"
#include "process.h"
#include "list.h"
#include "gapDB.h"
#include "misc.h"

#ifdef NOMEMMOVE
#    define memmove(a,b,c) bcopy(b,a,c)
#endif

#ifndef WINKLE
#define WINKLE 30000
#endif

/* complements a sequence */
char *compl(char *cp) {
    int i, l = strlen(cp)-1, l2 = (l+1)/2;
    char swap[256];

    for (i=0; i<256; i++)
	swap[i] = i;

    swap['A'] = 'T';
    swap['C'] = 'G';
    swap['G'] = 'C';
    swap['T'] = 'A';
    swap['a'] = 't';
    swap['c'] = 'g';
    swap['g'] = 'c';
    swap['t'] = 'a';

    for (i=0; i<l2; i++) {
	int t;
	
	t = cp[i];
	cp[i] = swap[(unsigned)cp[l-i]];
	cp[l-i] = swap[t];
    }

    if (i*2 == l)
	cp[i] = swap[(unsigned)cp[i]];

    return cp;
}

static GapIO *io;

void gap_open_for_read(List *l) {
}

List *gap_read_header(void) {
    return NULL;
}

void gap_close(List *l) {
    GClones c;
    GVectors v;

    /* Allocate a default clone */
    io->db.Nclones=1;
    ArrayRef(io->clones, 0);
    arr(GCardinal, io->clones, 0) = allocate(io, GT_Clones);
    memset(&c, 0, sizeof(c));
    c.name = allocate(io, GT_Text);
    TextWrite(io, c.name, "unknown", strlen("unknown"));
    c.vector = 1;
    GT_Write(io, arr(GCardinal, io->clones, 0), &c, sizeof(c), GT_Clones);

    /* Allocate a default vector */
    io->db.Nvectors=1;
    ArrayRef(io->vectors, 0);
    arr(GCardinal, io->vectors, 0) = allocate(io, GT_Vectors);
    memset(&v, 0, sizeof(v));
    v.name = allocate(io, GT_Text);
    TextWrite(io, v.name, "unknown", strlen("unknown"));
    v.level = 1;
    GT_Write(io, arr(GCardinal, io->vectors, 0), &v, sizeof(v), GT_Vectors);

    ArrayWrite(io, io->db.vectors, io->db.Nvectors, io->vectors);
    ArrayWrite(io, io->db.clones, io->db.Nclones, io->clones);
    ArrayWrite(io, io->db.templates, io->db.Ntemplates, io->templates);

    close_db(io);
}

List *gap_read_gel_data(void) {
    return NULL;
}

List *gap_read_contig_data(void) {
    return NULL;
}

void gap_open_for_write(List *l) {
    char *name;
    char *version;
    int status;

    name = assoc(l, db_name);
    if (!name) crash("No database name specified\n");

    version = assoc(l, db_version);
    if (!version) crash("No version specified\n");

    io = open_db(name, version, &status, 1, 0);
    if (!io) {
	puts("Failed to open Gap database - maybe database already exists.");
	exit(1);
    }
}

void gap_write_header(List *l) {
    char *cp;

    if (cp == assoc(l, db_data_class))
	io->db.data_class = atoi(cp);

    if (cp == assoc(l, db_max_gel_length))
	io->db.max_gel_len = atoi(cp);

    if (cp == assoc(l, db_max_db_size) && atoi(cp) > io->db.maximum_db_size) {
	io->db.maximum_db_size = io->db.actual_db_size = atoi(cp);
    }

    GT_Write(io, GR_Database, &io->db, sizeof(io->db), GT_Database);
}

void gap_write_gel_data(List *l) {
    GReadings r;
    GTemplates t;
    int length = 0, start = 0, end = 0, ostart = 1;
    int left, right, pos, comp, i, gel_num;
    char seq[WINKLE], *cp;
    int1 conf[WINKLE];
    int2 opos[WINKLE];
    List *notes, *edits;

    if (cp = assoc(l, gel_l_nbr))
	left = atoi(cp);

    if (cp = assoc(l, gel_r_nbr))
	right = atoi(cp);

    if (cp = assoc(l, gel_l_nbr))
	left = atoi(cp);

    if (cp = assoc(l, gel_pos))
	pos = atoi(cp);

    if (cp = assoc(l, gel_comp))
	comp = atoi(cp);

    seq[0] = '\0';
    if (comp) {
	if (cp = assoc(l, gel_r_cut_seq)) {
	    strcat(seq, compl(cp));
	    start = strlen(cp);
	} else
	    start = 0;

	if (NULL == (cp = assoc(l, gel_seq)))
	    crash("No sequence for gel!\n");
	else
	    strcat(seq, cp);

	end = strlen(seq) + 1;
	if (cp = assoc(l, gel_l_cut_seq)){
	    strcat(seq, compl(cp));
	    ostart = strlen(cp);
	}
    } else {
	if (cp = assoc(l, gel_l_cut_seq)) {
	    strcat(seq, cp);
	    ostart = start = strlen(cp);
	} else
	    start = 0;

	if (NULL == (cp = assoc(l, gel_seq)))
	    crash("No sequence for gel!\n");
	else
	    strcat(seq, cp);

	end = strlen(seq) + 1;
	if (cp = assoc(l, gel_r_cut_seq))
	    strcat(seq, cp);
    }

    /* gel_rd_length can't be trusted? */
    length = strlen(seq);

    /* init conf and opos */
    memset(conf, 100, length * sizeof(*conf));
    if (comp) {
	for (i = 0; i < length; i++) {
	    opos[i] = length - i;
	}
    } else {
	for (i = 0; i < length; i++) {
	    opos[i] = i+1;
	}
    }

    /* Write reading structure */
    gel_num = Nreadings(io)+1;

    io->db.Ntemplates++;
    io_write_seq(io, gel_num, &length, &start, &end, seq, conf, opos);

    /* Update other gel details */
    GT_Read(io, arr(GCardinal, io->readings, gel_num - 1),
	    &r, sizeof(r), GT_Readings);
    
    r.left = left;
    r.right = right;
    r.position = pos;
    r.sense = comp;
    r.start = start;
    r.end = end;
    r.sequence_length = r.end - r.start - 1; 
    r.template = gel_num;

    /* name */
    r.name = allocate(io, GT_Text);
    cp = assoc(l, gel_name);
    TextWrite(io, r.name, cp, strlen(cp));

    /* trace */
    if (cp = assoc(l, gel_rd_file)) {
	r.trace_name = allocate(io, GT_Text);
	TextWrite(io, r.trace_name, cp, strlen(cp));
    }

    if (cp = assoc(l, gel_rd_type)) {
	r.trace_type = allocate(io, GT_Text);
	TextWrite(io, r.trace_type, cp, strlen(cp));
    }

    /* Annotations */
    r.annotations = 0;
    notes = index_list_by_str(l, gel_annotation);
    if (!isNil(notes)) {

	notes = cdr(notes);
	if (!isNil(notes))
	    r.annotations = Nannotations(io) + 1;

	while(!isNil(notes)) {
	    List *n;
	    GAnnotations a;

	    n = car(notes);

	    if (cp = assoc(n, gel_an_pos))
		a.position = atoi(cp) + ostart;
	    else
		crash("No position for annotation");

	    /* Fix for corrupted databases */
	    if (a.position < 1) {
		puts("Fixing tag position\n");
		a.position = 1;
	    }

	    if (a.position > r.length) {
		puts("Fixing tag position\n");
		a.position = r.length;
	    }

	    if (cp = assoc(n, gel_an_len))
		a.length = atoi(cp);
	    else
		crash("No length for annotation");

	    /* Fix for corrupted databases */
	    if (a.length < 0) {
		puts("Fixing tag length\n");
		a.position = 0;
	    }

	    if (a.length + a.position > r.length + 1) {
		puts("Fixing tag length\n");
		a.length = r.length + 1 - a.position;
	    }

	    if (cp = assoc(n, gel_an_type))
		a.type = (cp[0] << 24) + (cp[1] << 16)
		    + (cp[2] << 8) + cp[3];
	    else
		crash("No type for annotation");

	    if (cp = assoc(n, gel_an_comment)) {
		a.annotation = allocate(io, GT_Text);
		TextWrite(io, a.annotation, cp, strlen(cp));
	    } else
		a.annotation = 0;

	    a.strand = 2; /* both */

	    notes = cdr(notes);
	    if (!isNil(notes))
		a.next = Nannotations(io) + 2;
	    else
		a.next = 0;
     
	    io_init_annotations(io, Nannotations(io)+1);
	    GT_Write(io, arr(GCardinal, io->annotations, Nannotations(io)-1),
		     &a, sizeof(a), GT_Annotations);
	}
    }

    /*
     * Original positions, already allocated by io_write_seq.
     * Scan through edits, shifting and changing opos as necessary.
     */
    edits = index_list_by_str(l, gel_edits);
    if (!isNil(edits)) {
	for (edits = cdr(edits); !isNil(edits); edits = cdr(edits)) {
	    List *e;
	    int pos;

	    e = car(edits);

	    if (comp) {
		pos = r.sequence_length - atoi(assoc(e, gel_ed_base_pos))
		    + start;

		if ((cp = assoc(e, gel_ed_op)) && pos <= r.length && pos > 0) {
		    if (strcmp(cp, gel_ed_insert) == 0) {
			memmove(&opos[0], &opos[1],
				(pos-1) * sizeof(*opos));
			opos[pos-1] = 0;

		    } else if (strcmp(cp, gel_ed_delete) == 0) {
			memmove(&opos[1], &opos[0],
				(pos) * sizeof(*opos));
			opos[0]++;
		    }
		}

	    } else {
		pos = atoi(assoc(e, gel_ed_base_pos)) + start;

		if ((cp = assoc(e, gel_ed_op)) && pos <= r.length && pos>= 0) {
		    if (strcmp(cp, gel_ed_insert) == 0) {
			memmove(&opos[pos], &opos[pos - 1],
				(r.length - pos) * sizeof(*opos));
			opos[pos-1] = 0;
			
		    } else if (strcmp(cp, gel_ed_delete) == 0) {
			memmove(&opos[pos-1], &opos[pos],
				(r.length - pos) * sizeof(*opos));
			opos[r.length-1]++;
		    }
		}
	    }
	}
    }

    DataWrite(io, r.orig_positions, opos, r.length * sizeof(int2),
	      sizeof(int2));

    /* rewrite gel reading structure */
    GT_Write(io, arr(GCardinal, io->readings, gel_num - 1),
	     &r, sizeof(r), GT_Readings);

    /* Create template structure */
    ArrayRef(io->templates, gel_num-1);
    arr(GCardinal, io->templates, gel_num-1) = allocate(io, GT_Templates);;
    memset(&t, 0, sizeof(t));

    /* Update template structure */
    t.name = allocate(io, GT_Text);
    cp = assoc(l, gel_name);
    TextWrite(io, t.name, cp, strlen(cp));
    t.strands = 1;
    t.vector = 1;
    t.clone = 1;
    t.insert_length_min = 1;
    t.insert_length_max = 99999;
    GT_Write(io, arr(GCardinal, io->templates, gel_num - 1),
	     &t, sizeof(t), GT_Templates);
}

void gap_write_contig_data(List *l) {
    GContigs c;
    char *cp;
    
    /* printf("io_init_contig %d\n", NumContigs(io)+1); */
    io_init_contig(io, NumContigs(io)+1);

    GT_Read(io, arr(GCardinal, io->contigs, NumContigs(io)-1),
	    &c, sizeof(c), GT_Contigs);

    if (cp = assoc(l, contig_left_end))
	c.left = atoi(cp);

    if (cp = assoc(l, contig_right_end))
	c.right = atoi(cp);

    if (cp = assoc(l, contig_length))
	c.length = atoi(cp);

    c.annotations = 0;

    GT_Write(io, arr(GCardinal, io->contigs, NumContigs(io)-1),
	     &c, sizeof(c), GT_Contigs);
}

/*
 * Stub functions to avoid link errors from IO.c.
 * This is the 'note execution' code, which we wish to avoid anyway.
 */
void execute_database_notes(GapIO *io, char *type) {}
void process_rawdata_note(GapIO *io) {}
void fix_notes(GapIO *io) {}
