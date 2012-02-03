#ifndef _TG_SEQUENCE_H_
#define _TG_SEQUENCE_H_

#define CONS_BIN_SIZE 100000

/*
 * 'get' functions - simply returns the structure member.
 *
 * <type> sequence_get_XXX(seq_t **s)
 */
#define sequence_get_pos(s)          ((*(s))->pos)
#define sequence_get_len(s)          ((*(s))->len)
#define sequence_get_left(s)         ((*(s))->left)
#define sequence_get_right(s)        ((*(s))->right)
#define sequence_get_mapping_qual(s) ((*(s))->mapping_qual)
#define sequence_get_name_len(s)     ((*(s))->name_len)
#define sequence_get_name(s)         ((*(s))->name)
#define sequence_get_trace_name(s)   ((*(s))->trace_name)
#define sequence_get_seq(s)          ((*(s))->seq)
#define sequence_get_conf(s)         ((*(s))->conf)
#define sequence_get_flags(s)        ((*(s))->flags)
#define sequence_get_bin_index(s)    ((*(s))->bin_index)
#define sequence_get_seq_tech(s)     ((*(s))->seq_tech)
#define sequence_get_parent_type(s)  ((*(s))->parent_type)
#define sequence_get_parent_rec(s)   ((*(s))->parent_rec)
#define sequence_get_aux_len(s)      ((*(s))->aux_len)
#define sequence_get_aux(s)          ((*(s))->sam_aux)


/*
 * 'set' functions all have prototype:
 *
 * int sequence_set_XXX(GapIO *io, seq_t **s, <type> new_value)
 *
 * Returns 0 for success, possibly also modifying *s pointer
 *        -1 for failure
 */
int sequence_set_pos         (GapIO *io, seq_t **s, int value);
int sequence_set_len         (GapIO *io, seq_t **s, int value);
int sequence_set_left        (GapIO *io, seq_t **s, int value);
int sequence_set_right       (GapIO *io, seq_t **s, int value);
int sequence_set_bin_index   (GapIO *io, seq_t **s, int value);
int sequence_set_parent_type (GapIO *io, seq_t **s, int value);
int sequence_set_parent_rec  (GapIO *io, seq_t **s, int value);
int sequence_set_flags       (GapIO *io, seq_t **s, int value);
int sequence_set_seq_tech    (GapIO *io, seq_t **s, int value);
int sequence_set_mapping_qual(GapIO *io, seq_t **s, uint8_t value);
int sequence_set_name        (GapIO *io, seq_t **s, char *name);
int sequence_set_trace_name  (GapIO *io, seq_t **s, char *trace_name);
int sequence_set_seq         (GapIO *io, seq_t **s, char *seq);
int sequence_set_conf        (GapIO *io, seq_t **s, char *conf);

int sequence_set_left_no_invalidate (GapIO *io, seq_t **s, int value);
int sequence_set_right_no_invalidate(GapIO *io, seq_t **s, int value);

tg_rec sequence_new_from(GapIO *io, seq_t *s);

tg_rec sequence_index_query(GapIO *io, char *name);
tg_rec *sequence_index_query_all(GapIO *io, char *name, int prefix,
				 int *nrecs);
int sequence_index_update(GapIO *io, char *name, int name_len, tg_rec rec);
btree_iter_t *sequence_index_iter(GapIO *io, char *name);

int sequence_get_position(GapIO *io, tg_rec snum, tg_rec *contig,
			  int *start, int *end, int *orient);
int sequence_get_position2(GapIO *io, tg_rec snum, tg_rec *contig,
			   int *start, int *end, int *orient,
			   tg_rec *brec, range_t *r_out, seq_t **s_out);
int sequence_get_clipped_position(GapIO *io, tg_rec snum, tg_rec *contig,
				  int *start, int *end,
				  int *clipped_start, int *clipped_end,
				  int *orient);

int sequence_get_orient(GapIO *io, tg_rec snum);
tg_rec sequence_get_contig(GapIO *io, tg_rec snum);
tg_rec sequence_get_pair(GapIO *io, seq_t *s);

/*
 * Trivial one-off sequence query functions
 */
int seq_pos(GapIO *io, tg_rec rec);
int seq_len(GapIO *io, tg_rec rec);
int seq_left(GapIO *io, tg_rec rec);
int seq_right(GapIO *io, tg_rec rec);
int seq_mapping_qual(GapIO *io, tg_rec rec);
char *seq_name(GapIO *io, tg_rec rec);
char *seq_trace_name(GapIO *io, tg_rec rec);
char *seq_seq(GapIO *io, tg_rec rec);
char *seq_conf(GapIO *io, tg_rec rec);

/*
 * Reverses and complements a piece of DNA
 */
void complement_seq_conf(char *seq, char *conf, int seq_len, int nconf);

seq_t *dup_seq(seq_t *s);
size_t sequence_extra_len(seq_t *s);

void complement_seq_t(seq_t *s);

int sequence_get_base(GapIO *io, seq_t **s, int pos, char *base, int *conf,
		      int *cutoff, int contig_orient);
int sequence_get_base4(GapIO *io, seq_t **s, int pos, char *base, double *conf,
		       int *cutoff, int contig_orient);
int sequence_replace_base(GapIO *io, seq_t **s, int pos, char base, int conf,
			  int contig_orient);
int sequence_insert_base(GapIO *io, seq_t **s, int pos, char base, char conf,
			 int contig_orient);
int sequence_delete_base(GapIO *io, seq_t **s, int pos,
			 int contig_orient);
int sequence_delete_base2(GapIO *io, seq_t **s, int pos, int contig_orient,
			  int check_base);
			  
/* Fix range_t length values after insert or delete */
int sequence_range_length(GapIO *io, seq_t **s);

int sequence_invalidate_consensus(GapIO *io, seq_t *s);

/*
 * Given a seq_t struct this updates the internal pointers to be valid offsets
 * into the s->data field. This is useful if the structure has been copied to
 * a new address.
 */
void sequence_reset_ptr(seq_t *s);

/*
 * Copies the 'f' seq_t struct to the 's' seq_t struct.
 * Assumes 's' has already been allocated to be large enough to hold 'f'.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int  sequence_copy(seq_t *s, seq_t *f);

/*
 * Moves all annotations attached to a sequence left or right by a certain
 * amount 'dist'. If dist is negative the annotations move left, otherwise
 * they move right.
 *
 * This is used when moving sequences within the editor, just prior to the
 * move itself as we only look for annotations covering the current
 * coordinates.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_move_annos(GapIO *io, seq_t **s, int dist);

/*
 * If we've done something that changed the bin a sequence is within without
 * also moving the tags to the new bin, eg by changing its size or location,
 * then this function will check annotations and ensure they're in the same
 * bin as the sequence they belong to.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sequence_fix_anno_bins(GapIO *io, seq_t **s);

#define TEMPLATE_ERR		-1 /* Unknown due to error */
#define TEMPLATE_SINGLE		0
#define TEMPLATE_PAIRED		1  /* Normal OK pair */
#define TEMPLATE_DISTANCE	2  /* Distance incorrect */
#define TEMPLATE_ORIENT		3  /* Invalid orientation */
#define TEMPLATE_SPANNING	4  /* Spans two contigs */

/*
 * Fetches information about a template - the size, status, type, etc.
 * Status is the primary return and all other returned fields passed as
 * arguments may be NULL.
 *
 * Returns status code on success
 *        -1 on failure
 */
int sequence_get_template_info(GapIO *io, seq_t *s,
			       tg_rec *library,
			       int *size);

range_t *sequence_get_range(GapIO *io, seq_t *s);

#endif /* _TG_SEQUENCE_H_ */
