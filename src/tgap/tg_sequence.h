#ifndef _TG_SEQUENCE_H_
#define _TG_SEQUENCE_H_

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

int sequence_new_from(GapIO *io, seq_t *s);

GRec sequence_index_query(GapIO *io, char *name);
int sequence_index_update(GapIO *io, char *name, int name_len, GRec rec);

int sequence_get_position(GapIO *io, GRec snum, int *contig,
			  int *start, int *end, int *orient);
int sequence_get_contig(GapIO *io, GRec snum);
int sequence_get_pair(GapIO *io, seq_t *s);

/*
 * Trivial one-off sequence query functions
 */
int seq_pos(GapIO *io, int rec);
int seq_len(GapIO *io, int rec);
int seq_left(GapIO *io, int rec);
int seq_right(GapIO *io, int rec);
int seq_mapping_qual(GapIO *io, int rec);
char *seq_name(GapIO *io, int rec);
char *seq_trace_name(GapIO *io, int rec);
char *seq_seq(GapIO *io, int rec);
char *seq_conf(GapIO *io, int rec);

/*
 * Reverses and complements a piece of DNA
 */
void complement_seq_conf(char *seq, char *conf, int seq_len, int nconf);

seq_t *dup_seq(seq_t *s);

void complement_seq_t(seq_t *s);

int sequence_get_base(GapIO *io, seq_t **s, int pos, char *base, int *conf,
		      int contig_orient);
int sequence_get_base4(GapIO *io, seq_t **s, int pos, char *base, double *conf,
		       int contig_orient);
int sequence_replace_base(GapIO *io, seq_t **s, int pos, char base, int conf,
			  int contig_orient);
int sequence_insert_base(GapIO *io, seq_t **s, int pos, char base, char conf,
			 int contig_orient);
int sequence_delete_base(GapIO *io, seq_t **s, int pos,
			 int contig_orient);


#endif /* _TG_SEQUENCE_H_ */
