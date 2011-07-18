#ifndef _TG_TRACK_H_
#define _TG_TRACK_H_

/*
 * 'get' functions - simply returns the structure member.
 *
 * <type> track_get_XXX(seq_t **s)
 */
#define track_get_type(t)          ((*(t))->type)
#define track_get_flag(t)          ((*(t))->flag)
#define track_get_rec(t)           ((*(t))->rec)
#define track_get_bin_size(t)      ((*(t))->bin_size)
#define track_get_item_size(t)     ((*(t))->item_size)
#define track_get_nitems(t)        ((*(t))->nitems)
#define track_get_data(t)          ((*(t))->data)

/*
 * 'set' functions all have prototype:
 *
 * int track_set_XXX(GapIO *io, track_t **t, <type> new_value)
 *
 * Returns 0 for success, possibly also modifying *s pointer
 *        -1 for failure
 */
int track_set_type         (GapIO *io, track_t **t, int value);
int track_set_flag         (GapIO *io, track_t **t, int value);
int track_set_bin_size     (GapIO *io, track_t **t, int value);
int track_set_item_size    (GapIO *io, track_t **t, int value);
int track_set_nitems       (GapIO *io, track_t **t, int value);
int track_set_data         (GapIO *io, track_t **t, Array value);

#endif /* _TG_TRACK_H_ */
