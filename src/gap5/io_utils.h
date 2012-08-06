#ifndef _IO_UTILS_H
#define _IO_UTILS_H

#include <errno.h>
#include <string.h>
#include "tg_gio.h"

/*
 * ---------------------------------------------------------------------------
 * Useful macros
 * ---------------------------------------------------------------------------
 */

#if 0
#define contig_read(io, cn, c) \
    GT_Read(io, arr(GCardinal, io->contigs, cn - 1), \
	    &c, sizeof(c), GT_Contigs)

#define contig_write(io, cn, c) \
    GT_Write(io, arr(GCardinal, io->contigs, cn - 1), \
	     &c, sizeof(c), GT_Contigs)


#define template_read(io, tpn, tp) \
    GT_Read(io, arr(GCardinal, io->templates, tpn - 1), \
	    &tp, sizeof(tp), GT_Templates)

#define template_write(io, tpn, tp) \
    GT_Write(io, arr(GCardinal, io->templates, tpn - 1), \
	     &tp, sizeof(tp), GT_Templates)


#if GAP_CACHE!=0
#    define gel_read(io, gn, g) \
        (gn > 0 ? (memcpy(&g, arrp(GReadings, (io)->reading, (gn)-1), sizeof(g)), 0) : -1)
#    define gel_write(io, gn, g) \
        GT_Write_cached(io, gn, &g)
#else
#    define gel_read(io, gn, g) \
        GT_Read(io, arr(GCardinal, io->readings, gn - 1), \
	        &g, sizeof(g), GT_Readings)
#    define gel_write(io, gn, g) \
        GT_Write(io, arr(GCardinal, io->readings, gn - 1), \
	        &g, sizeof(g), GT_Readings)
#endif

#if 0 /* gap4 bits */

#define tag_read(io, tn, t) \
    GT_Read(io, arr(GCardinal, io->annotations, tn - 1), \
	    &t, sizeof(t), GT_Annotations)

#define tag_write(io, tn, t) \
    GT_Write(io, arr(GCardinal, io->annotations, tn - 1), \
	     &t, sizeof(t), GT_Annotations)

#define note_read(io, nn, n) \
    GT_Read(io, arr(GCardinal, io->notes, nn - 1), \
	    &n, sizeof(n), GT_Notes)

#define note_write(io, nn, n) \
    GT_Write(io, arr(GCardinal, io->notes, nn - 1), \
	     &n, sizeof(n), GT_Notes)

#define vector_read(io, vn, v) \
    GT_Read(io, arr(GCardinal, io->vectors, vn - 1), \
	    &v, sizeof(v), GT_Vectors)

#define vector_write(io, vn, v) \
    GT_Write(io, arr(GCardinal, io->vectors, vn - 1), \
	     &v, sizeof(v), GT_Vectors)


#define clone_read(io, cn, c) \
    GT_Read(io, arr(GCardinal, io->clones, cn - 1), \
	    &c, sizeof(c), GT_Clones)

#define clone_write(io, cn, c) \
    GT_Write(io, arr(GCardinal, io->clones, cn - 1), \
	     &c, sizeof(c), GT_Clones)
#endif

#define gel_read(io, n, c)      x_gel_read(io, n, &c)
#define contig_read(io, n, c)   x_contig_read(io, n, &c)
#define tag_read(io, n, c)      x_tag_read(io, n, &c)
#define template_read(io, n, c) x_template_read(io, n, &c)
#define vector_read(io, n, c)   x_vector_read(io, n, &c)
#define clone_read(io, n, c)    x_clone_read(io, n, &c)

#define gel_write(io, n, c)      x_gel_write(io, n, &c)
#define contig_write(io, n, c)   x_contig_write(io, n, &c)
#define tag_write(io, n, c)      x_tag_write(io, n, &c)
#define template_write(io, n, c) x_template_write(io, n, &c)
#define vector_write(io, n, c)   x_vector_write(io, n, &c)
#define clone_write(io, n, c)    x_clone_write(io, n, &c)

#endif /* Gap4 bits */


/*
 * ---------------------------------------------------------------------------
 * Functions
 * ---------------------------------------------------------------------------
 */


#define GGN_NAME 1
#define GGN_ID 0

tg_rec get_gel_num(GapIO *io, char *gel_name, int is_name);
tg_rec chain_left(GapIO *io, tg_rec gel);
tg_rec rnumtocnum(GapIO *io, tg_rec gel);
tg_rec get_contig_num(GapIO *io, char *gel_name, int is_name);

int lget_gel_num(GapIO *io, int listArgc, char **listArgv,
		 int *rargc, tg_rec **rargv);
int lget_contig_num(GapIO *io, int listArgc, char **listArgv,
		    int *rargc, contig_list_t **rargv);
int lget_contig_num2(GapIO *io, int listArgc, char **listArgv,
		     int *rargc, contig_list_t **rargv);
int *to_contigs_only(int num_contigs, contig_list_t *cl);
int lget_scaffold_num(GapIO *io, int listArgc, char **listArgv,
		      int *rargc, tg_rec **rargv);
tg_rec scaffold_name_to_number(GapIO *io, char *scaf_name);

char *get_read_name(GapIO *io, tg_rec number);
char *get_contig_name(GapIO *io, tg_rec number);
char *get_vector_name(GapIO *io, tg_rec vector);
char *get_template_name(GapIO *io, tg_rec tmplate);
char *get_clone_name(GapIO *io, tg_rec clone);

void cache_template_name(GapIO *io, tg_rec number, char *name);
void cache_read_name(GapIO *io, tg_rec number, char *name);
void cache_delete_read_name(GapIO *io, tg_rec number);

/*
 * Converts a template name to a template number.
 *
 * Arguments:
 *     io	- GapIO *
 *     tname    - the string described above
 *
 * Returns:
 *    0 for failure, otherwise the template number
 */
tg_rec template_name_to_number(GapIO *io, char *tname);

#endif /* _IO_UTILS_H */
