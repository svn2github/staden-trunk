#ifndef _IO_UTILS_H
#define _IO_UTILS_H

#include <errno.h>
#include <string.h>
#include "IO.h"

typedef struct contig_list {
    int contig;
    int start;
    int end;
} contig_list_t;

/*
 * ---------------------------------------------------------------------------
 * Useful macros
 * ---------------------------------------------------------------------------
 */

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



/*
 * ---------------------------------------------------------------------------
 * Functions
 * ---------------------------------------------------------------------------
 */


#define GGN_NAME 1
#define GGN_ID 0

int get_gel_num(GapIO *io, char *gel_name, int is_name);
int chain_left(GapIO *io, int gel);
int rnumtocnum(GapIO *io, int gel);
int get_contig_num(GapIO *io, char *gel_name, int is_name);

int lget_gel_num(GapIO *io, int listArgc, char **listArgv,
		 int *rargc, int **rargv);
int lget_contig_num(GapIO *io, int listArgc, char **listArgv,
		    int *rargc, contig_list_t **rargv);
int *to_contigs_only(int num_contigs, contig_list_t *cl);

char *get_read_name(GapIO *io, int number);
char *get_contig_name(GapIO *io, int number);
char *get_vector_name(GapIO *io, int vector);
char *get_template_name(GapIO *io, int tmplate);
char *get_clone_name(GapIO *io, int clone);

void cache_template_name(GapIO *io, int number, char *name);
void cache_read_name(GapIO *io, int number, char *name);
void cache_delete_read_name(GapIO *io, int number);

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
int template_name_to_number(GapIO *io, char *tname);

void update_rnumtocnum(GapIO *io, int gel, int contig);
void invalidate_rnumtocnum(GapIO *io, int disable);

#endif /* _IO_UTILS_H */
