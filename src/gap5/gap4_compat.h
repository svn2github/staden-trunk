#include "tg_gio.h"

#define DB_NAMELEN 1024
#define GGN_NAME 1
#define GGN_ID 0

int io_clength(GapIO *io, tg_rec cnum);
int io_cclength(GapIO *io, tg_rec cnum);
tg_rec io_clnbr(GapIO *io, tg_rec cnum);
tg_rec io_crnbr(GapIO *io, tg_rec cnum);
int io_length(GapIO *io, tg_rec rnum);
int io_relpos(GapIO *io, tg_rec rnum);
tg_rec io_lnbr(GapIO *io, tg_rec rnum);
tg_rec io_rnbr(GapIO *io, tg_rec rnum);

tg_rec contig_name_to_number(GapIO *io, char *name);
tg_rec get_gel_num(GapIO *io, char *gel_name, int is_name);

#define NumContigs(io) ((io)->db->Ncontigs)
#define flush2t cache_flush

/* hack to make it compile until we create a struct */
/*
 * Annotations are singly linked lists, terminated with next == 0
 */
typedef struct {
    GCardinal type;
    GCardinal position;
    GCardinal length;
    GCardinal strand;
    GCardinal annotation;
    GCardinal next;
} GAnnotations;

int64_t CalcTotalContigLen(GapIO *io);
void bell(void);

/* FIXME */
#define io_rdonly(io) ((io)->read_only)
#define find_max_gel_len(io,a,b) 65536
#define io_name(io) ((io)->name)

/*
 * Complements an individual contig.
 * Returns 0 for success
 *        -1 for failure
 */
int complement_contig(GapIO *io, tg_rec crec);

/*
 * Converts a reading name to a reading number.
 *
 * Arguments:
 *     io	- GapIO *
 *     rname    - the string described above
 *
 * Returns:
 *    0 for failure, otherwise the gel number
 */
tg_rec read_name_to_number(GapIO *io, char *gel_name);

void bell(void);
