#include "tg_gio.h"

#define DB_NAMELEN 40
#define GGN_NAME 1
#define GGN_ID 0

int io_clength(GapIO *io, int cnum);
int io_clnbr(GapIO *io, int cnum);
int io_crnbr(GapIO *io, int cnum);
int io_length(GapIO *io, int rnum);
int io_relpos(GapIO *io, int rnum);
int io_lnbr(GapIO *io, int rnum);
int io_rnbr(GapIO *io, int rnum);

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

int CalcTotalContigLen(GapIO *io);
void bell(void);

/* FIXME */
#define io_rdonly(io) (0)
#define find_max_gel_len(io,a,b) 65536
#define io_name(io) "FIXME"
#define io_dbsize(io) 1
