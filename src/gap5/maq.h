#ifndef _MAQ_H_
#define _MAQ_H_

#include <tg_index.h>

#define PAIRFLAG_FF      0x01
#define PAIRFLAG_FR      0x02
#define PAIRFLAG_RF      0x04
#define PAIRFLAG_RR      0x08
#define PAIRFLAG_PAIRED  0x10
#define PAIRFLAG_DIFFCHR 0x20
#define PAIRFLAG_NOMATCH 0x40
#define PAIRFLAG_SW      0x80

int parse_maqmap(GapIO *io, char *dat_fn, tg_args *a);


#endif /* _MAQ_H_ */
