#include "IO.h"

/* ARGSUSED */
void busy_dialog(GapIO *io, int contig) {
    printf("The contig %d is busy (probably due to an editor in use for this "
	   "contig). Changes will not be made for this contig.",
	   contig);
}

/* ARGSUSED */
void update_results(GapIO *io) {
}
