#include "IO.h"
#include "misc.h"
#include "io-reg.h"
#include "fort.h"


int break_contig(GapIO *io, int r_num) {
    f_int ngels, nconts;
    f_int iok;
    int contig;

    ngels = NumReadings(io);
    nconts = NumContigs(io);

    if (r_num < 1 || r_num > ngels) {
	verror(ERR_FATAL, "break_contig", "Invalid reading number");
	return 0;
    }

    if (-1 == (contig = rnumtocnum(io, r_num))) {
	verror(ERR_FATAL, "break_contig", "reading is not in a contig!");
	return 0;
    }

    if (contig_lock_write(io, contig) == -1) {
	verror(ERR_WARN, "break_contig", "Contig is busy");
	return 0;
    }

    breakc_(&io_relpos(io,1), &io_length(io,1), &io_lnbr(io,1),
           &io_rnbr(io,1), &ngels, &nconts, &io_dbsize(io),
           handle_io(io), &r_num, &iok);

    return 0;
} /* end break_contig */

