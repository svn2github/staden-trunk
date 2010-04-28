#include <staden_config.h>

#include <tcl.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "os.h"
#include "capture.h"
#include "misc.h"
#include "tcl_utils.h"

int tcl_capture(ClientData clientData, Tcl_Interp *interp,
		int argc, char **argv) {
    int old_stdout;
    static int fd = 0;
    char *buf;
    struct stat statbuf;
    char *tmpfile;
    int result;

    if (argc != 2 && argc != 3) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], " command ?varName?\"", NULL);
	return TCL_ERROR;
    }

    /* File descriptor mangling */
    if (!fd) {
	tmpfile = tmpnam(NULL);
	fd = open(tmpfile, O_RDWR|O_CREAT|O_TRUNC, 0666);
    } else {
	lseek(fd, 0, SEEK_SET);
    }

    old_stdout = dup(1);
    close(1);
    dup2(fd, 1);

    /* Run the command */
    result = Tcl_Eval(interp, argv[1]);

    /* Reset file descriptors */
    dup2(old_stdout, 1);
    close(old_stdout);

    /* Reload the output */
    fstat(fd, &statbuf);
    if (NULL == (buf = (char *)xmalloc(statbuf.st_size+1)))
	return TCL_ERROR;
    lseek(fd, 0, SEEK_SET);
    read(fd, buf, statbuf.st_size);
    buf[statbuf.st_size]=0;

    /* Return it to Tcl */
    if (argc == 3) {
	Tcl_ResetResult(interp);
	vTcl_SetResult(interp, "%d", result);
	return Tcl_SetVar(interp, argv[2], buf, 0) ? TCL_OK : TCL_ERROR;
    } else {
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	free(buf);
    }

    return TCL_OK;
}
