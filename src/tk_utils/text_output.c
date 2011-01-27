/*
 * Duplicates the stdout and stderr streams to files for use with displaying
 * in a tk window.
 */

#include <staden_config.h>

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <tcl.h>
#include <tk.h>
#include <time.h>
#include <string.h>

#include "os.h"

#ifndef NOPIPE
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <sys/signal.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/wait.h>
#endif

#include "text_output.h"
#include "tcl_utils.h"
#include "vlen.h"
#include "os.h"
#include "getfile.h"
#include "FtoC.h"
#include "misc.h"
#include "getfile.h"

/* 7/1/99 johnt - added definitions for WINNT support */
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
static HANDLE eventLogHandle = NULL;
#undef WIN32_LEAN_AND_MEAN
/* EVEN_SOURCE MUST be same name as entry in Registry under
   KKEY_LOCAL_MACHINE\System\CurrentControlSet\Services\EventLog\Application */
#define EVENT_SOURCE "Staden Package" 
#endif


static char stdout_win[100], stderr_win[100];
static char cur_tag[100];
static int win_init = 0;
static int stdout_scroll = 1, stderr_scroll = 1;
static int header_outputted = 0;
static FILE *stdout_fp = NULL, *stderr_fp = NULL;
static Tcl_Interp *_interp = NULL;
static int noisy = 1;

static int info_win;
static Tcl_DString message;

static int logging = 1;
static int log_vmessage_st = 0;

void start_message(void)
{

    Tcl_DStringInit(&message);
    info_win = 1;

}

/*
 * Must be called for each and every start message in order to free the
 * DString allocated.
 * "parent" is a tk window path for the dialogue. Specify NULL to indicate
 * a bail-out state (where we free the DString but decide not to bring up
 * a message box afterall.
 */
void end_message(const char *parent)
{
    int argc = 1;
    char *argv[1], *merged;

    argv[0] = Tcl_DStringValue(&message);

    if (NULL == (merged = Tcl_Merge(argc, argv))) {
	info_win = 0;
	Tcl_DStringFree(&message);
	return;
    }

    /* display message box */
    if (parent && _interp) {
	Tcl_VarEval(_interp, "messagebox ", parent, " ", merged, NULL);
    }
    info_win = 0;
    Tcl_DStringFree(&message);
    Tcl_Free(merged);
}


/*
 * Outputs data to a persistent log file. The log file is always appended to,
 * but only contains headers and error messages.
 *
 * 'fn' is the filename to log to, or NULL if we use whichever one is currently
 * open. Not initialising the function will not cause any problems.
 *
 * Giving a filename of "" will close the log file.
 */
void log_file(const char *fn, const char *message) {
    static FILE *fp = NULL;
    char tbuf[100];
    static char hname[256];
    static int hname_done = 0;
    time_t t = time(NULL);

    strftime(tbuf, sizeof(tbuf)-1, "%a %d %b %H:%M:%S %Y", localtime(&t));

    if (!logging)
	return;

    if (!hname_done) {
#ifdef _WIN32
	sprintf(hname,"?");
#else
	if (gethostname(hname, 256))
	    sprintf(hname, "?");
#endif
	hname_done = 1;
    }

    if (fn) {
	if (fn && *fn == 0) {
	    if (fp) {
		if (message) {
		    fseeko(fp, 0, SEEK_END);
		    fprintf(fp, "%s [%d@%s] %s\n",
			    tbuf, (int)getpid(), hname, message);
		}
		fclose(fp);
		fp = NULL;
	    }
	} else {
	    if (fp)
		fclose(fp);
	    fp = fopen(fn, "a");
	}
    }

    if (fp && message) {
	fseeko(fp, 0, SEEK_END);
	fprintf(fp, "%s [%d@%s] %s\n", tbuf, (int)getpid(), hname, message);
	fflush(fp);
    }
}

/*
 * Controls whether vmessage output should also be written to the log file
 * (in addition to vfuncheader and verror messages).
 * A value of 0 means do not log. Any other values implies logging.
 */
void log_vmessage(int log) {
    log_vmessage_st = log;
}

static void tout_update_stream(int fd, const char *buf, int header,
			       const char *tag) {
    char * win;
    char tag_list[1024];

    if (!win_init) {
#ifdef _WIN32
	/* WINNT will not have stdout/err defined unless running in console mode
	 * so use a message box
	 */
	if( fileno(stdout) == -1 || fileno(stderr) == -1 ){
	    MessageBox(NULL,buf,"Error",MB_OK|MB_ICONERROR|MB_TASKMODAL);
	    return;
	}
#endif
	fprintf(fd == 1 ? stdout : stderr, "%s", buf);
	fflush(fd == 1 ? stdout : stderr);
	return;
    }

    win = fd == 1 ? stdout_win : stderr_win;

    /* Add to the redirection streams */
    if (fd == 1 && stdout_fp) {
	fprintf(stdout_fp, "%s", buf);
	fflush(stdout_fp);
    } else if (fd == 2 && stderr_fp) {
	fprintf(stderr_fp, "%s", buf);
	fflush(stderr_fp);
    }

    if (info_win) {
	Tcl_DStringAppend(&message, buf, strlen(buf));
    }

    if (tag) {
	sprintf(tag_list, "{%s%s %s}",
		cur_tag, header ? "_h" : "_t",
		tag);
    } else {
	sprintf(tag_list, "%s%s", cur_tag, header ? "_h" : "_t");
    }

    /* Add to the text widget */
    if (win_init) {
	Tcl_SetVar(_interp, "TEMP", buf, 0);

	Tcl_VarEval(_interp, win, " insert end ", "\"$TEMP\" ",
		    tag_list, NULL);

	if (fd == 1 ? stdout_scroll : stderr_scroll) {
	    /* scroll to bottom of output window */
	    Tcl_VarEval(_interp, win, " see end", NULL);
	}
    }
}

void bell(void) {
    if (_interp)
	Tcl_GlobalEval(_interp, "bell");
}

void funcparams(char *params) {

     if (win_init) {
	 Tcl_VarEval(_interp, "tout_tag_params ",
		     stdout_win,
		     " ", cur_tag,
		     " {", params, "}", NULL);
     }
 }

static void funcheader(const char *name) {
    char tbuf[100], buf[100+8192];
    time_t t = time(NULL);

    /* set 'header_outputted' to notify the next vfuncgroup() */
    header_outputted = 1;

    sprintf(cur_tag, "%d", atoi(cur_tag)+1);

    if (win_init) {
	Tcl_VarEval(_interp, "tout_new_header ",
		    stdout_win,
		    " ", cur_tag,
		    " {", name, "}", NULL);

	/* FIXME - add time stamp here */
	strftime(tbuf, sizeof(tbuf)-1, "%a %d %b %H:%M:%S %Y", localtime(&t));
	tout_update_stream(1, "============================================================\n", 1, NULL);
	sprintf(buf, "%s: %s\n", tbuf, name);
	log_file(NULL, name);
	tout_update_stream(1, buf, 1, NULL);
	tout_update_stream(1, "------------------------------------------------------------\n", 1, NULL);
    }
}


static void funcgroup(int group, const char *name) {
    static int group_num = 0;

    if (header_outputted || group != group_num) {
	funcheader(name);

	header_outputted = 0;
	group_num = ABS(group);
    }
}

#ifndef NOPIPE
/*
 *-----------------------------------------------------------------------------
 * A module for piping data to a program.
 *
 * This chunk is highly UNIX dependant. Goodness knows what it'll do on
 * creaky systems with micky mouse multitasking.
 *
 * As users are known to be awkward beasts at best, we're assuming they'll
 * try using shell style actions. Let's cheat, like system(), by using good
 * 'ol "sh -c".
 *-----------------------------------------------------------------------------
 */

/* Total time to take in fetching the output, in ms */
#define TIME_OUT 5000000 /* 5 seconds */

int pipe_mania(char *data, int len, char *command, int forever) {
    int pid, off = 0, count = 0, ret = -1;
    int fdp[3][2];
    char buf[8192+1];

    /*
     * Make the connections:
     *
     * fdp[0] is stdin for the child
     * fdp[1] is stdout for the child
     * fdp[2] is stderr for the child
     * fdp[x][0] is the read end, and fdp[x][1] is the write end.
     * Hence:
     *     fdp[0][0] = child's stdin
     *     fdp[1][1] = child's stdout
     *     fdp[2][1] = child's stderr
     *     fdp[0][1] = parent's output
     *     fdp[1][0] = parent's input of data
     *     fdp[2][0] = parent's input of errors
     */
    if (-1 == pipe(fdp[0]))
	return -1;

    if (-1 == pipe(fdp[1])) {
	close(fdp[0][0]);
	close(fdp[0][1]);
	return -1;
    }

    if (-1 == pipe(fdp[2])) {
	close(fdp[0][0]);
	close(fdp[0][1]);
	close(fdp[1][0]);
	close(fdp[1][1]);
	return -1;
    }

    switch(pid = fork()) {
    case 0: /* child */
	dup2(fdp[0][0], 0);
	dup2(fdp[1][1], 1);
	dup2(fdp[2][1], 2);
	close(fdp[0][1]);
	close(fdp[1][0]);
	close(fdp[2][0]);

	execlp("sh", "sh", "-c", command, NULL);
	exit(1);

    default: /* parent */
	close(fdp[0][0]);
	close(fdp[1][1]);
	close(fdp[2][1]);
	break;

    case -1: /* error */
	goto error;
    }

    /*
     * Set both parent ends to be non blocking. Deadlock can not be
     * completely avoided in a double pipe, so if something's going to
     * break we want to make sure it'll be the child, not the parent.
     */
    (void)fcntl(fdp[0][1], F_SETFL, O_NONBLOCK);
    (void)fcntl(fdp[1][0], F_SETFL, O_NONBLOCK);
    (void)fcntl(fdp[2][0], F_SETFL, O_NONBLOCK);

    do {
	int tmp, done_something;

	done_something = 0;

	/*
	 * Send our output.
	 */
	if (len) {
	    while (len > 0 && (tmp = write(fdp[0][1], &data[off], len)) >= 0) {
		off += tmp;
		len -= tmp;
		done_something = 1;
	    }
	    
	    /*
	     * If we've finished the write stage of things then we close the
	     * write end of the pipe which should send an EOF to the other
	     * end, which also helps to prevent deadlock.
	     */
	    if (0 == len) {
		close(fdp[0][1]);
	    }

	    /*
	     * We've got an error. If it's because of blocking then we could
	     * try to read data instead (until that blocks too)
	     */
	    if (-1 == len && errno != EWOULDBLOCK)
		goto error;
	}

	/*
	 * Read our input
	 */
	while ((tmp = read(fdp[1][0], buf, 8192)) > 0) {
	    buf[tmp] = 0;
	    tout_update_stream(1, buf, 0, NULL);
	    done_something = 1;
	}

	if (-1 == tmp && errno != EWOULDBLOCK)
	    goto error;


	/*
	 * Assume that reading zero bytes implies we've finished. Not ideal.
	 * We should check whether the process is still running amongst
	 * other things. However it's not such a loss as this module isn't
	 * going to be used to rot13 dictionaries (ok, ok, so that _is_ what
	 * I'm using for my testing!).
	 */
	if (tmp == 0)
	    break;

	/*
	 * We've written, and we've read back.
	 * If we've blocked both times then assume we're chasing a lost
	 * cause. Otherwise loop around again.
	 */
	if (!done_something) {
	    /* usleep(100000); */
	    sleep(1);
	    count += 1000000;
	}
    } while (count < TIME_OUT || forever);
    
    ret = (count < TIME_OUT || forever) ? 0 : -2;

    /*
     * Check we haven't had any stderr from the child - we'll redirect this
     * to the error window.
     *
     * NB: IO here is much sloppier with regards to robustness, but errs on
     * the side of caution for the parent.
     */
    if ((len = read(fdp[2][0], buf, 8192)) > 0) {
	char *p,  *o = buf;
	buf[len-1]=0;
	while (NULL != (p = strchr(o, '\n'))) {
	    *p = 0;
	    verror(ERR_WARN, "pipe", "stderr=%s", o);
	    o = p+1;
	}
	if (*o)
	    verror(ERR_WARN, "pipe", "stderr=%s", o);
    }

 error:
    kill(pid, SIGKILL);
    close(fdp[2][0]);
    close(fdp[1][0]);
    close(fdp[0][1]);
    waitpid(pid, &count, WNOHANG);

    return ret;
}
#endif

/*
 *-----------------------------------------------------------------------------
 * Tcl interface to these routines
 *-----------------------------------------------------------------------------
 */

int TextOutput_Init(Tcl_Interp *interp) {
    _interp = interp;

    Tcl_CreateCommand(interp, "tout_init", tcl_tout_init,
                      (ClientData) NULL,
                      NULL);
    Tcl_CreateCommand(interp, "tout_set_scroll", tcl_tout_set_scroll,
                      (ClientData) NULL,
                      NULL);
    Tcl_CreateCommand(interp, "tout_set_redir", tcl_tout_set_redir,
                      (ClientData) NULL,
                      NULL);
#ifndef NOPIPE
    Tcl_CreateCommand(interp, "tout_pipe", tcl_tout_pipe,
                      (ClientData) NULL,
                      NULL);
#endif
    Tcl_CreateCommand(interp, "vmessage", tcl_vmessage,
                      (ClientData) NULL,
                      NULL);
    Tcl_CreateCommand(interp, "vmessage_tagged", tcl_vmessage_tagged,
                      (ClientData) NULL,
                      NULL);
    Tcl_CreateCommand(interp, "verror", tcl_verror,
                      (ClientData) NULL,
                      NULL);
    Tcl_CreateCommand(interp, "vfuncheader", tcl_vfuncheader,
                      (ClientData) NULL,
                      NULL);
    Tcl_CreateCommand(interp, "vfuncgroup", tcl_vfuncgroup,
                      (ClientData) NULL,
                      NULL);
    Tcl_CreateCommand(interp, "error_bell", tcl_error_bell,
                      (ClientData) NULL,
                      NULL);

    Tcl_LinkVar(interp, "logging", (char *)&logging, TCL_LINK_INT);
    
    return TCL_OK;
}

int tcl_tout_init(ClientData clientData, Tcl_Interp *interp,
		  int argc, char **argv) {
    if (argc != 3)
	return TCL_ERROR;

    strcpy(stdout_win, argv[1]);
    strcpy(stderr_win, argv[2]);
    strcpy(cur_tag, "0");

    win_init++;

    return TCL_OK;
}

int tcl_tout_set_scroll(ClientData clientData, Tcl_Interp *interp,
			int argc, char **argv) {
    if (argc != 3)
	return TCL_ERROR;

    if (strcmp(argv[1], "stdout") == 0) {
	stdout_scroll = atoi(argv[2]);
    } else if (strcmp(argv[1], "stderr") == 0) {
	stderr_scroll = atoi(argv[2]);
    } else {
	return TCL_ERROR;
    }

    return TCL_OK;
}

int tcl_tout_set_redir(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv) {
    FILE **fp;

    if (argc != 3)
	return TCL_ERROR;

    if (strcmp(argv[1], "stdout") == 0) {
	fp = &stdout_fp;
    } else if (strcmp(argv[1], "stderr") == 0) {
	fp = &stderr_fp;
    } else {
	return TCL_ERROR;
    }

    /* Close the old stream if it's still open */
    if (*fp != NULL) {
	fclose(*fp);
	*fp = NULL;
    }
    
    /* Open the new stream */
    if (*argv[2] && NULL == (*fp = fopen(argv[2], "w"))) {
	Tcl_SetResult(interp, "0", TCL_STATIC);
    } else {
	Tcl_SetResult(interp, "1", TCL_STATIC);
    }

    return TCL_OK;
}

#ifndef NOPIPE
/*
 * Sends some text to a command and adds the command's stdout and stderr to
 * the output and error windows.
 *
 * Usage: tout_pipe command text
 */
int tcl_tout_pipe(ClientData clientData, Tcl_Interp *interp,
		  int argc, char **argv) {
    int ret;

    if (argc != 4)
	return TCL_ERROR;

    vfuncheader("Output from command '%s'", argv[1]);

    ret = pipe_mania(argv[2], strlen(argv[2]), argv[1], atoi(argv[3]));
    if (-1 == ret) {
	verror(ERR_WARN, "pipe", "command '%s' failed", argv[1]);
    } else if (-2 == ret) {
	verror(ERR_WARN, "pipe", "timeout - output from command truncated");
    }

    vTcl_SetResult(interp, "%d", ret);

    return TCL_OK;
}
#endif

/*
 * Tcl interface to vmessage.
 */
int tcl_vmessage(ClientData clientData, Tcl_Interp *interp,
		 int argc, char **argv) {
    int i, l;
    char buf[8192], *p2 = buf, *p = buf, *z;
    int newline = 1;
    int start_argc = 1;

    if (strcmp(argv[1], "-nonewline") == 0) {
	newline = 0;
	start_argc++;
    }

    for (i = start_argc, l = 0; i < argc; i++) {
	l += 1 + strlen(argv[i]);
    }
    l+=2;

    if (l >= 8192) {
	p = p2 = (char *)xmalloc(l);
    }
    *p = 0;
    for (i = start_argc; i < argc; i++, *p++ = ' ') {
	z = argv[i];
	while (*z)
	    *p++=*z++;
    }
    p--;
    if (newline)
	strcpy(p, "\n");
    else
	strcpy(p, "");

    tout_update_stream(1, p2, 0, NULL);

    if (p2 != buf)
	xfree(p2);

    return TCL_OK;
}

/*
 * Tcl interface to vmessage.
 */
int tcl_vmessage_tagged(ClientData clientData, Tcl_Interp *interp,
			int argc, char **argv) {
    int i;
    char *nl = "\n";
    int newline = 1;
    int start_argc = 1;

    if (strcmp(argv[1], "-nonewline") == 0) {
	newline = 0;
	start_argc++;
    }

    for (i = start_argc; i < argc-1; i+=2) {
	tout_update_stream(1, argv[i], 0, argv[i+1]);
    }
    if (newline)
	tout_update_stream(1, nl, 0, NULL);

    return TCL_OK;
}

/*
 * Tcl interface to verror
 * Restrictions - cannot output more than 8K.
 */
int tcl_verror(ClientData clientData, Tcl_Interp *interp,
	       int argc, char **argv) {
    int i, level, len;
    char buf[8192], tbuf[100], *p, *bufp = buf;
    time_t t = time(NULL);

    if (argc < 3) {
	return TCL_ERROR;
    }

    if (strcmp(argv[1], "ERR_WARN") == 0)
	level = ERR_WARN;
    else
	level = ERR_FATAL;

    len = 0;
    for (i = 2; i < argc; i++) {
	len += strlen(argv[i]);
    }
    len += 100; /* allow for time, etc */
    if (len > 8192) {
	if (NULL == (bufp = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "verror", "out of memory");
	    return TCL_OK;
	}
    }

    strftime(tbuf, sizeof(tbuf)-1, "%a %d %b %H:%M:%S %Y", localtime(&t));
    sprintf(bufp, "%s %.7500s: ", tbuf, argv[2]);
    p = bufp + strlen(bufp);

    for (i = 3; i < argc; i++, *p++ = ' ') {
	strcpy(p, argv[i]);
	p += strlen(p);
    }
    *(p-1) = '\n';
    *p = 0;

    if (level == ERR_FATAL && win_init)
	fprintf(stderr, "%s\n", bufp);

#ifdef _WIN32
    if (level == ERR_FATAL) {
    /* 7/1/99 johnt - log the messages to the Event Viewer on windows, as we don't always have stderr */
    char *a[] = { bufp };
    if( !eventLogHandle){
	eventLogHandle = RegisterEventSource(NULL,EVENT_SOURCE); /* get default application handle */
    }
    ReportEvent(eventLogHandle,
		level==ERR_FATAL?EVENTLOG_ERROR_TYPE:EVENTLOG_WARNING_TYPE,
		0,
		0,
		NULL,
		1,
		0,
		(LPCTSTR *)a,
		NULL);
    }
#endif

    tout_update_stream(2, bufp, 0, NULL);

    if (bufp != buf) {
	xfree(bufp);
    }

    return TCL_OK;
}

    
/*
 * Tcl interface to vfuncheader
 * Restrictions - cannot output more than 8K.
 */
int tcl_vfuncheader(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {

    if (argc != 2) {
	return TCL_ERROR;
    }

    funcheader(argv[1]);

    return TCL_OK;
}


/*
 * Tcl interface to vgroupheader
 * Restrictions - cannot output more than 8K.
 */
int tcl_vfuncgroup(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    if (argc != 3) {
	return TCL_ERROR;
    }

    funcgroup(atoi(argv[1]), argv[2]);

    return TCL_OK;
}

/*
 * Enables and disables the ringing of the bell when errors are displayed.
 * Essential for scripts!
 */
int tcl_error_bell(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " error_bell 0/1\"", (char *)NULL);
        return TCL_ERROR;
    }

    noisy = atoi(argv[1]);
    return TCL_OK;
}

/*
 *-----------------------------------------------------------------------------
 * C callable output routines
 *-----------------------------------------------------------------------------
 */

void UpdateTextOutput() {
    while (Tcl_DoOneEvent(TCL_WINDOW_EVENTS | TCL_IDLE_EVENTS | TCL_DONT_WAIT)
	   != 0)
	;
}

/*
 * Usage: verror(priority, name, format, args...);
 * NB: don't pass more than 8K per call
 */
__PRINTF_FORMAT__(3,4)
void verror(int priority, const char *name, const char *fmt, ...) {
    char buf[8192], tbuf[100], *bufp = buf;
    va_list args;
    time_t t = time(NULL);
    size_t l;
    static time_t last_time = 0;

    /* To improve error reporting */
    if (priority == ERR_FATAL && t - last_time > 10 && _interp)
	dump_tcl_stack();
    last_time = t;

    if (noisy) bell();
    fflush(stdout);

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((l = vflen(fmt, args)) > 8192 - sizeof(tbuf)+2) {
	if (NULL == (bufp = (char *)xmalloc(l + sizeof(tbuf)+2))) {
	    verror(ERR_FATAL, "verror", "out of memory");
	    return;
	}
    }

    strftime(tbuf, sizeof(tbuf)-1, "%a %d %b %H:%M:%S %Y", localtime(&t));
    sprintf(bufp, "%s %s: ", tbuf, name);

    if (priority == ERR_FATAL && win_init) {
	fputs(bufp, stderr);
	vfprintf(stderr, fmt, args);
	fputc('\n', stderr);
    }

    l = strlen(bufp) - strlen(name) - 2; /* "%s: ",name */
    vsprintf(&bufp[l], fmt, args);
    log_file(NULL, &bufp[l]);
    strcat(&bufp[l], "\n");
#ifdef _WIN32
    if (priority == ERR_FATAL) {
    /* 7/1/99 johnt - log the messages to the Event Viewer on windows, as we don't always have stderr */
    char *a[] = {bufp};
    if( !eventLogHandle){
	eventLogHandle = RegisterEventSource(NULL,EVENT_SOURCE); /* get default application handle */
    }
    ReportEvent(eventLogHandle,
		priority==ERR_FATAL?EVENTLOG_ERROR_TYPE:EVENTLOG_WARNING_TYPE,
		0,
		0,
		NULL,
		1,
		0,
		(LPCTSTR *)a,
		NULL);
    }
#endif
    tout_update_stream(2, bufp, 0, NULL);

    if (bufp != buf) {
	xfree(bufp);
    }

    va_end(args);
}

/*
 * Usage: vmessage(format, args...);
 */
__PRINTF_FORMAT__(1,2)
void vmessage(const char *fmt, ...) {
    char buf[8192], *bufp = buf;
    int len;
    va_list args;

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (bufp = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vmessage", "out of memory");
	    return;
	}
    }

    vsprintf(bufp, fmt, args);
    if (log_vmessage_st)
	log_file(NULL, bufp);
    tout_update_stream(1, bufp, 0, NULL);
    va_end(args);

    if (bufp != buf) {
	xfree(bufp);
    }
}

/*
 * Usage: vmessage_tagged(format, tag, args...);
 */
__PRINTF_FORMAT__(2,3)
void vmessage_tagged(const char *tag, const char *fmt, ...) {
    char buf[8192], *bufp = buf;
    int len;
    va_list args;

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (bufp = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vmessage", "out of memory");
	    return;
	}
    }

    vsprintf(bufp, fmt, args);
    if (log_vmessage_st)
	log_file(NULL, bufp);
    tout_update_stream(1, bufp, 0, tag);
    va_end(args);

    if (bufp != buf) {
	xfree(bufp);
    }
}

/*
 * Adds a new header to the text output window.
 */
__PRINTF_FORMAT__(1,2)
void vfuncheader(const char *fmt, ...) {
    char name[8192], *namep = name;
    va_list args;
    int len;

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (namep = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vfuncheader", "out of memory");
	    return;
	}
    }

    vsprintf(namep, fmt, args);
    funcheader(namep);
    va_end(args);

    if (namep != name)
	xfree(namep);
}


/*
 * Used for grouping outputs together (such as the 2D plot results).
 * Basically we don't output a header if the last output was from this
 * group and there haven't been function headers outputted since.
 *
 * group numbers:
 * 1	2D plot matches
 * 2	Information from template display
 */
__PRINTF_FORMAT__(2,3)
void vfuncgroup(int group, const char *fmt, ...) {
    char name[8192], *namep = name;
    va_list args;
    int len;

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (namep = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vfuncheader", "out of memory");
	    return;
	}
    }

    vsprintf(namep, fmt, args);	
    funcgroup(group, namep);
    va_end(args);

    if (namep != name)
	xfree(namep);
}

/*
 * Usage: vparams(format, args...);
 */
__PRINTF_FORMAT__(1,2)
void vfuncparams(const char *fmt, ...) {
    char params[8192], *paramsp = params;
    va_list args;
    int len;

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (paramsp = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vfuncheader", "out of memory");
	    return;
	}
    }

    vsprintf(paramsp, fmt, args);
    funcparams(paramsp);
    va_end(args);

    if (paramsp != params)
	xfree(paramsp);
}

/*
 *-----------------------------------------------------------------------------
 * Fortran callable output routines
 *-----------------------------------------------------------------------------
 */

f_proc_ret fverr_(f_int *priority, char *fname, char *fmess,
		  f_implicit name_l, f_implicit mess_l) {
    char mess[1024];
    char name[1024];

    Fstr2Cstr(fmess, mess_l, mess, 1024);
    Fstr2Cstr(fname, name_l, name, 1024);
    verror(*priority, name, "%s", mess);

    f_proc_return();
}

void updout_(void) {
    UpdateTextOutput();
}

/*
 * This doesn't really belong here, but this and the updout_ function used
 * to be together in nxspec.f.
 */
void initrs_(void) {
}

