/* 
 * tkAppInit.c --
 *
 *	Provides a default version of the Tcl_AppInit procedure for
 *	use in wish and similar Tk-based applications.
 *
 * Copyright (c) 1993 The Regents of the University of California.
 * Copyright (c) 1994 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include <stdlib.h>		/* For getenv(), _putenv(), malloc(), free(), ... */
#include <string.h>		/* For strlen(), strcat(), strcpy() */
#include <fcntl.h>
#include "tk.h"
#include "tclXkeylist.h"


#ifdef _WIN32
#include <windows.h>
#include <malloc.h>
#include <locale.h>
#endif

#ifdef TRAP_SIGNALS
#include <signal.h>
#include "gap-error.h"
#endif

/*
 * The following variable is a special hack that is needed in order for
 * Sun shared libraries to be used for Tcl.
 */

#ifdef NEED_MATHERR
extern int matherr();
int *tclDummyMathPtr = (int *) matherr;
#endif


void Tk_Utils_Main( int, char**, Tcl_AppInitProc*, Tcl_Interp* );

#ifdef _WIN32
static int needConsole = FALSE;	/* so Tcl_AppInit knows to Init the console */


/*
 * The following declarations refer to internal Tk routines.  These
 * interfaces are available for use, but are not supported.
 */

EXTERN void Tk_InitConsoleChannels( Tcl_Interp *interp );
EXTERN int  Tk_CreateConsoleWindow( Tcl_Interp *interp );

/*
 * Forward declarations for procedures defined later in this file:
 */

static void  setargv _ANSI_ARGS_((int *argcPtr, char ***argvPtr));
static void  WishPanic _ANSI_ARGS_(TCL_VARARGS(char *,format));
static void  WinEnvToUnixEnv( void );
static void  WinSlashToUnixSlash( char* s );
static void  WinEnvVarToUnixEnvVar( const char* pVarName );



/*
 *----------------------------------------------------------------------
 *
 * WinMain --
 *
 *	Main entry point from Windows.
 *
 * Results:
 *	Returns false if initialization fails, otherwise it never
 *	returns. 
 *
 * Side effects:
 *	Just about anything, since from here we call arbitrary Tcl code.
 *
 *----------------------------------------------------------------------
 */

int APIENTRY
WinMain(hInstance, hPrevInstance, lpszCmdLine, nCmdShow)
    HINSTANCE hInstance;
    HINSTANCE hPrevInstance;
    LPSTR lpszCmdLine;
    int nCmdShow;
{
    int   n;
    char  **argv, *p; 
    int   argc;
    char  buffer[MAX_PATH];
    Tcl_SetPanicProc(WishPanic);


    /*
     * Set up the default locale to be standard "C" locale so parsing
     * is performed correctly.
     */
    setlocale(LC_ALL, "C");



    /*
     * Increase the application queue size from default value of 8.
     * At the default value, cross application SendMessage of WM_KILLFOCUS
     * will fail because the handler will not be able to do a PostMessage!
     * This is only needed for Windows 3.x, since NT dynamically expands
     * the queue.
     */
    SetMessageQueue(64);



    /*
     *  On windows, we don't use any switches or anything fancy, so
     *  here we just convert all slashes in one hit to unix style.
     *  If we don't do this, GAP4 will complain about invalid file
     *  names etc.
     */
    setargv(&argc, &argv);
    n = 1;
    while( n<argc )
    	WinSlashToUnixSlash( argv[n++] );



    /*
     * Via this route, we want a console.
     */
    needConsole = 1;


 
    /*
     * Replace argv[0] with full pathname of executable, and forward
     * slashes substituted for backslashes.
     */
    GetModuleFileName(NULL, buffer, sizeof(buffer));
    argv[0] = buffer;
    WinSlashToUnixSlash( buffer );
    WinEnvToUnixEnv();


          
    /*
     * Invoke TCL interpreter
     */
    Tk_Utils_Main(argc, argv, Tcl_AppInit, Tcl_CreateInterp());
    return 1;
}

int main(int argc, char** argv)
{
    char *p;
    char buffer[MAX_PATH];
    Tcl_SetPanicProc(WishPanic);


    /*
     * Set up the default locale to be standard "C" locale so parsing
     * is performed correctly.
     */
    setlocale(LC_ALL, "C");


    /*
     * Increase the application queue size from default value of 8.
     * At the default value, cross application SendMessage of WM_KILLFOCUS
     * will fail because the handler will not be able to do a PostMessage!
     * This is only needed for Windows 3.x, since NT dynamically expands
     * the queue.
     */
    SetMessageQueue(64);


    /*
     * Replace argv[0] with full pathname of executable, and forward
     * slashes substituted for backslashes.
     */
    GetModuleFileName(NULL, buffer, sizeof(buffer));
    argv[0] = buffer;
    WinSlashToUnixSlash( buffer );
    WinEnvToUnixEnv();


    Tk_Utils_Main(argc, argv, Tcl_AppInit, Tcl_CreateInterp());
    return 1;
}

#else

#define WishPanic(x) 

/*
 *----------------------------------------------------------------------
 *
 * main --
 *
 *	This is the main program for the application.
 *
 * Results:
 *	None: Tk_Main never returns here, so this procedure never
 *	returns either.
 *
 * Side effects:
 *	Whatever the application does.
 *
 *----------------------------------------------------------------------
 */

int
main(int argc, char **argv)
{
#if defined(__sun__) && !defined(__svr4__)
    /*
     * SunOS 4 dlopen call is broken. This is fixed either by applying patch
     * 101783 or by opening /dev/zero as file descriptor 4.
     * We use the latter method as the patch only fixes it for the local
     * system, meaning that everyone using the programs would have to
     * apply the patch.
     * The bug is in SunOS4.1.1, 4.1.2, 4.1.3, but not 4.1.3U1.
     */
    {
	int fd = open("/dev/zero", O_RDWR);
	close(4); /* Just in case - although if this works it'll break
		     something else anyway. */
	dup2(fd, 4);
	close(fd);
    }
#endif

    Tk_Utils_Main(argc, argv, Tcl_AppInit, Tcl_CreateInterp());
    return 0;			/* Needed only to prevent compiler warning. */
}

#endif

/*
 *----------------------------------------------------------------------
 *
 * Tcl_AppInit --
 *
 *	This procedure performs application-specific initialization.
 *	Most applications, especially those that incorporate additional
 *	packages, will have their own version of this procedure.
 *
 * Results:
 *	Returns a standard Tcl completion code, and leaves an error
 *	message in interp->result if an error occurs.
 *
 * Side effects:
 *	Depends on the startup script.
 *
 *----------------------------------------------------------------------
 */


extern int Tk_utils_Init(Tcl_Interp *interp);

int tkinit(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
#ifdef _WIN32
    if( needConsole)
	/* Tk_Init will already have been called by Tcl_AppInit */
	return TCL_OK;
    else if (Tk_Init(interp) == TCL_ERROR) {
	WishPanic(interp->result);
	return TCL_ERROR;
    }
    return TCL_OK;
#else
    return Tk_Init(interp);
#endif
}

#ifdef TRAP_SIGNALS
static void error_sig(int sig) {
	char *err = "Program terminated unexpectedly with signal %d.\nThis is probably a bug.\nPlease email all bug reports to staden-package@mrc-lmb.cam.ac.uk.\n";
	char message[200];
	sprintf(message,err,sig);
#ifdef _WIN32
    MessageBeep(MB_ICONEXCLAMATION);
    MessageBox(NULL, buf, "Fatal Error in stash",
	    MB_ICONSTOP | MB_OK | MB_TASKMODAL | MB_SETFOREGROUND);
    DebugBreak();	
    ExitProcess(1);
#else
    fprintf(stderr,message);
    abort();
#endif

}
#endif

int
Tcl_AppInit(Tcl_Interp *interp)
{
   char *lib, buf[1025];

    /*
     * Redefine TCL_LIBRARY and TK_LIBRARY to the staden package versions
     * in $STADLIB/{tcl,tk}.
     */
#ifndef _WIN32 /* 11/1/98 johnt - not required for WIN32 */
    if (NULL != (lib = getenv("STADLIB"))) {
        sprintf(buf, "TCL_LIBRARY=%s/tcl", lib);
        Tcl_PutEnv(buf);
        sprintf(buf, "TK_LIBRARY=%s/tk", lib);
        Tcl_PutEnv(buf);
    }
#endif

#ifdef TRAP_SIGNALS
    /*
     * Use the BSD signal() command to trap probable program crashes.
     * This then adds a debug message to make sure that we tell people
     * to email us.
     */
#if defined(SIGBUS) /* 11/1/99 johnt - SIGBUS not defined under WINNT */
    signal(SIGBUS,  error_sig);
#endif
    signal(SIGSEGV, error_sig);
    signal(SIGILL,  error_sig);
    signal(SIGFPE,  error_sig);
#if defined(SIGSYS)
    signal(SIGSYS,  error_sig);
#endif /* SIGSYS */
#endif /* TRAP_SIGNALS */

    if (Tcl_Init(interp) == TCL_ERROR) {
	WishPanic(interp->result);
	return TCL_ERROR;
    }

#ifdef _WIN32

    if( needConsole ){
	/* running in Windows mode, so initialise TK, and Init console */
	if (Tk_Init(interp) == TCL_ERROR) {
	    WishPanic(interp->result);
	    return TCL_ERROR;
	}
	Tcl_StaticPackage(interp, "Tk", Tk_Init, Tk_SafeInit);

	Tk_InitConsoleChannels(interp);

        if (Tk_CreateConsoleWindow(interp) == TCL_ERROR) {
	    WishPanic(interp->result);
	    return TCL_ERROR;
	}
    }

#endif


    Tcl_CreateCommand(interp, "tkinit", tkinit,
		      (ClientData)NULL, NULL);

    /*
     * Call the init procedures for included packages.  Each call should
     * look like this:
     *
     * if (Mod_Init(interp) == TCL_ERROR) {
     *     return TCL_ERROR;
     * }
     *
     * where "Mod" is the name of the module.
     */

    /* Library init routines */
    if (Tk_utils_Init(interp) == TCL_ERROR) {
	WishPanic(interp->result);
	return TCL_ERROR;
    }

    /*
     * Call Tcl_CreateCommand for application-specific commands, if
     * they weren't already created by the init procedures called above.
     */

    /*
     * Specify a user-specific startup file to invoke if the application
     * is run interactively.  Typically the startup file is "~/.apprc"
     * where "app" is the name of the application.  If this line is deleted
     * then no user-specific startup file will be run under any conditions.
     */

    Tcl_SetVar(interp, "tcl_rcFileName", "~/.stashrc", TCL_GLOBAL_ONLY);
    return TCL_OK;
}

#ifdef _WIN32

void
WishPanic TCL_VARARGS_DEF(char *,arg1)
{
    va_list argList;
    char buf[1024];
    char *format;
    
    format = TCL_VARARGS_START(char *,arg1,argList);
    vsprintf(buf, format, argList);

    MessageBeep(MB_ICONEXCLAMATION);
    MessageBox(NULL, buf, "Fatal Error in Wish",
	    MB_ICONSTOP | MB_OK | MB_TASKMODAL | MB_SETFOREGROUND);
    DebugBreak();
    ExitProcess(1);
}



/*
 *-------------------------------------------------------------------------
 *
 * setargv --
 *
 *	Parse the Windows command line string into argc/argv.  Done here
 *	because we don't trust the builtin argument parser in crt0.  
 *	Windows applications are responsible for breaking their command
 *	line into arguments.
 *
 *	2N backslashes + quote -> N backslashes + begin quoted string
 *	2N + 1 backslashes + quote -> literal
 *	N backslashes + non-quote -> literal
 *	quote + quote in a quoted string -> single quote
 *	quote + quote not in quoted string -> empty string
 *	quote -> begin quoted string
 *
 * Results:
 *	Fills argcPtr with the number of arguments and argvPtr with the
 *	array of arguments.
 *
 * Side effects:
 *	Memory allocated.
 *
 *--------------------------------------------------------------------------
 */

static void
setargv(argcPtr, argvPtr)
    int *argcPtr;		/* Filled with number of argument strings. */
    char ***argvPtr;		/* Filled with argument strings (malloc'd). */
{
    char *cmdLine, *p, *arg, *argSpace;
    char **argv;
    int argc, size, inquote, copy, slashes;
    
    cmdLine = GetCommandLine();


    /*
     * Precompute an overly pessimistic guess at the number of arguments
     * in the command line by counting non-space spans.
     */

    size = 2;
    for (p = cmdLine; *p != '\0'; p++) {
	if (isspace(*p)) {
	    size++;
	    while (isspace(*p)) {
		p++;
	    }
	    if (*p == '\0') {
		break;
	    }
	}
    }
    argSpace = (char *) ckalloc((unsigned) (size * sizeof(char *) 
	    + strlen(cmdLine) + 1));
    argv = (char **) argSpace;
    argSpace += size * sizeof(char *);
    size--;

    p = cmdLine;
    for (argc = 0; argc < size; argc++) {
	argv[argc] = arg = argSpace;
	while (isspace(*p)) {
	    p++;
	}
	if (*p == '\0') {
	    break;
	}

	inquote = 0;
	slashes = 0;
	while (1) {
	    copy = 1;
	    while (*p == '\\') {
		slashes++;
		p++;
	    }
	    if (*p == '"') {
		if ((slashes & 1) == 0) {
		    copy = 0;
		    if ((inquote) && (p[1] == '"')) {
			p++;
			copy = 1;
		    } else {
			inquote = !inquote;
		    }
                }
                slashes >>= 1;
            }

            while (slashes) {
		*arg = '\\';
		arg++;
		slashes--;
	    }

	    if ((*p == '\0') || (!inquote && isspace(*p))) {
		break;
	    }
	    if (copy != 0) {
		*arg = *p;
		arg++;
	    }
	    p++;
        }
	*arg = '\0';
	argSpace = arg + 1;
    }
    argv[argc] = NULL;

    *argcPtr = argc;
    *argvPtr = argv;
}



/*
 * WinSlashToUnixSlash()
 *
 * Converts all back slashes to forward slashes.
 * 
 */

static void WinSlashToUnixSlash( char* s )
{
    if( s )
    {
	while( *s )
	{
	    if( *s == '\\' )
		*s = '/';
	    s++;
	}
    }
}



/*
 * WinEnvVarToUnixEnvVar()
 *
 * Converts a windows style environment variable to a unix style one with
 * back slashes replaced with forward slashes.
 *
 */

static void WinEnvVarToUnixEnvVar( const char* pVarName )
{
    char* pBuffer;
    char* pVariable;


    /* Get existing variable */
    pVariable = getenv( pVarName );
    if( pVariable )
    {
	/* Allocate some storage, +2 for NULL and '=' characters */
	pBuffer = (char*) malloc( (size_t) (strlen(pVariable)+strlen(pVarName)+2) );
        if( pBuffer )
	{
	    /* Assemble the environment variable */
	    strcpy( pBuffer, pVarName );
	    strcat( pBuffer, "=" );
	    strcat( pBuffer, pVariable );


	    /* Do slash conversion */
	    WinSlashToUnixSlash( pBuffer );


	    /* Replace old environment variable */
	    _putenv( pBuffer );


	    /* Cleanup */
	    free( pBuffer );
	}
    }
}



/*
 * WinEnvToUnixEnv()
 *
 * Converts all the Staden package environment variables to the unix style.
 *
 */

static void WinEnvToUnixEnv( void )
{
    WinEnvVarToUnixEnvVar( "STADLIB" );
    WinEnvVarToUnixEnvVar( "STADTABL" );
    WinEnvVarToUnixEnvVar( "STADENROOT" );
    WinEnvVarToUnixEnvVar( "TK_LIBRARY" );
    WinEnvVarToUnixEnvVar( "TCL_LIBRARY" );
    WinEnvVarToUnixEnvVar( "TEMP" );
}


#endif
