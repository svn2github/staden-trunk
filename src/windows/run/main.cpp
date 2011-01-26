 
/*
 * This program is designed to be compiled and put in to
 * $STADENROOT/bin.
 *
 * It sets up the environment and then runs wish.exe on
 * share/staden/tcl/PROG/PROG.tcl where PROG is the name of this
 * executable (copy it around to the various names it has to stand
 * for). This means that it can be used from both a cmd32.exe (DOS)
 * environment, the start menu or in a file association without
 * requiring any arguments. (It helps to give a feel more similar to
 * the unix systems too.)
 *
 * With the -console argument it instead invokes a command shell with the
 * environment set up. This allows for scripting and manual use of the command
 * line tools (vector_clip, extract_seq, convert_trace, etc...).
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include <errno.h>
#include <direct.h>
#include <string.h>
#include <ctype.h>
#include <windows.h>

// Buffers
const int MAXARG = 10;
const int MAXARGLEN = 1024;
static char argv[MAXARG][MAXARGLEN];

// Other data
HINSTANCE g_hInstance = 0;

//-------------------
// Parser for WinMain
//-------------------

unsigned int ParseCommandLine( char* pCmdLine, char argv[MAXARG][MAXARGLEN] )
{
   enum { STATE_INIT_ARGUMENTS, STATE_CONSUME_WHITESPACE, STATE_SCAN_CHARS,
          STATE_SCAN_QUOTED_CHARS, STATE_ADD_ARGUMENT };
   char   pBuffer[MAX_PATH];
   unsigned int j;
   char   c            = 1;
   unsigned int i      = 0;
   unsigned int n      = 1;
   unsigned int nState = STATE_INIT_ARGUMENTS;

   // Command Line Parser State Machine

   pBuffer[0] = 0;
   while(1)
   {
      switch( nState )
      {
         case STATE_INIT_ARGUMENTS:
            // Initialise argv[0] with filename, zero remainder
            ::GetModuleFileName( g_hInstance, pBuffer, MAX_PATH );
            strcpy( &argv[0][0], pBuffer );
            for( j=1; j<MAXARG; j++ )
               argv[j][0] = 0;
            nState = STATE_CONSUME_WHITESPACE;
            break;


         case STATE_CONSUME_WHITESPACE:
            // Eat through unwanted whitespace
            while( isspace(c=*pCmdLine++) );
            nState = STATE_SCAN_CHARS;
            pCmdLine--;
            break;


         case STATE_SCAN_CHARS:
            // Scan characters for whitespace/quotes
            c = *pCmdLine++;
            pBuffer[i++] = c;
            if( isspace(c) ) {
               nState = STATE_ADD_ARGUMENT;
	       i--;
	    }
            else if( c == '"' )
               nState = STATE_SCAN_QUOTED_CHARS;
            break;


         case STATE_SCAN_QUOTED_CHARS:
            // Fast forward until next set of quotes
            c = *pCmdLine++;
            pBuffer[i++] = c;
            if( c == '"' )
                nState = STATE_ADD_ARGUMENT;
            break;


         case STATE_ADD_ARGUMENT:
            // Add another argument to argv[]
	    if (i) {
                pBuffer[i] = 0;
                strcpy( argv[n++], pBuffer );
	    }
            nState = STATE_CONSUME_WHITESPACE;
            pBuffer[0] = 0;
            i = 0;
            break;
      }

      // Exit test
      if( !c )
      {
         if( strlen(pBuffer) == 0 )
            break;
         nState = STATE_ADD_ARGUMENT;
      }
   }
   return n;
}



//------
// Main
//------

int Main( int argc, char argv[MAXARG][MAXARGLEN] )
{
    int   k;
    int   n;
    char* p;
    int   status;
    char* argp[MAXARG+1];
    char  buffer[8192];
    char  rootdir[MAX_PATH];
    char  our_cmd[MAX_PATH];
    char  console[MAX_PATH];
    char  unix_rootdir[MAX_PATH];
    char  basename[MAX_PATH];
    char  arg1[MAX_PATH];

    // Get the location of the executable name
    if( !::GetModuleFileName(0,rootdir,MAX_PATH) )
    {
        ::MessageBox(0,"Unable to get Staden Package directory!", "SPRUN.EXE Message", MB_OK );
        return -3;
    }
    _strlwr( rootdir ); // lowercase

    // Lop off the .exe and the bin to get the Staden Package
    // root directory.
    n = strlen(rootdir)-1;
    while (n && rootdir[n] != '\\' && rootdir[n] != '/')
    	n--;
    strcpy(basename, &rootdir[n+1]);
    if (p = strstr(basename, ".exe"))
	*p = 0;
    if (n) n--;
    while (n && rootdir[n] != '\\' & rootdir[n] != '/')
    	n--;
    rootdir[n] = 0;

   // Create a copy of root directory in unix format (fwd slashes)
   strcpy( unix_rootdir, rootdir );
   p = unix_rootdir;
   while( *p )
   {
      if( *p == '\\' )
        *p = '/';
      p++;
   }


    // Add Staden Package environment variables to current environment
    sprintf( buffer, "TK_LIBRARY=%s/lib/tk8.4", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "TCL_LIBRARY=%s/lib/tcl8.4", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "STADLIB=%s/lib/staden", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "STADTCL=%s/share/staden/tcl", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "STADTABL=%s/share/staden/etc", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "STADENROOT=%s", unix_rootdir );
    _putenv( buffer );
    //    _putenv( "MACHINE=windows" );



    // Prepend our paths to PATH environment variable. Normally this is done
    // automatically by the App Path registry entry for sprun.exe, however if
    // we want to run more a remote copy of gap4 then the App Path will be
    // invalid, so we need to do it here again just to be sure we pick up the
    // appropriate dlls.
    sprintf( buffer, "PATH=%s\\bin;%s/lib/staden;", rootdir, rootdir );
    p = getenv( "PATH" );
    assert(p);
    strcat( buffer, p );
    _putenv( buffer );
    


    // Set default command to be winstash, don't use quotes otherwise execve
    // won't work properly, but quotes are required in argp[0]!
    strcpy( our_cmd, rootdir );
    strcat( our_cmd, "\\bin\\wish84.exe" );


    // Get command console
    p = getenv( "COMSPEC" );
    if(!p)
    {
        // Assume windows NT
        ::GetWindowsDirectory( console, MAX_PATH );
        strcat( console, "\\system32\\cmd.exe" );
    }
    else
    {
        // Console is specified by environment
        strcpy( console, p );
    }
    


    // Construct suitable argp array for execve
    k = 1;
    n = 1;
    // -console option runs the console with the staden environment intact
    // on windows NT platforms only. Command.com on win9x is crippled
    // in that it ignores the environment given to execve()
    if (n < argc && strcmp(argv[n], "-console") == 0) {
	strcpy(our_cmd, console);
	n++;
    }
    if (strcmp(basename, "sprun") == 0) {
	// Just wish by itself
        argp[1] = 0;
    } else {
	// Otherwise start with a tcl file based on argv[0] and append
	// all other arguments
        sprintf(arg1, "\"%s\\share\\staden\\tcl\\%s\\%s.tcl\"", rootdir, basename, basename);
	argp[k++] = arg1;

        while( n<argc )
            argp[k++] = &argv[n++][0];
        argp[k] = 0;
        assert(k<(MAXARG+1));
    }

    // Assemble argv[0], we must put quotes around this
    strcpy( buffer, "\"" );
    strcat( buffer, our_cmd );
    strcat( buffer, "\"" );
    argp[0] = buffer;

    for (n = 0; argp[n]; n++) {
	printf("'%s' ", argp[n]);
    }
    printf("\n");

    // Execute the program
    status = _execve( our_cmd, argp, _environ );
    if( status == -1 )
    {
        sprintf( buffer, "Sorry: Unable to execute %s - (%d).", our_cmd, errno );
        ::MessageBox(0,buffer,"SPRUN.EXE Message", MB_OK );
    }
    return 0;
}




//---------------------
// Windows Entry Point
//---------------------

int WINAPI WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow )
{
   // Global variable initialisation
   g_hInstance = hInstance;

   // Process the command line arguments
   int argc = ParseCommandLine( lpCmdLine, argv );

   // Call user entry point
   return Main( argc, argv );
}

