#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include <errno.h>
#include <direct.h>
#include <string.h>
#include <windows.h>



// Buffers
const int MAXARG = 10;
const int MAXARGLEN = 1024;
static char argv[MAXARG][MAXARGLEN];


// Other data
HINSTANCE g_hInstance = 0;



#ifdef RUNTIME_FONT_INSTALL

// Fonts
typedef struct
{
    int  yn;
    char FaceName[LF_FACESIZE];

}FONTINFO;



//-------------------
// Font Installation
//-------------------

int CALLBACK EnumFontProc( ENUMLOGFONTEX* lpelfe, NEWTEXTMETRICEX* lpntme, DWORD FontType, LPARAM lParam )
{
    FONTINFO* fi;
    if( FontType & TRUETYPE_FONTTYPE )
    {
        fi = (FONTINFO*) lParam;
        if( strcmp((const char*)lpelfe->elfFullName,fi->FaceName) == 0 )
        {
            fi->yn = 1;
            return 0;
        }
    }
    return 1;
}

bool IsFontInstalled( const char* pFaceName )
{
    assert(pFaceName);
    assert(strlen(pFaceName)<LF_FACESIZE);
    LOGFONT  lf;
    FONTINFO fi;
    fi.yn = 0;
    strcpy( fi.FaceName, pFaceName );
    lf.lfCharSet        = DEFAULT_CHARSET;
    lf.lfPitchAndFamily = 0;
    strcpy( lf.lfFaceName, pFaceName );
    EnumFontFamiliesEx( GetDC(0), &lf, (FONTENUMPROC)EnumFontProc, (LPARAM)&fi, 0 );
    return fi.yn ? true : false;
}

#endif /* RUNTIME_FONT_INSTALL */




//--------
// Parser
//--------

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
            if( isspace(c) )
               nState = STATE_ADD_ARGUMENT;
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
            pBuffer[i] = 0;
            strcpy( argv[n++], pBuffer );
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




    // Get the Staden Package directory
    if( !::GetModuleFileName(0,rootdir,MAX_PATH) )
    {
        ::MessageBox(0,"Unable to get Staden Package directory!", "SPRUN.EXE Message", MB_OK );
        return -3;
    }
    _strlwr( rootdir );
    p = strstr( rootdir, "\\sprun.exe" );
   *p = 0;



   // Create a copy of root directory in unix format
   strcpy( unix_rootdir, rootdir );
   p = unix_rootdir;
   while( *p )
   {
      if( *p == '\\' )
        *p = '/';
      p++;
   }



#ifdef RUNTIME_FONT_INSTALL
    /* Install the pregap font if it isn't already installed */
    if( IsFontInstalled("Pregap") == false )
    {
        //MessageBox(0,"Installing pregap.ttf","Debug",MB_OK);
        strcpy( buffer, rootdir );
        strcat( buffer, "\\windows-bin\\pregap.ttf" );
        ::GetShortPathName( buffer, buffer, MAX_PATH );
        _strlwr( buffer );
        if( !AddFontResource(buffer) )
            ::MessageBox(0,"Warning: Failed to install font pregap.ttf.","SPRUN.EXE Message", MB_OK );
        else
            ::SendMessage( HWND_BROADCAST, WM_FONTCHANGE, 0, 0 );
    }
#endif



    // Add Staden Package environment variables to current environment
    sprintf( buffer, "TK_LIBRARY=%s/lib/tk", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "TCL_LIBRARY=%s/lib/tcl", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "STADLIB=%s/lib", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "STADTABL=%s/tables", unix_rootdir );
    _putenv( buffer );
    sprintf( buffer, "STADENROOT=%s", unix_rootdir );
    _putenv( buffer );
    _putenv( "MACHINE=windows" );



    // Prepend our paths to PATH environment variable. Normally this is done
    // automatically by the App Path registry entry for sprun.exe, however if
    // we want to run more a remote copy of gap4 then the App Path will be
    // invalid, so we need to do it here again just to be sure we pick up the
    // appropriate dlls.
    sprintf( buffer, "PATH=%s\\windows-bin;%s\\lib\\windows-binaries;", rootdir, rootdir );
    p = getenv( "PATH" );
    assert(p);
    strcat( buffer, p );
    _putenv( buffer );
    


    // Set default command to be winstash, don't use quotes otherwise execve
    // won't work properly, but quotes are required in argp[0]!
    strcpy( our_cmd, rootdir );
    strcat( our_cmd, "\\windows-bin\\wish.exe" );



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
    while( n<argc )
    {
        // -c option runs the console with the staden environment intact
        // on windows NT platforms only. Command.com on win9x is crippled
        // in that it ignores the environment given to execve()
        if( (argv[n][0]=='-') && (argv[n][1]=='c') && !argv[n][2] )
            strcpy( our_cmd, console );
        else
        {
            argp[n] = &argv[n][0];
            k++;
        }
        n++;
    }
    argp[k] = 0;
    assert(k<(MAXARG+1));



    // Assemble argv[0], we must put quotes around this
    strcpy( buffer, "\"" );
    strcat( buffer, our_cmd );
    strcat( buffer, "\"" );
    argp[0] = buffer;



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

