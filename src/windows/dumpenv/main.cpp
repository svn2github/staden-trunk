#include <stdio.h>
#include <stdlib.h>
#include <windows.h>


char Buffer[32768];


int WINAPI WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance,
                    LPSTR lpCmdLine, int nCmdShow )
{
	int n = 0;
	strcpy( Buffer, "Short Command Line: " );
   strcat( Buffer, lpCmdLine );
   strcat( Buffer, "\n" );
   strcat( Buffer, "Long  Command Line: " );
   strcat( Buffer, GetCommandLine() );
	strcat( Buffer, "\n");
   while(1)
   {
      if( !_environ[n] )
         break;
      else
		{
         strcat( Buffer, _environ[n] );
			strcat( Buffer, "\n" );
		}
      n++;
   }
	::MessageBox( 0, Buffer, "Environment Dump", MB_OK );
   return 0;
}
 