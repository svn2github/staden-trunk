#include <stdio.h>
#include <ctype.h>
/*
** SD 2 August 1991
**   Changes way of interpreting uncertainty codes so that
**   we now only generate C A G T and -
*/
main (int argc, char **argv)
{
    char c ;
    int i = 0 ;

    if (argc != 2 || *argv[1] != '-') {
	fprintf (stderr,"Usage: '2rs -form' : form is abi or alf; a filter\n") ;
	exit (1) ;
    }

    if (!strcmp (&argv[1][1],"abi"))
	while ((c = getc (stdin)) != EOF) {
	    switch (c) {
		case 'N' : c = '-' ; break ;
	    }
	    putc (c,stdout) ;
	    if (!(++i%50))
		putc ('\n',stdout) ;
	}
      else if (! strcmp (&argv[1][1],"alf"))
    	/* ALF lower case uncertainty codes mean that the base may be missing.
	   RS uncertain length codes mean that there may be an extra base.
	   So I have to delay output by one character. Use oldc for this.
	*/
	while ((c = getc (stdin)) != EOF) { 
	    switch (c) {
	    case 'A':
	    case 'C':
	    case 'G':
	    case 'T':
		break;
	    default:
		if (isupper(c))
		    c = '-';
		else
		    c = '\0';
		break;
	    }
	    if (c) {
		putc (c,stdout);
		if (!(++i%50))
		    putc ('\n',stdout) ;
	    }
	}

    putc ('\n',stdout) ;
}
