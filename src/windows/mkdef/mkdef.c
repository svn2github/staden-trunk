/*
   mkdef.c
   *******

*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#define MAXLINE 256

int Usage( char* message)
{
   fprintf(stderr,"mkdef: %s\nUsage:\n\n",message);

   fprintf(stderr,"mkdef -c objfile.o\n\n");

   fprintf(stderr,"mkdef -l output.dll in1.o in2.o in3.o in4.o ...\n\n");

   fprintf(stderr,"With -c flag: use bindump to extract public symbols from specified object file and create a partial .def file containing the definitions.\n\n");

   fprintf(stderr,"With -l flag: create a full .def file by contatenating def files for supplied .o files, and adding a header. Uses basename on the output dll for the library name, and creates output.def.\n\n");

   return 1;
}

/*
 * partdef
 *********
 * produce a partial module definition file, by
 * running dumpbin on an object file
 */

int partdef( int argc, char *argv[]){
	char *objname = argv[2];
	char *defname; 
	char *dot ;
	char *command;
	char *bindump = "dumpbin /symbols ";
	char line[MAXLINE];
	char *ptr;
	FILE *fpipe,*fdef;

	/* build the def filename */
	if(!(dot=strrchr(objname,'.')))
		dot=objname+strlen(objname);
        defname = malloc( (dot-objname)+5);
	strncpy(defname,objname,dot-objname);
	strcpy(defname+(dot-objname),".def");

	/*
	 * create the pipe to get the public symbols
	 */
	/* build the pipe command line */
	command = malloc(strlen(bindump)+strlen(objname)+2);
	strcpy(command,bindump);
	strcat(command,objname);

	/* create the pipe */
	if( !(fpipe=_popen(command,"r")) ){
		perror("Unable to create pipe to bindump");
		return 1;
	}

	/* open the output file */
	if( !(fdef=fopen(defname,"w")) ){
		sprintf(line,"Unable to create output file \"%s\"",defname);
		perror(line);
		fclose(fpipe);
		return 1;
	}
	printf("%s:\n",defname);

	/* process each line of the input, and produce the DEF file */
	while( !(feof(fpipe)||ferror(fpipe)))
	{
		if( fgets(line,MAXLINE-1,fpipe))
		{
			if(! strstr(line,"UNDEF") )
			{
				if( (ptr=strstr(line,"External")))
				{
					if( (ptr=strstr(line,"| _")))
					{
						if( !strstr(ptr,"@") && !strstr(ptr,"?") )
						{
							ptr+=3;
							printf("exporting %s",ptr);
							fprintf(fdef,"\t%s",ptr);
						}
					}
				}
			}
		}
	}

	fclose(fdef);
	fclose(fpipe);

	return 0;
}
	
/*
 * fulldef
 *********
 * produce a full module definition file
 * by concatenating the partial def files belonging
 * to the specified object files
 */

int fulldef(int argc, char *argv[]){

	char *outname = argv[2];
	char *outdefname;
	char *libraryname;
	char *dot;
	char *start;
	FILE *fout, *fin;
	char line[MAXLINE];
	int stat=0;
	int i;

	/* build the output filename and libraryname */
	if(!(dot=strrchr(outname,'.')))
		dot=outname+strlen(outname);
	for(start=dot;start!=outname;start--)
		if( (*start == '/') || (*start == '\\') || (*start == ':' )){
			start++;
			break;
		}
	outdefname = malloc( (dot-outname)+5);
	libraryname = malloc(dot-start);
	strncpy(outdefname,outname,dot-outname);
	strcpy(outdefname+(dot-outname),".def");
	strncpy(libraryname,start,dot-start);
	*(libraryname+(dot-start))='\0';

	/* create the output file */
	if( !(fout=fopen(outdefname,"w"))){
		sprintf(line,"Unable to open \"%s\"",outdefname);
		perror(line);
		return 1;
	}

	printf("%s:\n",outdefname);

	/* write header */
	fprintf(fout,"LIBRARY\t%s\n",libraryname);
	fprintf(fout,"EXPORTS\n");

	/* get the exports from the def file for each input file */
	for(i=3;i<argc;i++){
		char *objname = argv[i];
		char *defname;

		/* build the input def filename from the object filename */
		if(!(dot=strrchr(objname,'.')))
			dot=objname+strlen(objname);
		defname = malloc( (dot-objname)+5);
		strncpy(defname,objname,dot-objname);
		strcpy(defname+(dot-objname),".def");

		/* open the file */
		if( !(fin=fopen(defname,"r"))){
			sprintf(line,"Warning: Unable to open input file \"%s\"",defname);
			perror(line);
		}
		else{
			while(!(ferror(fin) || feof(fin))){
				if( fgets(line,MAXLINE-1,fin)){
					printf("%s",line);
					fputs(line,fout);
				}
			}
			fclose(fin);
		}
		free(defname);
	}

	fclose(fout);

	return stat;

}

int main(int argc,char *argv[])
{

	if( argc < 3 )
		return Usage("Not enough parameters");

	if( !stricmp(argv[1],"-c")){
		return partdef(argc,argv);
	}
	if( !stricmp(argv[1],"-l")){
		return fulldef(argc,argv);
	}
	return Usage("Must specify -c or -l");
}