/*
 * A simple program to generate a config.h file with machine
 * specifics, ie MACHINE, ENDIAN, etc
 */

#include <stdio.h>
#include <inttypes.h>
#include <sys/utsname.h>

int main(void) {
    struct utsname uts;
    union {
	int32_t i;
	char c[4];
    } u;

    printf("#ifndef _CONFIG_H_\n");
    printf("#define _CONFIG_H_\n\n");

    /* Endianness */
    u.i = 1;
    if (u.c[0] == 1) {
	printf("#define SP_LITTLE_ENDIAN\n\n");
    } else {
	printf("#define SP_BIG_ENDIAN\n\n");
    }
    
    /* System bits */
    if (-1 == uname(&uts))
	return -1;
    printf("#ifdef MACHINE\n");
    printf("#   undef MACHINE\n");
    printf("#endif\n");
    printf("#define MACHINE %s\n",   uts.machine);
    printf("#ifdef OSTYPE\n");
    printf("#   undef OSTYPE\n");
    printf("#endif\n");
    printf("#define OSTYPE  %s\n\n", uts.sysname);

    printf("#endif /* _CONFIG_H_ */\n");

    return 0;
}
