#if !defined(SPLIT_H)
#define SPLIT_H

#include <stdio.h>

char **split(char *s, char *ct);
void split_xfree(char **token);

#endif

