#include <stdlib.h>
#include <string.h>
#include "misc.h"
#include "split.h"

char **split(char *s, char *ct) {
    char	*copy, *temp;
    char	**token;
    int		i;

    copy = strdup(s);
    if(NULL == (token = (char **) xmalloc(strlen(copy) * sizeof(char *)))) {
	xfree(copy);
	return NULL;
    }

    i = 0;
    temp = (char *) strtok(copy, ct);
    while(temp) {
	token[i] = strdup(temp);
	i++;
	temp = (char *) strtok(NULL, ct);
    }
    if(NULL == (token = (char **) xrealloc(token, (i + 1) * sizeof(char *) + 1))) {
	xfree(copy);
	return NULL;
    }
    token[i] = NULL;

    xfree(copy);
    return token;
}

void split_xfree(char **token) {
    int		i;

    i = 0;
    while(token[i]) {
	xfree(token[i]);
	token[i] = NULL;
	i++;
    }
    xfree(token);
}
