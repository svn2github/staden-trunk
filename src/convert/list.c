#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include "list.h"

static char *mystrdup(char *s)
/*
 * A quick implementation of strdup()
 */
{
    char *copy;
    if ( (copy = (char *)malloc(strlen(s)+1)) != NULL ) strcpy(copy,s);
    return copy;
}

Node *nodeList = nil;

void destroy_node_list(void)
{
    Node *n,*m;

    for(n=nodeList;!isNil(n);) {
	m = cdr(n);
	free(n);
	n = m;
    }

}

static Node *create_node (int type)
{
    Node *n;

    if (isNil(nodeList)) {
	n = (Node *) malloc(sizeof(Node));
    } else {
	n = nodeList;
	nodeList = cdr(nodeList);
    }
    if (! isNil(n)){
	n->type = type;
	n->val.list.head = nil;
	n->val.list.tail = nil;
    }

    return n;
}


static void free_node (Node *n)
{

    if (isNil(n))
	;
    else {
	n->type = Node_List; /* entries on the nodeList are ALWAYS lists */
	n->val.list.head = nil;
	n->val.list.tail = nodeList;
	nodeList = n;
    }

}


#define create_node_list() create_node(Node_List)
#define create_node_atom() create_node(Node_Atom)

void destroy_list (Node *n)
{
    
    if (isNil(n))
	;
    else if (isAtom(n)){
	free(atomVal(n)); /* free the atom here */
        free_node(n);
    } else {
	destroy_list(car(n));
	destroy_list(cdr(n));
	free_node(n);
    }
}



Node *copy_list(Node *n)
{
    Node *new;

    if (isNil(n))
	new = nil;
    else if (isAtom(n)){
	new = create_node_atom();
	new->val.atom.ptr = mystrdup(atomVal(n));
    } else {
	new = create_node_list();
	new->val.list.head = copy_list(car(n));
	new->val.list.tail = copy_list(cdr(n));
    }

    return new;
}


Node *atom_str(char *s)
{
    Node *new;

    new = create_node_atom();
    new->val.atom.ptr = mystrdup(s);

    return new;
}


Node *atom_int(int i)
{
    char buffer[200];

    sprintf(buffer,"%d",i);
    return atom_str(buffer);
}



Node *build_list(Node *n1, ...)
{
    Node *root;
    Node *next;
    Node *ROOT;

    va_list ap;

    va_start(ap, n1);

    root = ROOT = create_node_list();
    
    /*
     * no action if n1 is nil
     * we are trying to construct '()'
     */

    if (! isNil(n1)) {
	root->val.list.head = n1;
	next = va_arg(ap, Node *);
	while (! isNil(next)) {
	    Node *new;
	    new = create_node_list();
	    new->val.list.head = next;
	    root->val.list.tail = new;
	    root = new;
	    next = va_arg(ap, Node *);
	}
    }

    va_end(ap);

    return ROOT;
}



Node *tail_list(Node *n)
{
    Node *tail;

    if (isNil(n))
	tail = nil;
    else if (isAtom(n))
	tail = nil;
    else {
	for ( tail = n; !isNil(cdr(tail)); tail = cdr(tail)) ;
    }

    return tail;
}


Node *join_list(Node *n1,...)
{
    Node *tail;
    Node *ROOT;
    Node *next;

    va_list ap;

    va_start(ap, n1);

    ROOT = n1;
    tail = tail_list(ROOT);

    next = va_arg(ap, Node *);
    while (! isNil(next)) {
	Node *newtail;
	newtail = tail_list(next);
	tail->val.list.tail = next;
	tail = newtail;
	next = va_arg(ap, Node *);
    }

    va_end(ap);

    return ROOT;
}



static void _print_list (Node *n)
{

    if ( isNil(n) ) {
	printf("nil ");
    } else if (isAtom(n))
	printf("%s ", atomVal(n));
    else {
	Node *next;
	printf("( ");
	for (next = n; ! isNil(next); next = cdr(next))
	    _print_list(car(next));
	printf(") ");
    }

}

void print_list (Node *n)
{
    _print_list(n);
    nl;
}



static void printAtom(FILE *f, Atom *n)
{
    char *val = atomVal(n);
    char *q;

    for (q=val; *q && !(isspace(*q)) && *q!='"' && *q!='(' && *q!=')'; q++);

    if (*q) putc('"',f);

    for (; *val; val++) {
	if (*val=='"') putc('"',f);
	putc(*val,f);
    }

    if (*q) putc('"',f);

}


static void _print_list_f (FILE *f,Node *n,int level)
{

    if ( isNil(n) ) {
	fprintf(f,"nil ");
    } else if (isAtom(n)) {
	printAtom(f,n);
	fprintf(f," ");
    } else {
	Node *next;
	int i;
	if (level) fprintf(f,"\n");
	for (i=0;i<level;i++) fprintf(f,"  ");
	fprintf(f,"( ");
	for (next = n; ! isNil(next); next = cdr(next))
	    _print_list_f(f,car(next),level+1);
	fprintf(f,") ");
    }

}

void print_list_f (Node *n)
{
    _print_list_f(stdout, n,0);
    nl;
}


void fprint_list_f (FILE *f, Node *n)
{
    _print_list_f(f, n,0);
    fnl(f);
}


Node *index_list(Node *n, Atom *a)
{
    Node *index;

    if (! isList(n) || ! isAtom(a))
	index = nil;
    else {
	Node *i;
	int found = 0;
	for (i = n; ! isNil(i) && !found; i=cdr(i)) {
	    index = car(i);
	    if (isList(index) && isAtom(car(index))) {
		found = (strcmp(atomVal(car(index)),atomVal(a))==0);
	    }
	}
	if (! found) index = nil;
    }

    return index;
}


    

Node *index_list_by_str(Node *n, char *a)
{
    Atom *atom;
    List *node;

    atom = atom_str(a);

    node = index_list(n,atom);

    destroy_list(atom);

    return node;
}
    


char *assoc(Node *n, char *a)
{
    return atomVal(car(cdr(index_list_by_str(n,a))));

}




int read_string(FILE *fp, char **s)
{
    char buff[4096];
    int c;
    int len;
    int ret;

    /* skip over while space */
    for (c = getc(fp);c != EOF && isspace(c); c= getc(fp));

    len = 0;
    ret = 0;
    if (c == EOF)
	*s = NULL;
    else if (c == '"') {
	int l;
	for (c = getc(fp), l = getc(fp);
	     c != EOF && !(c=='"' && (l == EOF ||  l!='"') );
	     c=l,l = getc(fp)) {
	    buff[len++] = c;
	    if (c=='"' && l=='"') {
		c = '\0';
		l = getc(fp);
	    }
	}
	ret = 1;
	buff[len] = '\0';
	*s = (char *)mystrdup(buff);
    } else {
	buff[len++] = c;
	for (c = getc(fp); c != EOF && !isspace(c); c = getc(fp))
	    buff[len++] = c;
	buff[len] = '\0';
	*s = (char *)mystrdup(buff);
    }

    return ret;
}



List *read_list(FILE *fp)
{
    List *l;
    char *s;
    int q;

    q = read_string(fp,&s);

    if (s == NULL)
	l = nil;
    else if ( strcmp(s,")") == 0 && !q) {
	free(s);
	l = nil;
    } else if ( strcmp(s,"(") == 0 && !q) {
	List *item, *t;
	free(s);
	l = nil;
	t = nil;
	for (item = read_list(fp); ! isNil(item); item = read_list(fp)) {
	    List *m;
	    m = create_node_list();
	    m->val.list.head = item;
	    if (isNil(l)) l = m;
	    if (isNil(t))
		t = m;
	    else {
		t->val.list.tail = m;
		t = m;
	    }
	}
	if (isNil(l)) l = create_node_list();
    } else {
	l = atom_str(s);
    }

    return l;

}

