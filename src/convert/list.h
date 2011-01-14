#ifndef _list_h
#define _list_h
#include <stdio.h>
#include <stdarg.h>

#define Node_Nil  0
#define Node_List 1
#define Node_Atom 2
#define nl putchar('\n')
#define fnl(F) putc('\n',F)


#define nil NULL
#define isAtom(N) ((N)==nil?0:(N)->type==Node_Atom)
#define isList(N) ((N)==nil?0:(N)->type==Node_List)
#define isNil(N)  ((N)==nil)
#define car(N) ( isNil(N)||isAtom(N)?nil:(N)->val.list.head )
#define cdr(N) ( isNil(N)||isAtom(N)?nil:(N)->val.list.tail )
#define atomVal(N) ( isAtom(N)?(N)->val.atom.ptr:"" )

typedef struct _node{
    int type;
    union {
	struct {
	    struct _node *head;
	    struct _node *tail;
	} list;
	struct {
	    char *ptr;
	    char *spare;
	} atom;
    }val;
} Node,List,Atom;

extern void destroy_list(Node *n);
extern void destroy_node_list(void);
extern Node *copy_list(Node *n);
extern Node *atom_str(char *s);
extern Node *atom_int(int i);
extern Node *build_list(Node *n1, ...);
extern Node *tail_list(Node *n);
extern Node *join_list(Node *n1,...);
extern Node *index_list(Node *n, Atom *a);
extern void print_list (Node *n);
extern Node *index_list_by_str(Node *n, char *s);
extern char *assoc(Node *n, char *s);
extern void print_list_l (Node *n);
extern void fprint_list_l (FILE *f, Node *n);
extern List *read_list(FILE *f);
extern void fprint_list_f (FILE *f, Node *n);

#endif /* _list_h */
