#ifndef _HASH_H
#define _HASH_H

#define HASHMODULUS 100


typedef int ItemType;

typedef struct tableitem *TablePtr;
typedef struct tableitem {
    int key;
    ItemType info;
    TablePtr next;
} Table;

void
InitialiseTable(TablePtr T[]);

void
ChainInsert(TablePtr T[],
	    int newkey,
	    ItemType newinfo);

void
ChainSearch(TablePtr T[],
	    int search_key,
	    int *found,
	    ItemType *search_info);

void 
ChainDelete(TablePtr T[], 
	    int search_key);
#endif /* _HASH_H */
