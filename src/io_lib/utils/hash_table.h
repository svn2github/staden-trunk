#ifndef _HASH_TABLE_H_
#define _HASH_TABLE_H_

#include <inttypes.h>

/* The data referenced by the hash table */
typedef union {
    uint64_t i;
    void *p;
    struct {
	uint64_t pos;
	uint32_t len;
    } idx;
} HashData;

/* A hash item with "next" pointer to use in a linked list */
typedef struct HashItemStruct {
    HashData data;        /* user defined data attached to this key */
    char    *key;         /* key we hashed on */
    int      key_len;     /* and its length */
    struct HashItemStruct *next;
} HashItem;

/* The main hash table structure itself */
typedef struct {
    int       options;  /* HASH_FUNC & HASH_OPT macros */
    uint32_t  nbuckets; /* Number of hash buckets; power of 2 */
    uint32_t  mask;	/* bit-mask equiv of nbuckets */
    int       nused;    /* How many hash entries we're storing */
    HashItem **bucket;  /* The bucket "list heads" themselves */
} HashTable;

/* File format: the hash table header */
typedef struct {
    char magic[4];
    char vers[3];
    char hfunc;
    uint32_t nbuckets;
    uint32_t size;
} HashFileHeader;

typedef struct {
    HashFileHeader h;
    FILE *hfp;		/* hash FILE */
    FILE *afp;		/* archive FILE */
    int header_size;	/* size of header + archive filename */
} HashFile;

/* Functions to to use HashTable.options */
#define HASH_FUNC_HSIEH       0
#define HASH_FUNC_TCL         1
#define HASH_FUNC_JENKINS     2
#define HASH_FUNC_MASK        7

/* Other HashTable.options values */
#define HASH_NONVOLATILE_KEYS (1<<3)
#define HASH_ALLOW_DUP_KEYS   (1<<4)
#define HASH_DYNAMIC_SIZE     (1<<5)

/* Hashing prototypes */
uint32_t hash(int func, uint8_t *key, int key_len);
uint32_t HashJenkins(uint8_t *k, int length);
uint32_t HashTcl(uint8_t *data, int len);
uint32_t HashHsieh(uint8_t *k, int length);

/* HashTable management prototypes */
HashTable *HashTableCreate(int size, int options);
void HashTableDestroy(HashTable *h);
HashItem *HashTableAdd(HashTable *h, char *key, int key_len,
		       HashData data, int *new);
HashItem *HashTableSearch(HashTable *h, char *key, int key_len);
void HashTableStats(HashTable *h, FILE *fp);

/* HashFile prototypes */
void HashFileSave(HashTable *h, FILE *fp, char *archive);
HashTable *HashFileLoad(FILE *fp);
int HashFileQuery(HashFile *hf, uint8_t *key, int key_len,
		  uint64_t *r_pos, uint32_t *r_size);

HashFile *HashFileOpen(char *fname);
void HashFileClose(HashFile *hf);

#endif /* _HASH_TABLE_H_ */
