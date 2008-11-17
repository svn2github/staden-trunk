#ifndef _BAF_H_
#define _BAF_H_

#include <hache_table.h>

#define CC2(a,b) ((((unsigned char)a)<<8) | ((unsigned char)b))

enum line_type {
    XX=0, /* Blank line */
    CO=CC2('C','O'), /* Contig */
    LN=CC2('L','N'), /*    Length */
    SO=CC2('S','O'), /* DNA Source */
    ST=CC2('S','T'), /*    Source type */
    SI=CC2('S','I'), /*    Insert size mean */
    SS=CC2('S','S'), /*    Insert size standard deviation */
    SV=CC2('S','V'), /*    vector */
    PA=CC2('P','A'), /*    parent source */
    RD=CC2('R','D'), /* Reading */
    SQ=CC2('S','Q'), /*    Sequence */
    FQ=CC2('F','Q'), /*    Fastq quality */
    AP=CC2('A','P'), /*    Contig position */
    QL=CC2('Q','L'), /*    Left quality clip */
    QR=CC2('Q','R'), /*    Right quality clip */
    TN=CC2('T','N'), /*    Template name */
    DR=CC2('D','R'), /*    Direction, 1=>uncomp, -1=>complemented */
    TR=CC2('T','R'), /*    Trace name */
    MQ=CC2('M','Q'), /*    Mapping quality */
    AL=CC2('A','L'), /*    Alignment */

    /* Regexp versions of the above */
    ln=CC2('l','n'),
    st=CC2('s','t'),
    si=CC2('s','i'),
    ss=CC2('s','s'),
    sv=CC2('s','v'),
    pa=CC2('p','a'),
    sq=CC2('s','q'),
    fq=CC2('f','q'),
    ap=CC2('a','p'),
    ql=CC2('q','l'),
    qr=CC2('q','r'),
    tn=CC2('t','n'),
    dr=CC2('d','r'),
    tr=CC2('t','r'),
    mq=CC2('m','q'),
    al=CC2('a','l'),
};

typedef struct {
    /* Allocated memory and size */
    char *str;
    size_t len;

    /* Key, value and assignment type, eg CO=contig#1 */
    char *value;
    enum line_type type;
    int assign;  /* = or : */

    /* Order so we can reconstruct the file exactly */
    int order;
} line_t;

typedef struct {
    int type;
    HacheTable *h;
} baf_block;

/* The data attached to the hache */
typedef struct {
    char *value;
    int assign;
} baf_data;



void free_line(line_t *l);
char *linetype2str(int lt);
line_t *get_line(FILE *fp, line_t *in);
baf_block *baf_next_block(FILE *fp);
void baf_block_destroy(baf_block *b);
line_t *baf_line_for_type(baf_block *b, int type);
int parse_baf(GapIO *io, char *fn, int no_tree, int pair_reads,
	      int merge_contigs);

#endif /* _BAF_H_ */
