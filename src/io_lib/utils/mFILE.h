#ifndef _MFILE_H_
#define _MFILE_H_

#include <stdio.h>

typedef struct {
    char *fname;
    FILE *fp;
    char *data;
    size_t alloced;
    int eof;
    size_t size;
    size_t offset;
    size_t flush_pos;
} mFILE;

mFILE *mfreopen(const char *path, const char *mode, FILE *fp);
mFILE *mfopen(const char *path, const char *mode);
int mfclose(mFILE *mf);
int mfdestroy(mFILE *mf);
int mfseek(mFILE *mf, long offset, int whence);
long mftell(mFILE *mf);
void mrewind(mFILE *mf);
void mftruncate(mFILE *mf, long offset);
int mfeof(mFILE *mf);
size_t mfread(void *ptr, size_t size, size_t nmemb, mFILE *mf);
size_t mfwrite(void *ptr, size_t size, size_t nmemb, mFILE *mf);
int mfgetc(mFILE *mf);
int mungetc(int c, mFILE *mf);
mFILE *mfcreate(char *data, int size);
void mfrecreate(mFILE *mf, char *data, int size);
char *mfgets(char *s, int size, mFILE *mf);
int mfflush(mFILE *mf);
int mfprintf(mFILE *mf, char *fmt, ...);
mFILE *mstdin(void);
mFILE *mstdout(void);
mFILE *mstderr(void);

#endif /* _MFILE_H_ */
