#include <stdio.h>
#include "bitmap.h"
#include "misc.h"

int main() {
    int i, freerec;
    Bitmap b;
    
    b = BitmapCreate(0, 0);
    for (i=0; i<1000; i++) {
	freerec = BitmapFree(b);
	BIT_SET(b, freerec);
	printf(">>>> %d:%d <<<<\n", i, freerec);
	BitmapPrint(stdout, b);
	while (random()%3 == 0) {
	    int bit = random()%b->Nbits;
	    printf("clear %d\n", bit);
	    BIT_CLR(b, bit);
	    BitmapPrint(stdout, b);
	}
    }

    BitmapPrint(stdout, b);
}
