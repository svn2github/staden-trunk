#include <stdio.h>

double delta_g[5][5] = {
    /* A A */ -1.2,
    /* A C */ -1.5,
    /* A G */ -1.5,
    /* A T */ -0.9,
    /* A - */ 0,
    /* C A */ -1.7,
    /* C C */ -2.1,
    /* C G */ -2.8,
    /* C T */ -1.5,
    /* C - */ 0,
    /* G A */ -1.5,
    /* G C */ -2.3,
    /* G G */ -2.1,
    /* G T */ -1.5,
    /* G - */ 0,
    /* T A */ -0.9,
    /* T C */ -1.5,
    /* T G */ -1.7,
    /* T T */ -1.2,
    /* T - */ 0,
    /* - A */ 0,
    /* - C */ 0,
    /* - G */ 0,
    /* - T */ 0,
    /* - - */ 0
};

double init_g = 3.4;

int lookup[256];
init_lookup(void) {
    int i;
    for (i = 0; i < 256; i++)
	lookup[i] = 5;

    lookup['A'] = 0;
    lookup['C'] = 1;
    lookup['G'] = 2;
    lookup['T'] = 3;
}

int main(int argc, char **argv) {
    char *str;
    double score;

    if (argc != 2) {
	return 1;
    }
    str = argv[1];
    init_lookup();

    score = init_g;
    while (str[0] && str[1]) {
	int a = lookup[str[0]];
	int b = lookup[str[1]];

	printf("%c%c => %f\n", str[0], str[1], delta_g[a][b]);

	score += delta_g[a][b];
	str++;
    }

    printf("Score=%f\n", score);
}
