#include <math.h>

/*
 * Distribution of template sizes. Probability of being this size => sums to 1
 */
double tsize_dist[5000];
int dist_len = 5000;

void set_tsize(void) {
    int i, sum;

    for (i = 0; i < 5000; i++)
	tsize_dist[i] = 0;

    sum = 0;
    for (i = 2000; i <= 4000; i++) {
	sum += (tsize_dist[i] = 2000 - abs(3000 - i));
    }
    for (i = 0; i < 5000; i++) {
	tsize_dist[i] /= sum;
    }
}

/*
 * Returns the probability of a sequence of length 'length' from either
 * the start (end = 0) or end (end = 1) of the template covering position
 * 'pos'.
 *
 * FIXME: This in itself should take note that the sequences are of varying
 * lengths so we should multiply the tsize_dist[] element with a rsize_dist[]
 * element too.
 */
double pcover(int pos, int length, int end) {
    int i;
    double prob;

    if (end) {
	for (prob = i = 0; i < length && pos + i < dist_len; i++) {
	    prob += tsize_dist[pos + i];
	}
    } else {
	for (prob = i = 0; i < length && pos - i >= 0; i++) {
	    prob += tsize_dist[pos - i];
	}
    }

    return prob;
}

int main(void) {
    set_tsize();
    printf("cover=%f\n", pcover(2500, 500, 0));
    printf("cover=%f\n", pcover(2500, 500, 1));
    return 0;
}
