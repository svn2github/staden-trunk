#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <dna_utils.h>

#define TRUE   1
#define FALSE  0

static int word = 3; 
static int window = 48; 
static int window2 = 24; 
static int level = 20;

static int mv, iv, jv;

void set_dust_level(int value)
{
	level = value;
}

void set_dust_window(int value)
{
	window = value;
	window2 = window / 2;
}

void set_dust_word(int value)
{
	word = value;
}

static void wo1(int len, char *s, int ivv)
{
	int i, ii, j, v, t, n, n1, sum;
	static int counts[32*32*32];
	static int iis[32*32*32];
	int js, nis;

	n = 32 * 32 * 32;
	n1 = n - 1;
	nis = 0;
	i = 0;
	ii = 0;
	sum = 0;
	v = 0;
	for (j=0; j < len; j++, s++) {
		ii <<= 5;
		if (isalpha(*s)) {
			if (islower(*s)) {
				ii |= *s - 'a';
			} else {
				ii |= *s - 'A';
			}
		} else {
			i = 0;
			continue;
		}
		ii &= n1;
		i++;
		if (i >= word) {
			for (js=0; js < nis && iis[js] != ii; js++) ;
			if (js == nis) {
				iis[nis] = ii;
				counts[ii] = 0;
				nis++;
			}
			if ((t = counts[ii]) > 0) {
				sum += t;
				v = 10 * sum / j;
				if (mv < v) {
					mv = v;
					iv = ivv;
					jv = j;
				}
			}
			counts[ii]++;
		}
	}
}

static int wo(int len, char *s, int *beg, int *end)
{
	int i, l1;

	l1 = len - word + 1;
	if (l1 < 0) {
		*beg = 0;
		*end = len - 1;
		return 0;
	}
	mv = 0;
	iv = 0;
	jv = 0;
	for (i=0; i < l1; i++) {
		wo1(len-i, s+i, i);
	}
	*beg = iv;
	*end = iv + jv;
	return mv;
}

void dust(int len, char *s)
{
	int i, j, l, from, to, a, b, v;
	char *depadded = (char *)malloc(len);
	int *depad_to_pad = (int *)calloc(len, sizeof(int));
	int depadded_len;

	if (!depadded || !depad_to_pad)
	    return;

	memcpy(depadded, s, len);
	depadded_len = len;
	depad_seq(depadded, &depadded_len, depad_to_pad);

	from = 0;
	to = -1;
	for (i=0; i < depadded_len; i += window2) {
		from -= window2;
		to -= window2;
		l = (depadded_len > i+window) ? window : depadded_len-i;
		v = wo(l, depadded+i, &a, &b);
		for (j = from; j <= to; j++) {
			if (isalpha(s[depad_to_pad[i+j]]))
				s[depad_to_pad[i+j]] = '#';
		}
		if (v > level) {
			for (j = a; j <= b && j < window2; j++) {
				if (isalpha(s[depad_to_pad[i+j]]))
					s[depad_to_pad[i+j]] = '#';
			}
			from = j;
			to = b;
		} else {
			from = 0;
			to = -1;
		}
	}

	free(depadded);
	free(depad_to_pad);
}
