/* sorting algorithms */

#pragma global noalias
#include <stdio.h>
#include "vtc.h"
#include "vtclocal.h"

static void print_array(int *a, int cnt);
static void sort_sanity_check(Mortonkey *key, int n);
static void insertionsort(Body *a, Mortonkey *key, int cnt);
static void insertionsort_reverse(Body *a, Mortonkey *key, int cnt);
static void bitonicmerge_hybrid(Body *a, Mortonkey *key, int lo, int cnt, int direction);
static void bitonicmergeS(Body *a, Mortonkey *key, int lo, int cnt, int direction);
static void bitonicmerge(Body *a, Mortonkey *key, int lo, int cnt, int direction);
static void bitonicsort_hybrid(Body *a, Mortonkey *key, int lo, int cnt, int direction);
static void bitonicsort(Body *a, Mortonkey *key, int lo, int cnt, int direction);
static void myqsort(Body *a, Mortonkey *key, int lo, int up);
static void myqsort_hybrid(Body *a, Mortonkey *key, int lo, int up, int crit);
static int ispow2(int n);

#define ASCENDING (1)
#define DESCENDING (0)
#define BITONIC_MIN_LEN (8)

static int narray;

/*
 * global function(s)
 */
void
sort_body(Body *b, Mortonkey *key, int n)
{
    narray = n;
#if 0
    myqsort(b, key, 0, n-1);
#elif 1
    myqsort_hybrid(b, key, 0, n-1, 16);
#elif 0
    if (!ispow2(n)) {
	fprintf(stderr, "%s sort_body: n (%d) is not a power of 2\n", __FILE__, n);
	exit(1);
    }
    bitonicsort(b, key, 0, n, ASCENDING);
#else
    if (!ispow2()) {
	fprintf(stderr, "%s sort_body: n (%d) is not a power of 2\n", __FILE__, n);
	exit(1);
    }
    bitonicsort_hybrid(b, key, 0, n, ASCENDING);
#endif
    /*
    sort_sanity_check(key, n);
    */
}

/*
 * sorting algorithm(s)
 */
static void
insertionsort(Body *a, Mortonkey *key, int cnt)
{
    int i, j, k;
    Body tmpa;
    Mortonkey tmpk;

#pragma loop scalar
    for (i = 1; i < cnt; i++) {
#pragma loop scalar
	for (j = 0; j < i; j++) {
	    if (key[j] > key[i]) {
		tmpa = a[i];
		tmpk = key[i];
#pragma loop scalar
		for (k = i; k > j; k--) {
		    a[k] = a[k-1];
		    key[k] = key[k-1];
		}
		a[j] = tmpa;
		key[j] = tmpk;
		break;
	    }
	}
    }
}

static void
insertionsort_reverse(Body *a, Mortonkey *key, int cnt)
{
    int i, j, k;
    Body tmpa;
    Mortonkey tmpk;

#pragma loop scalar
    for (i = 1; i < cnt; i++) {
#pragma loop scalar
	for (j = 0; j < i; j++) {
	    if (key[j] < key[i]) {
		tmpa = a[i];
		tmpk = key[i];
#pragma loop scalar
		for (k = i; k > j; k--) {
		    a[k] = a[k-1];
		    key[k] = key[k-1];
		}
		a[j] = tmpa;
		key[j] = tmpk;
		break;
	    }
	}
    }
}

static void
bitonicmerge_hybrid(Body *a, Mortonkey *key, int lo, int cnt, int direction)
{
    int k, i, tmpa;
    Mortonkey tmpk;

    if (cnt <= 1) {
	return;
    }
    if (cnt <= BITONIC_MIN_LEN) {
	if (direction == ASCENDING) {
	    insertionsort(a+lo, key+lo, cnt);
	}
	else {
	    insertionsort_reverse(a+lo, key+lo, cnt);
	}
	return;
    }
    k = cnt / 2;
    for (i = lo; i < lo+k; i++) {
	if (direction == (key[i] > key[i+k])) {
	    tmpa = a[i];
	    a[i] = a[i+k];
	    a[i+k] = tmpa;

	    tmpk = key[i];
	    key[i] = key[i+k];
	    key[i+k] = tmpk;
	}
    }
    bitonicmerge_hybrid(a, key, lo, k, direction);
    bitonicmerge_hybrid(a, key, lo+k, k, direction);
}

static void
bitonicmergeS(Body *a, Mortonkey *key, int lo, int cnt, int direction)
{
#pragma procedure scalar
    int k, i, tmpa;
    Mortonkey tmpk;

    if (cnt <= 1) {
	return;
    }
    k = cnt / 2;
    for (i = lo; i < lo+k; i++) {
	if (direction == (key[i] > key[i+k])) {
	    tmpa = a[i];
	    a[i] = a[i+k];
	    a[i+k] = tmpa;

	    tmpk = key[i];
	    key[i] = key[i+k];
	    key[i+k] = tmpk;
	}
    }
    bitonicmergeS(a, key, lo, k, direction);
    bitonicmergeS(a, key, lo+k, k, direction);
}

static void
bitonicmerge(Body *a, Mortonkey *key, int lo, int cnt, int direction)
{
    int k, i, tmpa;
    Mortonkey tmpk;

    if (cnt <= 1) {
	return;
    }
    k = cnt / 2;
    for (i = lo; i < lo+k; i++) {
	if (direction == (key[i] > key[i+k])) {
	    tmpa = a[i];
	    a[i] = a[i+k];
	    a[i+k] = tmpa;

	    tmpk = key[i];
	    key[i] = key[i+k];
	    key[i+k] = tmpk;
	}
    }
    if (k > 64) {
	bitonicmerge(a, key, lo, k, direction);
	bitonicmerge(a, key, lo+k, k, direction);
    }
    else {
	bitonicmergeS(a, key, lo, k, direction);
	bitonicmergeS(a, key, lo+k, k, direction);
    }
}

static void
bitonicsort_hybrid(Body *a, Mortonkey *key, int lo, int cnt, int direction)
{
    int k;
    if (cnt <= 1) {
	return;
    }
    k = cnt / 2;
    /*
    print_array(a, narray);
    */
    bitonicsort_hybrid(a, key, lo, k, ASCENDING);
    bitonicsort_hybrid(a, key, lo+k, k, DESCENDING);
    bitonicmerge_hybrid(a, key, lo, cnt, direction);
}

static void
bitonicsort(Body *a, Mortonkey *key, int lo, int cnt, int direction)
{
    int k;
    if (cnt <= 1) {
	return;
    }
    k = cnt / 2;
    /*
    print_array(a, narray);
    */
    bitonicsort(a, key, lo, k, ASCENDING);
    bitonicsort(a, key, lo+k, k, DESCENDING);
    bitonicmerge(a, key, lo, cnt, direction);
}

static void
myqsort(Body *a, Mortonkey *key, int lo, int up)
{
    int i, j;
    Body tempa;
    Mortonkey tempk;

    while (up > lo) {
	i = lo;
	j = up;
	tempa = a[lo];
	tempk = key[lo];
	while (i < j) {
	    for (; key[j] > tempk; j--);
	    for (a[i] = a[j], key[i] = key[j]; i < j && key[i] <= tempk; i++);
	    a[j] = a[i];
	    key[j] = key[i];
	}
	a[i] = tempa;
	key[i] = tempk;

	if (i-lo < up-i) {
	    myqsort(a, key, lo, i-1);
	    lo = i+1;
	}
	else {
	    myqsort(a, key, i+1, up);
	    up = i-1;
	}
    }
}

static void
myqsort_hybrid(Body *a, Mortonkey *key, int lo, int up, int crit)
{
    int i, j;
    Body tempa;
    Mortonkey tempk;

    if (up - lo < crit) {
	insertionsort(a+lo, key+lo, up-lo+1);
    }
    else {
	while (up > lo) {
	    i = lo;
	    j = up;
	    tempa = a[lo];
	    tempk = key[lo];
	    while (i < j) {
		for (; key[j] > tempk; j--);
		for (a[i] = a[j], key[i] = key[j]; i < j && key[i] <= tempk; i++);
		a[j] = a[i];
		key[j] = key[i];
	    }
	    a[i] = tempa;
	    key[i] = tempk;

	    if (i-lo < up-i) {
		myqsort_hybrid(a, key, lo, i-1, crit);
		lo = i+1;
	    }
	    else {
		myqsort_hybrid(a, key, i+1, up, crit);
		up = i-1;
	    }
	}
    }
}

static int
ispow2(int n)
{
    int i;
    int bitset = 0;

    for (i = 0; i < 32; i++) {
	if (n & (1<<i)) {
	    bitset++;
	}
	if (bitset > 1) {
	    return (FALSE);
	}
    }
    return (TRUE);
}


/*
 * for debugging
 */
static void
print_array(int *a, int cnt)
{
    int i;

    for (i = 0; i < cnt; i++) {
	fprintf(stdout, "\n%04d", a[i]);
    }
    fprintf(stdout, "\n");
}

static void
sort_sanity_check(Mortonkey *key, int n)
{
    int i;

    for (i = 0; i < n-1; i++) {
	if (key[i] > key[i+1]) {
	    fprintf(stderr, "sort failed. key[%d] = %d key[%d] = %d\n",
		    i, key[i], i+1, key[i+1]);
	    exit(1);
	}
    }
}

