//
// procedures to create all different cyclic permutations of the relators
//

#define RCPR__C
#include "lowx.h"

#ifdef XREF

#define MAX_NO_RELATORS 64

typedef int rxref[MAX_NO_GENS][MAX_NO_RELATORS];

rxref * make_rxref_table(struct Parameters *p, int no_relators, int *relators)
{
    rxref * x;
    int i,j;

    if ((x = (rxref *) malloc(sizeof(rxref))) == NULL) {
        Error(ABORT, errno, "Can't allocate memory");
    }

    for (i = 0; i < MAX_NO_GENS; i++) {
        for (j = 0; j < MAX_NO_RELATORS; j++) {
            (*x)[i][j] = EMPTY;
        }
    }

    return x;
}

#else

static int compare_words(struct Word *p, struct Word *q)
// compare two words lexicographically
// return
//      -1 if *p < *q
//      0  if *p == *q
//      +1 if *p > *q
{
int i;
int shorter;

    if (p->length < q->length) {
        shorter = p->length;
    } else {
        shorter = q->length;
    }
    for (i=0; i<shorter; i++) {
        if (p->element[i] < q->element[i]) {
            return -1;
        } else if (p->element[i] > q->element[i]) {
            return +1;
        }
    }
    if (p->length == q->length) {
        return 0;
    } else if (p->length < q->length) {
        return -1;
    } else {
        return +1;
    }
}


static void sort_words(struct Warray *wap)
// using Shell sort
{
int gap, i, j;
struct Word temp;

    for(gap = wap->no_words / 2; gap > 0; gap /= 2) {
        for (i = gap; i < wap->no_words; i++) {
            for (j = i - gap; j >= 0 &&
                 compare_words(&(wap->word[j]), &(wap->word[j+gap])) > 0;
                 j -= gap) {
                     temp.length = wap->word[j].length;
                     temp.element = wap->word[j].element;
                     wap->word[j].length = wap->word[j+gap].length;
                     wap->word[j].element = wap->word[j+gap].element;
                     wap->word[j+gap].length = temp.length;
                     wap->word[j+gap].element = temp.element;
            }
        }
    }
}


static void make_reduced_cyclic_permutations(struct Warray *wap, int no_relators, int *relators)
{
int i,j;
int n = 0;
struct Warray w1;

    // calculate the number of all cyclic permutations of all relators
    for (i=j=0; i<no_relators; i++, j+=relators[j]+1) {
        n += relators[j];
    }

    // allocate an array of Words to hold all cyclic permutations
    w1.no_words = n;
    if ((w1.word = (struct Word *) calloc((size_t) n, sizeof(struct Word)))
        == NULL) {
        Error(ABORT, errno, "Can't allocate memory");
    }
 
    for (i=j=n=0; i<no_relators; i++, j+=relators[j]+1) {
        int k,l;
        int first, last, length, tmp;

        first = j + 1;
        last = j + relators[j];
        length = relators[j];

        for (k=0; k<length; k++) {
            // copy the current cyclic permutation of the j-th relator
            // into the n-th word of w1
            w1.word[n].length = length;
            if ((w1.word[n].element
                = (int *) calloc((size_t) length, sizeof(int))) == NULL) {
                Error(ABORT, errno, "Can't allocate memory");
            }
            for (l=0; l<length; l++) {
                w1.word[n].element[l] = relators[first +l];
            }

            // cyclic left-shift of the current (j-th) relator
            tmp = relators[first];
            for (l=first; l<last; l++) {
                relators[l] = relators[l+1];
            }
            relators[last] = tmp;

            // go to the next word of w1
            ++n;
        }
    }

    sort_words(&w1);

    // count the number of different words in w1
    for (i = 0, j = n = 1; j < w1.no_words; j++) {
        if (compare_words(&w1.word[i], &w1.word[j]) != 0) {
            ++n;
            i = j;
        }
    }

    // allocate an array of words big enough to hold the different words
    wap->no_words = n;
    if ((wap->word = (struct Word *) calloc((size_t) n, sizeof(struct Word)))
        == NULL) {
        Error(ABORT, errno, "Can't allocate memory");
    }
 
    // copy the different words into the new array of words
    wap->word[0].length = w1.word[0].length;
    wap->word[0].element = w1.word[0].element;
    for (i = 0, j = 1; j < w1.no_words; j++) {
        if (compare_words(&(wap->word[i]), &w1.word[j]) != 0) {
            ++i;
            wap->word[i].length = w1.word[j].length;
            wap->word[i].element = w1.word[j].element;
        } else {
            // we already have this word in the new array
            // so we can 'delete' (actually 'free') it
            free((void *) w1.word[j].element);
        }
    }

    // free the array of word allocated to w1
    free((void *) w1.word);
}


struct Warray * make_rcpr_table(int no_columns, int no_relators, int *relators)
{
// make a table of the Reduced Cyclic Permutations of the Relators
// the i-th row contains the array of RCPRs which starts with generator i
// when this function returns wap disappears
//
int i;
int offset;
struct Warray wap;
struct Warray *rcpr_table;

    make_reduced_cyclic_permutations(&wap, no_relators, relators);

    // allocate an RCPR table
    if ((rcpr_table =
        (struct Warray *) malloc(no_columns * sizeof(struct Warray))) == NULL) {
        Error(ABORT, errno, "Can't allocate memory");
    }

    // calculates the lengths af the word arrays
    // starting with the same generator
    for (i = 0; i < no_columns; i++) {
        rcpr_table[i].no_words = 0;
    }
    for (i = 0; i < wap.no_words; i++) {
        ++rcpr_table[wap.word[i].element[0]].no_words;
    }

    // "copy" the words into the table
    // we use that *wap is sorted in increasing order
    rcpr_table[0].word = wap.word;
    offset = rcpr_table[0].no_words;
    for (i = 1; i < no_columns; i++) {
        rcpr_table[i].word = &(wap.word[offset]);
        offset += rcpr_table[i].no_words;
    }

    return rcpr_table;
}


void free_rcpr_table(struct Warray *rcpr_t, int no_columns)
{
int i, j;
struct Warray *wap;
struct Word *wp;

    for (i = 0; i < no_columns; i++) {
        wap = &(rcpr_t[i]);
        for (j = 0; j < wap->no_words; j++) {
            wp = &(wap->word[j]);
            free(wp->element);
        }
    }
    // all Wordarrays for the table was allocated in one chunk
    // that is why we can free them this way
    wap = &(rcpr_t[0]);
    free(wap->word);
    free(rcpr_t);
}

#endif
