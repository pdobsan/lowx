
/*
This module contains the subroutines used by the parser for the language
describing group presentations and subgroups.

Here is an example for input accepted by the parser:

    #
    # sample input for 'cosetenum'
    #

    group = <a, b | b * a * b * a^-2, a * b * a * b^-2>

    subgroup = <a^2>


The computer data representations of the algebraic structures:

Let
E = {r1,r2,...rN,r1_,r2_,...,rN_} generators and their inverses
G = <r1,r2,...,rN | R1(e1,e2,..) R2(e1,e2,...), ... RK(e1,e2,...)>
H = <w1,w2,...,wM> is a subgroup of G generated by w-s
where e-s are elements of E  and  w-s and R-s are words on E

mapping of the generators and their inverses:

E={r1,r2,...rN,r1_,r2_,...,rN_} --> [0,1,...,N-1,N,N+1,...2N-1]

Because we allow only one letter generator symbols the symbol table is a
vector of the corresponding characters. The index of the symbol in the
array is simply the code. The code of the inverse of rI is the code of
rI plus N.

The words on E are represented as a vector of integers whose first
element is the number of codes and the rest is the codes of the elements
in the word. e^k is expressed by placing k consecutive codes of e in the
array.

for example if G = <a, b, | ...>
    a     --> 0
    b     --> 1
    a^-1  --> 2
    b^-1  --> 3

    a^3 * b^-1  --> [4,0,0,0,3]

The relators of the group are represented as a one-dimensional array of
integers created by packing the computer representations of the
individual relators after each other.

The computer representation of the subgroup generators is created in a
similar manner.

*/

#define PARSER_AUX__C
#include "lowx.h"

extern struct Parameters params;

int *counter;   // points to no_XXX depending what is currently parsed
int *vector;    // point to the vector currently filled
int length_x;   // the index of the current word's lenght in "vector"
int code_x;     // current index of the next free code-entry in "vector"

// a frame is a segment of the current word satisfying the following relation
//
//      vector[i] where frame[fx] <= i < code_x
//
int fx = -1;            // index of the current frame pointer
int frame[MAX_VECTOR];  // array of frame pointers (i.e. indices in 'vector')


void new_frame()
{
    ++fx;
    frame[fx] = code_x;
}


void remove_frame()
{
    --fx;
}


void invert_word(int *w, int from, int to)
{
int i, j, tmp;

    // reverse
    for (i = from, j = to - 1; i < (from + to) / 2; i++, j--) {
        tmp = w[i];
        w[i] = w[j];
        w[j] = tmp;
    }

    // invert
    for (i = from; i < to; i++) {
        w[i] = (w[i] + params.no_generators) % (2 * params.no_generators);
    }
}


void append_frame(int power)
{
int i, j;
int stop;

    stop = code_x;      // mark the end of the current frame


    if (power < 0) {
        invert_word(vector,frame[fx], stop);
    }

    // the first copy of the frame is already there
    // so append it abs(power)-1 times more
    for (i = 0; i < abs(power) - 1; i++) {
        for (j = frame[fx]; j < stop; j++) {
            vector[code_x++] = vector[j];
            if (code_x >= MAX_VECTOR) {
                Error(ABORT, 0, "Insufficient space to hold input data.");
            }
        }
    }

    vector[length_x] += (abs(power) - 1) * (stop - frame[fx]);
}


int lookup_code(char sym, int power)
{
int i;

    for (i = 0; i< params.no_generators; i++) {
        if (params.generators[i] == sym) {
            if (power < 0) {
                return(i + params.no_generators);
            } else {
                return i;
            }
        }
    }

    fprintf(stderr, "unknown generator: %c^%d\n", sym,power);
    exit(1);
}


void put_code(char sym, int power)
{
int i;
int code;

    code = lookup_code(sym, power);

    for (i=0; i<abs(power); i++) {
        vector[code_x++] = code;
        if (code_x >= MAX_VECTOR) {
            Error(ABORT, 0, "Insufficient space to hold input data.");
        }
    }
    vector[length_x] += abs(power);
}


void init_vector(int *counter_p, int *vector_p)
{
    counter = counter_p;
    vector = vector_p;
    *counter = 0;
    length_x = 0;
    vector[length_x] = 0;
    code_x = 1;
}


void next_word()
{
    ++(*counter);
    length_x = code_x;
    if (++code_x >= MAX_VECTOR) {
        Error(ABORT, 0, "Insufficient space to hold input data.");
    }
    vector[length_x] = 0;
}


void commutator() {
int i,j;
int w1len, w2len;
int w1[MAX_WORD];
int w1inv[MAX_WORD];
int w2[MAX_WORD];
int w2inv[MAX_WORD];

    w1len =  frame[fx] - frame[fx-1];
    if (w1len > MAX_WORD) {
        Error(ABORT, 0, "Insufficient space to hold word.");
    }
    w2len = code_x - frame[fx];
    if (w2len > MAX_WORD) {
        Error(ABORT, 0, "Insufficient space to hold word.");
    }

    // copy the last two frames and produce their inverses
    for (i = 0; i < w1len; i++) {
        w1[i] = w1inv[i] = vector[frame[fx-1] + i];
    }
    invert_word(w1inv, 0, i);

    for (i = 0; i < w2len; i++) {
        w2[i] = w2inv[i] = vector[frame[fx] + i];
    }
    invert_word(w2inv, 0, i);

    // replace the last two frames with the expanded commutator
    i = frame[fx-1] ;
    for (j = 0; j < w1len; j++) {
        vector[i++] = w1inv[j];
    }
    for (j = 0; j < w2len; j++) {
        vector[i++] = w2inv[j];
    }
    for (j = 0; j <w1len; j++) {
        vector[i++] = w1[j];
    }
    for (j = 0; j <w2len; j++) {
        vector[i++] = w2[j];
    }
    // reset parameters in vector[] accordingly
    vector[length_x] += w1len + w2len;
    code_x = i;
}
