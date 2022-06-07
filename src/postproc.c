
/*
 * This module contains various auxiliary routines for postprocessing tasks.
 * The procedures assume some particular fixed order of generators (columns)
 * in the coset-tables.
 *
 * The general idea is to perform some transformation on the coset-table
 * renumber the cosets in normal order then compare the new and
 * the original coset-tables.
 *
 * regular maps problem:
 *      check for duality:
 *              swap 'u' and 'v'
 *              if original < new  --> distinct
 *              if original == new --> self-dual
 *              if original > new  --> reject
 *
 * trivalent symmetric graphs (mixed):
 *      filter out some duplications caused by outer automorphisms:
 *              make a mix specific transformation and the renumbering
 *              if original <= new  --> OK
 *              if original > new  --> reject
 *
 */

#define POSTPROC__C
#include "lowx.h"

#define NO_ROWS         2048
#define NO_COLUMNS      32

struct coset_table {
    int no_rows;
    int no_columns;
    int ct[NO_ROWS][NO_COLUMNS];
} ;


/***  common routines  *******************************************************/

static void
copy_ct(struct Enum *p, struct coset_table *t1, struct coset_table *t2)
// transform/duplicate the coset-table of 'p' into two 'coset-table'
// data structure: 't1' and 't2'
// all the other routines deal with these new coset-tables
{
int i, j, k;

    if (p->no_rows > NO_ROWS || p->no_columns > NO_COLUMNS) {
        Error(ABORT, 0, "%dx%d coset-table does not fit in %dx%d.",
                p->no_rows, p->no_columns, NO_ROWS, NO_COLUMNS);
    }

    t1->no_rows = t2->no_rows = p->no_rows;
    t1->no_columns = t2->no_columns = p->no_columns;

    for (i = 0; i < p->no_rows; i++) {
        for (j = 0; j < p->no_columns; j++) {
            k = p->row[p->index[i]].column[j];
            k = p->order[k];
            t1->ct[i][j] = k;
            t2->ct[i][j] = k;
        }
    }
}


static void swap_columns(struct coset_table *t, int i, int j)
{
int n;
int tmp;

    for (n = 0; n < t->no_rows; n++) {
        tmp = t->ct[n][i];
        t->ct[n][i] = t->ct[n][j];
        t->ct[n][j] = tmp;
    }
}


static int multiply(struct coset_table *o, int coset, int f, int g)
// multiply 'coset' with the word 'f * g' from the right
// where 'f' and 'g' generators
// return the coset resulted
{
int x;

    x = o->ct[coset][f];
    return o->ct[x][g];
}


static void find_new_order(struct coset_table *t, int *order)
{
int i, j, k;
int flag[NO_ROWS];
int coset;

    for (i = 0; i < t->no_rows; i++) {
        flag[i] = 1;
    }
    flag[0] = 0;        // the first coset (0) is fixed
    order[0] = 0;
    k = 1;

    for (i = 0; i < t->no_rows; i++) {
        for (j = 0; j < t->no_columns; j++) {
            coset = t->ct[order[i]][j];
            if (flag[coset]) {
                order[k++] = coset;
                flag[coset] = 0;
            }
        }
    }
}


static void
rearrange_rows(struct coset_table *t, struct coset_table *v, int *order)
{
int i, j;

    v->no_rows = t->no_rows;
    v->no_columns = t->no_columns;

    for (i = 0; i < t->no_rows; i++) {
        for (j = 0; j < t->no_columns; j++) {
            v->ct[i][j] = t->ct[order[i]][j];
        }
    }
}


static void print_order(int *order, int n)
{
int i;

    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%d:%d ", i+1, order[i]+1);
    }
    printf("\n");
}


static void renumber(struct coset_table *t, int *order)
{
int i, j;
int new[NO_ROWS];
int coset;

    for (i = 0; i < t->no_rows; i++) {
        new[order[i]] = i;
    }

    for (i = 0; i < t->no_rows; i++) {
        for (j = 0; j < t->no_columns; j++) {
            coset = t->ct[i][j];
            t->ct[i][j] = new[coset];
        }
    }
}


static void print_coset_table(struct coset_table *t)
{
int i, j;

    for (i = 0; i < t->no_rows; i++) {
        printf("%4d: ", i + 1);
        for (j = 0; j < t->no_columns; j++) {
            printf("%6d", t->ct[i][j] + 1);
        }
        printf("\n");
    }
    printf("\n");
}


static int compare(struct coset_table *t, struct coset_table *v)
{
int i, j;

    for (i = 0; i < t->no_rows; i++) {
        for (j = 0; j < t->no_columns; j++) {
            if (t->ct[i][j] == v->ct[i][j]) {
                continue;
            } else if (t->ct[i][j] < v->ct[i][j]) {
                return -1;
            } else {
                return 1;
            }
        }
    }

    return 0;
}


/***  regular maps duality check  ********************************************/

static void swap_u_v(struct coset_table *t)
// swap 'u' and 'v'
// assume that the columns of the coset-table are in the order below
// t  u  v  u_  v_      <-- generators
// 0  1  2  3   4       <-- column indices
{
int n;
int tmp;

    for (n = 0; n < t->no_rows; n++) {
        tmp = t->ct[n][1];
        t->ct[n][1] = t->ct[n][2];
        t->ct[n][2] = tmp;
        tmp = t->ct[n][3];
        t->ct[n][3] = t->ct[n][4];
        t->ct[n][4] = tmp;
    }
}


/***  reflecting transformation for chiral rotary maps  **********************/

static void chiral_reflect(struct coset_table *t)
// tranformation:
//      u --> u^-1
//      v --> v^-1
// assume that the columns of the coset-table are in the order below
// u  v  u_  v_         <-- generators
// 0  1  2   3          <-- column indices
{
int n;
int tmp;

    for (n = 0; n < t->no_rows; n++) {
        tmp = t->ct[n][0];
        t->ct[n][0] = t->ct[n][2];
        t->ct[n][2] = tmp;
        tmp = t->ct[n][1];
        t->ct[n][1] = t->ct[n][3];
        t->ct[n][3] = tmp;
    }
}


/*** chiral duality test  ****************************************************/

static void chiral_swap_u_v(struct coset_table *t)
// swap 'u' and 'v'
// assume that the columns of the coset-table are in the order below
// u  v  u_  v_         <-- generators
// 0  1  2   3          <-- column indices
{
int n;
int tmp;

    for (n = 0; n < t->no_rows; n++) {
        tmp = t->ct[n][0];
        t->ct[n][0] = t->ct[n][1];
        t->ct[n][1] = tmp;
        tmp = t->ct[n][2];
        t->ct[n][2] = t->ct[n][3];
        t->ct[n][3] = tmp;
    }
}


/***  symmetric trivalent graphs -- mixed adaptation  ************************/

/***  mix1  -----------------------------------------------------------------*/
// transformation: h <--> h^-1
// we assume that the coset-table has the following form
//  h  a  h_            <-- generators
//  0  1  2             <-- column indices

static void mix1_transform(struct coset_table *t)
{
    swap_columns(t, 0, 2);
}


/***  mix2  -----------------------------------------------------------------*/
// transformation: a <--> a*p
// we assume that the coset-table has the following form
//  h  p  a  h_         <-- generators
//  0  1  2  3          <-- column indices

static void mix2_transform(struct coset_table *o, struct coset_table *t)
{
int i;

    for (i = 0; i < t->no_rows; i++) {
        t->ct[i][2] = multiply(o, i, 2, 1);
    }
}


/***  mix3  -----------------------------------------------------------------*/
// transformations:
//      a   --> a*p
//      a_  --> p*a_
// and these can more simply be expressed as:
//      a <--> a_       (just interchanging the two columns)
//
// since from relations a*a*p=1 and p*p=1 we get a*p = a_  and also  a*a = p
// so a --> a*p == a --> a_
// and p*a_ = a*a*a_ = a  so  a_ --> p*a_ == a_ --> a
// 
// we assume that the coset-table has the following form
//  h  p  a  h_ a_      <-- generators
//  0  1  2  3  4       <-- column indices

static void mix3_transform(struct coset_table *o, struct coset_table *t)
{
    swap_columns(t, 2, 4);
}


/***  mix5  -----------------------------------------------------------------*/
// transformation: a <--> a*p
// we assume that the coset-table has the following form
//  h  p  q  r  a  h_           <-- generators
//  0  1  2  3  4  5            <-- column indices

static void mix5_transform(struct coset_table *o, struct coset_table *t)
{
int i;

    for (i = 0; i < t->no_rows; i++) {
        t->ct[i][4] = multiply(o, i, 4, 1);
    }
}


/***  mix6  -----------------------------------------------------------------*/
// transformations (as in mix3): a <--> a_
// we assume that the coset-table has the following form
//  h  p  q  r  a  h_ a_        <-- generators
//  0  1  2  3  4  5  6         <-- column indices

static void mix6_transform(struct coset_table *o, struct coset_table *t)
{
    swap_columns(t, 4, 6);
}


/***  P O S T P R O C  *******************************************************/

void postproc(struct Enum *p)
{
char code = p->params->postproc_code;
struct coset_table original;
struct coset_table t;
int order[NO_ROWS];
struct coset_table new;
int cmp;

    copy_ct(p, &original, &t);
#ifdef DEBUG
    print_lx(p->params);
    printf("original ");
    print_cosetenum_result(p);
    // printf("pp-original:\n");
    // print_coset_table(&original);
#endif

    if (code == 'd') {
        swap_u_v(&t);
    } else if (code == 'D') {
        chiral_swap_u_v(&t);
    } else if (code == 'c') {
        chiral_reflect(&t);
    } else if (code == '1') {
        mix1_transform(&t);
    } else if (code == '2') {
        mix2_transform(&original, &t);
    } else if (code == '3') {
        mix3_transform(&original, &t);
    } else if (code == '5') {
        mix5_transform(&original, &t);
    } else if (code == '6') {
        mix6_transform(&original, &t);
    } else {
        Error(ABORT, 0, "Invalid post-processing mode: %c", code);
    }

#ifdef DEBUG
    printf("\ntransformed:\n");
    print_coset_table(&t);
#endif

    find_new_order(&t, order);
    rearrange_rows(&t, &new, order);
    renumber(&new, order);

#ifdef DEBUG
    printf("renumbered in normal order:\n");
    print_coset_table(&new);
    printf("comparing the original and the new coset-table: ");
#endif

    cmp = compare(&original, &new);
    if (cmp < 0) {
        printf("less\n");
    } else if (cmp == 0) {
        printf("equal\n");
    } else {
        printf("greater\n");
    }

    exit(0);
}
