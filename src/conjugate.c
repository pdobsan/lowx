
#define CONJUGATE__C
#include "lowx.h"

#ifdef OLD_CONJUGATE

int conjugate(struct Enum *p)
{
int i,j,k,m,n;
int *original;  // conjugate coset-number --> original coset-number map;
int *conjugate; // original coset-number --> conjugate coset-number map;
int *po, *pn;
int x, y;
int return_value = 0;

    // allocate the maps
    if ((original = (int *) malloc(p->no_rows * sizeof(int))) == NULL) {
        Error(ABORT, errno, "can't allocate original[]");
    }
    if ((conjugate = (int *) malloc(p->no_rows * sizeof(int))) == NULL) {
        Error(ABORT, errno, "can't allocate conjugate[]");
    }

    for (i = 1; i < p->no_rows; i++) {
        // clear the maps
        for (j = 0; j < p->no_rows; j++) {
            original[j] = conjugate[j] = EMPTY;
        }

        // row i ==> conjugate row 0
        original[0] = i;
        conjugate[i] = 0;

        // build the maps
        for (j = 1; j < p->no_rows; j++) {
            // find the next coset-number to be conjugate coset-number j
            // these loops (k,m) scan the original rows in the conjugate order
            for (k = 0; k < j; k++) {
                for (m = 0; m < no_columns; m++) {
                    pn = &(p->row[p->index[original[k]]].column[m]);
                    if (*pn > EMPTY) {
                        // check if we already renumbered *pn
                        // i.e. if it is in original[]
                        y = p->order[*pn];
                        for (n = 0; n < j; n++) {
                            if (y == original[n]) {
                                // renumbered already
                                goto label_1;
                            }
                        }
                        // found the next coset-number
                        original[j] = y;
                        conjugate[y] = j;
                        goto label_2;
                    } else {
                        // we have encountered a gap, can't build maps further
                        goto lexico_test;
                    }
                    label_1: ;
                }
            }
            label_2: ;
        }

        lexico_test:

        // compare the original and the conjugate renumbered coset-table
        // in lexicographical order
        for (j = 0; j < p->no_rows; j++) {
            for (k = 0; k < no_columns; k++) {
                if (original[j] == EMPTY) {
                    // because the j-th row could not have been renumbered
                    // we have to abandon the conjugacy check
                    return_value = 0;
                    goto EXIT;
                }
                po = &(p->row[p->index[j]].column[k]);
                pn = &(p->row[p->index[original[j]]].column[k]);
                // not sure about these EMPTY checks ???
                if (*po != EMPTY && *pn != EMPTY) {
                    x = p->order[*po];
                    y = p->order[*pn];
                    if (conjugate[y] != EMPTY) {
                        if (conjugate[y] == x) {
                            continue;
                        } else if (conjugate[y] < x) {
                            // the conjugate table precedes the original one
                            // so the current subgroup is a conjugate of one
                            // found earlier
                            return_value = 1;
                            goto EXIT;
                        } else {
                            // the original table precedes the conjugate one
                            goto next;
                        }
                    }
                }
                // the coset-tables are not comparable
                goto next;
            }
        }
        next: ;
    }

EXIT:
    free(original);
    free(conjugate);
    return return_value;
}

#elif defined CTTR_CONJUGATE
//----------------------------------------------------------------------------

//@?
// This was an attempt to clean up the conjugate test
// using routines from postproc.c.
// -- unfinished --

#define NO_ROWS         2048
#define NO_COLUMNS      32

struct coset_table {
    int no_rows;
    int no_columns;
    int ct[NO_ROWS][NO_COLUMNS];
} ;

#define CT_INCOMPARABLE
#define CT_LESS
#define CT_EQUAL
#define CT_GREATER


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


static int compare(struct coset_table *t, struct coset_table *v)
{
int i, j;

    for (i = 0; i < t->no_rows; i++) {
        for (j = 0; j < t->no_columns; j++) {
            if (t->ct[i][j] == EMPTY || v->ct[i][j] == EMPTY) {
                // some times partial coset-tables are incomparable.
                return CT_INCOMPARABLE;
            } else if (t->ct[i][j] == v->ct[i][j]) {
                continue;
            } else if (t->ct[i][j] < v->ct[i][j]) {
                return CT_LESS;
            } else {
                return CT_GREATER;
            }
        }
    }

    return CT_EQUAL;
}


int conjugate(struct Enum *p) {
struct coset_table original;
struct coset_table t;
int order[NO_ROWS];
struct coset_table new;
int i;
int cmp;

    copy_ct(p, &original, &t);

    for (i = 1; i < p->no_rows; i++) {
        find_new_order(&t, order);
        rearrange_rows(&t, &new, order);
        renumber(&new, order);

        cmp = compare(&original, &new);

        if (cmp == CT_INCOMPARABLE) {
            return 0;
        } else if (cmp == CT_GREATER) {
            return 1;
        }
    }

    return 0;
}

//----------------------------------------------------------------------------
#else

int conjugate(struct Enum *p)
{
int i, j, k, m;
int compare_from_row, compare_from_column;
int new_from_row, new_from_column;
int *orig;      // conjugate coset-number --> original coset-number map;
int *conj;      // original coset-number --> conjugate coset-number map;
int *po, *pn;
int x, y;
int retval = 0;

    // allocate the maps
    if ((orig = (int *) malloc(p->no_rows * sizeof(int))) == NULL) {
        Error(ABORT, errno, "Can't allocate orig[]");
    }
    if ((conj = (int *) malloc(p->no_rows * sizeof(int))) == NULL) {
        Error(ABORT, errno, "Can't allocate conj[]");
    }

    for (i = 1; i < p->no_rows; i++) {
        // clear the maps
        for (j = 0; j < p->no_rows; j++) {
            orig[j] = conj[j] = EMPTY;
        }

        // row i ==> conjugate row 0
        orig[0] = i;
        conj[i] = 0;

        compare_from_row = 0;
        compare_from_column = 0;

        new_from_row = 0;
        new_from_column = 0;

        pn = &(p->row[p->index[i]].column[0]);
        if (*pn == EMPTY) {
            // gap --> coset-tables are not comparable
            goto next_coset_table;
        }

        // some times it happens that after taking the initial row
        // the first element in the new row can be renumbered immediatelly
        // so we can start comparing the coset-tables right here
        //
        if ((y = p->order[*pn]) == i) {         // i.e. conj[y] == 0
            x = p->order[p->row[p->index[0]].column[0]];
            // no need to check, x (the first entry) can never be EMPTY
            if (0 < x) {
                // the conjugate table precedes the original one
                retval = 1;
                goto EXIT;
            }
        }

        // build the maps
        for (j = 1; j < p->no_rows; j++) {
            // find the next coset-number to be conjugate coset-number j
            // these loops (k,m) scan the original rows in the conjugate order
            for (k = new_from_row; k < j; k++) {
                for (m = new_from_column; m < p->no_columns; m++) {
                    pn = &(p->row[p->index[orig[k]]].column[m]);
                    if (*pn > EMPTY) {
                        // check if we already renumbered *pn that is
                        // if the corresponding entry in conj[] is set
                        y = p->order[*pn];
                        if (conj[y] != EMPTY) {
                            continue;
                        }
                        // found the next coset-number
                        orig[j] = y;
                        conj[y] = j;
                        new_from_column = m + 1;
                        if (new_from_column >= p->no_columns) {
                            new_from_row = k + 1;
                            new_from_column = 0;
                        } else {
                            new_from_row = k;
                        }
                        goto next_coset;
                    } else {
                        // we have encountered a gap, can't build maps further
                        goto lexico_test;
                    }
                }
                new_from_column = 0;
            }
            next_coset: ;
        }

        lexico_test:

        // compare the original and the conjugate renumbered coset-table
        // in lexicographical order
        for (j = 0; j < p->no_rows; j++) {
            for (k = 0; k < p->no_columns; k++) {
                if (orig[j] == EMPTY) {
                    printf("===>\n");
                    // because the j-th row could not have been renumbered
                    // we have to abandon the conjugacy check ???
                    goto next_coset_table;
                }
                po = &(p->row[p->index[j]].column[k]);
                pn = &(p->row[p->index[orig[j]]].column[k]);
                // not sure about these EMPTY checks ???
                if (*po != EMPTY && *pn != EMPTY) {
                    x = p->order[*po];
                    y = p->order[*pn];
                    if (conj[y] != EMPTY) {
                        if (conj[y] == x) {
                            continue;
                        } else if (conj[y] < x) {
                            // the conjugate table precedes the original one
                            // so the current subgroup is a conjugate of one
                            // found earlier
                            retval = 1;
                            goto EXIT;
                        } else {
                            // the original table precedes the conjugate one
                            goto next_coset_table;
                        }
                    }
                }
                // the coset-tables are not comparable
                goto next_coset_table;
            }
        }
        next_coset_table: ;
    }

  EXIT:
    free(conj);
    free(orig);
    return retval;
}

#endif
