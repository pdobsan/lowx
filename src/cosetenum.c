
#define COSETENUM__C
#include "lowx.h"


struct Enum * make_enum(int size, struct Parameters *params)
//
// dynamically allocate an 'Enum' structure and all of its components
//
// 'size' is the maximal number of rows in the coset-table
//
{
struct Enum *p;
int i;

    // struct Enum
    if ((p=(struct Enum *)malloc(sizeof(struct Enum))) == NULL) {
        Error(ABORT, errno, "Can't allocate memory for Enum.");
    }

    // coset-table
    // allocate array of pointers to the rows
    if ((p->row = (struct Row *) malloc(size * sizeof(struct Row))) == NULL) {
        Error(ABORT, errno, "Can't allocate memory for coset-table.");
    }
    p->max_no_rows = size;
    p->params = params;
    p->no_columns = params->no_columns;

    // allocate the rows for the coset-table
    for (i = 0; i < size; i++) {
        if ((p->row[i].column = (int *) malloc(p->no_columns * sizeof(int)))
            == NULL) {
            Error(ABORT, errno, "Can't allocate memory for rows");
        }
    }

    // free[p->max_no_rows]
    if ((p->free = (int *) malloc(size * sizeof(int))) == NULL) {
        Error(ABORT, errno, "Can't allocate free[]");
    }
 
    // order[p->max_no_rows]
    if ((p->order = (int *) malloc(size * sizeof(int))) == NULL) {
            Error(ABORT, errno, "Can't allocate order[]");
    }

    // index[p->max_no_rows + 1]
    if ((p->index = (int *) malloc((size + 1) * sizeof(int))) == NULL) {
            Error(ABORT, errno, "Can't allocate index[]");
    }

    p->tests = NULL;
    p->coincs = NULL;
 
    return p;
}


void free_enum(struct Enum *p)
    // NOTE:
    // free_enum() DO NOT free the 'coincs' and 'tests' data structures
    // associated with the 'Enum'
{
int i;

    for (i = 0; i < p->max_no_rows; i++) {
        free(p->row[i].column);
    }
    free(p->row);

    free(p->free);
    free(p->order);
    free(p->index);

    if (p->params->normal_flag) {
        free_rcpr_table(p->new_rcpr_table, p->no_columns);
    }

    free(p);
}


void init_enum(struct Enum *p)
// initialize an 'Enum'
{
int i,j;
int *q;

    // stack of free rows
    for (i = 0, j = p->max_no_rows - 1; j > 0; i++, j--) {
        p->free[i] = j;
    }
    p->freep = i;

    /* the first row is always active */
    p->no_rows = 1;
    p->last_def = 0;
    q = p->row[0].column;
    for (j = 0; j < p->no_columns; j++) {
        *q++ = EMPTY;
    }

    // stack of free rows
    for (i = 0, j = p->max_no_rows - 1; j > 0; i++, j--) {
        p->free[i] = j;
    }
    p->freep = i;

    // order[] & index[]
    for (i = 0; i < p->max_no_rows; i++) {
        p->order[i] = p->index[i] = EMPTY;
    }
    p->index[p->max_no_rows] = EMPTY;
    p->order[0] = p->index[0] = 0;

    p->coincp = 0;

    // node labels
    p->parent = NULL;
    p->marker = 1;
    p->node_i = -1;
    p->node_k = -1;
    p->branch_n = 0;

    p->no_subgens = 0;
    p->length_x = 0;
    p->code_x = 1;

    p->rcpr_table = NULL;
    p->new_rcpr_table = NULL;
}


struct Enum * copy_enum(struct Enum *p)
    // NOTE:
    // after copy_enum() the new and old 'Enum'
    // SHARE the same 'tests set' data structure
{
struct Enum *q;
int i, j;

    q = make_enum(p->max_no_rows, p->params);

    // coset-table
    q->no_rows = p->no_rows;
    q->no_columns = p->no_columns;
    q->last_def = p->last_def;

    for (i = 0; i < p->max_no_rows; i++) {
        for (j = 0; j < p->no_columns; j++) {
            q->row[i].column[j] = p->row[i].column[j];
        }
    }

    q->freep = p->freep;
    for (i = 0; i < p->freep; i++) {
        q->free[i] = p->free[i];
    }

    // order[] & index[]
    for (i = 0; i < p->max_no_rows; i++) {
        q->order[i] = p->order[i];
        q->index[i] = p->index[i];
    }
    q->index[p->max_no_rows] = p->index[p->max_no_rows];

    // coincs
    q->coincs = p->coincs;
    q->coincp = 0;

    q->tests = p->tests;

    // node labels
    q->parent = p;
    q->marker = p->marker;
    q->node_i = p->node_i;
    q->node_k = p->node_k;
    q->branch_n = p->branch_n;

    // subgroup generators;
    q->no_subgens = p->no_subgens;
    q->length_x = p->length_x;
    q->code_x = p->code_x;
    for (i = 0; i < p->length_x; i++) {         // ???
        q->subgens[i] = p->subgens[i];
    }

    q->rcpr_table = p->rcpr_table;
    return q;
}


int (* make_coincs(int max_no_coincs)) [2]
{
int (*p)[2];

    // coincs[max_no_coincs][2]
    if ((p = (int(*)[2]) malloc(max_no_coincs * sizeof(int [2]))) == NULL) {
        Error(ABORT, errno, "Can't allocate coincs[][2]");
    }
    return p;
}


//
// procedures dealing with the coset table
//

static int add_row(struct Enum *p)
{
int new_row, i;
int *q;

    // get a free row
    if (p->freep <= 0) {
        Error(ABORT, 0, "No more free rows.");
    }
    new_row = p->free[--(p->freep)];

    // update order <--> index maps
    p->order[new_row] = p->no_rows;
    p->index[p->no_rows] = new_row;

    p->no_rows++;

    // clear the new row
    q = p->row[new_row].column;
    for (i = 0; i < p->no_columns; i++) {
        *q++ = EMPTY;
    }

    return new_row;
}


static void set_ct_entry(struct Enum *p, int a, int g, int b)
{
    p->row[a].column[g] = b;
    p->row[b].column[p->params->inverse[g]] = a;
    push_test(p->tests, a, g);
    push_test(p->tests, b, p->params->inverse[g]);
}


static Symbol def_new_coset(struct Enum *p)
//
// make a new definition in the first (row, column order)
// empty place in the coset table
// the definition has the form:
//      Ag = B
//      where B is the new cosetnumber
//
// return
//      SUCCESS if a new definition has been made
//      FAILURE if the coset table is closed
{
int a, b;
int g;
int n;

    for (n = p->last_def; n < p->no_rows; n++) {
        a = p->index[n];
        for (g = 0; g < p->no_columns; g++) {
            if (p->row[a].column[g] == EMPTY) {
                p->last_def = n;
                b = add_row(p);
                set_ct_entry(p, a, g, b);
                return SUCCESS;
            }
        }
    }
    return FAILURE;
}


//
// COINCIDENCE processing
//

static void push_coinc(struct Enum *p, int A, int B)
{
    if (p->coincp >= p->params->max_no_coincs) {
        Error(ABORT, 0, "The coinc stack is full. p = %d  x = %d\n",
                p->coincp, p->params->max_no_coincs);
    }
    p->coincs[p->coincp][0] = A;
    p->coincs[p->coincp++][1] = B;
}


static Symbol pop_coinc(int *A, int *B, struct Enum *p)
{
    // skip the 'dead' coincidences already processed
    do {
        if (p->coincp <= 0) {
            return FAILURE;
        }
        p->coincp--;
    } while (p->coincs[p->coincp][0] == p->coincs[p->coincp][1]);

    *A = p->coincs[p->coincp][0];
    *B = p->coincs[p->coincp][1];
    return SUCCESS;
}


static void replace_in_ct(struct Enum *p, int a, int b)
//
// replace all 'b' with 'a' in the coset-table
// the input parameters 'a' and 'b' are indices
//
{
int i,j;

#ifdef DEBUG
    if (p->params->debug_flag) {
        printf("-D coincidence: %d = %d\n", a+1, b+1);
    }
#endif

    for (i = 0; i < p->no_rows; i++) {
        for (j = 0; j < p->no_columns; j++) {
            if (p->row[p->index[i]].column[j] == b) {
                p->row[p->index[i]].column[j] = a;
                push_test(p->tests, p->index[i], j);
            }
        }
    }
}


static void delete_row(struct Enum *p, int b)
// delete the 'b'-th row,
{
int i;
int n;

    if (p->max_no_rows <= b) {
        Error(ABORT, 0, "Can't delete nonexisting row.");
    } else if (p->no_rows == 1) {
        Error(ABORT, 0, "The last remaining row cannot be deleted.");
    } else if (p->order[b] == EMPTY) {
        Error(ABORT, 0, "Can't delete unlinked row.");
    }

    p->no_rows--;

    n = p->order[b];
    if (n == p->no_rows) {      // i.e. we are deleting the last row
        if (n == p->last_def) {
            p->last_def--;
        }
    } else {
        if (n < p->last_def) {
            p->last_def--;
        }
    }

    // adjust index <--> order maps
    for (i = p->order[b]; p->index[i] != EMPTY; i++) {
        (p->order[ p->index[i] ])--;
        p->index[i] = p->index[i+1];
    }
    p->order[b] = EMPTY;

    // add row 'b' to the stack of free rows
    // ??? !!! that was wrong: using 'max_no_rows' only
    if (p->freep >= p->max_no_rows) {
        Error(ABORT, 0, "Stack of free rows is full.");
    }
    p->free[p->freep++] = b;
}


static void replace_in_coincs(struct Enum *p, int A, int B)
{
int i;

    for (i = 0; i < p->coincp; i++) {
        if (p->coincs[i][0] == B) {
            p->coincs[i][0] = A;
        }
        if (p->coincs[i][1] == B) {
            p->coincs[i][1] = A;
        }
    }
}


Symbol coincidence(struct Enum *p, int a, int b)
{
int x, y;
int X, Y;
int A, B;
int i;

    while(1) {
        A = p->order[a];
        B = p->order[b];

        if (p->params->lowindex_flag) {
            // during adjusting the coincidence-set
            // a coincidence can get such what we should reject
            // that is why we need this check here too
            if (A < B) {
                if (B < p->marker) {
                    return ABANDON;
                }
            } else {
                if (A < p->marker) {
                    return ABANDON;
                }
            }
        }

        if (B < A) {
            // swap the physical rows if it is necessary
            i = b;
            b = a;
            a = i;
        }

        replace_in_ct(p, a, b);

        // scan row 'a' & 'b'
        for (i = 0; i < p->no_columns; i++) {
            x = p->row[a].column[i];
            y = p->row[b].column[i];
            if (y == EMPTY) {
                // do nothing
                continue;
            } else if (x == EMPTY) {
                // can make a deduction
                set_ct_entry(p, a, i, y);
            } else { // the i-th positions are filled in both rows
                if (x == y) {
                    // there is no new information
                    continue;
                } else {
                    // we have just found another coincidence
                    //
                    // we do not check here if we already have
                    // this coincidence, just store the new coincidence
                    //
                    X = p->order[x];
                    Y = p->order[y];
                    if (p->params->lowindex_flag) {
                        if (X < Y) {
                            if (Y < p->marker) {
                                return ABANDON;
                            }
                        } else {
                            if (X < p->marker) {
                                return ABANDON;
                            }
                        }
                    }
                    push_coinc(p, x, y);
                }
            }
        }

        delete_row(p, b);

        replace_in_coincs(p, a, b);
        replace_in_tests(p->tests, a, b);
        if (p->coset == b) {
            p->coset = a;
        }

        if (pop_coinc(&a, &b, p) == FAILURE) {
            return COINCIDENCE;
        }
    }
}


static Symbol coset_x_word(struct Enum *p, int coset, struct Word *w)
//
// multiplying a coset by a word of generators and their inverses on the right
//
{
int a,b;
int A,B;
int i,j;

    /* scanning forward */
    a = coset;
    for (i = 0; i < w->length; i++) {
        if (p->row[a].column[w->element[i]] == EMPTY) {
            goto backward;
        } else {
            a = p->row[a].column[w->element[i]];
        }
    }

    if (a == coset) {
        return CLOSED;
    } else {
        A = p->order[a];
        B = p->order[coset];
        if (A < B) {
            if (B < p->marker) return ABANDON;
        } else {
            if (A < p->marker) return ABANDON;
        }
        return coincidence(p, a, coset);
    }

  backward: /* scanning bacward */
    b = coset;
    for (j = w->length-1; i <= j; --j) {
        if (p->row[b].column[p->params->inverse[w->element[j]]] == EMPTY) {
            goto deduction_or_gap;
        } else {
            b = p->row[b].column[p->params->inverse[w->element[j]]];
        }
    }

    if (a == b) {
        return COMPLETE;
    } else {
        A = p->order[a];
        B = p->order[b];
        if (A < B) {
            if (B < p->marker) return ABANDON;
        } else {
            if (A < p->marker) return ABANDON;
        }
        return coincidence(p, a, b);
    }

  deduction_or_gap:
    if (i == j) {
        set_ct_entry(p, a, w->element[i], b);
        return DEDUCTION;
    } else {
        return GAP;
    }
}


static Symbol scan_relators(struct Enum *p)
{
int coset, gen;
int i;
int result;

    while (pop_test(&coset, &gen, p->tests) == SUCCESS) {
        p->coset = coset;
        p->gen = gen;

        // scanning the original relators
        for (i = 0; i < p->rcpr_table[gen].no_words; i++) {
            result = coset_x_word(p, coset, &(p->rcpr_table[gen].word[i]));
            if (result == COINCIDENCE) {
                // the cosettable has changed
                push_test(p->tests, p->coset, p->gen);
                goto coincidence_break;
            } else if (result == ABANDON) {
                return ABANDON;
            }
        }

        // scanning the 'new relators' if there is any
        if (p->new_rcpr_table != NULL) {
            for (i = 0; i < p->new_rcpr_table[gen].no_words; i++) {
                result = coset_x_word(p, coset,
                                      &(p->new_rcpr_table[gen].word[i]));
                if (result == COINCIDENCE) {
                    // the cosettable has changed
                    push_test(p->tests, p->coset, p->gen);
                    goto coincidence_break;
                } else if (result == ABANDON) {
                    return ABANDON;
                }
            }
        }

        coincidence_break: ;
    }

    return COMPLETE;
}


static void scan_subgens(struct Enum *p)
{
int i,j;
struct Word gen;

    for (i = j = 0; i < p->params->no_subgens;
                                        i++, j += p->params->subgens[j] + 1) {
        gen.length = p->params->subgens[j];
        gen.element = &(p->params->subgens[j + 1]);
        coset_x_word(p, 0, &gen);
    }
}


Symbol cosetenum(struct Enum *p)
{
    do {
        if (tests_empty(p->tests)) {
            if (def_new_coset(p) == FAILURE) {
                return CLOSED;
            }
        }
        if (scan_relators(p) == ABANDON) {
            reset_tests(p->tests);
            return ABANDON;
        }
#ifdef DEBUG
        if (p->params->debug_flag) {
            printf("-D ");
            print_cosetenum_result(p);
        }
#endif
    } while (p->no_rows < p->max_no_rows);

    return COMPLETE;
}


static Symbol mixed_scan_subgens(struct Enum *p)
{
int i,j;
struct Word gen;

    for (i = j = 0; i < p->params->no_subgens;
                                        i++, j += p->params->subgens[j] + 1) {
        gen.length = p->params->subgens[j];
        gen.element = &(p->params->subgens[j + 1]);
        if (coset_x_word(p, 0, &gen) == ABANDON) return ABANDON;
    }

    return COMPLETE;
}


Symbol mixed_initial_cosetenum(struct Enum *p)
{
    do {
        if (tests_empty(p->tests)) {
            if (def_new_coset(p) == FAILURE) {
                return CLOSED;
            }
        }
        if (mixed_scan_subgens(p) == ABANDON) {
            reset_tests(p->tests);
            return ABANDON;
        }
        if (scan_relators(p) == ABANDON) {
            reset_tests(p->tests);
            return ABANDON;
        }
    } while (p->no_rows < p->max_no_rows);

    return COMPLETE;
}


void plain_cosetenum(struct Parameters *params)
{
struct Enum *p;

    p = make_enum(params->max_no_rows, params);
    init_enum(p);
    p->rcpr_table = make_rcpr_table
        (params->no_columns, params->no_relators, params->relators);
    p->tests = make_tests(params->max_no_rows, params->no_columns);
    p->coincs = make_coincs(params->max_no_coincs);

    while(1) {
        if (def_new_coset(p) == FAILURE) {
#ifdef POSTPROC
            postproc(p);
#else
            print_cosetenum_result(p);
#ifdef DEBUG
            printf("-------------------------------------------------\n\n");
#endif
#endif
            break;
        } else {
            // this check is to print out partial coset-tables
            // which contain only restricted number of rows
            if (p->no_rows >= params->max_no_rows) {
                printf("\n# Warning: coset-table has reached the maximum size.\n");
                print_cosetenum_result(p);
                break;
            }
        }
#ifdef DEBUG
        // if (params->debug_flag) {
        if (p->params->debug_flag) {
            printf("-D ");
            print_cosetenum_result(p);
        }
#endif
        scan_subgens(p);
        scan_relators(p);
    }
}
