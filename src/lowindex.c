
#define LOWINDEX__C
#include "lowx.h"


static void find_definition(int *b, int *g, struct Enum *p, int a)
//
// to a given 'a' (>0) coset find the defining equation
// in the coset table in the form of "a = bg"
{
int i,j;
int n;

    for (n = 0; n < p->no_rows; n++) {
        i = p->index[n];
        for (j = 0; j < p->no_columns; j++) {
            if (p->row[i].column[j] == a) {
                *b = i;
                *g = j;
                return;
            }
        }
    }

    *b = *g = EMPTY;
}


static void extend_tests(struct Enum *p)
{
int i;
int gen;

    gen = p->subgens[p->length_x - 1];

    for (i = 0; i < p->no_rows; i++) {
        push_test(p->tests, p->index[i], gen);
    }
}


void add_subgen(struct Enum *p, int a, int b)
{
int x, g;
int i, j, tmp, half;

    p->subgens[p->length_x] = 0;

    // left side
    x = p->index[a];
    while (x > 0) {
        find_definition(&x, &g, p, x);
        p->subgens[p->code_x++] = g;
        if (p->code_x >= MAX_VECTOR) {
            Error(ABORT, 0, "Insufficient space to hold subgens.");
        }
        p->subgens[p->length_x] ++;
    }
    // reverse
    half = (p->length_x + 1 + p->code_x) / 2;   // who knows how smart CC is :-)
    for (i = p->length_x + 1, j = p->code_x - 1; i < half; i++, j--) {
        tmp = p->subgens[i];
        p->subgens[i] = p->subgens[j];
        p->subgens[j] = tmp;
    }

    // right side
    x = p->index[b];
    while (x > 0) {
        find_definition(&x, &g, p, x);
        p->subgens[p->code_x++] = p->params->inverse[g];
        if (p->code_x >= MAX_VECTOR) {
            Error(ABORT, 0, "Insufficient space to hold subgens.");
        }
        p->subgens[p->length_x] ++;
    }

    p->no_subgens++;
    p->length_x = p->code_x++;
    if (p->code_x >= MAX_VECTOR) {
        Error(ABORT, 0, "Insufficient space to hold subgens.");
    }

    extend_tests(p);
}


void add_mixed_relator(struct Enum *p, int a, int b, int wlen, int *w)
{
int x, g;
int i, j, tmp, half;
int subgens[MAX_VECTOR];
int length_x;
int code_x;

    length_x = 0;
    code_x = 1;
    subgens[length_x] = 0;

    // left side
    x = p->index[a];
    while (x > 0) {
        find_definition(&x, &g, p, x);
        subgens[code_x++] = g;
        if (code_x >= MAX_VECTOR) {
            Error(ABORT, 0, "Insufficient space to hold subgens.");
        }
        subgens[length_x] ++;
    }
    // reverse
    half = (length_x + 1 + code_x) / 2; // who knows how smart CC is :-)
    for (i = length_x + 1, j = code_x - 1; i < half; i++, j--) {
        tmp = subgens[i];
        subgens[i] = subgens[j];
        subgens[j] = tmp;
    }

    // right side
    x = p->index[b];
    while (x > 0) {
        find_definition(&x, &g, p, x);
        subgens[code_x++] = p->params->inverse[g];
        if (code_x >= MAX_VECTOR) {
            Error(ABORT, 0, "Insufficient space to hold subgens.");
        }
        subgens[length_x] ++;
    }

    for (i = 0; i < wlen; i++) {
        p->subgens[p->code_x++] = w[i];
    }
    for (i = 1; i <= subgens[0]; i++) {
        p->subgens[p->code_x++] = subgens[i];
    }

    p->subgens[p->length_x] = wlen + subgens[0];
    p->length_x += wlen + subgens[0] + 1;
    p->code_x = p->length_x + 1;

    if (p->code_x >= MAX_VECTOR) {
        Error(ABORT, 0, "Insufficient space to hold subgens.");
    }
    p->no_subgens++;

    extend_tests(p);
}


#define MAX_TREE 100

static void print_node_labels(struct Enum *p) {
    struct Enum * stack[MAX_TREE];
    int n = -1;

    while (p->parent != NULL) {
        if (++n >= MAX_TREE) {
            Error(ABORT, 0, "Tree is too deep: %dn", n);
        }
        stack[n] = p;
        p = p->parent;
    }
    // stack[++n] = p;  // @ root node

    printf("#");
    while (n >= 0) {
        if (stack[n]->params->mixed_flag) {
            printf(" (%d %d %d)", stack[n]->marker, stack[n]->node_i, stack[n]->node_k);
        } else {
            printf(" (%d %d)", stack[n]->marker, stack[n]->node_i);
        }
        --n;
    }
    printf("\n");
}


void mixed_descend(struct Enum *p) {
int i, j, boundary;
struct Enum *q;
int cosetenum_result;
int k, n;
int *h;
int no_h;
#ifndef DISTRIBUTED
    static int first_flag = 1;
#endif

    h = p->params->subgroup;
    no_h = p->params->no_subgroup;

#ifndef DISTRIBUTED
    if (first_flag) {
        if((cosetenum_result = mixed_initial_cosetenum(p)) == ABANDON) {
            return;
        }
        first_flag = 0;
    } else {
#endif
        if((cosetenum_result = cosetenum(p)) == ABANDON) {
            return;
        }
#ifndef DISTRIBUTED
    }
#endif

    if (!(p->params->normal_flag)) {
        if (conjugate(p)) {
            return;
        }
    }

    if (cosetenum_result == CLOSED && p->no_rows <= p->params->max_index) {
        print_result(p);
        // print_node_labels(p);
#ifdef DEBUG
        if (p->params->debug_flag) {
            printf("-------------------------------------------------\n\n");
        }
#endif
    }

    boundary = (p->params->max_index + 1 < p->no_rows) ?
                                        p->params->max_index + 1 : p->no_rows;

    for (j = p->marker; j < boundary; j++) {
        for (i = 0; i < j; ++i) {
            // iterate over the subgroup elements
            for (k = n = 0; k < no_h; k++, n += h[n] + 1) {
                q = copy_enum(p);
                // after copy_enum() q->tests and p->tests point
                // to the same 'tests set' data structure
                q->marker = j;
                q->node_i = i;
                q->node_k = k;
                if (h[n] > 0) {
                    add_mixed_relator(q, i, j, h[n], &(h[n+1]));
                } else {
                    add_subgen(q, i, j);
                }
                if (p->params->normal_flag) {
                    q->new_rcpr_table =
                    make_rcpr_table(q->no_columns, q->no_subgens, q->subgens);
                }
                if (coincidence(q, q->index[i], q->index[j]) == COINCIDENCE) {
                    mixed_descend(q);
                } else {
                    reset_tests(q->tests);
                }
                free_enum(q);
            }
        }
    }
}


void descend(struct Enum *p)
{
int i, j, boundary;
struct Enum *q;
int cosetenum_result;

        // print_node_labels(p);
    if((cosetenum_result = cosetenum(p)) == ABANDON) {
        return;
    }

    if (!(p->params->normal_flag)) {
        if (conjugate(p)) {
            return;
        }
    }

    if (cosetenum_result == CLOSED && p->no_rows <= p->params->max_index) {
        print_result(p);
#ifdef DEBUG
        if (p->params->debug_flag) {
            printf("-------------------------------------------------\n\n");
        }
#endif
    }

    boundary = (p->params->max_index + 1 < p->no_rows) ?
                                        p->params->max_index + 1 : p->no_rows;

    for (j = p->marker; j < boundary; j++) {
        for (i = 0; i < j; ++i) {
            q = copy_enum(p);
            // after copy_enum() q->tests and p->tests point
            // to the same 'tests set' data structure
            q->marker = j;
            q->node_i = i;
            add_subgen(q, i, j);
            if (p->params->normal_flag) {
                q->new_rcpr_table =
                    make_rcpr_table(q->no_columns, q->no_subgens, q->subgens);
            }
            if (coincidence(q, q->index[i], q->index[j]) == COINCIDENCE) {
                descend(q);
            } else {
                reset_tests(q->tests);
            }
            free_enum(q);
        }
    }
}


void mixed_descend_root(struct Enum *p)
{
int i, j, n, boundary;
struct Enum *q;
int cosetenum_result;
int k, m;
int *h;
int no_h;
#ifndef DISTRIBUTED
    static int first_flag = 1;
#endif

    h = p->params->subgroup;
    no_h = p->params->no_subgroup;

#ifndef DISTRIBUTED
    if (first_flag) {
        if((cosetenum_result = mixed_initial_cosetenum(p)) == ABANDON) {
            return;
        }
        first_flag = 0;
    } else {
#endif
        if((cosetenum_result = cosetenum(p)) == ABANDON) {
            return;
        }
#ifndef DISTRIBUTED
    }
#endif

    if (!(p->params->normal_flag)) {
        if (conjugate(p)) {
            return;
        }
    }

    if (cosetenum_result == CLOSED && p->no_rows <= p->params->max_index) {
        print_result(p);
#ifdef DEBUG
        if (p->params->debug_flag) {
            printf("-------------------------------------------------\n\n");
        }
#endif
    }

    boundary = (p->params->max_index + 1 < p->no_rows) ?
                                        p->params->max_index + 1 : p->no_rows;

    if (!p->params->slave_flag &&
        (p->params->checkpoint_flag || p->params->no_skips > 0)) {
        init_checkpoint(boundary, p->params->no_skips, p->params->skips);
    }

    n = -1;
    for (j = p->marker; j < boundary; j++) {
        for (i = 0; i < j; ++i) {
            ++n;
            if (!p->params->slave_flag &&
                (p->params->checkpoint_flag || p->params->no_skips > 0)) {
                if (checkpoint_done(n)) continue;
            }

            for (k = m = 0; k < no_h; k++, m += h[m] + 1) {

                q = copy_enum(p);
                reset_tests(q->tests);  // @? 
                q->marker = j;
                
                if (h[m] > 0) {
                    add_mixed_relator(q, i, j, h[m], &(h[m+1]));
                } else {
                    add_subgen(q, i, j);
                }

                if (p->params->normal_flag) {
                    q->new_rcpr_table =
                        make_rcpr_table(q->no_columns, q->no_subgens, q->subgens);
                }

                if (coincidence(q, q->index[i], q->index[j]) == COINCIDENCE) {
                    q->branch_n = n;
                    mixed_descend(q);
                } else {
                    reset_tests(q->tests);
                }
                if (p->params->checkpoint_flag) {
                    set_checkpoint(n);
                }
                free_enum(q);
            }
        }
    }
}


void descend_root(struct Enum *p)
{
int i, j, n, boundary;
struct Enum *q;
int cosetenum_result;

    if((cosetenum_result = cosetenum(p)) == ABANDON) {
        return;
    }

    if (!(p->params->normal_flag)) {
        if (conjugate(p)) {
            return;
        }
    }

    if (cosetenum_result == CLOSED && p->no_rows <= p->params->max_index) {
        print_result(p);
#ifdef DEBUG
        if (p->params->debug_flag) {
            printf("-------------------------------------------------\n\n");
        }
#endif
    }

    boundary = (p->params->max_index + 1 < p->no_rows) ?
                                        p->params->max_index + 1 : p->no_rows;

    if (!p->params->slave_flag &&
        (p->params->checkpoint_flag || p->params->no_skips > 0)) {
        init_checkpoint(boundary, p->params->no_skips, p->params->skips);
    }

    n = -1;
    for (j = p->marker; j < boundary; j++) {
        for (i = 0; i < j; ++i) {
            ++n;
            if (!p->params->slave_flag &&
                (p->params->checkpoint_flag || p->params->no_skips > 0)) {
                if (checkpoint_done(n)) continue;
            }
            // after copy_enum() q->tests and p->tests will point
            // to the same 'tests set' data structure
            q = copy_enum(p);
            // entering into the loops first time we may loose a few
            // tests here inherited from the initial cosetenum above. ???
            // this can be more significant when descend_root() is called
            // descned_slave() in the distributed multi-threaded case. ???
            reset_tests(q->tests);
            q->marker = j;
            add_subgen(q, i, j);
            if (p->params->normal_flag) {
                q->new_rcpr_table =
                    make_rcpr_table(q->no_columns, q->no_subgens, q->subgens);
            }

            if (coincidence(q, q->index[i], q->index[j]) == COINCIDENCE) {
                q->branch_n = n;
                if (p->params->mixed_flag) {
                    mixed_descend(q);
                } else {
                    descend(q);
                }
            } else {
                reset_tests(q->tests);
            }
            if (p->params->checkpoint_flag) {
                set_checkpoint(n);
            }
            free_enum(q);
        }
    }
}


void lowindex(struct Parameters *params)
{
struct Enum *enump;

    enump = make_enum(params->max_no_rows, params);
    init_enum(enump);
    enump->rcpr_table = make_rcpr_table(params->no_columns,
                                        params->no_relators, params->relators);
    enump->tests = make_tests(params->max_no_rows, params->no_columns);
    enump->coincs = make_coincs(params->max_no_coincs);

#ifdef DISTRIBUTED

    if (params->master_flag) {
        descend_master(enump);
    } else {
        if (params->mixed_flag) {
            mixed_descend_slave(enump);
        } else {
            descend_slave(enump);
        }
    }
 
#else

    if (params->peek_flag) {
        peek(enump, params->min_nob, params->max_nob);
        printf("\n");
    } else {
        if (params->checkpoint_flag || params->no_skips > 0) {
            // in the original MIXED case descend_root was never called
            // since neither checkpoint nor skips were allowed.
            if (params->mixed_flag) {
                mixed_descend_root(enump);
            } else {
                descend_root(enump);
            }
        } else {
            if (params->mixed_flag) {
                mixed_descend(enump);
            } else {
                descend(enump);
            }
        }
    }

#endif

    printf("# The End\n");
}


#ifndef DISTRIBUTED

void peek(struct Enum *p, int min_nob, int max_nob)
{
int i, j, n;
struct Enum *q;
int cosetenum_result;
char outbuf[MAXOUTLEN];

    if((cosetenum_result = cosetenum(p)) == ABANDON) {
        return;
    }

    if (!(p->params->normal_flag)) {
        if (conjugate(p)) {
            return;
        }
    }

    if (cosetenum_result == CLOSED) {
        print_result(p);
#ifdef DEBUG
        if (p->params->debug_flag) {
            printf("-------------------------------------------------\n\n");
        }
#endif
    }

    printf("total number of branches on level 1: %d\n\n",
            (p->params->max_index * (p->params->max_index + 1)) / 2);

    n = 1;
    for (j = p->marker; j < p->params->max_index + 1; j++) {
        for (i = 0; i < j; ++i) {
            if (min_nob <= n) {
                q = copy_enum(p);
                // after copy_enum() q->tests and p->tests point
                // to the same 'tests set' data structure
                q->marker = j;
                add_subgen(q, i, j);
                if (q->params->normal_flag) {
                    q->new_rcpr_table = make_rcpr_table(q->no_columns,
                                                    q->no_subgens, q->subgens);
                }
                if (coincidence(q, q->index[i], q->index[j]) == COINCIDENCE) {
                    int len = 0;
                    // print subgens
                    printf("%7d:  %d = %d  -->  ", n, j+1, i+1);
                    print_words(p->params,outbuf,&len,q->no_subgens,q->subgens);
                    printf("%s\n", outbuf);
                    fflush(stdout);
                }
                reset_tests(q->tests);
                free_enum(q);
            }
            if (++n > max_nob) return;
        }
    }
}

#endif
