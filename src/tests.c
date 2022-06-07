
#define TESTS__C
#include "lowx.h"

struct Tests * make_tests(int no_rows, int no_columns)
{
struct Tests *t;
int i;

    // struct Tests
    if ((t = (struct Tests *) malloc(sizeof(struct Tests))) == NULL) {
        Error(ABORT, errno, "Can't allocate tests.");
    }
    t->no_rows = no_rows;
    t->no_columns = no_columns;
    t->no_tests = 0;
    t->first = 0;

    // allocate array of pointers to the rows
    if ((t->row = (struct Row *) malloc(no_rows * sizeof(struct Row)))
                                                                    == NULL) {
        Error(ABORT, errno, "Can't allocate tests.");
    }

    // allocate the rows
    for (i = 0; i < no_rows; i++) {
        if ((t->row[i].column = (int *) malloc(no_columns * sizeof(int)))
                                                                    == NULL) {
            Error(ABORT, errno, "Can't allocate memory for rows");
        }
    }
 
    reset_tests(t);
    return t;
}


void reset_tests(struct Tests *t)
{
int i,j;

    for (i = 0; i < t->no_rows; i++) {
        for (j = 0; j < t->no_columns; j++) {
            t->row[i].column[j] = 0;
        }
    }

    t->first = 0;
    t->no_tests = 0;
}


void delete_tests(struct Tests *t)
{
int i;

    for (i = 0; i < t->no_rows; i++) {
        free(t->row[i].column);
    }
    free(t->row);
    free(t);
}


void copy_tests(struct Tests * old, struct Tests * new)
{
int i,j;

    for (i = old->first; i < old->no_rows; i++) {
        for (j = 0; j < old->no_columns; j++) {
            new->row[i].column[j] = old->row[i].column[j];
        }
    }

    new->first = old->first;
    new->no_tests = old->no_tests;
}


struct Tests * dup_tests(struct Tests * old)
{
struct Tests * new;

    new = make_tests(old->no_rows, old->no_columns);
    copy_tests(old, new);
    return new;
}


int tests_empty(struct Tests *t)
{
    return (t->no_tests <= 0);
}


void push_test(struct Tests *t, int a, int g)
{
    if (t->row[a].column[g] == 0) {
        if (t->no_tests == 0 || a < t->first) {
            t->first = a;
        }
        t->no_tests++;
        t->row[a].column[g] = 1;
    }
}


Symbol pop_test(int *a, int *g, struct Tests *t)
{
int i,j;

    if (t->no_tests > 0) {
        for (i = t->first; i < t->no_rows; i++) {
            for (j = 0; j < t->no_columns; j++) {
                if (t->row[i].column[j] == 1) {
                    t->row[i].column[j] = 0;
                    t->no_tests--;
                    if (t->no_tests == 0) {
                        t->first = 0;
                    } else {
                        t->first = i;
                    }
                    *a = i;
                    *g = j;
                    return SUCCESS;
                }
            }
        }
    }

    return FAILURE;
}


void replace_in_tests(struct Tests *t, int a, int b)
{
int i;
int change_flag = 0;

    for (i = 0; i < t->no_columns; i++) {
        if (t->row[b].column[i] == 1) {
            if (t->row[a].column[i] == 1) {
                t->no_tests--;
            } else {
                t->row[a].column[i] = 1;
                change_flag = 1;
            }
            t->row[b].column[i] = 0;
        }
    }

    if (change_flag && a < t->first) {
        t->first = a;
    }
}
