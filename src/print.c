
#define PRINT__C
#include "lowx.h"

int printing_in_progress = 0;

// @!
// my_printf should be replaced with a more safe procedure
//
#define my_printf sprintf


static void send_print(struct Parameters *params, char *buf, int len)
{
#ifdef DISTRIBUTED
char tmp[MAXOUTLEN+20];
    if (params->print_server == NULL) {
        printf("%s\n", buf);
        fflush(stdout);
    } else {
        int s;
        s = connect_to_server(params->print_server, params->print_port);
        len = sprintf(tmp, "output %s", buf);
        if (write(s, tmp, len) < len) {
            Error(ABORT, 0, "Can't send print: %s", tmp);
        }
        close(s);
    }
#else
    printf("%s", buf);
    fflush(stdout);
#endif
}


void print_lx(struct Parameters *params)
{
char buf[MAXOUTLEN];
int i,len;

    len = sprintf(buf, "# group=<");
    for (i = 0; i < params->no_generators; i++) {
        if (i > 0) {
            len += my_printf(&buf[len], ",");
        }
        len += my_printf(&buf[len], "%c", params->generators[i]);
    }
    len += my_printf(&buf[len], "|");

    print_words(params, buf, &len, params->no_relators, params->relators);
    len += my_printf(&buf[len], ">;\n");

    if (params->no_subgens > 0) {
        len += sprintf(&buf[len], "# subgens=<");
        print_words(params, buf, &len, params->no_subgens, params->subgens);
        len += my_printf(&buf[len], ">;\n");
    } 

    // this for the mixed adaptation
    if (params->no_subgroup > 0) {
        len += sprintf(&buf[len], "# subgroup=<");
        print_words(params, buf, &len, params->no_subgroup, params->subgroup);
        len += my_printf(&buf[len], ">;\n");
    } 
    
    if (params->normal_flag) {
        len += my_printf(&buf[len], "# normal_subgroups;\n");
    }
    if (params->mixed_flag) {
        len += my_printf(&buf[len], "# mixed;\n");
    }
    len += sprintf(&buf[len], "# max_index=%d;\n", params->max_index);
    len += my_printf(&buf[len], "# max_no_rows=%d;\n", params->max_no_rows);
    if (params->checkpoint_flag) {
        len += my_printf(&buf[len], "# checkpoint;\n");
    }

    // this should always be printed on 'stdout'
    printf("%s\n", buf);
    fflush(stdout);
}


static void print_power(struct Parameters *params, char *buf, int *len, int gen, int n)
{
char gensym;

    if (gen < params->no_generators) {
        gensym = params->generators[gen];
    } else {
        gensym = params->generators[params->inverse[gen]];
        n *= -1;
    }
    if (n == 1) {
        *len += my_printf(&buf[*len], "%c", gensym);
    } else {
        *len += my_printf(&buf[*len], "%c^%d", gensym, n);
    }
}


void print_words(struct Parameters *params, char *buf, int *len, int no_words, int * words)
{
int i,j,k;
int gen,n;

    for (i = j = 0; i < no_words; i++, j += words[j] + 1) {
        if (i > 0) {
            *len += my_printf(&buf[*len], " ");
        }
        if (words[j] == 0) {
            // the empty word stands for the identity
            *len += my_printf(&buf[*len], "1");
            continue;
        }
        for (k = j + 1; k <= j + words[j]; ) {
            if (k > j + 1) {
                *len += my_printf(&buf[*len], "*");
            }
            // collect a power
            gen = words[k];
            n = 0;
            while (words[k] == gen && k <= j + words[j]) {
                ++n;
                ++k;
            }
            print_power(params, buf, len, gen, n);
        }
    }
}


#ifdef REGMAPS
// this is for the "regular maps" problem only
// calculate and print the type of the regular map
// i.e. the orders of generator 'u' and 'v'

static void print_maptype(char *buf, int *len, struct Enum *p)
{
int *column;
int i,j,k,tmp;
int cycle_length;

    if ((column = (int *) malloc(p->no_rows * sizeof(int))) == NULL) {
            Error(ABORT, errno, "Can't allocate memory");
    }

    for (i = 0; i < p->params->no_generators; i++) {
        if (p->params->generators[i] == 't') continue;
        cycle_length = 0;
        for (j = 0; j < p->no_rows; j++) {
            column[j] = p->order[p->row[p->index[j]].column[i]];
        }
        for (j = 0; j < p->no_rows; j++) {
            while (column[j] != j) {
                ++cycle_length;
                k = column[j];
                column[j] = j;
                while ( k != j) {
                    ++cycle_length;
                    tmp = column[k];
                    column[k] = k;
                    k = tmp;
                }
                // in this case all cyles of u or v are of the same length
                // so we need to calculate the length of the first cycle only
                goto next;
            }
        }
next:
        if (p->params->generators[i] == 'u') {
            *len += my_printf(&buf[*len], "p");
        } else {
            *len += my_printf(&buf[*len], "q");
        }
        *len += my_printf(&buf[*len], "=%d ", cycle_length);
    }

    free(column);
}

#endif


void print_result(struct Enum *p)
{
static int counter = 0;
char buf[MAXOUTLEN];
int len;

    printing_in_progress = 1;

    ++counter;
    // index
    len = my_printf(buf, "%d ", p->no_rows);
#ifdef REGMAPS
    print_maptype(buf, &len, p);
#endif
    // subgroup generators or addititonal relators
    print_words(p->params, buf, &len, p->no_subgens, p->subgens);

    len += my_printf(&buf[len], "\n");
    send_print(p->params, buf, len);

#ifndef DISTRIBUTED
    if (p->params->print_mode != PR_INDEX) {
        print_cosetenum_result(p);
    }
#endif

    printing_in_progress = 0;
}


/***  coset-table functions  ***/

void print_cosetenum_result(struct Enum *p)
{
    switch (p->params->print_mode) {
        case PR_COSET_TABLE:
            print_ct(p);
            break;
        case PR_ALL:
            print_permrep(p);
            print_ct(p);
            break;
        case PR_PERM_REP:
            print_permrep(p);
        default:
            if (!p->params->lowindex_flag) {
                printf("index = %d\n", p->no_rows);
            }
            break;
    }
    printf("\n");
}


static void print_ct(struct Enum *p)
{
int i, j, n;
int coset;

    printf("coset-table:\n       ");
    for (i = 0; i < p->params->no_columns; i++) {
        if (i < p->params->no_generators) {
            printf("%5c ", p->params->generators[i]);
        } else {
            printf("%5c_", p->params->generators[p->params->inverse[i]]);
        }
    }
    printf("\n");
    for (n = 0; n < p->no_rows; n++) {
        printf("%4d: ", n + 1);
        for (j = 0; j < p->params->no_columns; j++) {
            coset = p->row[p->index[n]].column[j];
            if (coset == EMPTY) {
                printf("%6s", "-");
            } else {
                printf("%6d", p->order[coset] + 1);
            }
        }
        printf("\n");
    }
}


static void print_permrep(struct Enum *p)
{
int *column;
int i,j,k,tmp;
int len, ident;

    if ((column = (int *) malloc(p->no_rows * sizeof(int))) == NULL) {
            Error(ABORT, errno, "Can't allocate memory");
    }

    printf("cosets permuted by the generators:\n");
    for (i = 0; i < p->params->no_generators; i++) {
        ident = len = printf("  %c --> ", p->params->generators[i]);
        for (j = 0; j < p->no_rows; j++) {
            column[j] = p->order[p->row[p->index[j]].column[i]];
        }
        for (j = 0; j < p->no_rows; j++) {
            while (column[j] != j) {
                len += printf("(%d", j+1);
                k = column[j];
                column[j] = j;
                while ( k != j) {
                    len += printf(",");
                    len += printf("%d", k+1);
                    tmp = column[k];
                    column[k] = k;
                    k = tmp;
                }
                len += printf(")");
            }
        }
        printf("\n");
    }

    free(column);
}


/****  utility functions  ***/

void print_times(FILE *fp)
{
static long clktck = 0;
static struct tms Start;
static clock_t start = 0;
struct tms End;
clock_t end = 0;
char buf[512];
int len;

    if (fp == NULL) {
        // start measuring times
        if ((start = times(&Start)) == -1) Error(ABORT, errno, "times() failed");
        return;
    }

    if (clktck == 0) {
        if ((clktck = sysconf(_SC_CLK_TCK)) < 0)
            Error(ABORT, errno, "sysconf failed");
    }
    
    if ((end = times(&End)) == -1) Error(ABORT, errno, "times() failed");

    len = sprintf(buf, "# times:  real=%.2f  user=%.2f  system=%.2f",
                (end - start) / (double) clktck,
                (End.tms_utime - Start.tms_utime) / (double) clktck,
                (End.tms_stime - Start.tms_stime) / (double) clktck);

    // send_print(p->params, buf, len);
    printf("%s\n", buf);
    fflush(stdout);
}


void print_bye(int sig_number)
// called from signal_handler() at termination
{
    printf("\n# interrupted by signal %d ...\n", sig_number);
    printf("... bye.\n");
}



#ifdef DEBUG

void dump_words(char *words_name, int no_words, int *words)
{
int i, j, k;

    printf("\n%s[%d]:\n", words_name, no_words);
    for (i = j = 0; i < no_words; i++, j += words[j] + 1) {
        printf("%3d:", i);
        for (k = j + 1; k <= j + words[j]; k++) {
            printf(" %d", words[k]);
        }
        printf("\n");
    }
}


void print_rcpr_table(struct Parameters *params, struct Warray *rcpr_table)
{
int i, j, k;
struct Warray *wap;
struct Word *wp;

    printf("\nrcpr_table[%d]:\n", params->no_columns);
    for (i = 0; i < params->no_columns; i++) {
        wap = &(rcpr_table[i]);
        printf("%2d[%d]:\n", i, wap->no_words);
        for (j = 0; j < wap->no_words; j++) {
            printf("    ");
            wp = &(wap->word[j]);
            for (k = 0; k < wp->length; k++) {
                printf(" %d", wp->element[k]);
            }
            printf("\n");
        }
    }
}

#endif
