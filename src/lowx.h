
#ifndef LOWX__H
#define LOWX__H

#include <print_error.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <signal.h>
#ifdef DISTRIBUTED
#  include <sys/types.h>
#  include <sys/socket.h>
#  include <netdb.h>
#  include <netinet/in.h>
#  include <strings.h>
#endif

#define MAX_WORD                2048
#define MAX_VECTOR              8192
#define MAX_NO_GENS             50
#define MAX_NO_ROWS             20000
#define MAX_NO_COINCS           5000
#define MAX_NO_SKIPS            100

typedef enum {
    PR_INDEX, PR_COSET_TABLE, PR_PERM_REP, PR_ALL
} PrintMode;

struct Parameters {
    // lx input parameters
    int no_generators;
    int no_columns;
    char generators[MAX_NO_GENS];
    int inverse[2 * MAX_NO_GENS];
    int no_relators;
    int relators[MAX_VECTOR];
    int no_subgens;
    int subgens[MAX_VECTOR];
    int no_subgroup;
    int subgroup[MAX_VECTOR];
    int max_no_rows;
    int max_no_tests;
    int max_no_coincs;
    int max_index;
    // mode parameters
    int lowindex_flag;
    int normal_flag;
    int mixed_flag;
    int peek_flag;
    int min_nob;
    int max_nob;
    int checkpoint_flag;
    int skips[MAX_NO_SKIPS][2];
    int no_skips;
    PrintMode print_mode;
    // distributed
    int master_flag;
    int slave_flag;
    char *server;
    int port;
    char *print_server;
    int print_port;
#ifdef POSTPROC
    char postproc_code;
#endif
#ifdef DEBUG
    int debug_flag;
#endif
};

typedef enum {
    EMPTY = -1,
    FAILURE, SUCCESS,
    COMPLETE, GAP, CLOSED, DEDUCTION, COINCIDENCE, ABANDON
} Symbol;

// a row of the coset table
struct Row {
    int* column;        // pointer to the actual row which is an array of ints
};


// a two dimensional (no_rows x no_columns) array of Test entries
// the array is a mirror of the cosettable and as such avoids generating
// multiple tests

struct Tests {
    int no_rows;
    int no_columns;
    int no_tests;
    int first;          // first row containing tests
    struct Row *row;
};


// data structures representing the relators and
// the subgroup generators

// a word of the generators and their inverses
struct Word {
    int length;
    int *element;
};

// an array of words of the generators and their inverses
struct Warray {
    int no_words;
    struct Word *word;
};


// this structure incorporates all data specific to the
// enumeration of a given subgroup
struct Enum {

    // input parameters
    struct Parameters * params;

    // the coset table
    int max_no_rows;    // the maximal number of rows
    int no_rows;        // the current number of active rows
    int no_columns;     //
    int last_def; // order of the row where the last coset definition was made
    struct Row *row;    // pointer to the array of rows
    int *free;
    int freep;

    // dynamically allocated arrays holding the order <--> index maps
    //
    // order[max_no_rows]
    // order[i] gives the logical order of the i-th physical row in ct
    // init: order[x] = EMPTY for all x
    int *order;

    // index[max_no_rows + 1]
    // index[o] gives the physical row-index in ct of the o-th logical row
    // init: index[x] = EMPTY for all x, the last entry is a sentinel
    int *index;

    struct Tests *tests;

    // set of <order1, order2> pairs representing equivalent cosetnumbers
    // implemented as a stack
    int (*coincs) [2];
    int coincp;         // stack pointer

    // for labeling the nodes of the search tree
    struct Enum *parent; // points to the parent 'Enum'
    int marker;    // (outer loop) marker for the 'DoneAlready' test ('node_j')
    int node_i;    // inner loop iterating over possible forced coincidences
    int node_k;    // iterate over subgroup elements to be "mixed"
    int branch_n;  // ordinal number of the level 1 branch (0,1,...)

    // the subgroup generators
    int no_subgens;
    int subgens[MAX_VECTOR];
    int length_x;
    int code_x;

    // the RCPR table
    struct Warray *rcpr_table;

    // the RCPR table of the 'new relators'
    // this is used in searching for normal subgroups
    struct Warray *new_rcpr_table;

    // the <coset,generator> pair currently being tested
    int coset;
    int gen;
};


//---------------------------------------------  checkpoint.c  ----------------

void make_checkpoint(int checkpoint_flag) ;
void init_checkpoint(int no_rows, int no_skips, int (* skips)[2]) ;
void set_checkpoint(int branch) ;
int checkpoint_done(int branch) ;

#ifdef CHECKPOINT__C
    typedef enum {
        RAW, IN_PROGRESS, DONE
    } BranchStatus;

    typedef struct {
        char    *done;          // will be allocated as char done[nob]
        int     nob;            // number of level 1 branches
        int     changed;        // progress indicator
    } checkpoint_t;

    static void signal_handler(int sig_number) ;
    static void open_chkp_file(void) ;
    static void checkpoint(FILE *fp) ;
#endif


//---------------------------------------------  conjugate.c  -----------------

int conjugate(struct Enum *p) ;


//---------------------------------------------  cosetenum.c  -----------------

struct Enum * make_enum(int size, struct Parameters *params) ;
void free_enum(struct Enum *p) ;
void init_enum(struct Enum *p) ;
struct Enum * copy_enum(struct Enum *p) ;
int (* make_coincs(int max_no_coincs)) [2] ;
Symbol coincidence(struct Enum *p, int a, int b) ;
Symbol cosetenum(struct Enum *p) ;
Symbol mixed_initial_cosetenum(struct Enum *p) ;
void plain_cosetenum(struct Parameters *params) ;

#ifdef COSETENUM__C
    static int add_row(struct Enum *p) ;
    static void set_ct_entry(struct Enum *p, int a, int g, int b) ;
    static Symbol def_new_coset(struct Enum *p) ;
    static void push_coinc(struct Enum *p, int A, int B) ;
    static Symbol pop_coinc(int *A, int *B, struct Enum *p) ;
    static void replace_in_ct(struct Enum *p, int a, int b) ;
    static void delete_row(struct Enum *p, int b) ;
    static void replace_in_coincs(struct Enum *p, int A, int B) ;
    static Symbol coset_x_word(struct Enum *p, int coset, struct Word *w) ;
    static Symbol scan_relators(struct Enum *p) ;
    static void scan_subgens(struct Enum *p) ;
#endif


//---------------------------------------------  gens_init.c  -----------------

void gens_init(struct Parameters *params) ;

#ifdef GENS_INIT__C
#endif


//---------------------------------------------  lowindex.c  ------------------

void add_subgen(struct Enum *p, int a, int b) ;
void add_mixed_relator(struct Enum *p, int a, int b, int wlen, int *w) ;
void descend(struct Enum *p) ;
void mixed_descend(struct Enum *p) ;
void descend_root(struct Enum *p) ;
void lowindex(struct Parameters *params) ;
#ifndef DISTRIBUTED
    void peek(struct Enum *p, int min_nob, int max_nob) ;
#endif

#ifdef LOWINDEX__C
    static void find_definition(int *b, int *g, struct Enum *p, int a) ;
    static void print_node_labels(struct Enum *p) ;
    static void extend_tests(struct Enum *p) ;
#endif


//---------------------------------------------  distributed.c  ---------------

void descend_master(struct Enum *p) ;
void descend_slave(struct Enum *p) ;
void mixed_descend_slave(struct Enum *p) ;
int connect_to_server(char *server, int port) ;

#ifdef DISTRIBUTED__C
    typedef enum {
        JOB, NO_MORE_JOB, INCOMPLETE
    } Message;
    static Message get_job(char *server, int port, int *n, int *j, int *i) ;
    static void job_done(char *server, int port, int n, int j, int i) ;
#endif

//---------------------------------------------  postproc.c  ------------------

void postproc(struct Enum *p) ;

#ifdef POSTPROC__C
#endif

//---------------------------------------------  parser_aux.c  ----------------

void new_frame(void);
void remove_frame(void);
void invert_word(int *w, int from, int to);
void append_frame(int power);
int lookup_code(char sym, int power);
void put_code(char sym, int power);
void init_vector(int *counter_p, int *vector_p);
void next_word(void);
void commutator(void);


//---------------------------------------------  print.c  ---------------------

// the size of the buffer in which we build an output line
#define MAXOUTLEN 5000

void print_lx(struct Parameters *params) ;
void print_words(struct Parameters *params, char *buf, int *len, int no_words, int * words) ;
void print_result(struct Enum *p) ;
void print_cosetenum_result(struct Enum *p) ;
void print_times(FILE *fp) ;
void print_bye(int sig_number) ;

#ifndef PTHREAD
    extern int printing_in_progress;
#endif

#ifdef PRINT__C
    static void send_print(struct Parameters *params, char *buf, int len) ;
    static void print_power(struct Parameters *params, char *buf, int *len, int gen, int n) ;
#   ifdef REGMAPS
        static void print_maptype(char *buf, int *len, struct Enum *p) ;
#   endif
    static void print_ct(struct Enum *p);
    static void print_permrep(struct Enum *p);
#endif

#ifdef DEBUG
    void dump_words(char *words_name, int no_words, int *words) ;
    void print_rcpr_table(struct Parameters *params, struct Warray *rcpr_table);
#endif


//---------------------------------------------  rcpr.c  ----------------------

struct Warray * make_rcpr_table(int no_columns, int no_relators, int *relators);
void free_rcpr_table(struct Warray *rcpr_t, int no_columns) ;

#ifdef RCPR__C
    static int compare_words(struct Word *p, struct Word *q) ;
    static void sort_words(struct Warray *wap) ;
    static void make_reduced_cyclic_permutations
        (struct Warray *wap, int no_relators, int *relators) ;
#endif


//---------------------------------------------  tests.c  -----------------

struct Tests * make_tests(int no_rows, int no_columns) ;
void reset_tests(struct Tests *t) ;
void delete_tests(struct Tests *t) ;
void copy_tests(struct Tests * old, struct Tests * new) ;
struct Tests * dup_tests(struct Tests * old) ;
int tests_empty(struct Tests *t) ;
void push_test(struct Tests *t, int a, int g) ;
Symbol pop_test(int *a, int *g, struct Tests *t) ;
void replace_in_tests(struct Tests *t, int a, int b) ;

#ifdef TESTS__C
#endif

//-----------------------

#endif
