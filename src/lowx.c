
char program_version[] = "1.1";
char copyright [] = "Copyright (©) 2022 Peter Dobsan\n\
This program comes with absolutely no warranty.";

#define SHORT_VERSION   1
#define LONG_VERSION    2

#define LOWX__C
#include "lowx.h"

#ifdef DISTRIBUTED

char usage_distributed[] ="\n\
distributed lowx\n\n\
    usage: lowx [-v] {-s|-m} lx_file\n\n\
        -m 'master' producing the first level branches as 'jobs'\n\
        -s specify the 'server host' for distributed lowx\n\
";

#elif defined POSTPROC

char usage_postproc[] ="\n\
lowx for post-processing\n\n\
    usage: lowx [-v] -xN  lx_file\n\n\
        N=d duality check for regular maps\n\
        N=D duality check for chiral maps\n\
        N=c reflection check for chiral rotary maps\n\
        N=1 mix1 check for symmetric trivalent graphs\n\
        N=2 mix2 check for symmetric trivalent graphs\n\
        N=3 mix3 check for symmetric trivalent graphs\n\
        N=5 mix5 check for symmetric trivalent graphs\n\
        N=6 mix6 check for symmetric trivalent graphs\n\
";

#else

char usage[] ="\n\
usage: lowx [-v -bn,m -p{t|p|a}] lx_file\n\n\
        -pt print the full coset-table\n\
        -pp print permutations of cosets induced by the generators\n\
        -pa print both coset-table and permutation representation\n\
            default: print only the index and the extra relators/generators\n\
        -bn,m dump the first level branches from branch n to branch m\n\
        -v display version information and exit\n\
";

#endif


extern FILE *yyin;
int yyparse(void);

struct Parameters params;

void init_params(struct Parameters *p)
{
    p->no_generators            = 0;
    p->no_columns               = 0;
    p->no_relators              = 0;
    p->no_subgens               = 0;
    p->no_subgroup              = 0;
    p->max_no_rows              = 0;
    p->max_no_coincs            = MAX_NO_COINCS;
    p->max_index                = MAX_NO_ROWS;
    p->print_mode               = PR_INDEX;
    p->lowindex_flag            = 1;
    p->normal_flag              = 0;
    p->mixed_flag               = 0;
    p->peek_flag                = 0;
    p->min_nob                  = 0;
    p->max_nob                  = 0;
    p->checkpoint_flag          = 0;
    p->no_skips                 = 0;
    p->master_flag              = 0;
    p->slave_flag               = 0;
    p->port                     = 0;
    p->server                   = NULL;
    p->print_port               = 0;
    p->print_server             = NULL;
#ifdef POSTPROC
    p->postproc_code            = '\0';
#endif
#ifdef DEBUG
    p->debug_flag               = 0;
#endif
}



static void Usage();
static void Version(int quit);


int main(int argc, char *argv[])
{
int c;
extern char *optarg;
extern int optind;
int error_flag = 0;

#ifdef PTMALLOC
    ptmalloc_init();
#endif

    // to help remote start
    fclose(stdin);

    init_params(&params);

    while ((c = getopt(argc, argv, "b:p:vVms:x:D")) != EOF) {
        switch (c) {
            case 'p':
                switch (optarg[0]) {
                    case 't':
                        params.print_mode = PR_COSET_TABLE;
                        break;
                    case 'p':
                        params.print_mode = PR_PERM_REP;
                        break;
                    case 'a':
                        params.print_mode = PR_ALL;
                        break;
                    default:
                        params.print_mode = PR_INDEX;
                        break;
                }
                break;
            case 'x':
#ifdef POSTPROC
                params.postproc_code = optarg[0];
#else
                Error(ABORT,0, "You can't use '-x, this lowx is not 'postproc'.");
#endif
                break;
            case 'b':
#ifdef DISTRIBUTED
                Error(ABORT, 0, "You can't use '-b', this lowx is 'distributed'.");
#endif
#ifdef POSTPROC
                Error(ABORT, 0, "You can't use '-b', this lowx is 'postproc'.");
#endif
                params.peek_flag = 1;
                sscanf(optarg,"%d,%d", &params.min_nob, &params.max_nob);
                break;
            case 'm':
#ifdef DISTRIBUTED
                params.master_flag = 1;
#else
                Error(ABORT, 0,
                    "You can't use '-m', this lowx is not 'distributed'.");
#endif
                break;
            case 's':
#ifdef DISTRIBUTED
                params.server = optarg;
                params.print_server = optarg;
                params.slave_flag = 1;
#else
                Error(ABORT, 0,
                    "You can't use '-s', this lowx is not 'distributed'.");
#endif
                break;
            case 'v':
                Version(SHORT_VERSION);
                break;
            case 'V':
                Version(LONG_VERSION);
                break;
            case 'D':
#ifdef DEBUG
                params.debug_flag = 1;
#else
                Error(ABORT, 0,
                "You can't use '-D', lowx was not built with DEBUG defined.");
#endif
                break;
            case '?':
                error_flag++;
        }
    }

    if (error_flag || argc < 2) {
        Usage();
    }

    if (optind == argc) {
        Error(ABORT, 0, "Missing input file");
    }

#ifdef DISTRIBUTED
    if (!params.master_flag && !params.slave_flag) {
        Error(ABORT,0,"Distributed 'lowx' should be either 'master' or 'slave'.");
    }
#endif

    // start measuring times
    print_times(NULL);

    if ((yyin = fopen(argv[optind], "r")) == NULL) {
        Error(ABORT, errno, "Can't open file");
    }

    if (yyparse()) {
        exit(1);
    }
    fclose(yyin);

    if (params.mixed_flag) {
        if (params.no_subgens <= 0) {
            Error(ABORT, 0, "'mixed' is given but 'subgens' is missing.");
        }

        if (params.no_subgroup <= 0) {
            Error(ABORT, 0, "'mixed' is given but 'subgroup' is missing.");
        }
    }

    gens_init(&params);

    if (params.max_no_rows <= 0) {
        params.max_no_rows = (int) (params.max_index * 1.3 + 1);
    } else if (params.max_no_rows <= params.max_index) {
        Error(ABORT, 0, "'max_no_rows' should be greater than 'max_index'.");
    }

#ifndef POSTPROC
    print_lx(&params);
#endif

    if (params.checkpoint_flag || params.no_skips > 0) {
        // create the checkpoint data structure
        make_checkpoint(params.checkpoint_flag);
    }

    if (params.lowindex_flag) {
        lowindex(&params);
    } else {
        plain_cosetenum(&params);
    }

    print_times(stdout);

    exit(0);
}


static void Usage()
{
#ifdef DISTRIBUTED
    fprintf(stderr, "%s\n", usage_distributed);
#elif defined POSTPROC
    fprintf(stderr, "%s\n", usage_postproc);
#else
    fprintf(stderr, "%s\n", usage);
#endif
    exit(1);
}


static void Version(int mode)
{
    fprintf(stderr,"lowx version %s\n", program_version);
    fprintf(stderr,"%s\n", copyright);
    if (mode == SHORT_VERSION) exit(0);
    fprintf(stderr,"compiler: ");
#if defined __GNUC__
    fprintf(stderr,"gcc %d.%d\n", __GNUC__, __GNUC_MINOR__);
#elif defined __SUNPRO_C
    fprintf(stderr,"Sun WorkShop Compilers C 4.2\n");
#else
    fprintf(stderr,"unknown\n");
#endif
    exit(0);
}
