
#define CHECKPOINT__C
#include "lowx.h"

// make a checkpoint in every ALARM_PERIOD seconds
#define ALARM_PERIOD    600

static checkpoint_t *chkp;
static FILE    *fp;

static void signal_handler(int sig_number)
{
    switch (sig_number) {
        case SIGHUP:
            if (signal(SIGHUP,signal_handler) == SIG_ERR) {
                Error(ABORT, errno, "signal(SIGHUP) error");
            }
            checkpoint(fp);
            break;
        case SIGALRM:
            if (signal(SIGALRM,signal_handler) == SIG_ERR) {
                Error(ABORT, errno, "signal(SIGALRM) error");
            }
            alarm(ALARM_PERIOD);
            checkpoint(fp);
            break;
        case SIGUSR1:
            break;
        case SIGUSR2:
            break;
        case SIGINT:
        case SIGTERM:
        default:
            fprintf(fp, "\n# interrupted by signal %d ...\n", sig_number);
            checkpoint(fp);
            fprintf(fp, "# use the last skip list to restart\n");
            if (printing_in_progress) {
                fprintf(fp, "# WARNING!!!\n");
                fprintf(fp, "# printing was in progress at termination.\n");
            }
            fprintf(fp, "# ... bye.\n");
            print_bye(sig_number);
            exit(1);
    }
}


static void open_chkp_file()
{
pid_t   pid;
char    buf[512];

    pid = getpid();
    sprintf(buf, "%d.chkp", (int) pid);
    fp = fopen(buf, "w");
    if (fp == NULL) Error(ABORT, errno, "Can't open log file.");
}


void make_checkpoint(int checkpoint_flag)
// called from main()
// create the checkpoint data structure
// then open the checkpoint file
{
    // we need this part for both skip lists and checkpointing
    chkp = (checkpoint_t *) malloc(sizeof(checkpoint_t));
    if (chkp == NULL) Error(ABORT, errno, "Can't allocate memory for chkp");

    chkp->done = NULL;
    chkp->nob = 0;
    chkp->changed = 0;

    if (checkpoint_flag) {
        // the rest is needed only for checkpointing
        if (signal(SIGHUP, signal_handler) == SIG_ERR) {
            Error(ABORT, errno, "signal(SIGHUP) error");
        }
        if (signal(SIGINT, signal_handler) == SIG_ERR) {
            Error(ABORT, errno, "signal(SIGINT) error");
        }
        if (signal(SIGTERM, signal_handler) == SIG_ERR) {
            Error(ABORT, errno, "signal(SIGTERM) error");
        }
        if (signal(SIGALRM, signal_handler) == SIG_ERR) {
            Error(ABORT, errno, "signal(SIGALRM) error");
        } else {
            alarm(ALARM_PERIOD);
        }
        open_chkp_file();
    }
}


void init_checkpoint(int no_rows, int no_skips, int (* skips)[2])
// allocate and initilaize the chkp->done[] array
{
int i, j;

    if (chkp == NULL)
        Error(ABORT, 0, "Can't initialize checkpoint, for 'chkp' is NULL");

    chkp->nob = (no_rows - 1) * no_rows;
    chkp->nob /= 2;

    chkp->done = (char *) malloc(chkp->nob * sizeof(char));
    if (chkp->done == NULL)
        Error(ABORT, errno, "Can't allocate memory for chkp->done[]");

    for (i = 0; i < chkp->nob; i++) chkp->done[i] = RAW;

    if (no_skips > 0) {
        for (i = 0; i < no_skips; i++) {
            if (skips[i][0] > chkp->nob || skips[i][1] > chkp->nob ||
                skips[i][0] > skips[i][1]) Error(ABORT, 0, "Invalid skip list.");
            for (j = skips[i][0]; j < skips[i][1]; j++) {
                chkp->done[j] = DONE;
            }
        }
        chkp->changed = 1;
        printf("# initial skip list\n");
        checkpoint(stdout);
    }
}


void set_checkpoint(int branch)
{
    chkp->done[branch] = DONE;
    chkp->changed = 1;
}


int checkpoint_done(int branch)
{
    return (chkp->done[branch] == DONE);
}


static void checkpoint(FILE *fp)
{
int i;
int from, to;
int raw_flag;
int done;

    if (chkp->changed) {
        fprintf(fp, "# skip ");
        done = 0;
        from = to = 0;
        raw_flag = 1;
        for (i = 0; i < chkp->nob; i++) {
            if (raw_flag) {
                from = to = i;
            }
            if (chkp->done[i] == DONE) {
                ++to;
                ++done;
                raw_flag = 0;
            } else {
                raw_flag = 1;
                if (from != to) {
                    if (from+1 == to) {
                        fprintf(fp, "%d ", to);
                    } else {
                        fprintf(fp, "%d-%d ", from+1, to);
                    }
                }
            }
        }
        if (!raw_flag) {
            if (from+1 == to) {
                fprintf(fp, "%d ", to);
            } else {
                fprintf(fp, "%d-%d ", from+1, to);
            }
        }
        fprintf(fp, ";\n");

        fprintf(fp, "# skipped %d (%.2f%%) out of %d\n",
            done, (double) done * 100.0 / (double) chkp->nob, chkp->nob);
        fflush(fp);

        chkp->changed = 0;
    }
}
