
#define DISTRIBUTED__C
#include "lowx.h"


void descend_master(struct Enum *p)
{
int i, j, n, boundary;
struct Enum *q;
int cosetenum_result;

    if (p->params->mixed_flag) {
        cosetenum_result = mixed_initial_cosetenum(p);
    } else {
        cosetenum_result = cosetenum(p);
    }
    if(cosetenum_result == ABANDON) {
        return;
    }

    if (!(p->params->normal_flag)) {
        if (conjugate(p)) {
            return;
        }
    }

    if (cosetenum_result == CLOSED && p->no_rows <= p->params->max_index) {
        print_result(p);
    }

    boundary = (p->params->max_index + 1 < p->no_rows) ?
                                        p->params->max_index + 1 : p->no_rows;

    if (p->params->checkpoint_flag || p->params->no_skips > 0) {
        init_checkpoint(boundary, p->params->no_skips, p->params->skips);
    }

    n = -1;
    for (j = p->marker; j < boundary; j++) {
        for (i = 0; i < j; ++i) {
            ++n;
            if (p->params->checkpoint_flag || p->params->no_skips > 0) {
                if (checkpoint_done(n)) continue;
            }
            q = copy_enum(p);
            // after copy_enum() q->tests and p->tests point
            // to the same 'tests set' data structure
            q->marker = j;
            add_subgen(q, i, j);
            if (p->params->normal_flag) {
                q->new_rcpr_table =
                    make_rcpr_table(q->no_columns, q->no_subgens, q->subgens);
            }
            if (coincidence(q, q->index[i], q->index[j]) == COINCIDENCE) {
                q->branch_n = n;
                printf("%d %d %d\n", n, j, i);
            } else {
                // why to call only here ???
                reset_tests(q->tests);
            }
            if (p->params->checkpoint_flag) {
                set_checkpoint(n);
            }
            free_enum(q);
        }
    }
}


int connect_to_server(char *server, int port)
{
int s;
struct hostent *h;
struct sockaddr_in a;

    // prepare address
    bzero(&a, sizeof(a));
    // port number
    a.sin_port = htons((u_short) port);
    // map host name to IP address
    if ((h = gethostbyname(server)) == NULL) {
        Error(ABORT, errno, "Can't get host entries for %s.", server);
    }
    bcopy(h->h_addr, &a.sin_addr, h->h_length);
    a.sin_family = h->h_addrtype;

    // allocate socket
    if ((s = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        Error(ABORT, errno, "Can't create socket.");
    }

    // connect to server
    while (connect(s, (struct sockaddr *) &a, sizeof(a)) < 0) {
        if (errno != ETIMEDOUT) {
            Error(ABORT, errno, "Can't connect to %s.", server);
        }
    }

    return s;
}


static Message get_job(char *server, int port, int *np, int *jp, int *ip)
{
int s;
char buf[100];
int len;
static pid_t pid = (pid_t) 0;
static char pid_string[50];

    if (!pid) {
        pid = getpid();
        sprintf(pid_string, "%d", (int) pid);
    }

    s = connect_to_server(server, port);

    len = sprintf(buf, "get %s\n", pid_string);
    if (write(s, buf,len) < len) {
        Error(ABORT, 0, "Can't send 'get' message.");
    }

    bzero(buf, sizeof(buf));
    if ((len = read(s, buf, sizeof(buf))) < 0) {
        Error(ABORT, errno, "Reading from socket failed.");
    }
    close(s);

    // messages are sent with a trailing '\n'
    if (buf[len-1] == '\n') {
        // we have got a complete message
        if (buf[0] == '-') {
            // no more jobs
            return NO_MORE_JOB;
        } else {
            sscanf(buf, "%d %d %d", np, jp, ip);
            return JOB;
        }
    } else {
        fprintf(stderr, "PID <%s> received incomplete job: %s\n",
                pid_string, buf);
        return INCOMPLETE;
    }

    return 0;
}


static void job_done(char *server, int port, int n, int j, int i)
{
int s;
char buf[100];
int len;

    s = connect_to_server(server, port);
    len = sprintf(buf, "done %d %d %d\n", n, j, i);
    if (write(s, buf, len) < len) {
        Error(ABORT, 0, "Can't send 'done' message: %s", buf);
    }
    close(s);
}


void mixed_descend_slave(struct Enum *p)
{
int i, j, branch;
struct Enum *q;
int cosetenum_result;
Message msg;
char *server = p->params->server;
int port = p->params->port;
int k, n;
int *h;
int no_h;

    h = p->params->subgroup;
    no_h = p->params->no_subgroup;

    if((cosetenum_result = mixed_initial_cosetenum(p)) == ABANDON) {
        return;
    }

    if (!(p->params->normal_flag)) {
        if (conjugate(p)) {
            return;
        }
    }

    if (cosetenum_result == CLOSED && p->no_rows <= p->params->max_index) {
        print_result(p);
    }

    while ((msg = get_job(server, port, &branch, &j, &i)) != NO_MORE_JOB) {
        if (msg == INCOMPLETE) {
            continue;
        }
        // iterate over the subgroup elements
        for (k = n = 0; k < no_h; k++, n += h[n] + 1) {
            q = copy_enum(p);
            // after copy_enum() q->tests and p->tests point
            // to the same 'tests set' data structure
            q->marker = j;
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
                q->branch_n = branch;
                mixed_descend(q);
            } else {
                reset_tests(q->tests);
            }
            free_enum(q);
        }
        job_done(server, port, branch, j, i);
    }

}


void descend_slave(struct Enum *p)
{
int i, j, branch;
struct Enum *q;
int cosetenum_result;
Message msg;
char *server = p->params->server;
int port = p->params->port;

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
    }

    while ((msg = get_job(server, port, &branch, &j, &i)) != NO_MORE_JOB) {
        if (msg == INCOMPLETE) {
            continue;
        }
        q = copy_enum(p);
        // after copy_enum() q->tests and p->tests point
        // to the same 'tests set' data structure
        q->marker = j;
        add_subgen(q, i, j);
        if (p->params->normal_flag) {
            q->new_rcpr_table =
                make_rcpr_table(q->no_columns, q->no_subgens, q->subgens);
        }
        if (coincidence(q, q->index[i], q->index[j]) == COINCIDENCE) {
            q->branch_n = branch;
            descend(q);
        } else {
            reset_tests(q->tests);
        }
        job_done(server, port, branch, j, i);
        free_enum(q);
    }

}
