
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>

//
//  This is the procedure to be called for printing any error
//  messages to stderr. It can optionally abort the program.
//  It expects a fixed argument list and optionally excepts
//  a variable part made up of a valid argument list for printf().
//
// Fixed input parameters:
//      abort_flag      indicates if the program should be terminated
//      file            pointer to a string containig the name of the source
//                      file where Error() was called from.
//      number          line number position in the above file
//      The first three arguments should be supplied by using one of the
//      following preprocessor macros defined in error.h:
//                      ABORT   --> terminate program
//                      WARNING --> return from Error()
//      err_code        either a valid error code or zero
//
// Variable arg list:
//      Either NULL or an argument list for printf().
//

void Error(int abort_flag, char *file, int line, int err_code, ...)
{
va_list args;
char *fmt;

    fprintf (stderr, "\n%s[%d]: ", file, line);
    if (err_code) {
        fprintf (stderr, "%s (%d)\n", strerror(err_code), err_code);
    }

    va_start(args, err_code);

    fmt = va_arg(args, char *);
    if (fmt) {
        vfprintf(stderr, fmt, args);
    }
    fprintf(stderr, "\n");

    if (abort_flag) {
        abort();
    }

    va_end(args);
}
