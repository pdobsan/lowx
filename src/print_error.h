#ifndef __error_h
#define __error_h

#include <errno.h>

#define ABORT   1, __FILE__, __LINE__
#define WARNING 0, __FILE__, __LINE__

void Error(int abort_flag, char *file, int line, int err_code, ...);

#ifdef DEBUG
# define debug_printf(arg) printf arg
#else
# define debug_printf(arg)
#endif

#endif
