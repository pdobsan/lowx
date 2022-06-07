/*
    lexer module for the 'cosetenum' program
    here is a simple input file:

        #
        # sample input
        #

        group = <a, b | b * a * b * a^-2, a * b * a * b^-2>;
        subgens = <a^2>;
*/

%{

#include "parser.tab.h"
#include <stdlib.h>

int lineno = 1;

%}

%%

[\t ]+                  /* ignore white space */ ;
\n              ++lineno;
#.*                     /* ignore comment */ ;

group                   { return GROUP; }
subgroup                { return SUBGROUP; }
subgens                 { return SUBGENS; }
normal_subgroups        { return NORMAL; }
mixed                   { return MIXED_; }
max_no_coincs           { return MAX_COINCS; }
max_no_rows             { return MAX_ROWS; }
max_index               { return MAX_INDEX; }
checkpoint              { return CHECKPOINT; }
max_queue               { return MAX_QUEUE; }
skip                    { return SKIP; }
server                  { return SERVER; }
port                    { return PORT; }

[a-zA-Z]                { yylval.sym = (int) yytext[0]; return SYMBOL; }
[0-9]+                  { yylval.num = atoi(yytext); return NUMBER; }

.                       { return yytext[0]; }

%%

void yyerror(char *msg)
{
    fprintf(stderr, "%d: %s at '%s'\n", lineno, msg, yytext);
}
