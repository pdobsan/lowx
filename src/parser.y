/*
    This module contains the parser definition for a simple language
    describing group presentations and subgroups.

    Here is an example for input accepted by the parser:

        #
        # sample input for 'lowx'
        #

        group = <a, b | b * a * b * a^-2, a * b * a * b^-2> ;

        max_index = 50;
*/

%{

#include "lowx.h"

extern struct Parameters params;

int yylex();
void yyerror(char *msg);

%}

%union {
    char sym;
    int num;
}

%token GROUP SUBGENS SUBGROUP MAX_TESTS MAX_COINCS MAX_ROWS MAX_INDEX SKIP
%token NORMAL MIXED_ CHECKPOINT MAX_QUEUE PORT SERVER
%token <sym> SYMBOL
%token <num> NUMBER

%%

problem: parameters ';'
    ;

parameters: expression 
    | parameters ';' expression

expression: group
    | subgens
    | subgroup
    | MAX_INDEX '=' NUMBER
        {
            params.max_index = $3;
        }
    | NORMAL
        {
            params.normal_flag = 1;
        }
    | MIXED_
        {
            params.mixed_flag = 1;
            params.lowindex_flag = 1;
        }
    | CHECKPOINT
        {
            params.checkpoint_flag = 1;
        }
    | SKIP skip_list
        {
            /* printf("\n") */
            ;
        }
    | MAX_COINCS '=' NUMBER
        {
            params.max_no_coincs = $3;
        }
    | MAX_ROWS '=' NUMBER
        {
            params.max_no_rows = $3;
        }
    | MAX_QUEUE '=' NUMBER
        {
            fprintf(stderr, "#-- Warning: max_queue is ignored.\n");
        }
    | PORT '=' NUMBER
        {
#ifdef DISTRIBUTED
            params.port = $3;
            params.print_port = $3 + 1;
#else
        fprintf(stderr, "#-- Warning: port number is ignored.\n");
#endif
        }
    ;

group: GROUP '=' presentation
    ;

presentation: '<' group_generators bar_sign word_list '>'
    ;

group_generators:
    SYMBOL
        {
            params.generators[params.no_generators++] = $1;
        }
    | group_generators ',' SYMBOL
        {
            params.generators[params.no_generators++] = $3;
        }
    ;

bar_sign: '|'
        {
            init_vector(&(params.no_relators), params.relators);
        }
    ;

subgens: SUBGENS '=' subgens_begin word_list '>'
        {
            params.lowindex_flag = 0;
        }
    ;

subgens_begin: '<'
        {
            init_vector(&(params.no_subgens), params.subgens);
        }
    ;

subgroup: SUBGROUP '=' subgroup_begin word_list '>'
    ;

subgroup_begin: '<'
        {
            init_vector(&(params.no_subgroup), params.subgroup);
        }
    ;

word_list:
    word
        {
            next_word();
        }
    | word_list ',' word
        {
            next_word();
        }
    ;

word:
    factor
    | word '*' factor
    ;

factor:
    SYMBOL
        {
            put_code($1, 1);
        }
    | SYMBOL '^' NUMBER
        {
            put_code($1, $3);
        }
    | SYMBOL '^' '-' NUMBER
        {
            put_code($1, -$4);
        }
    | left_paren word ')' '^' NUMBER
        {
            append_frame($5);
            remove_frame();
        }
    | left_paren word ')' '^' '-' NUMBER
        {
            append_frame(-$6);
            remove_frame();
        }
    | commexp
    ;

commexp:
    commutator
        {
            remove_frame();
        }
    | commutator '^' NUMBER
        {
            append_frame($3);
            remove_frame();
        }
    | commutator '^' '-' NUMBER
        {
            append_frame(-$4);
            remove_frame();
        }
    ;

commutator:
    | left_comm ',' word ']'
        {
            commutator();
            remove_frame();
        }
    ;

left_comm: left_bracket word
        {
            new_frame();
        }
    ;

left_bracket: '['
        {
            new_frame();
        }
    ;

left_paren: '('
        {
            new_frame();
        }
    ;

skip_list: skip_entry
    | skip_list skip_entry
    ;

skip_entry:
    NUMBER
        {
            params.skips[params.no_skips][0] = $1 - 1;
            params.skips[params.no_skips][1] = $1;
            ++(params.no_skips);
        }
    | NUMBER '-' NUMBER
        {
            params.skips[params.no_skips][0] = $1 - 1;
            params.skips[params.no_skips][1] = $3;
            ++(params.no_skips);
        }
    ;

%%
