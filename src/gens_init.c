
#define GENS_INIT__C
#include "lowx.h"

void gens_init(struct Parameters *params)
{
int i,j;
int k,m,n,g;
int involution[2 * MAX_NO_GENS];

    params->no_columns = 2 * params->no_generators;

    for (i = 0; i < params->no_columns; i++) {
        involution[i] = 0;
    }

    // find involutions (g^2 == 1)
    for (i = j = 0; i < params->no_relators; i++, j += params->relators[j] + 1) {
        if (params->relators[j] == 2) {
            if (params->relators[j+1] == params->relators[j+2]) {
                involution[params->relators[j+1]] = 1;
                --(params->no_columns);
            }
        } else if (params->relators[j] == -2) {
            // for the case if an involution is given in the form of g^-2
            if (params->relators[j+1] == params->relators[j+2]) {
                involution[params->relators[j+1] - params->no_generators] = 1;
                --(params->no_columns);
            }
        }
    }

    // initialize inverse table
    for (i = 0, j = params->no_generators; i < params->no_generators; i++) {
        if (involution[i]) {
            params->inverse[i] = i;
        } else {
            params->inverse[i] = j;
            params->inverse[j] = i;
            ++j;
        }
    }

    // modify params->relators[] and params->subgens[] because the codes of inverses
    // of generators have been changed
    //
    // if generator g is an involution then g and g^-1 are represented
    // by the same code
    //
    // if g is not an involution then its code should be decreased
    // by the number of involutions which precede g in generators[]
    //
    for (i = j = 0; i < params->no_relators; i++, j += params->relators[j] + 1) {
        for (k = j + 1; k <= j + params->relators[j]; k++) {
            if (params->relators[k] >= params->no_generators) {
                g = params->relators[k] - params->no_generators;
                if (involution[g]) {
                    params->relators[k] = g;
                } else {
                    for (m = n = 0; m < g; m++) {
                        if (involution[m]) {
                            ++n;
                        }
                    }
                    params->relators[k] -= n;
                }
            }
        }
    }

    for (i = j = 0; i < params->no_subgens; i++, j += params->subgens[j] + 1) {
        for (k = j + 1; k <= j + params->subgens[j]; k++) {
            if (params->subgens[k] >= params->no_generators) {
                g = params->subgens[k] - params->no_generators;
                if (involution[g]) {
                    params->subgens[k] = g;
                } else {
                    for (m = n = 0; m < g; m++) {
                        if (involution[m]) {
                            ++n;
                        }
                    }
                    params->subgens[k] -= n;
                }
            }
        }
    }
}
