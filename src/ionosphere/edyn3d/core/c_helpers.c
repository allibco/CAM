#include "superlu_ddefs.h"

typedef long long int fptr;  /* match superlu_c2f_wrap.c */

int get_refine_steps(SuperLUStat_t *stat) {
    return stat->RefineSteps;
}
void f_set_superlu_diagpivotthresh_(fptr *opt, double *thresh)
{
    superlu_dist_options_t *options = (superlu_dist_options_t *)(*opt);
    options->DiagPivotThresh = *thresh;
}
