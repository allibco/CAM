#include "superlu_ddefs.h"

int get_refine_steps(SuperLUStat_t *stat) {
    return stat->RefineSteps;
}
