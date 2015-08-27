#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// These are problem definitions.

#define DIM 2
#define DG_ORDER 0
#define INIT_REF_NUM 5
#define COMPONENT_COUNT 1

#define T_FINAL 1.0
#define DELTA_T 0.1

const bool PRINT_ALGEBRA = true;

#define FLUX (1 - x) * x * nx + (0.5 - x) * ny

#endif