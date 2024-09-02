#ifndef HEATSOLVERLIB_H
#define HEATSOLVERLIB_H

#include "structs.h"  // Including necessary structures like cart_str, master_str, etc.
#include <stdbool.h>

slc_str get_serial_dims(master_str *master);

slc_str get_parallel_dims(master_str *master);

void initialize_cell_buffers(**double values, **double levels, dim_str dimensions);

void initialize_dt(double *dt, double dx, double dy);

void solve_heat_equation(**double values, double dx, double dy, double dt);

void refine_mesh(**double levels, **double values, double dx, double dy);

void update_neighbors(int i, int j, **double levels);

double laplacian(**double values, int i, int j, double dx, double dy);


#endif // HEATSOLVERLIB_H