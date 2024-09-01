#ifndef HEATSOLVERLIB_H
#define HEATSOLVERLIB_H

#include "structs.h"  // Including necessary structures like cart_str, master_str, etc.
#include <stdbool.h>


void solve_heat_equation(cell **local_grid, double dx, double dy, double dt);
void refine_mesh(cell **local_grid, double dx, double dy);
void update_neighbors(int i, int j, cell **local_grid);
double laplacian(cell **local_grid, int i, int j, double dx, double dy);
void initialize_grid(cell  **global_grid);
void initialize_dt(double *dt, double dx, double dy);

#endif // HEATSOLVERLIB_H