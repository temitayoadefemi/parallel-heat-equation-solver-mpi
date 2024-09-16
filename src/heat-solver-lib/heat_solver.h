#ifndef HEATSOLVERLIB_H
#define HEATSOLVERLIB_H

#include "structs.h"  // Including necessary structures like cart_str, master_str, etc.
#include <stdbool.h>
#include <math.h> // Ensure math.h is included for math functions

// Function prototypes
slc_str get_serial_dims(int landscape);
slc_str get_parallel_dims(cart_str cart, int landscape);
void initialize_cell_buffers(double **values, double **levels, dim_str dimensions, cart_str cart, int landscape, int rho, int seed);
void update_neighbors(int i, int j, double **levels, slc_str slice);
double laplacian(double **values, int i, int j, double dx, double dy);
void initialize_dt(double *dt, double dx, double dy);
void solve_heat_equation(double **values, double dx, double dy, double dt, slc_str slice);
void refine_mesh(double **levels, double **values, double dx, double dy, slc_str slice, int MAX_LEVEL);
void distribute_cells(double **local, double **global, cart_str cart, slc_str slice);
void zerotmpcell(double **values, double **levels, dim_str dimensions);
void copy_buff_to_mini(double **mini, double **local, slc_str slice);
void initialise_edges(double **values, double **levels, dim_str slice);
void copy_buff_to_local(double **mini_values, double **mini_levels, double **local_values, double **local_levels, slc_str slice);

#endif // HEAT_SOLVER_H