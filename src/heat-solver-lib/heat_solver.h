#ifndef HEATSOLVERLIB_H
#define HEATSOLVERLIB_H

#include "structs.h"  // Including necessary structures like cart_str, master_str, etc.
#include <stdbool.h>
#include <math.h> // Ensure math.h is included for math functions

// Function to get serial dimensions
slc_str get_serial_dims(master_str *master);

// Function to get parallel dimensions
slc_str get_parallel_dims(master_str *master);

// Function to initialize cell buffers
void initialize_cell_buffers(double **values, double **levels, dim_str dimensions);

// Function to initialize dt
void initialize_dt(double *dt, double dx, double dy);

// Function to solve the heat equation
void solve_heat_equation(double **values, double dx, double dy, double dt);

// Function to refine the mesh
void refine_mesh(double **levels, double **values, double dx, double dy, slc_str slice, int MAX_LEVEL);

// Function to distribute cells
void distribute_cells(double **local, double **global, cart_str cart, slc_str slice, master_str *master);

// Function to zero out temporary cell buffer
void zerotmpcell(double **temp, dim_str dimensions);

// Function to copy buffer to mini
void copy_buff_to_mini(double **mini, double **local, slc_str slice);

// Function to update neighbors
void update_neighbors(int i, int j, double **levels, int LX, int LY);

// Function to calculate the laplacian
double laplacian(double **values, int i, int j, double dx, double dy);

#endif // HEATSOLVERLIB_H