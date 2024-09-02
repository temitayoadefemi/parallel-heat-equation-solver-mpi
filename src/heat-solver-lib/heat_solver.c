#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "misc.h"
#include <stdbool.h>
#include <math.h> // Ensure math.h is included for math functions

// Function to get serial dimensions
slc_str get_serial_dims(master_str *master) {
    slc_str slice;

    slice.actual.width = master->params.landscape;
    slice.actual.height = master->params.landscape;

    slice.rem.width = 0;
    slice.rem.height = 0;

    slice.padded.width  = slice.actual.width  + slice.rem.width;
    slice.padded.height = slice.actual.height + slice.rem.height;

    slice.halo.width = slice.actual.width + 2;
    slice.halo.height = slice.actual.height + 2;

    return slice;
}

// Function to get parallel dimensions
slc_str get_parallel_dims(master_str *master) {
    slc_str slice;

    // Ensure cart is defined or passed as a parameter
    slice.actual.width  = (int) floor((double)master->params.landscape  / (double)cart.dims[0]);
    slice.actual.height = (int) floor((double)master->params.landscape / (double)cart.dims[1]);

    slice.rem.width  = master->params.landscape - slice.actual.width  * cart.dims[0];
    slice.rem.height = master->params.landscape - slice.actual.height * cart.dims[1];

    slice.padded.width  = slice.actual.width  + slice.rem.width;
    slice.padded.height = slice.actual.height + slice.rem.height;

    if(cart.coords[0] == cart.dims[0] - 1)
        slice.actual.width = slice.padded.width;

    if(cart.coords[1] == cart.dims[1] - 1)
        slice.actual.height = slice.padded.height;

    slice.halo.height = slice.actual.height + 2;
    slice.halo.width  = slice.actual.width  + 2;

    return slice;
}

// Function to initialize cell buffers
void initialize_cell_buffers(double **values, double **levels, dim_str dimensions) {
    // Add bounds checking in loops
    for (int i = 0; i < dimensions.width; i++) {
        for (int j = 0; j < dimensions.height; j++) {
            double r = uni();
            if (r < 0.51) {
                values[i][j] = 1.5;
                levels[i][j] = 3;
            } else {
                values[i][j] = 8.5;
                levels[i][j] = 10;
            }
        }
    }
}

// Function to initialize dt
void initialize_dt(double *dt, double dx, double dy) {
    *dt = 0.25 * (dx * dx + dy * dy); 
}

// Function to solve the heat equation
void solve_heat_equation(double **values, double dx, double dy, double dt, slc_str slice)  {
    // Ensure local_lx and local_ly are defined or passed as parameters
    for (int i = 1; i < slice.actual.width - 1; i++) {
        for (int j = 1; j < slice.actual.height - 1; j++) {
            double laplacian_value = laplacian(values, i, j, dx, dy);
            // Correct array access
            values[i][j] += dt * laplacian_value;
        }
    }
}

// Function to refine the mesh
void refine_mesh(double **levels, double **values, double dx, double dy, slc_str slice) {
    // Ensure local_lx, local_ly, ca, and MAX_LEVEL are defined or passed as parameters
    for (int i = 1; i < slice.actual.width - 1; i++) {
        for (int j = 1; j < slice.actual.height - 1; j++) {
            if (ca[i][j].level < MAX_LEVEL && fabs(laplacian(values, i, j, dx, dy)) > 0.1) {
                levels[i][j]++;
                update_neighbors(i, j, levels);
            }
        }
    }
}

// Function to distribute cells
void distribute_cells(double **local, double **global, cart_str cart, slc_str slice) {
    for (int i = 0; i < slice.actual.width; i++) {
        for (int j = 0; j < slice.actual.height; j++) {
            // Ensure master and its dimensions are defined or passed as parameters
            int global_row = cart.coords[0] * master->dimensions.rows + i;
            int global_col = cart.coords[1] * master->dimensions.cols + j;
            local[i][j] = global[global_row][global_col];
        }
    }
}

// Function to zero out temporary cell buffer
void zerotmpcell(double **temp, dim_str dimensions) {
    for (int i = 0; i < dimensions.width; i++) {
        for (int j = 0; j < dimensions.height; j++) {
            temp[i][j] = 0;
        }
    }
}

// Function to copy buffer to mini
void copy_buff_to_mini(double **mini, double **local, slc_str slice) {
    for (int i = 1; i <= slice.actual.width; i++) { // Removed extra comma
        for (int j = 1; j <= slice.actual.height; j++) {
            mini[i - 1][j - 1] = local[i][j];
        }
    }
}

// Function to update neighbors
void update_neighbors(int i, int j, double **levels) {
    int neighbors[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int k = 0; k < 4; k++) {
        int ni = i + neighbors[k][0];
        int nj = j + neighbors[k][1];
        // Ensure LX and LY are defined or passed as parameters
        if (ni >= 0 && ni < LX && nj >= 0 && nj < LY) {
           levels[ni][nj] = fmax(levels[ni][nj], levels[i][j] - 1);
        }
    }
}

// Function to calculate the laplacian
double laplacian(double **values, int i, int j, double dx, double dy) {
    return (values[i+1][j] + values[i-1][j] - 2*values[i][j]) / (dx*dx) +
           (values[i][j+1] + values[i][j-1] - 2*values[i][j]) / (dy*dy);
}

