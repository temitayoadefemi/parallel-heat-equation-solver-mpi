#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "misc.h"
#include <stdbool.h>
#include <math.h> // Ensure math.h is included for math functions


slc_str get_serial_dims(int landscape) {
    slc_str slice;

    slice.actual.width = landscape;
    slice.actual.height = landscape;

    slice.rem.width = 0;
    slice.rem.height = 0;

    slice.padded.width  = slice.actual.width  + slice.rem.width;
    slice.padded.height = slice.actual.height + slice.rem.height;

    slice.halo.width = slice.actual.width + 2;
    slice.halo.height = slice.actual.height + 2;

    return slice;
}

slc_str get_parallel_dims(cart_str cart, int landscape) {
    slc_str slice;

    slice.actual.width  = (int) floor((double)landscape  / (double)cart.dims[0]);
    slice.actual.height = (int) floor((double)landscape / (double)cart.dims[1]);

    slice.rem.width  = landscape - slice.actual.width  * cart.dims[0];
    slice.rem.height = landscape - slice.actual.height * cart.dims[1];

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

void initialize_cell_buffers(double **values, double **levels, dim_str dimensions, cart_str cart, int landscape, int rho, int seed) {
    rinit(seed);
    int width =  cart.dims[0] * dimensions.width;
    int height = cart.dims[1] * dimensions.height;
    
    printf("Initializing buffer with dimensions: %d x %d\n", width, height);
    printf("Cart dimensions: %d x %d\n", cart.dims[0], cart.dims[1]);
    printf("Input dimensions: %d x %d\n", dimensions.width, dimensions.height);

    for (int i = 0; i < landscape; i++) {
        for (int j = 0; j < landscape; j++) {
            double r = uni();
            if (r < rho) {
                values[i][j] = 1.5;
                levels[i][j] = 3;
            } else {
                values[i][j] = 8.5;
                levels[i][j] = 10;
            }
            
            if (i == 0 && j < 10) {
                printf("Sample value at [%d][%d]: %f\n", i, j, values[i][j]);
            }
        }
    }
    
    printf("Buffer initialization complete\n");
}

void update_neighbors(int i, int j, double **levels, slc_str slice) {
    int neighbors[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int k = 0; k < 4; k++) {
        int ni = i + neighbors[k][0];
        int nj = j + neighbors[k][1];
        if (ni >= 0 && ni < slice.actual.width && nj >= 0 && nj < slice.actual.height) {
           levels[ni][nj] = fmax(levels[ni][nj], levels[i][j] - 1);
        }
    }
}

double laplacian(double **values, int i, int j, double dx, double dy) {
    return (values[i+1][j] + values[i-1][j] - 2*values[i][j]) / (dx*dx) +
           (values[i][j+1] + values[i][j-1] - 2*values[i][j]) / (dy*dy);
}

void initialize_dt(double *dt, double dx, double dy) {
    *dt = 0.25 * (dx * dx + dy * dy); 
}

void solve_heat_equation(double **values, double dx, double dy, double dt, slc_str slice)  {
    for (int i = 1; i < slice.actual.width - 1; i++) {
        for (int j = 1; j < slice.actual.height - 1; j++) {
            double laplacian_value = laplacian(values, i, j, dx, dy);
            values[i][j] += dt * laplacian_value;
        }
    }
}

void refine_mesh(double **levels, double **values, double dx, double dy, slc_str slice, int MAX_LEVEL) {
    for (int i = 1; i < slice.actual.width - 1; i++) {
        for (int j = 1; j < slice.actual.height - 1; j++) {
            if (levels[i][j] < MAX_LEVEL && fabs(laplacian(values, i, j, dx, dy)) > 0.1) {
                levels[i][j]++;
                update_neighbors(i, j, levels, slice);
            }
        }
    }
}

void distribute_cells(double **local, double **global, cart_str cart, slc_str slice) {
    for (int i = 0; i < slice.actual.width; i++) {
        for (int j = 0; j < slice.actual.height; j++) {
            int global_row = cart.coords[0] * slice.actual.width + i;
            int global_col = cart.coords[1] * slice.actual.height + j;
            local[i][j] = global[global_row][global_col];
        }
    }
}

void distribute_to_temp(double **temp, double **mini, cart_str cart, slc_str slice) {
    for (int i = 0; i < slice.actual.width; i++) {
        for (int j = 0; j < slice.actual.height; j++) {
            int global_row = cart.coords[0] * slice.actual.width + i;
            int global_col = cart.coords[1] * slice.actual.height + j;
            temp[global_row][global_col] = mini[i][j];
        }
    }
}

void zerotmpcell(double **values, double **levels, dim_str dimensions) {
    for (int i = 0; i < dimensions.width; i++) {
        for (int j = 0; j < dimensions.height; j++) {
            values[i][j] = 0;
            levels[i][j] = 0;
            
        }
    }
}

void copy_buff_to_mini(double **mini, double **local, slc_str slice) {
    for (int i = 1; i <= slice.actual.width; i++) {
        for (int j = 1; j <= slice.actual.height; j++) {
            mini[i - 1][j - 1] = local[i][j];
        }
    }
}

void initialise_edges(double **values, double **levels, dim_str slice){
    for (int i = 0; i <= slice.width+1; i++) {
        for (int j = 0; j <= slice.height+1; j++) {
            if (i == 0 || i == slice.width+1 || j == 0 || j == slice.height+1) {
                values[i][j] = 0;
                levels[i][j] = 0;
            }
        }
    }
}

void copy_buff_to_local(double **mini_values, double **mini_levels, double **local_values, double **local_levels, slc_str slice) {
    for (int i = 1; i <= slice.actual.width; i++) {
        for (int j = 1; j <= slice.actual.height; j++) {
            mini_values[i][j] = local_values[i - 1][j - 1];
            mini_levels[i][j] = local_levels[i - 1][j - 1];
        }
    }
}