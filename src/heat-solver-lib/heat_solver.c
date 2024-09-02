#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "misc.h"
#include <stdbool.h>
#include <math.h>

slc_str get_serial_dims(master_str *master)
{
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

slc_str get_parallel_dims(master_str *master)
{
    slc_str slice;

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

void initialize_cell_buffers(**double values, **double levels, dim_str dimensions ) {
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

void initialize_dt(double *dt, double dx, double dy) {
    *dt = 0.25 * (dx * dx + dy * dy); 
}



void solve_heat_equation(**double values, double dx, double dy, double dt)  {
    for (int i = 1; i < local_lx - 1; i++) {
        for (int j = 1; j < local_ly - 1; j++) {
            double laplacian_value = laplacian(values, i, j, dx, dy);
            values[i][j] += dt * laplacian_value;
        }
    }
}


void refine_mesh(**double levels, **double values, double dx, double dy) {
    for (int i = 1; i < local_lx - 1; i++) {
        for (int j = 1; j < local_ly - 1; j++) {
            if (ca[i][j].level < MAX_LEVEL && fabs(laplacian(values, i, j, dx, dy)) > 0.1) {
                levels[i][j]++;
                update_neighbors(i, j, local_grid);
            }
        }
    }
}


void update_neighbors(int i, int j, **double levels) {
    int neighbors[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int k = 0; k < 4; k++) {
        int ni = i + neighbors[k][0];
        int nj = j + neighbors[k][1];
        if (ni >= 0 && ni < LX && nj >= 0 && nj < LY) {
           levels[ni][nj] = fmax(levels[ni][nj], levels[i][j] - 1);
        }
    }
}


double laplacian(**double values, int i, int j, double dx, double dy) {
    return (values[i+1][j] + values[i-1][j] - 2*values[i][j]) / (dx*dx) +
           (values[i][j+1] + values[i][j-1] - 2*values[i][j]) / (dy*dy);
}

