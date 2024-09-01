#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "arralloc.h"
#include "misc.h"
#include <stdbool.h>


void initialize_grid(cell **global_grid, master *master) {
    printf("Entering initialize_grid\n");
    
    // Check allocation sizes
    printf("Grid dimensions: %d x %d\n", master->dimensions.landscape, master->dimensions.landscape);
    
    // Add checks before accessing arrays
    if (grid != NULL) {
        printf("Grid allocated successfully\n");
    } else {
        printf("Error: Grid not allocated\n");
        return;
    }
    
    // Add bounds checking in loops
    for (int i = 0; i < master->dimensions.landscape; i++) {
        for (int j = 0; j < master->dimensions.landscape; j++) {
            double r = uni();
            if (r < 0.51) {
                grid[i][j].value = 1.5;
                grid[i][j].level = 3;
            } else {
                grid[i][j].value = 8.5;
                grid[i][j].level = 10;
            }
        }
    }
    
    printf("Grid initialization complete\n");
}

void initialize_dt(double *dt, double dx, double dy) {
    *dt = 0.25 * (dx * dx + dy * dy); 
}


void solve_heat_equation(cell **local_grid, double dx, double dy, double dt)  {
    for (int i = 1; i < local_lx - 1; i++) {
        for (int j = 1; j < local_ly - 1; j++) {
            double laplacian_value = laplacian(local_grid, i, j, dx, dy);
            local_grid[i][j].value += dt * laplacian_value;
        }
    }
}


void refine_mesh(cell **local_grid, double dx, double dy) {
    for (int i = 1; i < local_lx - 1; i++) {
        for (int j = 1; j < local_ly - 1; j++) {
            if (local_grid[i][j].level < MAX_LEVEL && fabs(laplacian(local_grid, i, j, dx, dy)) > 0.1) {
                local_grid[i][j].level++;
                update_neighbors(i, j, local_grid);
            }
        }
    }
}


void update_neighbors(int i, int j, cell **local_grid) {
    int neighbors[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int k = 0; k < 4; k++) {
        int ni = i + neighbors[k][0];
        int nj = j + neighbors[k][1];
        if (ni >= 0 && ni < LX && nj >= 0 && nj < LY) {
            local_grid[ni][nj].level = fmax(local_grid[ni][nj].level, local_grid[i][j].level - 1);
        }
    }
}

double laplacian(cell **local_grid, int i, int j, double dx, double dy) {
    return (local_grid[i+1][j].value + local_grid[i-1][j].value - 2*local_grid[i][j].value) / (dx*dx) +
           (local_grid[i][j+1].value + local_grid[i][j-1].value - 2*local_grid[i][j].value) / (dy*dy);
}

