#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"
#include <mpi.h>


#define LANDSCAPE 1152 
#define MAX_LEVEL 3  // Maximum refinement level

// Structure to represent a grid cell
typedef struct {
    double value;
    int level;
} Cell;



// Updated function prototypes
void initialize_grid(int LX, int LY);
void initialize_dt();
void solve_heat_equation(int LX, int LY);
void refine_mesh(int LX, int LY);
// ... other function prototypes ...

int main(int argc, char **argv) {
    MPI_Comm comm;
    MPI_Status status;
    int size, rank;
    int dims[2], periods[2], coords[2];
    int LX, LY;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create a 2D Cartesian topology
    dims[0] = dims[1] = 0;  // Let MPI decide the dimensions
    MPI_Dims_create(size, 2, dims);
    
    LX = LANDSCAPE / dims[0];
    LY = LANDSCAPE / dims[1];


    Cell local_grid[LX + 2][LY + 2];
    Cell global_grid[LANDSCAPE][LANDSCAPE];

    Cell tempcell[LANDSCAPE][LANDSCAPE];
    Cell smallcell[LX][LY];

    double dx = 1.0 / LX;
    double dy = 1.0 / LY;
    double dt;  // Declare dt without initialization

    periods[0] = periods[1] = 0;  // Non-periodic
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Cart_coords(comm, rank, 2, coords);

    if (rank == 0) {
        initialize_dt();
        initialize_grid(LX, LY);
    }

    MPI_Bcast(&global_grid[0][0], LX*LY, MPI_INT, 0, cart.comm2d);

    for (i = 0; i <= LX; i++) {
        for (j = 0; j <= LY; j++) {
            smallcell[i][j] = allcell;
        }
    }

    for (i = 1; i <= LX; i++) {
        for (j = 1; j <= LY; j++) {
            cell[i][j] = smallcell;
        }
    }

    MPI_Type_vector(LX, 1, LY + 2, MPI_INT, column_type);
    MPI_Type_commit(column_type); // Commit the type to use it for MPI operations.

    // Create a contiguous type for transferring rows.
    MPI_Type_contiguous(LY, MPI_INT, row_type);
    MPI_Type_commit(row_type); // Commit the type to use it for MPI operations.
    
    for (int step = 0; step < 1000; step++) {
        MPI_Isend(&cell_grid[LX][1], 1, row_type, cart.down.val, 1, cart.comm2d, &reqs[0]); // Send bottom row.
        MPI_Isend(&cell_grid[1][1], 1, row_type, cart.up.val, 2, cart.comm2d, &reqs[2]); // Send top row.
        MPI_Isend(&cell_grid[1][LY], 1, column_type, cart.right.val, 3, cart.comm2d, &reqs[4]); // Send right column.
        MPI_Isend(&cell_grid[1][1], 1, column_type, cart.left.val, 4, cart.comm2d, &reqs[6]); // Send left column.

        MPI_Irecv(&cell_grid[0][1], 1, row_type, cart.up.val, 1, cart.comm2d, &reqs[1]); // Receive top row.
        MPI_Irecv(&cell_grid[LX+1][1], 1, row_type, cart.down.val, 2, cart.comm2d, &reqs[3]); // Receive bottom row.
        MPI_Irecv(&cell_grid[1][0], 1, column_type, cart.left.val, 3, cart.comm2d, &reqs[5]); // Receive left column.
        MPI_Irecv(&cell_grid[1][LY+1], 1, column_type, cart.right.val, 4, cart.comm2d, &reqs[7]); // Receive right column.


        solve_heat_equation(LX, LY);
        if (step % 10 == 0) {
            refine_mesh(LX, LY);
        }
    }

    save_to_pbm("heat_equation", grid);
    
    MPI_Finalize();
    return 0;
}


void initialize_grid(int L) {
    double r;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            r = uni();
            if (r < 0.51) {
                global_grid[i][j].value = 3.5;
                global_grid[i][j].level = 1;
            } else {
                global_grid[i][j].value = 7.0;
                global_grid[i][j].level = 2;
            }
        }
    }
}

void initialize_dt() {
    dt = 0.25 * (dx * dx + dy * dy); 
}


void solve_heat_equation(int LX, int LY) {
    for (int i = 1; i < LX - 1; i++) {
        for (int j = 1; j < LY - 1; j++) {
            double laplacian_value = laplacian(grid, i, j, dx, dy);
            local_grid[i][j].value += dt * laplacian_value;
        }
    }
}


void refine_mesh(int LX, int LY) {
    for (int i = 1; i < LX - 1; i++) {
        for (int j = 1; j < LY - 1; j++) {
            if (local_grid[i][j].level < MAX_LEVEL && fabs(laplacian(local_grid, i, j, dx, dy)) > 0.1) {
                local_grid[i][j].level++;
                update_neighbors(i, j);
            }
        }
    }
}


void update_neighbors(int i, int j) {
    int neighbors[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int k = 0; k < 4; k++) {
        int ni = i + neighbors[k][0];
        int nj = j + neighbors[k][1];
        if (ni >= 0 && ni < LX && nj >= 0 && nj < LY) {
            local_grid[ni][nj].level = fmax(local_grid[ni][nj].level, local_grid[i][j].level - 1);
        }
    }
}


double laplacian(Cell grid[NX][NY], int i, int j, double dx, double dy) {
    return (local_grid[i+1][j].value + local_grid[i-1][j].value - 2*local_grid[i][j].value) / (dx*dx) +
           (local_grid[i][j+1].value + local_grid[i][j-1].value - 2*local_grid[i][j].value) / (dy*dy);
}



void save_to_pbm(const char* filename, Cell (*grid)[NY]) {
    FILE *fp;
    char full_filename[100];
    snprintf(full_filename, sizeof(full_filename), "%s.pbm", filename);
    
    fp = fopen(full_filename, "w");
    if (fp == NULL) {
        printf("Error opening file %s\n", full_filename);
        return;
    }

    // Write PBM header
    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", LX, LY);
    fprintf(fp, "255\n");  // Maximum gray value

    // Write grid values
    for (int j = 0; j < LY; j++) {
        for (int i = 0; i < LX; i++) {
            // Scale the value to 0-255 range
            int gray_value = (int)((grid[i][j].value - 3.5) * 255 / 3.5);
            // Clamp the value to 0-255
            gray_value = (gray_value < 0) ? 0 : (gray_value > 255) ? 255 : gray_value;
            fprintf(fp, "%d ", gray_value);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Saved grid state to %s\n", full_filename);
}