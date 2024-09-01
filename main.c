#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"
#include <mpi.h>
#include <time.h>

#define LANDSCAPE 1152
#define MAX_LEVEL 3  // Maximum refinement level
#define NPROC 4
#define LX 576  // Replace with your desired value
#define LY 576

// Structure to represent a grid cell
typedef struct {
    double value;
    int level;
} Cell;



// Updated function prototypes
// Update the function prototype (near the top of the file)
// Update these function prototypes
void solve_heat_equation(int local_lx, int local_ly, Cell local_grid[LX+2][LY+2], double dx, double dy, double dt);
void refine_mesh(int local_lx, int local_ly, Cell local_grid[LX+2][LY+2], double dx, double dy);
void update_neighbors(int i, int j, Cell local_grid[LX+2][LY+2]);
double laplacian(Cell local_grid[LX+2][LY+2], int i, int j, double dx, double dy);
void initialize_grid(int l, Cell grid[LANDSCAPE][LANDSCAPE]);
void initialize_dt(double *dt, double dx, double dy);
void save_to_pbm(const char* filename, Cell grid[LANDSCAPE][LANDSCAPE]);
void cell_sum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);
void print_grid(int L, Cell grid[L][L]);
// ... other function prototypes ...


int main(int argc, char *argv[]) {

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    int size, rank;
    int local_lx, local_ly;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int i, j;
    int dims[2] = {0, 0};  
    int periods[2] = {0, 0};
    int coords[2];
    int reorder = 0; 

    double r; 

    MPI_Dims_create(size, 2, dims);

    local_lx = LANDSCAPE / dims[0];
    local_ly = LANDSCAPE / dims[1];

    MPI_Request reqs[8];
    int neighbors[4];
    
    Cell local_grid[LX + 2][LY + 2];
    Cell global_grid[LANDSCAPE][LANDSCAPE];
    Cell tempgrid[LANDSCAPE][LANDSCAPE];
    Cell smallgrid[LX][LY];

    double dx = 1.0 / LX;
    double dy = 1.0 / LY;
    double dt;
    initialize_dt(&dt, dx, dy);

    periods[0] = periods[1] = 0;  // Non-periodic
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Cart_coords(comm, rank, 2, coords);

    MPI_Cart_shift(comm, 0, 1, &neighbors[0], &neighbors[1]);
    MPI_Cart_shift(comm, 1, 1, &neighbors[2], &neighbors[3]);

    rinit(1);


    if (rank == 0) {
        initialize_grid(LANDSCAPE, global_grid);
    }

    // Create MPI datatype for Cell
    MPI_Datatype MPI_CELL;
    int blocklengths[2] = {1, 1};
    MPI_Aint offsets[2];
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
    offsets[0] = offsetof(Cell, value);
    offsets[1] = offsetof(Cell, level);
    MPI_Type_create_struct(2, blocklengths, offsets, types, &MPI_CELL);
    MPI_Type_commit(&MPI_CELL);

    MPI_Op MPI_CELL_SUM;
    MPI_Op_create((MPI_User_function *)cell_sum, 1, &MPI_CELL_SUM);


    MPI_Bcast(&global_grid[0][0], LANDSCAPE*LANDSCAPE, MPI_CELL, 0, comm);

    for (int i = 0; i < LX; i++) {
        for (int j = 0; j < LY; j++) {
            int global_row = coords[0] * LX + i;
            int global_col = coords[1] * LY + j;
            smallgrid[i][j] = global_grid[global_row][global_col];
        }
    }

    for (int i = 1; i <= LX; i++) {
        for (int j = 1; j <= LY; j++) {
            local_grid[i][j] = smallgrid[i - 1][j - 1];
        }
    }

    for (i=0; i <= LX+1; i++) {
        local_grid[i][0].value    = 0;
        local_grid[i][LY+1].value = 0;
    }

    for (j=0; j <= LY+1; j++) {
        local_grid[0][j].value    = 0;
        local_grid[LX+1][j].value = 0;
    }


    MPI_Datatype column_type, row_type;
    MPI_Type_vector(LX, 1, LY + 2, MPI_CELL, &column_type);
    MPI_Type_commit(&column_type);

    MPI_Type_contiguous(LY, MPI_CELL, &row_type);
    MPI_Type_commit(&row_type);
    
    int bsize = (LX + LY) * sizeof(Cell) + MPI_BSEND_OVERHEAD;
    void* buffer = malloc(bsize);
    if (buffer == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Buffer_attach(buffer, bsize);

    for (int step = 0; step < 1000; step++) {
        MPI_Isend(&local_grid[LX][1], 1, row_type, neighbors[1], 0, comm, &reqs[0]);
        MPI_Irecv(&local_grid[0][1], 1, row_type, neighbors[0], 0, comm, &reqs[1]);
        
        MPI_Isend(&local_grid[1][1], 1, row_type, neighbors[0], 1, comm, &reqs[2]);
        MPI_Irecv(&local_grid[LX+1][1], 1, row_type, neighbors[1], 1, comm, &reqs[3]);
        
        MPI_Isend(&local_grid[1][LY], 1, column_type, neighbors[3], 2, comm, &reqs[4]);
        MPI_Irecv(&local_grid[1][0], 1, column_type, neighbors[2], 2, comm, &reqs[5]);
        
        MPI_Isend(&local_grid[1][1], 1, column_type, neighbors[2], 3, comm, &reqs[6]);
        MPI_Irecv(&local_grid[1][LY+1], 1, column_type, neighbors[3], 3, comm, &reqs[7]);

        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

        solve_heat_equation(local_lx, local_ly, local_grid, dx, dy, dt);

        if (step % 10 == 0) {
            refine_mesh(local_lx, local_ly, local_grid, dx, dy);
        }
    }

    for (int i = 1; i <= LX; i++) {
        for (int j = 1; j <= LY; j++) {
            smallgrid[i-1][j-1].value = local_grid[i][j].value;
        }
    }

    for (int i = 0; i < LANDSCAPE; i++) {
        for (int j = 0; j < LANDSCAPE; j++) {
            tempgrid[i][j].value = 0;
            tempgrid[i][j].level = 0;
        }
    }

    for (int i = 0; i < LX; i++) {
        for (int j = 0; j < LY; j++) {
            int global_row = coords[0] * LX + i;
            int global_col = coords[1] * LY + j;
            tempgrid[global_row][global_col].value = smallgrid[i][j].value;
        }
    }

    MPI_Reduce(&tempgrid[0][0], &global_grid[0][0], LANDSCAPE*LANDSCAPE, MPI_CELL, MPI_CELL_SUM, 0, comm);

    if (rank == 0) {
        save_to_pbm("heat_equation", global_grid);
    }

    MPI_Op_free(&MPI_CELL_SUM);
    MPI_Type_free(&MPI_CELL);
    MPI_Type_free(&column_type);
    MPI_Type_free(&row_type);
    MPI_Buffer_detach(&buffer, &bsize);
    free(buffer);
    MPI_Finalize();

    return 0;
}



void initialize_grid(int L, Cell grid[L][L]) {
    printf("Entering initialize_grid\n");
    
    // Check allocation sizes
    printf("Grid dimensions: %d x %d\n", L, L);
    
    // Add checks before accessing arrays
    if (grid != NULL) {
        printf("Grid allocated successfully\n");
    } else {
        printf("Error: Grid not allocated\n");
        return;
    }
    
    // Add bounds checking in loops
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            double r = uni();
            if (r < 0.51) {
                grid[i][j].value = 1.5;
                grid[i][j].level = 3;
            } else {
                grid[i][j].value = 8.5;
                grid[i][j].level = ;
            }
        }
    }
    
    printf("Grid initialization complete\n");
}

void initialize_dt(double *dt, double dx, double dy) {
    *dt = 0.25 * (dx * dx + dy * dy); 
}

void solve_heat_equation(int local_lx, int local_ly, Cell local_grid[LX+2][LY+2], double dx, double dy, double dt)  {
    for (int i = 1; i < local_lx - 1; i++) {
        for (int j = 1; j < local_ly - 1; j++) {
            double laplacian_value = laplacian(local_grid, i, j, dx, dy);
            local_grid[i][j].value += dt * laplacian_value;
        }
    }
}


void refine_mesh(int local_lx, int local_ly, Cell local_grid[LX+2][LY+2], double dx, double dy) {
    for (int i = 1; i < local_lx - 1; i++) {
        for (int j = 1; j < local_ly - 1; j++) {
            if (local_grid[i][j].level < MAX_LEVEL && fabs(laplacian(local_grid, i, j, dx, dy)) > 0.1) {
                local_grid[i][j].level++;
                update_neighbors(i, j, local_grid);
            }
        }
    }
}


void update_neighbors(int i, int j, Cell local_grid[LX+2][LY+2]) {
    int neighbors[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int k = 0; k < 4; k++) {
        int ni = i + neighbors[k][0];
        int nj = j + neighbors[k][1];
        if (ni >= 0 && ni < LX && nj >= 0 && nj < LY) {
            local_grid[ni][nj].level = fmax(local_grid[ni][nj].level, local_grid[i][j].level - 1);
        }
    }
}

double laplacian(Cell local_grid[LX+2][LY+2], int i, int j, double dx, double dy) {
    return (local_grid[i+1][j].value + local_grid[i-1][j].value - 2*local_grid[i][j].value) / (dx*dx) +
           (local_grid[i][j+1].value + local_grid[i][j-1].value - 2*local_grid[i][j].value) / (dy*dy);
}



void save_to_pbm(const char* filename, Cell grid[LANDSCAPE][LANDSCAPE]) {
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
    fprintf(fp, "%d %d\n", LANDSCAPE, LANDSCAPE);
    fprintf(fp, "255\n");  // Maximum gray value

    // Write grid values
    for (int j = 0; j < LANDSCAPE; j++) {
        for (int i = 0; i < LANDSCAPE; i++) {
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

void cell_sum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
    Cell *in = (Cell *)invec;
    Cell *inout = (Cell *)inoutvec;
    for (int i = 0; i < *len; i++) {
        inout[i].value += in[i].value;  // Summing the 'value' field
        inout[i].level += in[i].level;  // Summing the 'level' field (or use a different logic if needed)
    }
}

void print_grid(int L, Cell grid[L][L]) {
    printf("Printing grid:\n");
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            printf("%.2f ", grid[i][j].value);
        }
        printf("\n");
    }
    printf("\n");
}
