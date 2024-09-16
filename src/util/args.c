#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"

#define RHO 0.51
#define PRINTFREQ 500
#define LANDSCAPE 1152
#define STEP_MULTIPLIER 10

// Reads parameters from command-line arguments and initializes them into the master structure
int read_parameters(master_str *master, int argc, char **argv) {
    // Check for minimum number of arguments
    if (argc < 2) {
        // Only the master node outputs the usage message
        if (master->comm.rank == 0) {
            printf("Usage: automaton <seed> [-rho value] [-printfreq value] [-landscape value] [-maxstep value]\n");
        }
        return 1;  // Return 1 to indicate failure due to insufficient arguments
    }

    // Initialize parameters with default values
    master->params.seed = atoi(argv[1]); // Seed is mandatory and taken from the first command-line argument
    master->params.rho = RHO;           // Default density
    master->params.printfreq = PRINTFREQ;      // Default print frequency
    master->params.landscape = LANDSCAPE;     // Default landscape size
    master->params.maxstep = STEP_MULTIPLIER * master->params.landscape;  // Default number of steps

    // Determine the version based on the number of processes
    if (master->comm.size > 1) {
        master->params.version = par2D;  // Parallel version
    } else if (master->comm.size == 1) {
        master->params.version = serial; // Serial version
    }

    // Loop over remaining arguments to set optional parameters
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-rho") == 0 && i + 1 < argc) {
            master->params.rho = atof(argv[++i]);  // Set density
        } else if (strcmp(argv[i], "-printfreq") == 0 && i + 1 < argc) {
            master->params.printfreq = atoi(argv[++i]);  // Set print frequency
        } else if (strcmp(argv[i], "-landscape") == 0 && i + 1 < argc) {
            master->params.landscape = atoi(argv[++i]);  // Set landscape size
        } else if (strcmp(argv[i], "-maxstep") == 0 && i + 1 < argc) {
            master->params.maxstep = atoi(argv[++i]);  // Set maximum steps
        }
    }

    return 0;  // Return 0 to indicate successful completion
}