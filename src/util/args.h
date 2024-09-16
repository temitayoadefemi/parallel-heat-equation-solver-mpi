#ifndef ARGS_H
#define ARGS_H

#include "structs.h"  // Include structs.h if it defines master_str

// Function declaration for reading and parsing command-line parameters
int read_parameters(master_str *master, int argc, char **argv);

#endif // ARGS_H