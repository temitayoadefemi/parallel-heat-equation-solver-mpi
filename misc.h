#ifndef MISC_H
#define MISC_H

// Function declarations for miscellaneous utilities:

// Writes the dynamic state of a cell array to a specified file
// Parameters:
//    cellfile - the file path where the cell data will be written
//    cell - pointer to the 2D array of cells
//    l - the length of the array (assumed square for simplicity)
void writecelldynamic(char *cellfile, int **cell, int l);

// Seeds the random number generator with a specific integer
// Parameter:
//    ijkl - the seed value
void rinit(int ijkl);

// Returns a uniformly distributed random number between 0.0 and 1.0
float uni(void);

#endif // MISC_H