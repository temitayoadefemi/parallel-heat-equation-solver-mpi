#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "arralloc.h"


#define FREE(ptr) \
({\
	if(ptr != NULL)\
	{\
		free(ptr);\
		ptr = NULL;\
	}\
})

double** allocate_2d_array(int rows, int cols) {
    double **array = (double**) arralloc(sizeof(double), 2, rows, cols);
    if (array == NULL) {
        handle_allocation_failure();
    }
    return array;
}



buf_str allocate_serial_buffers()
{
	buf_str buffer;

	buffer.global = NULL;
	buffer.local = NULL;
    buffer.temp = NULL;
    buffer.mini = NULL;

	buffer.global= allocate_2d_array(local_grid.halo.width, local_grid.halo.width);
	buffer.local = allocate_2d_array(local_grid.halo.width, local_grid.halo.width);
	buffer.temp = allocate_2d_array(local_grid.halo.width, local_grid.halo.width);
    buffer.mini = allocate_2d_array(local_grid.halo.width, local_grid.halo.width);


	return buffer;
}

buf_str allocate_parallel_buffers(comm_str comm, slc_str slice, master_str *master )
{
	buf_str buffer;

	buffer.master = NULL;
	buffer.local = NULL;
	buffer.temp = NULL;
	buffer.mini = NULL;

	/* Allocate arrays dynamically using special malloc routine. */
	buffer.master =  allocate_2d_array(master->dimensions.width, master->dimensions.height);
	buffer.local  = allocate_2d_array(slice.padded.width, slice.padded.height);
	buffer.temp  = allocate_2d_array(slice.halo.width, slice.halo.height);
	buffer.mini  = allocate_2d_array(slice.halo.width, slice.halo.height);

	return buffer;
}

void dealocate_buffers(buf_str *buffers) {

	FREE(buffers->master);
	FREE(buffers->local);
	FREE(buffers->temp);
	FREE(buffers->mini);
}

// Handle memory allocation failure.
void handle_allocation_failure() {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
}

// Define a function to allocate memory for a 2D array.
