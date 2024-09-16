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

// Handle memory allocation failure.
void handle_allocation_failure() {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
}


double** allocate_2d_array(int rows, int cols) {
    double **array = (double**) arralloc(sizeof(double), 2, rows, cols);
    if (array == NULL) {
        handle_allocation_failure();
    }
    return array;
}



buf_str allocate_serial_buffers(slc_str slice) {

	buf_str buffer;

	buffer.global.values = NULL;
	buffer.local.values = NULL;
	buffer.temp.values = NULL;
    buffer.mini.values = NULL;

    buffer.global.levels = NULL;
	buffer.local.levels = NULL;
	buffer.temp.levels = NULL;
    buffer.mini.levels = NULL;

	buffer.global.values = allocate_2d_array(slice.actual.width, slice.actual.height );
	buffer.local.values = allocate_2d_array(slice.halo.width, slice.halo.height);
    buffer.mini.values = allocate_2d_array(slice.actual.width, slice.actual.height);

    buffer.global.levels = allocate_2d_array(slice.actual.width, slice.actual.height);
	buffer.local.levels = allocate_2d_array(slice.halo.width, slice.halo.height);
    buffer.mini.levels = allocate_2d_array(slice.actual.width, slice.actual.height);


    return buffer;
}

buf_str allocate_parallel_buffers(slc_str slice, cart_str cart) {

	buf_str buffer;

	buffer.global.values = NULL;
	buffer.local.values = NULL;
	buffer.temp.values = NULL;
    buffer.mini.values = NULL;

    buffer.global.levels = NULL;
	buffer.local.levels = NULL;
	buffer.temp.levels = NULL;
    buffer.mini.levels = NULL;

	/* Allocate arrays dynamically using special malloc routine. */
	buffer.global.values = allocate_2d_array(cart.dims[0] * slice.padded.width, cart.dims[1] * slice.padded.height);
	buffer.local.values = allocate_2d_array(slice.halo.width, slice.halo.height);
	buffer.temp.values = allocate_2d_array(cart.dims[0] * slice.padded.width, cart.dims[1] * slice.padded.height);
    buffer.mini.values = allocate_2d_array(slice.actual.width, slice.padded.width);

    buffer.global.levels = allocate_2d_array(cart.dims[0] * slice.padded.width, cart.dims[1] * slice.padded.height);
	buffer.local.levels = allocate_2d_array(slice.halo.width, slice.halo.height);
	buffer.temp.levels = allocate_2d_array(cart.dims[0] * slice.padded.width, cart.dims[1] * slice.padded.height);
    buffer.mini.levels = allocate_2d_array(slice.actual.width, slice.padded.width);


	return buffer;
}

void dealocate_buffers(buf_str *buffers) {

	FREE(buffers->global.values);
	FREE(buffers->local.values);
	FREE(buffers->temp.values);
	FREE(buffers->mini.values);

    FREE(buffers->global.levels);
	FREE(buffers->local.levels);
	FREE(buffers->temp.levels);
	FREE(buffers->mini.levels);
}


// Define a function to allocate memory for a 2D array.
