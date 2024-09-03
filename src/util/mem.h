#ifndef MEM_H
#define MEM_H

#include <stdio.h>
#include <stdlib.h>
#include "structs.h"

// Macro for safely freeing memory
#define FREE(ptr) \
({\
    if(ptr != NULL)\
    {\
        free(ptr);\
        ptr = NULL;\
    }\
})

// Function prototypes
double** allocate_2d_array(int rows, int cols);
buf_str allocate_serial_buffers();
buf_str allocate_parallel_buffers(comm_str comm, slc_str slice, master_str *master);
void dealocate_buffers(buf_str *buffers);
void handle_allocation_failure();

#endif // MEM_H