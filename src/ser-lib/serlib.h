#ifndef SERLIB_H
#define SERLIB_H

#include "structs.h"  
#include <stdbool.h>
#include <math.h>


void serial_initialise_buffers(master_str *master);
void serial_process(master_str *master);
void serial_write_data(master_str *master);
void serial_clean_buffers(master_str *master);


#endif // SERLIB_H