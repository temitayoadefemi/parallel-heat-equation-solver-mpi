#ifndef PARLIB_H
#define PARLIB_H

#include "structs.h"  
#include <stdbool.h>
#include <math.h>


void par_initialise_comm(master_str *master);

void par_initialise_buffers(master_str *master);

void par_halo_exchange(master_str *master);

void par_process(master_str *master);

void par_gather_write_data(master_str *master);

void par_clean_buffers_stop_comm(master_str *master);

#endif // PARLIB_H