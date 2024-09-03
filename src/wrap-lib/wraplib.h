#ifndef __WRAPLIB_H__
#define __WRAPLIB_H__

#include "structs.h"


void setup_comm(master_str *master);
void setup_buffers(master_str *master);
void process(master_str *master);
void gather_write_data(master_str *master);
void clean_buffers_end_comm(master_str *master);

#endif //__WRAPLIB_H__