#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "wraplib.h"
#include "serlib.h"
#include "parlib.h"



void setup_comm(master_str *master)
{
	if(master->params.version == par2D)
	{
		par_initialise_comm(master);
	}
}


void setup_buffers(master_str *master)
{
	if(master->params.version == serial)
	{
		serial_initialise_buffers(master);

	}else if(master->params.version == par2D)
	{
		par_initialise_buffers(master);
	}
}


void process(master_str *master)
{
	if(master->params.version == serial)
	{
		serial_process(master);

	}else if(master->params.version == par2D)
	{
		par_process(master);
	}
}

void gather_write_data(master_str *master)
{
	if(master->params.version == serial)
	{
		serial_write_data(master);

	}else if(master->params.version == par2D)
	{
		par_gather_write_data(master);
	}
}

void clean_buffers_end_comm(master_str *master)
{
	if(master->params.version == serial)
	{
		serial_clean_buffers(master);

	}else if(master->params.version == par2D)
	{
		par_clean_buffers_stop_comm(master);
	}
}