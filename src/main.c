#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "wraplib.h"

int main(int argc, char **argv)
{
	master_str master;

	setup_comm(&master);

	setup_buffers(&master);

	process(&master);

	gather_write_data(&master);

	clean_buffers_end_comm(&master);

	return 0;
}