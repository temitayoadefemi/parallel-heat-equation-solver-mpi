#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "wraplib.h"
#include "args.h"
#include "mpilib.h"


int main(int argc, char **argv)
{

	master_str master;

	setup_comm(&master);

	if (read_parameters(&master, argc, argv) == FAILED) {

        mpstop(); // Stop the MPI environment

        return 0;
    }

	setup_buffers(&master);

	process(&master);

	gather_write_data(&master);

	clean_buffers_end_comm(&master);

	return 0;
}