#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "serlib.h"
#include "args.h"
#include "mem.h"
#include "pgmio.h"
#include "heat_solver.h"

#define PRINTFREQ 200
#define STOP_DELTA 0.1

void serial_initialise_buffers(master_str *master)
{
	/* Initialise edge slices to buf values 
	 * old and new image arrays to 255(white)
	 * and setup fixed sawtooth boundaries to left and right
	 */
	initialise_cell_buffers(master->img.buffers.old, master->img.buffers.new, 255.0, master->slice.halo);
    setup_sawtooth_boundaries(&master->img, master->slice);
}

void serial_process(master_str *master)
{


}

void serial_image_update(master_str *master)
{
	master->slice.ustart.width  = 1;
	master->slice.ustart.height = 1;
	master->slice.uend.width  = master->slice.actual.width  + 1;
	master->slice.uend.height = master->slice.actual.height + 1;

    /* Implement periodic boundary conditions on bottom and top sides */
    setup_periodic_boundaries(&master->img);

    update_image_slice(&master->img, master->slice.ustart, master->slice.uend);
}


double serial_gettime(void)
{ 
  	return MPI_Wtime(); 
} 

void serial_start_timing(master_str *master)
{
#ifdef TIME
    master->time.start = serial_gettime();
#endif
}

void serial_stop_timing(master_str *master)
{
#ifdef TIME
    master->time.average = serial_gettime() - master->time.start;
#endif	
}

void serial_print_timing(master_str master)
{
#ifdef TIME
    printf("Average Time for %d iterations = %f\n", master.params.citer, master.time.average);
#endif	
}


void serial_write_data(master_str *master)
{
    update_buffer(master->img.buffers.local, master->img.buffers.old, master->slice.actual);

    printf("\nWriting <%s>\n", master->img.ifilename); 
    pgmwrite(master->img.ifilename, &master->img.buffers.local[0][0], master->img.size.width, 
    		 master->img.size.height);
	serial_print_timing(*master);
}


void serial_clean_buffers(master_str *master)
{
	dealocate_buffers(&master->img.buffers);
}