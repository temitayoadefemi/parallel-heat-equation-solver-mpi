#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parlib.h"
#include "args.h"
#include "mem.h"
#include "pgmio.h"
#include "heat_solver.h"
#include "mplib.h"

#define PRINTFREQ 200
#define MASTER_PROCESS 0
#define STOP_DELTA 0.1
#define TRUE  1
#define FALSE 0

void par_initialise_comm(master_str *master)
{
	mpstart(&master->comm);
	setup_cartesian_topology(&master->comm, &master->cart);
}


void par_initialise_buffers(master_str *master)
{
	/* Initialise edge slices to buf values 
	 * old and new image arrays to 255(white)
	 * and setup fixed sawtooth boundaries to left and right
	 */
	initialise_cell_buffers(master->img.buffers.old, master->img.buffers.new, 255.0, master->slice.halo);
	setup_parallel_sawtooth_boundaries(master->cart, &master->img, master->slice);
}

void par_get_ustart_uend(master_str *master)
{
#ifdef OVERLAP
	master->slice.ustart.width  = 2;
	master->slice.ustart.height = 2;
	master->slice.uend.width  = master->slice.actual.width;
	master->slice.uend.height = master->slice.actual.height;
#else
	master->slice.ustart.width  = 1;
	master->slice.ustart.height = 1;
	master->slice.uend.width  = master->slice.actual.width  + 1;
	master->slice.uend.height = master->slice.actual.height + 1;
#endif	
}

void par_halo_exchange_image_update(master_str *master, MPI_Datatype hor_halo_type)
{
	par_get_ustart_uend(master);

    send_halos_2D(&master->cart, hor_halo_type, &master->img, master->slice);
    receive_halos_2D(&master->cart, hor_halo_type, &master->img, master->slice);

#ifdef OVERLAP
    /* update only the entries not on the boundaries until communication is complete*/
    update_image_slice(&master->img, master->slice.ustart, master->slice.uend);
#endif       

    complete_communication_2D(&master->cart);

#ifdef OVERLAP	
	/* update the boundaries received*/	
	update_image_slice_boundaries(&master->img, master->slice.actual);
#else
	update_image_slice(&master->img, master->slice.ustart, master->slice.uend);
#endif
}



void par_process(master_str *master) {

	 for (int step = 0; step < 1000; step++) {
        
        send_halos_2D(&master->cart, hor_halo_type, &master->img, master->slice);
        receive_halos_2D(&master->cart, hor_halo_type, &master->img, master->slice);
        complete_communication_2D(&master->cart);

        solve_heat_equation(cell **local_grid, double dx, double dy, double dt);

        if (step % 10 == 0) {
            
            refine_mesh(cell **local_grid, double dx, double dy);
        }
    }
}

void par_start_timing(master_str *master)
{
#ifdef TIME
    MPI_Barrier(master->cart.comm2d);
    master->time.start = gettime();
#endif
}

void par_stop_timing(master_str *master)
{
#ifdef TIME
    MPI_Barrier(master->cart.comm2d);
    master->time.local = gettime() - master->time.start;
    master->time.average = mpgsum(master->cart, &master->time.local) / master->comm.size;
#endif	
}

void par_print_timing(master_str master)
{
#ifdef TIME
    printf("Average Time for %d iterations = %f\n", master.params.citer, master.time.average);
#endif	
}


void par_gather_write_data(master_str *master)
{
    update_buffer(master->img.buffers.local, master->img.buffers.old, master->slice.actual);

  	if(master->comm.rank == MASTER_PROCESS)
  	{
  		printf("\nFinished %d iterations\nGathering data...", master->params.citer);
	}

	/* gather the results from all slices into master buffer */
    mpgather(master->cart, &master->img.buffers.local[0][0], &master->img.buffers.master[0][0], 
    		 master->slice.padded.width * master->slice.padded.height);

  	if(master->comm.rank == MASTER_PROCESS)
  	{
		printf("\t Done\nWriting <%s>\n", master->img.ifilename); 
	    pgmwrite_generalised_cascaded(master->img.ifilename, &master->img.buffers.master[0][0], master->cart, master->img.size, master->slice);
	    par_print_timing(*master);
	}
}

void par_clean_buffers_stop_comm(master_str *master)
{
	dealocate_buffers(&master->img.buffers);
	mpstop();
}