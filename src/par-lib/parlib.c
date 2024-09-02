#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parlib.h"
#include "args.h"
#include "mem.h"
#include "pgmio.h"
#include "heat_solver.h"
#include "mpilib.h"



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
    initialize_cell_buffers(&master->cell.buffers.global.values, &master->cell.buffers.global.levels, &master->dimensions )
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
        
        sen
        receive_halos_2D(&master->cart, hor_halo_type, &master->cell.buffers.local);
        complete_communication_2D(&master->cart);

        solve_heat_equation(&master->cell.buffers.local, double dx, double dy, double dt);

        if (step % 10 == 0) {

            refine_mesh(&master->cell.buffers.local, double dx, double dy);
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

}

void par_clean_buffers_stop_comm(master_str *master)
{
	dealocate_buffers(&master->img.buffers);
	mpstop();
}