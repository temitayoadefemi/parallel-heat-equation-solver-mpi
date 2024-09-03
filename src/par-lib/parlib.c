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
    initialize_cell_buffers(&master->cell.buffers.global.values, &master->cell.buffers.global.levels, &master->dimensions);
    distribute_cells(&master->cell.buffers.mini.values, &master->cell.buffers.global.values, &master->cart, &master->slice)
}

void par_halo_exchange(master_str *master)
{

    send_halos_2D(&master->cart, row_type, column_type);
    receive_halos_2D(&master->cart, hor_halo_type, &master->img, master->slice);     
    complete_communication_2D(&master->cart);

}

void par_process(master_str *master) {

	 for (int step = 0; step < 1000; step++) {
        
        par_halo_exchange(master);
        solve_heat_equation(&master->cell.buffers.local.values, dx, dy, dt, &master->slice);

        if (step % 10 == 0) {

            refine_mesh(&master->cell.buffers.local.levels, &master->cell.buffers.local.values, dx, dy)
        }
    }
}


void par_gather_write_data(master_str *master)
{
    copy_buff_to_mini(&master->cell.buffers.mini.values, &master->cell.buffers.local.values, &master->slice);
    zerotmpcell(&master->cell.buffers.temp.values, &master->dimensions);
    distribute_cells(&master->cell.buffers.local.values, &master->cell.buffers.global.values, &master->cart, &master->slice);
    mpireduce(&master->cart, &master->cell.buffers.temp.values, &master->cell.buffers.global.values, &master->dimensions.width*&master->dimensions.height);


    if (master->comm.rank == 0) {
        save_to_pbm("heat_equation", global_grid);
    }
}

void par_clean_buffers_stop_comm(master_str *master)
{
	dealocate_buffers(&master->img.buffers);
	mpstop();
}