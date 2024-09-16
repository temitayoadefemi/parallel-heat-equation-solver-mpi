#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parlib.h"
#include "mem.h"
#include "heat_solver.h"
#include "mpilib.h"
#include "misc.h"




void par_initialise_comm(master_str *master)
{
	mpstart(&master->comm);
	setup_cartesian_topology(&master->comm, &master->cart);

}


void par_initialise_buffers(master_str *master)
{

    master->slice = get_parallel_dims(master->cart, master->params.landscape);
    master->cell.buffers = allocate_parallel_buffers(master->slice, master->cart);
    initialize_cell_buffers(master->cell.buffers.global.values, master->cell.buffers.global.levels, master->slice.padded, master->cart, master->params.landscape, master->params.rho, master->params.seed);
    mpbcast(master->cart, master->cell.buffers.global.values, master->slice.padded);
    distribute_cells(master->cell.buffers.mini.values, master->cell.buffers.global.values, master->cart, master->slice);
    copy_buff_to_local(master->cell.buffers.local.values, master->cell.buffers.local.levels, master->cell.buffers.global.values, master->cell.buffers.global.levels, master->slice);
    initialise_edges(master->cell.buffers.local.values, master->cell.buffers.local.levels, master->slice.actual);

}

void par_halo_exchange(master_str *master) {
    MPI_Datatype column_type, row_type;
    void *buffer;
    int bsize;

    initialize_mpi_types(&column_type, &row_type, master->slice);
    initialize_mpi_buffer(&buffer, &bsize, master->slice);

    send_halos_2D(&master->cart, row_type, column_type, master->cell.buffers.local.values, master->slice);
    receive_halos_2D(&master->cart, row_type, column_type, master->cell.buffers.local.values, master->slice);     
    complete_communication_2D(&master->cart);

    send_halos_2D(&master->cart, row_type, column_type, master->cell.buffers.local.levels, master->slice);
    receive_halos_2D(&master->cart, row_type, column_type, master->cell.buffers.local.levels, master->slice); 
    complete_communication_2D(&master->cart);

}

void par_process(master_str *master) {

    double dt = 0.25 * (1.0 / master->slice.actual.width * 1.0 / master->slice.actual.width + 1.0 / master->slice.actual.height * 1.0 / master->slice.actual.height);

    for (int step = 0; step < 1000; step++) {
        
        solve_heat_equation(master->cell.buffers.local.values, 
                            1.0 / master->slice.actual.width, 
                            1.0 / master->slice.actual.height, 
                            dt, master->slice);

        if (step % 10 == 0) {
            refine_mesh(master->cell.buffers.local.levels, 
                        master->cell.buffers.local.values, 
                        1.0 / master->slice.actual.width, 
                        1.0 / master->slice.actual.height, master->slice, master->params.maxlevel);
        }

        par_halo_exchange(master);
    }

}


void par_gather_write_data(master_str *master)
{
    copy_buff_to_mini(master->cell.buffers.mini.values, master->cell.buffers.local.values, master->slice);
    zerotmpcell(master->cell.buffers.temp.values, master->cell.buffers.temp.levels , master->dimensions);
    distribute_cells(master->cell.buffers.local.values, master->cell.buffers.global.values, master->cart, master->slice);
    mpireduce(master->cart, master->cell.buffers.temp.values, master->cell.buffers.global.values, master->slice.padded);

    if (master->comm.rank == 0) {
        writecelldynamic("heat equation", master->cell.buffers.global.values, master->params.landscape);
    }
}

void par_clean_buffers_stop_comm(master_str *master)
{
	dealocate_buffers(&master->cell.buffers);
	mpstop();
}