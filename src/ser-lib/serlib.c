#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "serlib.h"
#include "mem.h"
#include "heat_solver.h"
#include "misc.h"

void serial_initialise_buffers(master_str *master)
{
	/* Initialise edge slices to buf values 
	 * old and new image arrays to 255(white)
	 * and setup fixed sawtooth boundaries to left and right
	 */
    master->slice = get_serial_dims(master->params.landscape);
    master->cell.buffers = allocate_serial_buffers(master->slice);
	initialize_cell_buffers(master->cell.buffers.global.values, master->cell.buffers.global.levels, master->slice.padded, master->cart, master->params.landscape, master->params.rho, master->params.seed);
    copy_buff_to_local(master->cell.buffers.local.values, master->cell.buffers.local.levels, master->cell.buffers.global.values, master->cell.buffers.global.levels, master->slice);
    initialise_edges(master->cell.buffers.local.values, master->cell.buffers.local.levels, master->slice.actual);

}

void serial_process(master_str *master)
{
    double dt = 0.25 * (1.0 / master->slice.actual.width * 1.0 / master->slice.actual.width + 1.0 / master->slice.actual.height * 1.0 / master->slice.actual.height);

    for (int step = 0; step < master->params.maxstep; step++) {
        
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
    }
}


void serial_write_data(master_str *master)
{

    copy_buff_to_mini(master->cell.buffers.global.values, master->cell.buffers.local.values, master->slice);
    writecelldynamic("heat equation", master->cell.buffers.global.values, master->params.landscape);
    
}



void serial_clean_buffers(master_str *master)
{
	dealocate_buffers(&master->cell.buffers);
}