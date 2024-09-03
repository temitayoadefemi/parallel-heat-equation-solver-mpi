#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "serlib.h"
#include "args.h"
#include "mem.h"
#include "pgmio.h"
#include "heat_solver.h"
#include "misc.h"

void serial_initialise_buffers(master_str *master)
{
	/* Initialise edge slices to buf values 
	 * old and new image arrays to 255(white)
	 * and setup fixed sawtooth boundaries to left and right
	 */
	  initialize_cell_buffers(&master->cell.buffers.global.values, &master->cell.buffers.global.levels, &master->slice.actual);
      initialise_edges(&master->cell.buffers.global.values, &master->slice.halo)

}

void serial_process(master_str *master)
{

	 for (int step = 0; step < 1000; step++) {
        
        solve_heat_equation(&master->cell.buffers.local.values, dx, dy, dt, &master->slice);

        if (step % 10 == 0) {

            refine_mesh(&master->cell.buffers.local.levels, &master->cell.buffers.local.values, dx, dy)
        }
    }
}


void serial_write_data(master_str *master)
{

    copy_buff_to_mini(&master->cell.buffers.mini.values, &master->cell.buffers.local.values, &master->slice);
    zerotmpcell(&master->cell.buffers.temp.values, &master->dimensions);
    writecelldynamic("heat equation", &master->cell.buffers.global.values)
    
}



void serial_clean_buffers(master_str *master)
{
	dealocate_buffers(&master->cell.buffers);
}