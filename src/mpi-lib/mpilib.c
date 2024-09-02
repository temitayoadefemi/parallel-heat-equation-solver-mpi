#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "structs.h"

#include "mpilib.h"

#define TRUE  1
#define FALSE 0
#define ndims 2   
#define vertical    1 
#define horizontal  0

MPI_Datatype hor_halo_type;

void mpstart(comm_str *comm)
{ 
  MPI_Init(NULL, NULL);
  comm->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm->comm, &comm->rank);
  MPI_Comm_size(comm->comm, &comm->size);
}
  
void mpstop(void)
{ 
  MPI_Finalize();
}

void setup_cartesian_topology(comm_str *comm, cart_str *cart)
{
	/* Setup cartesian topology */
	cart->dims[0] = 0;
	cart->dims[1] = 0;
	cart->period[0] = FALSE; 	// Acyclic to form a line in the horizontal
	cart->period[1] = FALSE; 	// Cyclic to form a ring in the vertical  
	cart->reorder = FALSE;    // Do not reorder
	cart->disp = 1;           // shift by 1

	/* Set the dimensions, create the 2D Cartesian topology 
	 * and get the rank of the current process
	 */
	MPI_Dims_create(comm->size, ndims, cart->dims);
	MPI_Cart_create(comm->comm, ndims, cart->dims, cart->period, cart->reorder, &cart->comm2d);

	MPI_Comm_rank(cart->comm2d, &comm->rank);
	MPI_Cart_coords(cart->comm2d, comm->rank, ndims, cart->coords);

	/* First compute the ranks of the neighbouring processes 
	 * along the horizontal axis (left to right) 
	 * and then along vertical axis (top to bottom)
	 */
	MPI_Cart_shift(cart->comm2d, horizontal, cart->disp, &cart->l.val, &cart->r.val);
	MPI_Cart_shift(cart->comm2d, vertical, cart->disp, &cart->b.val, &cart->t.val);
}

void send_halos_2D(cart_str *cart, MPI_Datatype row_type, MPI_Datatype column_type, **double values, slc_str slice)
{
	/*
	 * Use non-blocking communications.
	 * First use non-blocking asynchronous sent to sent the vertical halo values 
	 * from left to right and then from right to left
	 * Then sent the horisontal halo values 
	 * from top to bottom and then from bottom to top
	 */ 
    MPI_Issend(&values[slice.actual.width][1], 1, row_type, cart->down.val, 0, cart->comm2d, &cart->reqs[0]);
    MPI_Issend(&values[1][1], 1, row_type, cart->up.val, 1, cart->comm2d, &cart->reqs[1]);
    MPI_Issend(&values[1][1], 1, column_type, cart->left.val, 3, cart->comm2d, &cart->reqs[2]);
    MPI_Issend(&values[1][slice.actual.height], 1, column_type, cart->right.val, 2, cart->comm2d, &cart->reqs[3]);

}		

void receive_halos_2D(cart_str *cart, MPI_Datatype row_type, MPI_Datatype column_type, **double values, slc_str slice)
{

	MPI_Irecv(&values[0][1], 1, row_type, cart->up.val, 0, cart->comm2d, &cart->reqs[4]);
    MPI_Irecv(&values[slice.actual.width + 1][1], 1, row_type, cart->down.val, 1, cart->comm2d, &cart->reqs[5]);
    MPI_Irecv(&values[1][slice.actual.height + 1], 1, column_type, cart->right.val, 3, cart->comm2d, &cart->reqs[6]);
    MPI_Irecv(&local_grid[1][0], 1, column_type, cart->left.val, 2, cart->comm2d, &cart->reqs[7]);


}

void complete_communication_2D(cart_str *cart)
{
	/*
	 * Wait for left to right and right to left sends to complete
	 * Wait for bottom to top and bottom to top sends to complete
	 */ 
	MPI_Wait(&cart->l.sreq, &cart->r.sstat);
	MPI_Wait(&cart->r.sreq, &cart->l.sstat);
	MPI_Wait(&cart->t.sreq, &cart->b.sstat);
	MPI_Wait(&cart->b.sreq, &cart->t.sstat);
	MPI_Wait(&cart->l.rreq, &cart->r.rstat);
	MPI_Wait(&cart->r.rreq, &cart->l.rstat);
	MPI_Wait(&cart->t.rreq, &cart->b.rstat);
	MPI_Wait(&cart->b.rreq, &cart->t.rstat);

}

void initialize_mpi_types(MPI_Datatype *column_type, MPI_Datatype *row_type, master_str *master) {
    // Create a vector type for transferring columns.
    
	MPI_Datatype column_type, row_type;
    MPI_Type_vector(master->dimensions.rows, 1, master-dimensions.cols + 2, MPI_CELL, &column_type);
    MPI_Type_commit(&column_type);

    MPI_Type_contiguous(master->dimensions.cols, MPI_CELL, &row_type);
    MPI_Type_commit(&row_type);
    
}
void initialize_mpi_buffer(void **buffer, int *bsize, master_str *master) {
    int bsize = (master->dimensions.rows + master->dimensions.cols) * sizeof(Cell) + MPI_BSEND_OVERHEAD;
    void* buffer = malloc(bsize);
    if (buffer == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

     MPI_Buffer_attach(buffer, bsize);

}


void mpbcast(cart_str cart, cell **inbuff, int ndata)
{
	MPI_Bcast(&inbuff[0][0], size*size, MPI_INT, 0, cart.comm2d);
}


double gettime(void)
{ 
  return MPI_Wtime(); 
} 