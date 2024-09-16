#ifndef MPILIB_H
#define MPILIB_H


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "structs.h"

void mpstart(comm_str *comm);
void mpstop(void);
void setup_cartesian_topology(comm_str *comm, cart_str *cart);
void send_halos_2D(cart_str *cart, MPI_Datatype row_type, MPI_Datatype column_type, double **values, slc_str slice);
void receive_halos_2D(cart_str *cart, MPI_Datatype row_type, MPI_Datatype column_type, double **values, slc_str slice);
void complete_communication_2D(cart_str *cart);
void initialize_mpi_types(MPI_Datatype *column_type, MPI_Datatype *row_type, slc_str slice);
void initialize_mpi_buffer(void **buffer, int *bsize, slc_str slice);
void mpireduce(cart_str cart, double **inbuff, double **outbuff, dim_str dims);
void mpbcast(cart_str cart, double **inbuff, dim_str dims);
double gettime(void);

#endif // MPILIB_H