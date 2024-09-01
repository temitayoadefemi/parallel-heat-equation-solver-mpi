#ifndef __STRUCTS_H__
#define __STRUCTS_H__

#include <mpi.h>
#define ndims 2 

typedef enum version_enum
{
	serial,
	par2D
}version;

typedef struct dimensions_struct
{
	int width;
	int height;
} dim_str;


typedef struct cell_struct {
    double value;
    int level;
} cell;


/* Contains the inforamtion about a neighbour */
typedef struct dim_comm_struct
{
	int val;
	MPI_Status sstat,rstat;
	MPI_Request sreq,rreq;
}dir_str;

/* Variables for 2D Cartesian topology */
typedef struct cartesian_struct
{
	MPI_Comm comm2d;
	int disp;
	int dims[ndims];
	int period[ndims];
	int coords[ndims];
	int reorder;
	dir_str r, l, b, t;
} cart_str;

/* Communication variables */
typedef struct communication_struct
{
	MPI_Comm comm;
	int rank, size;
	dir_str r, l;
}comm_str;

typedef struct time_struct
{
	double start;
	double local;
	double average;
}time_str;


typedef struct master {

    params_str params;
    comm_str comm;
    cart_str cart;
    dim_str dimensions;
    int version;
	time_str time;
} master_str;

#endif	//__STRUCTS_H__