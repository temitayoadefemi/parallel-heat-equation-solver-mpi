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

typedef struct grid_struct 
{
	double **values;
	double **levels;

}grid_str;

typedef struct buf_struct
{
	grid_str global;
	grid_str local;
	grid_str temp;
	grid_str mini;

}buf_str;

/* structure containing all the neccessary 
 * information about the image
 */
typedef struct cell_struct
{
	buf_str buffers;
	dim_str size;

} cell_str;
/* structure containing all the neccessary 
 * information about the sliced image
 */
typedef struct slice_struct
{
	dim_str padded;
	dim_str actual;
	dim_str rem;
	dim_str halo;
	dim_str ustart,uend;
}slc_str;

typedef struct parameters
{
	  int seed;
      double rho;
	  int printfreq;
	  int landscape;
	  int maxstep;
	  double r;
	  version version;
} params_str;

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
	dir_str right, left, down, up;
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

typedef struct master_struct
{
	cell_str cell;
	slc_str slice;
	params_str params;
	dim_str dimensions;
	comm_str comm;
	cart_str cart;

	time_str time;
}master_str;

#endif	//__STRUCTS_H__