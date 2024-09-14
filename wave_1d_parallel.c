#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

// Option to change numerical precision.
typedef int64_t int_t;
typedef double real_t;


// TASK: T1b
int rank, size; 

// Arrays to hold the counts and displacements for each process
int *rank_partition_size = NULL;
int *rank_offset = NULL;
// END: T1b


// Simulation parameters: size, step count, and how often to save the state.
const int_t
    N = 65536,
    max_iteration = 100000,
    snapshot_freq = 500;

// Wave equation parameters, time step is derived from the space step.
const real_t
    c  = 1.0,
    dx = 1.0;
real_t
    dt;

// Buffers for three time steps, indexed with 2 ghost points for the boundary.
real_t
    *buffers[3] = { NULL, NULL, NULL };


#define U_prv(i) buffers[0][(i)+1]
#define U(i)     buffers[1][(i)+1]
#define U_nxt(i) buffers[2][(i)+1]


// Convert 'struct timeval' into seconds in double prec. floating point
#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)


// TASK: T8
// Save the present time step in a numbered file under 'data/'.
void domain_save ( int_t step )
{
// BEGIN: T8
    if ( rank != 0 )
    {
        return;
    }
    char filename[256];
    sprintf ( filename, "data/%.5ld.dat", step );
    FILE *out = fopen ( filename, "wb" );
    fwrite ( &U(0), sizeof(real_t), N, out );
    fclose ( out );
// END: T8
}


// TASK: T3
// Allocate space for each process' sub-grids
// Set up our three buffers, fill two with an initial cosine wave,
// and set the time step.
void domain_initialize ( void )
{
// BEGIN: T3
    
    // Allocate memory for the arrays containing the counts and displacements for each process
    rank_partition_size = malloc(size * sizeof(int));
    rank_offset = malloc(size * sizeof(int));

    // Calculate partition sizes and displacements for each process
    int base_partition_size = N / size;
    int remainder = N % size;

    // Adjust partition size for each process based on its rank
    int offset = 0;
    for (int i = 0; i < size; i++) {
        int current_partition_size = (i < remainder) ? base_partition_size + 1 : base_partition_size;

        // Each rank fills the arrays with partition sizes and displacements
        rank_partition_size[i] = current_partition_size;
        rank_offset[i] = offset;

        // Offset is updated for the next process
        offset += current_partition_size;
    }

    if (rank == 0) {
        // Allocate space for the entire domain on rank 0
        buffers[0] = malloc((N + 2) * sizeof(real_t));  // +2 for ghost points
        buffers[1] = malloc((N + 2) * sizeof(real_t));
        buffers[2] = malloc((N + 2) * sizeof(real_t));
    } else {
        // Allocate space for only the partition on other ranks
        buffers[0] = malloc((rank_partition_size[rank] + 2) * sizeof(real_t));  // +2 for ghost points
        buffers[1] = malloc((rank_partition_size[rank] + 2) * sizeof(real_t));
        buffers[2] = malloc((rank_partition_size[rank] + 2) * sizeof(real_t));
    }
    for ( int_t i = 0; i < rank_partition_size[rank]; i++ )
    {
        U_prv(i) = U(i) = cos ( (M_PI*(i+rank_offset[rank])) / (real_t)N );
    }
// END: T3

    // Set the time step for 1D case.
    dt = dx / c;
}


// Return the memory to the OS.
void domain_finalize ( void )
{
    free ( buffers[0] );
    free ( buffers[1] );
    free ( buffers[2] );
}


// Rotate the time step buffers.
void move_buffer_window ( void )
{
    real_t *temp = buffers[0];
    buffers[0] = buffers[1];
    buffers[1] = buffers[2];
    buffers[2] = temp;
}


// TASK: T4
// Derive step t+1 from steps t and t-1.
void time_step ( void )
{
// BEGIN: T4
    for ( int_t i=0; i<rank_partition_size[rank]; i++ )
    {
        U_nxt(i) = -U_prv(i) + 2.0*U(i)
                 + (dt*dt*c*c)/(dx*dx) * (U(i-1)+U(i+1)-2.0*U(i));
    }
// END: T4
}


// TASK: T6
// Neumann (reflective) boundary condition.
void boundary_condition ( void )
{
// BEGIN: T6

    if ( rank == 0 )
    {
        U(-1) = U(1);
    }
    if ( rank == size - 1 )
    {
        U(rank_partition_size[rank]) = U(rank_partition_size[rank]-2);
    }
// END: T6
}


// TASK: T5
// Communicate the border between processes.
void border_exchange( void )
{
// BEGIN: T5
    if ( rank > 0 )
    {
        MPI_Send ( &U(0), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD );
        MPI_Recv ( &U(-1), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }
    if ( rank < size - 1)
    {
        MPI_Send ( &U(rank_partition_size[rank]-1), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD );
        MPI_Recv ( &U(rank_partition_size[rank]), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }
// END: T5
}


// TASK: T7
// Every process needs to communicate its results
// to root and assemble it in the root buffer
void send_data_to_root()
{
    
    if (rank == 0) {
        MPI_Gatherv(MPI_IN_PLACE, rank_partition_size[rank], MPI_DOUBLE,
                    &U(0), rank_partition_size, rank_offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        // Other ranks send their data to rank 0
        MPI_Gatherv(&U(0), rank_partition_size[rank], MPI_DOUBLE,
                    NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

}



// Main time integration.
void simulate( void )
{
    // Go through each time step.
    for ( int_t iteration=0; iteration<=max_iteration; iteration++ )
    {
        if ( (iteration % snapshot_freq)==0 )
        {
            send_data_to_root();
            domain_save ( iteration / snapshot_freq );
        }

        // Derive step t+1 from steps t and t-1.
        border_exchange();
        boundary_condition();
        time_step();

        move_buffer_window();
    }
}


int main ( int argc, char **argv )
{
// TASK: T1c
// Initialise MPI

// BEGIN: T1c
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &size );
// END: T1c
    
    struct timeval t_start, t_end;

    domain_initialize();

// TASK: T2
// Time your code
// BEGIN: T2
    if ( rank == 0 )
    {
        gettimeofday ( &t_start, NULL );
    }
    simulate();
    if ( rank == 0 )
    {
        gettimeofday ( &t_end, NULL );
        printf ( "Elapsed time: %.6f seconds\n", WALLTIME(t_end)-WALLTIME(t_start) );
    }
// END: T2
   
    domain_finalize();

// TASK: T1d
// Finalise MPI
// BEGIN: T1d
    MPI_Finalize();
// END: T1d

    exit ( EXIT_SUCCESS );
}
