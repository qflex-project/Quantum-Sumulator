#include "dcpu.h"

#include <mpi.h>

void DCpuExecution1(float complex *state, PT **pts, int qubits, int procs, int it) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // We are assuming at least 2 processes for this task
    if (world_size < 2)
    {
        fprintf(stderr, "World size must be greater than 1 for this to work!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else {
        fprintf(stderr, "Its all right!\n");
    }

    MPI_Finalize();
}