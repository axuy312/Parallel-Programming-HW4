#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>

const int MAX = std::sqrt(RAND_MAX);

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: init MPI
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank > 0)
    {
        // TODO: handle workers
		pi_result = 0;
		unsigned int seed = world_rank;//time(NULL);
		for ( long long int toss = tosses / world_size; toss > 0; toss--) {
			int r = rand_r(&seed);
			double x = (double(r % MAX) / MAX - 0.5);
			double y = (double(r / MAX) / MAX - 0.5);
			
			double distance_squared = x * x + y * y;
			if ( distance_squared <= 0.250000)
				pi_result++;
		}
		MPI_Send(&pi_result, 1, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD);
    }
    else if (world_rank == 0)
    {
        // TODO: master
		pi_result = 0;
		unsigned int seed = world_rank;//time(NULL);
		for ( long long int toss = tosses / world_size; toss > 0; toss--) {
			int r = rand_r(&seed);
			double x = (double(r % MAX) / MAX - 0.5);
			double y = (double(r / MAX) / MAX - 0.5);
			
			double distance_squared = x * x + y * y;
			if ( distance_squared <= 0.250000)
				pi_result++;
		}
    }

    if (world_rank == 0)
    {
        // TODO: process PI result
		for ( int i = 1; i < world_size; i++)
		{
			double recv = 0;
			MPI_Recv(&recv, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			pi_result += recv;
		}
		
		pi_result = (pi_result * 4.0) / tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
