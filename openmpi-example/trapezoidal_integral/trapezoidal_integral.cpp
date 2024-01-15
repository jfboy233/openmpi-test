#include<iostream>
using namespace std;

#include"mpi.h"
#include<cmath>

double f(double x){
    return pow(x,3);
}

double Trap(double left_endpt, double right_endpt, double trap_count, double base_len){

    double estimate = ( f(left_endpt) + f(right_endpt) )/2, x = left_endpt;

    for(int i=1;i<=trap_count-1;i++){
        x += base_len;
        estimate += f(x);
    }
    estimate = estimate * base_len;
    return estimate;
}

int main(){

    int my_rank = 0, comm_sz = 0, n = 1E7, local_n = 0;
    double a = 0, b = 3, h = 0, local_a = 0, local_b = 0;
    double local_int = 0, total_int = 0;
    int source;
    clock_t t_start = clock();

    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD, & my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, & comm_sz );

    double start_time = MPI_Wtime();

    h = (b-a)/n; 
    local_n = n / comm_sz;
    if( my_rank == comm_sz - 1 ){
        local_n = n - local_n * (comm_sz-1);
    }// when n is not indivisible by comm_sz, the last chunk is larger than others.

    local_a = a + my_rank * ( n / comm_sz ) * h;
    local_b = local_a + local_n * h;
    local_int = Trap( local_a, local_b, local_n, h );
    //printf("%d of %d processes, local_int = %.15e\n", my_rank, comm_sz, local_int);

    if( my_rank != 0 ){
        MPI_Send( & local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
    }
    else{
        total_int = local_int;
        for( source = 1; source < comm_sz; source ++ ){
            MPI_Recv( & local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            total_int += local_int;
        }

        printf("With n = %d trapezoids, our estimate \n", n);
        printf("of the integral from %f to %f = %.15e\n", a, b, total_int );

        double end_time = MPI_Wtime();
        printf("That took %lf seconds.\n", end_time - start_time);    
    }

    MPI_Finalize();

    return 0;
}