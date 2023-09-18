#include "include/gol.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>



#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1


int   action = 0;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;

void update_data(const int size, const int num_proc, const int num_threads, const double time);
int get_time(double tbegin, double tend);

int main( int argc, char **argv )
{    
    int num_proc, rank, mpi_provided_thread_level;
    double tstart, tend;

    MPI_Init_thread(NULL,NULL, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);

    if(rank == 0) tstart = MPI_Wtime();

    if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) {
        printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n");
        MPI_Finalize();
        exit( 1 );
    }

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank == 0){
        get_args(argc, argv, &action, &k, &e, &n, &s, &fname);
    }

    MPI_Bcast(&action, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if( action == INIT && rank == 0){
        printf("generating image\n");
        void *board = random_board(k, fname);
        
        printf("saving image\n");
        write_pgm_image(board, 255, k, k, fname);
        
        free(board);
    }else{
        /* board stats */
        int     maxval = 0;
        int     size = 0;
        void    *board;
        
        if(rank == 0){
            printf("reading image\n");
            read_pgm_image( &board, &maxval, &size, &size, fname);
        }

        MPI_Bcast(&maxval, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        switch( maxval ){
            case -1: printf("I/O error in header\n"); break;
            case -2: printf("Memory not sufficient\n"); break;
            case -3: printf("I/O error in body\n"); break;
            default:
                if(rank == 0){
                    /* swap endianism */
                    if ( LITTLE_ENDIAN )swap_image( board, size, size, maxval);
                    printf("start evolution\n");
                }
                
                if( e == ORDERED ){
                    if(rank == 0) evolution_ordered(board, size, n, maxval, s);
                }else{ 
                    evolution_static(board, size, n, maxval, s, num_proc, rank);
                }
        }
        
        if(rank == 0) free(board);
    
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0){
            tend = MPI_Wtime();
            update_data(size, num_proc, omp_get_max_threads(), get_time(tstart, tend));
        }
    }
    
    if (rank == 0 && fname != NULL ) free ( fname );


    MPI_Finalize();
}

void update_data(const int size, const int num_proc, const int num_threads, const double time){
    FILE *data_file;

    data_file = fopen("./data.csv", "a");

    fprintf(data_file, "%d,%d,%d,%lf\n", size, num_proc, num_threads, time);

    fclose(data_file);
}

int get_time(double tbegin, double tend){
    return (int)((tend - tbegin)*1000000.0);
}

// for (int ii=0; ii<size; ii++) {
//                         for (int jj=0; jj<size; jj++) printf("%3d ",((unsigned char*)board)[ii*size+jj]);
//                         printf("\n");
//                     }