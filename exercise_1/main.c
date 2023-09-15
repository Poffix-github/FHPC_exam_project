#include "include/gol.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>



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



int main( int argc, char **argv )
{    
    int num_proc, rank;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank == 0){
        get_args(argc, argv, &action, &k, &e, &n, &s, &fname);
    }
    
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
        
        switch( maxval ){
            case -1: printf("I/O error in header\n"); break;
            case -2: printf("Memory not sufficient\n"); break;
            case -3: printf("I/O error in body\n"); break;
            default:
                if(rank == 0){
                    /* swap endianism */
                    if ( LITTLE_ENDIAN ) swap_image( board, size, size, maxval);
                    printf("start evolution\n");
                }
                
                if( e == ORDERED ){
                    if(rank == 0) evolution_ordered(board, size, n, maxval, s);
                }else evolution_static(board, size, n, maxval, s, num_proc, rank);

                printf("rank %d evolution done done", rank);
        }
        
        if(rank == 0) free(board);
    }
    
    if (rank == 0 && fname != NULL ) free ( fname );

    MPI_Finalize();
}
