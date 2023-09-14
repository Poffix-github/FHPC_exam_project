#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 10

int main(){
    int c_size, rank;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &c_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    /* partition grid into blocks and scatter them through the processes*/
    const int NDEC = 5; /* number of blocks in a coloumn in decomposition */
    const int BLOCKSIZE = SIZE/NDEC; /* number of rows and columns in a block */
    char board[SIZE*SIZE];

    if (rank == 0) {
        for (int ii=0; ii<SIZE*SIZE; ii++) {
            board[ii] = (char)ii;
        }
    }

    if(c_size != NDEC*NDEC){
        fprintf(stderr,"Error: number of PEs %d != %d x %d\n", c_size, NDEC, NDEC);
        MPI_Finalize();
        exit(-1);
    }

    char chunk[BLOCKSIZE*BLOCKSIZE];
    for (int ii=0; ii<BLOCKSIZE*BLOCKSIZE; ii++) chunk[ii] = 0;

    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    MPI_Type_vector(BLOCKSIZE, BLOCKSIZE, SIZE, MPI_CHAR, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(char), &blocktype);
    MPI_Type_commit(&blocktype);

    int disps[NDEC*NDEC];
    int counts[NDEC*NDEC];
    for (int ii=0; ii<NDEC; ii++) {
        for (int jj=0; jj<NDEC; jj++) {
            disps[ii*NDEC+jj] = ii*SIZE*BLOCKSIZE+jj*BLOCKSIZE;
            counts [ii*NDEC+jj] = 1;
        }
    }

    MPI_Scatterv(board, counts, disps, blocktype, chunk, BLOCKSIZE*BLOCKSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* printing, each process prints its "chunk" */
    for (int proc=0; proc<c_size; proc++) {
        if (proc == rank) {
            printf("Rank = %d\n", rank);
            if (rank == 0) {
                printf("Global matrix: \n");
                for (int ii=0; ii<SIZE; ii++) {
                    for (int jj=0; jj<SIZE; jj++) {
                        printf("%3d ",(int)board[ii*SIZE+jj]);
                    }
                    printf("\n");
                }
            }
            printf("Local Matrix:\n");
            for (int ii=0; ii<BLOCKSIZE; ii++) {
                for (int jj=0; jj<BLOCKSIZE; jj++) {
                    printf("%3d ",(int)chunk[ii*BLOCKSIZE+jj]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
}