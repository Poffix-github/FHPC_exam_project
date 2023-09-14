#include <gol.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

/* ==================================================
 * These functions evolve the board.
 * 
 * CELL NOMENCLATURE
 * 0:   dead
 * 255: alive
 * 127: is dead and will live next update
 * 128: is alive and will die next update
 *  
 * NOTE: to check if a cell is alive or dead at the current state, 
 * it is necessary to check if the cell is >=128 or <=127 respectively.
 * ==================================================
 */

#define MAX_NUMBER_OF_CELLS 8


/* Count the number of alive neighbours 
 */
char check_neighbours(const void* board, const int DIM, const int i, const int j){
    const int off_sets[MAX_NUMBER_OF_CELLS][2] = {  {-1, -1}, {-1, 0}, {-1, 1},
                                                     {0, -1},           {0, 1},
                                                     {1, -1},  {1, 0},  {1, 1}  };
    int count=0;    /* number of alive cells around the ij-th */
    
    for(int z=0; z<MAX_NUMBER_OF_CELLS; z++){
        /* compute coordinates couple */
        int k = i + off_sets[z][0] < 0 ? DIM-1 : (i + off_sets[z][0]) % DIM;
        int l = j + off_sets[z][1] < 0 ? DIM-1 : (j + off_sets[z][1]) % DIM;
        
        if( ((unsigned char*)board)[k*DIM + l] >= 128 ){
            count++;
        }
    }
    return count == 2 || count == 3 ? 1 : 0;    /* if 2 or 3 neighbour cells are alive, 
                                                 * the current cell has to become, or remain, alive; 
                                                 * otherwise it has to die 
                                                 */
}

/* Evolves the whole board once. 
 * The evolution is ordered, meaning by row and from the top left cell.
 */
void evolution_ordered(void* board, const int DIM, const int STEPS, const int maxval, const int SAVE){
    #pragma omp parallel for collapse(3)
    for(int s=0; s<STEPS; s++){
        for(int i=0; i<DIM; i++){
            for(int j=0; j<DIM; j++){
                if(check_neighbours(board, DIM, i, j) == 1){
                    *(((unsigned short int*)board) + i*DIM + j) = 255; /* board[i][j] = true, the cell is alive */
                }else{
                    *(((unsigned short int*)board) + i*DIM + j) = 0; /* board[i][j] = false, the cell is dead */
                }
            }
        }
        if(s % SAVE == 0){
            save_snap(board, DIM, maxval, s);
        }
    }
}

int top_block(const int rank, const int NDEC){
    return rank < NDEC ? (NDEC-1)*NDEC + rank : rank - NDEC;
}

int bottom_block(const int rank, const int NDEC){
    return rank >= NDEC*(NDEC-1) ? rank % NDEC : rank + NDEC;
}

int right_block(const int rank, const int NDEC){
    return rank%NDEC == NDEC-1 ? rank - NDEC + 1 : rank + 1;
}

int left_block(const int rank, const int NDEC){
    return rank%NDEC == 0 ? rank + NDEC - 1 : rank - 1;
}

/* Evolves the whole board once. 
 * The evolution is static, meaning the evaluation of the board is disentangled from the update. 
 */
void evolution_static(void* board, const int DIM, const int STEPS, const int maxval, const int SAVE){
    /* Options:
     * - save list of cells to modify;
     * - mark cells to be modified;
     * */
    
    int c_size, rank;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &c_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    /* partition grid into blocks and scatter them through the processes*/
    const int NDEC = 5; /* number of blocks in a coloumn in decomposition */
    const int BLOCKSIZE = DIM/NDEC; /* number of rows and columns in a block */

    if(c_size != NDEC*NDEC){
        fprintf(stderr,"Error: number of PEs %d != %d x %d\n", c_size, NDEC, NDEC);
        MPI_Finalize();
        exit(-1);
    }

    char block[BLOCKSIZE*BLOCKSIZE];
    for (int ii=0; ii<BLOCKSIZE*BLOCKSIZE; ii++) block[ii] = 0;

    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    MPI_Type_vector(BLOCKSIZE, BLOCKSIZE, DIM, MPI_CHAR, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(char), &blocktype);
    MPI_Type_commit(&blocktype);

    int disps[NDEC*NDEC];
    int counts[NDEC*NDEC];
    for (int ii=0; ii<NDEC; ii++) {
        for (int jj=0; jj<NDEC; jj++) {
            disps[ii*NDEC+jj] = ii*DIM*BLOCKSIZE+jj*BLOCKSIZE;
            counts [ii*NDEC+jj] = 1;
        }
    }

    MPI_Scatterv(board, counts, disps, blocktype, block, BLOCKSIZE*BLOCKSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);


    /* evolution */
    // #pragma omp parallel for collapse(3)
    for(int s=0; s<STEPS; s++){
        // TODO: block outer edge propagation
        MPI_status status;
        char bttm_row[BLOCKSIZE];
        char top_row[BLOCKSIZE];
        char left_clmn[BLOCKSIZE];
        char temp[BLOCKSIZE];
        for (int ii=0; ii<BLOCKSIZE; ii++) bttm_row[ii] = 0;
        for (int ii=0; ii<BLOCKSIZE; ii++) top_row[ii] = 0;
        for (int ii=0; ii<BLOCKSIZE; ii++) left_clmn[ii] = 0;

        /* send top row */
        if(rank == 0){
            // Send above
            MPI_Send(block, BLOCKSIZE, MPI_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD);
            // Receive from below
            MPI_Recv(bttm_row, BLOCKSIZE, MPI_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
        }else{
            // Receive from below
            MPI_Recv(bttm_row, BLOCKSIZE, MPI_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            // Send above
            MPI_Send(block, BLOCKSIZE, MPI_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD);
        }
        /* send bottom row */
        if(rank == 0){
            // Send below
            MPI_Sendv(block + (BLOCKSIZE*(BLOCKSIZE-1)), BLOCKSIZE, MPI_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD);
            // Receive from above
            MPI_Recv(bttm_row, BLOCKSIZE, MPI_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
        }else{
            // Receive from above
            MPI_Recv(bttm_row, BLOCKSIZE, MPI_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            // Send below
            MPI_Sendv(block + (BLOCKSIZE*(BLOCKSIZE-1)), BLOCKSIZE, MPI_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD);
        }
        /* send right coloumn */
        for(int i=1; i<=BLOCKSIZE; i++) temp[i-1] = block[i*BLOCKSIZE - 1];
        if(rank == 0){
            // Send right
            MPI_Send(temp, BLOCKSIZE, MPI_CHAR, right_block(rank, NDEC), 0, MPI_COMM_WORLD);
            // Receive from left
            MPI_Recv(left_clmn, BLOCKSIZE, MPI_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
        }else{
            MPI_Recv(left_clmn, BLOCKSIZE, MPI_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Send(temp, BLOCKSIZE, MPI_CHAR, right_block(rank, NDEC), 0, MPI_COMM_WORLD);
        }

        for(int i=0; i<BLOCKSIZE; i++){
            for(int j=0; j<BLOCKSIZE; j++){
                if(check_neighbours(board, BLOCKSIZE, i, j) == 1){ /* cell will be or remain alive */
                    if(*(((unsigned char*)block) + i*BLOCKSIZE + j) <= 127){  /* cell is currently dead */
                        *(((unsigned char*)block) + i*BLOCKSIZE + j) = 127;
                    }
                }else{ /* cell will be or remain dead */
                    if(*(((unsigned char*)block) + i*BLOCKSIZE + j) >= 128){ /* cell is currently alive */
                        *(((unsigned char*)block) + i*BLOCKSIZE + j) = 128;
                    }
                }
            }
        }
        
        for(int i=0; i<DIM; i++){
            for(int j=0; j<DIM; j++){
                if(*(((unsigned char*)board) + i*DIM + j) == 127){
                    *(((unsigned char*)board) + i*DIM + j) = 255;
                }else{
                    if(*(((unsigned char*)board) + i*DIM + j) == 128){
                        *(((unsigned char*)board) + i*DIM + j) = 0;
                    }
                }
            }
        }

        /* TODO: save board on one process */
        if(s % SAVE == 0){
            save_snap(board, DIM, maxval, s);
        }
    }

    MPI_Finalize();
}
