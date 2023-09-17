#include <gol.h>
#include <stdio.h>
#include <stdlib.h>
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
#define NUM_NEIGHBOURS 8

/* Count the number of alive neighbours 
 */
char check_neighbours(const void* board, const int DIM, const int i, const int j, 
                      const unsigned char t_l, const unsigned char* top, const unsigned char t_r, const unsigned char* left, 
                      const unsigned char* right, const unsigned char b_l, const unsigned char* bottom, const unsigned char b_r, int rank){
    const int off_sets[NUM_NEIGHBOURS][2] = {{-1, -1}, {-1, 0}, {-1, 1},
                                                  {0, -1},                       {0, 1},
                                                  {1, -1},  {1, 0},  {1, 1}};
    unsigned int count=0;    /* number of alive cells around the ij-th */

    // if(rank == 0){
    //     printf("neighbours of board[%d][%d] = %d\n", i, j, ((unsigned char*)board)[i*DIM + j]);
    // }

    /* inner cell */
    if(i != 0 && j != 0 && i != DIM-1 && j != DIM-1){
        // printf("rank %d, neighbours of (%d, %d):", rank, i, j); 
        for(int z=0; z<NUM_NEIGHBOURS; z++){
                // int k = i + off_sets[z][0];
                // int l = j + off_sets[z][1];
                // printf("[%d][%d] = %d\n", k, l, ((unsigned char*)board)[k*DIM + l]);
            if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
        }
    }else {     /* outer cell */
        if(i == 0){
            if(j == 0){
                count += t_l >= 128 ? 1 : 0;
                count += top[0] >= 128 ? 1 : 0;
                count += top[1] >= 128 ? 1 : 0;
                count += left[0] >= 128 ? 1 : 0;
                count += ((unsigned char*)board)[(i + off_sets[4][0])*DIM + (j + off_sets[4][1])] >= 128 ? 1 : 0;
                count += left[1] >= 128 ? 1 : 0;
                for(int z=6; z<NUM_NEIGHBOURS; z++)
                    if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
            }else{
                if(j == DIM-1){
                    count += top[DIM-2] >= 128 ? 1 : 0;
                    count += top[DIM-1] >= 128 ? 1 : 0;
                    count += t_r >= 128 ? 1 : 0;
                    count += ((unsigned char*)board)[(i + off_sets[3][0])*DIM + (j + off_sets[3][1])] >= 128 ? 1 : 0;
                    count += right[0] >= 128 ? 1 : 0;
                    for(int z=5; z<NUM_NEIGHBOURS-1; z++)
                        if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                    count += right[1] >= 128 ? 1 : 0;
                }else{
                    count += top[j-1] >= 128 ? 1 : 0;
                    count += top[j] >= 128 ? 1 : 0;
                    count += top[j+1] >= 128 ? 1 : 0;
                    for(int z=3; z<NUM_NEIGHBOURS; z++)
                        if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                }
            }
        }else{
            if(i == DIM-1){
                if(j == 0){
                    count += left[DIM-2] >= 128 ? 1 : 0;
                    for(int z=1; z<3; z++)
                        if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                    count += left[DIM-1] >= 128 ? 1 : 0;
                    count += ((unsigned char*)board)[(i + off_sets[4][0])*DIM + (j + off_sets[4][1])] >= 128 ? 1 : 0;
                    count += b_l >= 128 ? 1 : 0;
                    count += bottom[0] >= 128 ? 1 : 0;
                    count += bottom[1] >= 128 ? 1 : 0;
                }else{
                    if(j == DIM-1){
                        for(int z=0; z<2; z++)
                            if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                        count += right[DIM-2] >= 128 ? 1 : 0;
                        count += ((unsigned char*)board)[(i + off_sets[3][0])*DIM + (j + off_sets[3][1])] >= 128 ? 1 : 0;
                        count += right[DIM-1] >= 128 ? 1 : 0;
                        count += bottom[DIM-2] >= 128 ? 1 : 0;
                        count += bottom[DIM-1] >= 128 ? 1 : 0;
                        count += b_r >= 128 ? 1 : 0;
                    }else{
                        for(int z=0; z<5; z++)
                            if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                        count += bottom[j-1] >= 128 ? 1 : 0;
                        count += bottom[j] >= 128 ? 1 : 0;
                        count += bottom[j+1] >= 128 ? 1 : 0;
                    }
                }
            }else{
                if(j == 0){
                    count += left[i-1] >= 128 ? 1 : 0;
                    for(int z=1; z<3; z++)
                        if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                    count += left[i] >= 128 ? 1 : 0;
                    count += ((unsigned char*)board)[(i + off_sets[4][0])*DIM + (j + off_sets[4][1])] >= 128 ? 1 : 0;
                    count += left[i+1] >= 128 ? 1 : 0;
                    for(int z=6; z<NUM_NEIGHBOURS; z++)
                        if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                }else{
                    for(int z=0; z<2; z++)
                        if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                    count += right[i-1] >= 128 ? 1 : 0;
                    count += ((unsigned char*)board)[(i + off_sets[3][0])*DIM + (j + off_sets[3][1])] >= 128 ? 1 : 0;
                    count += right[i] >= 128 ? 1 : 0;
                    for(int z=5; z<7; z++)
                        if( ((unsigned char*)board)[(i + off_sets[z][0])*DIM + (j + off_sets[z][1])] >= 128 ) count++;
                    count += right[i+1] >= 128 ? 1 : 0;
                }
            }
        }
    }

    return count == 2 || count == 3 ? 1 : 0;    /* if 2 or 3 neighbour cells are alive, 
                                                 * the current cell has to become, or remain, alive; 
                                                 * otherwise it has to die 
                                                 */
}

char check_neighbours_ord(const void* board, const int dim, const int i, const int j){
    const int off_sets[NUM_NEIGHBOURS][2] = {  {-1, -1}, {-1, 0}, {-1, 1},
                                                     {0, -1},           {0, 1},
                                                     {1, -1},  {1, 0},  {1, 1}  };
    int count=0;    /* number of alive cells around the ij-th */

    for(int z=0; z<NUM_NEIGHBOURS; z++){
        /* compute coordinates couple */
        int k = i + off_sets[z][0] < 0 ? dim-1 : (i + off_sets[z][0]) % dim;
        int l = j + off_sets[z][1] < 0 ? dim-1 : (j + off_sets[z][1]) % dim;

        if( ((unsigned char*)board)[k*dim + l] >= 128 ){
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
                if(check_neighbours_ord(board, DIM, i, j) == 1){
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

int top_left_blk(const int rank, const int NDEC){
    return rank%NDEC == 0 ? top_block(rank, NDEC) + NDEC-1 : top_block(rank, NDEC) - 1;
}

int top_right_blk(const int rank, const int NDEC){
    return rank%NDEC == NDEC-1 ? top_block(rank, NDEC) - NDEC+1 : top_block(rank, NDEC) + 1;
}

int btm_left_blk(const int rank, const int NDEC){
    return rank%NDEC == 0 ? bottom_block(rank, NDEC) + NDEC-1 : bottom_block(rank, NDEC) - 1;
}

int btm_right_blk(const int rank, const int NDEC){
    return rank%NDEC == NDEC-1 ? bottom_block(rank, NDEC) - NDEC+1 : bottom_block(rank, NDEC) + 1;
}

/* Evolves the whole board once. 
 * The evolution is static, meaning the evaluation of the board is disentangled from the update. 
 */
void evolution_static(void* board, const int DIM, const int STEPS, const int maxval, const int SAVE, const int num_proc, const int rank){
    /* Options:
     * - save list of cells to modify;
     * - mark cells to be modified;
     * */
    
    

    /* partition grid into blocks and scatter them through the processes*/
    const int NDEC = 2; /* number of blocks in a coloumn in decomposition */
    const int BLOCKSIZE = DIM/NDEC; /* number of rows and columns in a block */

    if(num_proc != NDEC*NDEC){
        fprintf(stderr,"Error: number of PEs %d != %d x %d\n", num_proc, NDEC, NDEC);
        MPI_Finalize();
        exit(-1);
    }

    unsigned char block[BLOCKSIZE*BLOCKSIZE];
    for (int ii=0; ii<BLOCKSIZE*BLOCKSIZE; ii++) block[ii] = 0;

    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    MPI_Type_vector(BLOCKSIZE, BLOCKSIZE, DIM, MPI_UNSIGNED_CHAR, &blocktype2);
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

    MPI_Scatterv(board, counts, disps, blocktype, block, BLOCKSIZE*BLOCKSIZE, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    MPI_Status status;
    unsigned char btm_row[BLOCKSIZE];
    unsigned char top_row[BLOCKSIZE];
    unsigned char left_clmn[BLOCKSIZE];
    unsigned char right_clmn[BLOCKSIZE];
    unsigned char temp[BLOCKSIZE];
    unsigned char top_left, top_right, btm_left, btm_right;

    for (int ii=0; ii<BLOCKSIZE; ii++) btm_row[ii] = 0;
    for (int ii=0; ii<BLOCKSIZE; ii++) top_row[ii] = 0;
    for (int ii=0; ii<BLOCKSIZE; ii++) left_clmn[ii] = 0;
    for (int ii=0; ii<BLOCKSIZE; ii++) right_clmn[ii] = 0;

    for (int proc=0; proc<num_proc; proc++) {
        if (proc == rank) {
            printf("Rank = %d\n", rank);
            if (rank == 0) {
                printf("Global matrix: \n");
                for (int ii=0; ii<DIM; ii++) {
                    for (int jj=0; jj<DIM; jj++) {
                        printf("%3d ",((unsigned char*)board)[ii*DIM+jj]);
                    }
                    printf("\n");
                }
            }
            printf("Local Matrix:\n");
            for (int ii=0; ii<BLOCKSIZE; ii++) {
                for (int jj=0; jj<BLOCKSIZE; jj++) {
                    printf("%3d ", block[ii*BLOCKSIZE+jj]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* evolution */
    // #pragma omp parallel for collapse(3)
    for(int s=0; s<STEPS; s++){

        for(int p=0; p<num_proc; p++){
            if(rank == p){
                printf("rank %d\ntop row: ", rank);
                for(int i=0; i<BLOCKSIZE; i++){
                    printf("%d ", block[i]);
                }
                printf("\nbottom row: ");
                unsigned char *last_r = block + (BLOCKSIZE*(BLOCKSIZE-1));
                for(int i=0; i<BLOCKSIZE; i++){
                    printf("%d ", last_r[i]);
                }
                printf("\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        /* propagate rows */
        if((rank/NDEC)%2 == 0){

            // printf("rank: %d,    top: %d,    bottom: %d\n", rank, top_block(rank, NDEC), bottom_block(rank, NDEC));
            
            /* send top row */
            MPI_Send(block, BLOCKSIZE, MPI_UNSIGNED_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Recv(btm_row, BLOCKSIZE, MPI_UNSIGNED_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            /* send bottom row */
            MPI_Send(block + (BLOCKSIZE*(BLOCKSIZE-1)), BLOCKSIZE, MPI_UNSIGNED_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Recv(btm_row, BLOCKSIZE, MPI_UNSIGNED_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
        }else{
            // printf("rank: %d,    top: %d,    bottom: %d\n", rank, top_block(rank, NDEC), bottom_block(rank, NDEC));

            MPI_Recv(btm_row, BLOCKSIZE, MPI_UNSIGNED_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Send(block, BLOCKSIZE, MPI_UNSIGNED_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD);

            MPI_Recv(btm_row, BLOCKSIZE, MPI_UNSIGNED_CHAR, top_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Send(block + (BLOCKSIZE*(BLOCKSIZE-1)), BLOCKSIZE, MPI_UNSIGNED_CHAR, bottom_block(rank, NDEC), 0, MPI_COMM_WORLD);
        }
        /* propagate coloumns */
        if((rank%NDEC)%2 == 0){
            /* send right coloumn */
            for(int i=1; i<=BLOCKSIZE; i++) temp[i-1] = block[i*BLOCKSIZE - 1];
            MPI_Send(temp, BLOCKSIZE, MPI_UNSIGNED_CHAR, right_block(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Recv(left_clmn, BLOCKSIZE, MPI_UNSIGNED_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            /* send left coloumn */
            for(int i=0; i<BLOCKSIZE; i++) temp[i] = block[i*BLOCKSIZE];
            MPI_Send(temp, BLOCKSIZE, MPI_UNSIGNED_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Recv(right_clmn, BLOCKSIZE, MPI_UNSIGNED_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
        }else{
            MPI_Recv(left_clmn, BLOCKSIZE, MPI_UNSIGNED_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            for(int i=1; i<=BLOCKSIZE; i++) temp[i-1] = block[i*BLOCKSIZE - 1];
            MPI_Send(temp, BLOCKSIZE, MPI_UNSIGNED_CHAR, right_block(rank, NDEC), 0, MPI_COMM_WORLD);

            MPI_Recv(right_clmn, BLOCKSIZE, MPI_UNSIGNED_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            for(int i=0; i<BLOCKSIZE; i++) temp[i] = block[i*BLOCKSIZE];
            MPI_Send(temp, BLOCKSIZE, MPI_UNSIGNED_CHAR, left_block(rank, NDEC), 0, MPI_COMM_WORLD);
        }

        /* propagate corners in diagonal */
        if( (rank/NDEC)%2 == 0 ){   /* blocks in even rows send first */
            MPI_Send(block, 1, MPI_UNSIGNED_CHAR, top_left_blk(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Send(block + NDEC-1, 1, MPI_UNSIGNED_CHAR, top_right_blk(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Send(block + NDEC*(NDEC-1), 1, MPI_UNSIGNED_CHAR, btm_left_blk(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Send(block + NDEC*NDEC - 1, 1, MPI_UNSIGNED_CHAR, btm_right_blk(rank, NDEC), 0, MPI_COMM_WORLD);

            MPI_Recv(&btm_right, 1, MPI_UNSIGNED_CHAR, btm_right_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&btm_left, 1, MPI_UNSIGNED_CHAR, btm_left_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_right, 1, MPI_UNSIGNED_CHAR, top_right_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_left, 1, MPI_UNSIGNED_CHAR, top_left_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);
        }else{
            MPI_Recv(&btm_right, 1, MPI_UNSIGNED_CHAR, btm_right_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&btm_left, 1, MPI_UNSIGNED_CHAR, btm_left_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_right, 1, MPI_UNSIGNED_CHAR, top_right_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_left, 1, MPI_UNSIGNED_CHAR, top_left_blk(rank, NDEC), 0, MPI_COMM_WORLD, &status);

            MPI_Send(block, 1, MPI_UNSIGNED_CHAR, top_left_blk(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Send(block + NDEC-1, 1, MPI_UNSIGNED_CHAR, top_right_blk(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Send(block + NDEC*(NDEC-1), 1, MPI_UNSIGNED_CHAR, btm_left_blk(rank, NDEC), 0, MPI_COMM_WORLD);
            MPI_Send(block + NDEC*NDEC - 1, 1, MPI_UNSIGNED_CHAR, btm_right_blk(rank, NDEC), 0, MPI_COMM_WORLD);
        }

        // if(rank == 2){
        //     printf(" top_1: %d   top_2: %d    top_3: %d    top_4: %d\n", top_left, top_row[0], top_row[1], top_right);
        //     printf("left_1: %d                           right_1: %d\n", left_clmn[0], right_clmn[0]);
        //     printf("left_2: %d                           right_2: %d\n", left_clmn[1], right_clmn[1]);
        //     printf(" btm_1: %d   btm_2: %d    btm_3: %d    btm_4: %d\n", btm_left, btm_row[0], btm_row[1], btm_right);
        // }

	    // if(rank == 3){
        //     printf(" top_1: %d   top_2: %d    top_3: %d    top_4: %d\n", top_left, top_row[0], top_row[1], top_right);
        //     printf("left_1: %d                           right_1: %d\n", left_clmn[0], right_clmn[0]);
        //     printf("left_2: %d                           right_2: %d\n", left_clmn[1], right_clmn[1]);
        //     printf(" btm_1: %d   btm_2: %d    btm_3: %d    btm_4: %d\n", btm_left, btm_row[0], btm_row[1], btm_right);
        // }


        // if(rank == 0){
        //     printf("before fisrt pass: \n");
        //     print_board_minimal(block, BLOCKSIZE);
        // }

        for(int i=0; i<BLOCKSIZE; i++){
            for(int j=0; j<BLOCKSIZE; j++){
                if(check_neighbours(block, BLOCKSIZE, i, j, top_left, top_row, top_right, left_clmn, right_clmn, btm_left, btm_row, btm_right, rank) == 1){ /* cell will be or remain alive */
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

        // if(rank == 0){
        //     printf("after fisrt pass: \n");
        //     print_board_minimal(block, BLOCKSIZE);
        // }
        
        for(int i=0; i<BLOCKSIZE; i++){
            for(int j=0; j<BLOCKSIZE; j++){
                if(*(((unsigned char*)block) + i*BLOCKSIZE + j) == 127){
                    *(((unsigned char*)block) + i*BLOCKSIZE + j) = 255;
                }else{
                    if(*(((unsigned char*)block) + i*BLOCKSIZE + j) == 128){
                        *(((unsigned char*)block) + i*BLOCKSIZE + j) = 0;
                    }
                }
            }
        }

        // if(rank == 0){
        //     printf("after second pass: \n");
        //     print_board_minimal(block, BLOCKSIZE);
        // }
        
        if(s % SAVE == 0){
            MPI_Gatherv(block, BLOCKSIZE*BLOCKSIZE, MPI_UNSIGNED_CHAR, board, counts, disps, blocktype, 0, MPI_COMM_WORLD);
            
            // if(rank == 0){
            //     printf("after gather: \n");
            //     print_board_minimal(block, BLOCKSIZE);
            // }   
            
            if(rank == 0) save_snap(board, DIM, maxval, s);

            MPI_Scatterv(board, counts, disps, blocktype, block, BLOCKSIZE*BLOCKSIZE, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        }
    }
}
