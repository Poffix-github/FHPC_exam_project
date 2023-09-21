#include <gol.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

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
#define NPROWS 1
#define NPCOLS 1

/* Count the number of alive neighbours */
char check_neighbours(const void* board, const int BLROWS, const int BLCOLS, const int i, const int j, 
                      const unsigned char t_l, const unsigned char* top, const unsigned char t_r, const unsigned char* left, 
                      const unsigned char* right, const unsigned char b_l, const unsigned char* bottom, const unsigned char b_r, const int RANK){
    const int off_sets[NUM_NEIGHBOURS][2] = {{-1, -1}, {-1, 0}, {-1, 1},
                                             {0, -1},                       {0, 1},
                                             {1, -1},  {1, 0},  {1, 1}};
    unsigned int count=0;    /* number of alive cells around the ij-th */

    /* inner cell */
    if(i != 0 && j != 0 && i != BLROWS-1 && j != BLCOLS-1){
        for(int z=0; z<NUM_NEIGHBOURS; z++){
            if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
        }
    }else if(i == 0){     /* outer cell */
        if(j == 0){
            if(t_l >= 128) count++;
            if(top[0] >= 128) count++;
            if(top[1] >= 128) count++;
            if(left[0] >= 128) count++;
            if(((unsigned char*)board)[(i + off_sets[4][0])*BLCOLS + (j + off_sets[4][1])] >= 128) count++;
            if(left[1] >= 128) count++;
            for(int z=6; z<NUM_NEIGHBOURS; z++)
                if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
        }else if(j == BLCOLS-1){
            if(top[BLCOLS-2] >= 128) count++;
            if(top[BLCOLS-1] >= 128) count++;
            if(t_r >= 128) count++;
            if(((unsigned char*)board)[(i + off_sets[3][0])*BLCOLS + (j + off_sets[3][1])] >= 128) count++;
            if(right[0] >= 128) count++;
            for(int z=5; z<NUM_NEIGHBOURS-1; z++)
                if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
            if(right[1] >= 128) count++;
        }else{
            if(top[j-1] >= 128) count++;
            if(top[j] >= 128) count++;
            if(top[j+1] >= 128) count++;
            for(int z=3; z<NUM_NEIGHBOURS; z++)
                if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
        }
    }else if(i == BLROWS-1){
        if(j == 0){
            if(left[BLROWS-2] >= 128) count++;
            for(int z=1; z<3; z++)
                if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
            if(left[BLROWS-1] >= 128) count++;
            if(((unsigned char*)board)[(i + off_sets[4][0])*BLCOLS + (j + off_sets[4][1])] >= 128) count++;
            if(b_l >= 128) count++;
            if(bottom[0] >= 128) count++;
            if(bottom[1] >= 128) count++;
        }else if(j == BLCOLS-1){
            for(int z=0; z<2; z++)
                if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
            if(right[BLROWS-2] >= 128) count++;
            if(((unsigned char*)board)[(i + off_sets[3][0])*BLCOLS + (j + off_sets[3][1])] >= 128) count++;
            if(right[BLROWS-1] >= 128) count++;
            if(bottom[BLCOLS-2] >= 128) count++;
            if(bottom[BLCOLS-1] >= 128) count++;
            if(b_r >= 128) count++;
        }else{
            for(int z=0; z<5; z++)
                if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
            if(bottom[j-1] >= 128) count++;
            if(bottom[j] >= 128) count++;
            if(bottom[j+1] >= 128) count++;
        }
    }else if(j == 0){
        if(left[i-1] >= 128) count++;
        for(int z=1; z<3; z++)
            if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
    	if(left[i] >= 128) count++;
        if(((unsigned char*)board)[(i + off_sets[4][0])*BLCOLS + (j + off_sets[4][1])] >= 128) count++;
        if(left[i+1] >= 128) count++;
        for(int z=6; z<NUM_NEIGHBOURS; z++)
            if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
    }else{
        for(int z=0; z<2; z++)
            if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
        if(right[i-1] >= 128) count++;
        if(((unsigned char*)board)[(i + off_sets[3][0])*BLCOLS + (j + off_sets[3][1])] >= 128) count++;
        if(right[i] >= 128) count++;
        for(int z=5; z<7; z++)
            if( ((unsigned char*)board)[(i + off_sets[z][0])*BLCOLS + (j + off_sets[z][1])] >= 128 ) count++;
        if(right[i+1] >= 128) count++;
    }   

    // if(RANK == 0) printf("board[%d][%d] count = %d\n",i,j,count);

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
void evolution_ordered(void* board, const int DIM, const int STEPS, const int MAXVAL, const int SAVE){
    for(int s=0; s<STEPS; s++){
        // TODO: MPI
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
            save_snap(board, DIM, MAXVAL, s);
        }
    }
}

int top_block(const int RANK){
    return RANK < NPCOLS ? (NPROWS-1)*NPCOLS + RANK : RANK - NPCOLS;
}

int bottom_block(const int RANK){
    return RANK >= NPCOLS*(NPROWS-1) ? RANK % NPCOLS : RANK + NPCOLS;
}

int right_block(const int RANK){
    return RANK%NPCOLS == NPCOLS-1 ? RANK - NPCOLS + 1 : RANK + 1;
}

int left_block(const int RANK){
    return RANK%NPCOLS == 0 ? RANK + NPCOLS - 1 : RANK - 1;
}

int top_left_blk(const int RANK){
    return RANK%NPCOLS == 0 ? top_block(RANK) + NPCOLS-1 : top_block(RANK) - 1;
}

int top_right_blk(const int RANK){
    return RANK%NPCOLS == NPCOLS-1 ? top_block(RANK) - NPCOLS+1 : top_block(RANK) + 1;
}

int btm_left_blk(const int RANK){
    return RANK%NPCOLS == 0 ? bottom_block(RANK) + NPCOLS-1 : bottom_block(RANK) - 1;
}

int btm_right_blk(const int RANK){
    return RANK%NPCOLS == NPCOLS-1 ? bottom_block(RANK) - NPCOLS+1 : bottom_block(RANK) + 1;
}

/* Evolves the whole board once. 
 * The evolution is static, meaning the evaluation of the board is disentangled from the update. 
 */
void evolution_static(void* board, const int DIM, const int STEPS, const int MAXVAL, const int SAVE, const int NUM_PROC, const int RANK, int* evo_T, int* avg_propT){
    /* Options:
     * - save list of cells to modify;
     * - mark cells to be modified;
     * */
    
    /* partition grid into blocks and scatter them through the processes*/
    const int BLOCKROWS = DIM/NPROWS; /* number of rows and columns in a block */
    const int BLOCKCOLS = DIM/NPCOLS; /* number of rows and columns in a block */

    if(NUM_PROC != NPROWS*NPCOLS){
        fprintf(stderr,"Error: number of PEs %d != %d x %d\n", NUM_PROC, NPROWS, NPCOLS);
        MPI_Finalize();
        exit(-1);
    }

    unsigned char block[BLOCKROWS*BLOCKCOLS];
    #pragma omp parallel for schedule(static) shared(BLOCKROWS, BLOCKCOLS)
    for (int ii=0; ii<BLOCKROWS*BLOCKCOLS; ii++) block[ii] = 0;

    // if(RANK == 0){
    //     printf("block inizializzato\n");
    //     for (int ii=0; ii<BLOCKSIZE*BLOCKSIZE; ii++) printf("%d ",block[ii]);
    //     printf("\n");
    // }

    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    MPI_Type_vector(BLOCKROWS, BLOCKCOLS, DIM, MPI_UNSIGNED_CHAR, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(char), &blocktype);
    MPI_Type_commit(&blocktype);

    int disps[NPROWS*NPCOLS];
    int counts[NPROWS*NPCOLS];
    #pragma omp parallel for schedule(static) collapse(2) shared(BLOCKROWS, BLOCKCOLS, disps, counts)
    for (int ii=0; ii<NPROWS; ii++) {
        for (int jj=0; jj<NPCOLS; jj++) {
            disps[ii*NPCOLS+jj] = ii*DIM*BLOCKROWS+jj*BLOCKCOLS;
            counts [ii*NPCOLS+jj] = 1;
        }
    }

    // if(RANK == 0){
    //     printf("displacement initialized\n");
    //     for (int i=0; i<NPROWS*NPCOLS; i++) printf("%d ", disps[i]); /*0, 2, 8, 10*/
    //     printf("\n");

    //     printf("counts initialized\n");
    //     for (int i=0; i<NPROWS*NPCOLS; i++) printf("%d ", counts[i]);
    //     printf("\n");
    // }

    MPI_Scatterv(board, counts, disps, blocktype, block, BLOCKROWS*BLOCKCOLS, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    // printf("scatter done\n");

    MPI_Status status;
    unsigned char btm_row[BLOCKCOLS];
    unsigned char top_row[BLOCKCOLS];
    unsigned char left_clmn[BLOCKROWS];
    unsigned char right_clmn[BLOCKROWS];
    unsigned char temp[BLOCKROWS];
    unsigned char top_left, top_right, btm_left, btm_right;

    #pragma omp parallel
    {
        #pragma omp parallel for schedule(static) shared(BLOCKROWS, BLOCKCOLS, btm_row)
        for (int ii=0; ii<BLOCKCOLS; ii++) btm_row[ii] = 0;
        #pragma omp parallel for schedule(static) shared(BLOCKROWS, BLOCKCOLS, top_row)
        for (int ii=0; ii<BLOCKCOLS; ii++) top_row[ii] = 0;
        #pragma omp parallel for schedule(static) shared(BLOCKROWS, BLOCKCOLS, left_clmn)
        for (int ii=0; ii<BLOCKROWS; ii++) left_clmn[ii] = 0;
        #pragma omp parallel for schedule(static) shared(BLOCKROWS, BLOCKCOLS, right_clmn)
        for (int ii=0; ii<BLOCKROWS; ii++) right_clmn[ii] = 0;
    }

    // if(RANK == 0){
    //     printf("btm_row initialized\n");
    //     for (int i=0; i<BLOCKSIZE; i++) printf("%d ", btm_row[i]);
    //     printf("\n");

    //     printf("top_row initialized\n");
    //     for (int i=0; i<BLOCKSIZE; i++) printf("%d ", top_row[i]);
    //     printf("\n");

    //     printf("left_clmn initialized\n");
    //     for (int i=0; i<BLOCKSIZE; i++) printf("%d ", left_clmn[i]);
    //     printf("\n");
        
    //     printf("right_clmn initialized\n");
    //     for (int i=0; i<BLOCKSIZE; i++) printf("%d ", right_clmn[i]);
    //     printf("\n");   
    // }

    double tbegin, tend, tpstart, tpend;
    int tevo[2]={0,0}, tprop=0;
    if(RANK == 0) tbegin = MPI_Wtime();

    /* evolution */
    for(int s=0; s<STEPS; s++){

        if(RANK == 0) tpstart = MPI_Wtime();

        /* propagate rows */
        if((RANK/NPCOLS)%2 == 0){
            
            /* send top row */
            MPI_Send(block, BLOCKCOLS, MPI_UNSIGNED_CHAR, top_block(RANK), 0, MPI_COMM_WORLD);
            MPI_Recv(btm_row, BLOCKCOLS, MPI_UNSIGNED_CHAR, bottom_block(RANK), 0, MPI_COMM_WORLD, &status);
            /* send bottom row */
            MPI_Send(block + (BLOCKCOLS*(BLOCKROWS-1)), BLOCKCOLS, MPI_UNSIGNED_CHAR, bottom_block(RANK), 0, MPI_COMM_WORLD);
            MPI_Recv(top_row, BLOCKCOLS, MPI_UNSIGNED_CHAR, top_block(RANK), 0, MPI_COMM_WORLD, &status);
        }else{
            MPI_Recv(btm_row, BLOCKCOLS, MPI_UNSIGNED_CHAR, bottom_block(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Send(block, BLOCKCOLS, MPI_UNSIGNED_CHAR, top_block(RANK), 0, MPI_COMM_WORLD);

            MPI_Recv(top_row, BLOCKCOLS, MPI_UNSIGNED_CHAR, top_block(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Send(block + (BLOCKCOLS*(BLOCKROWS-1)), BLOCKCOLS, MPI_UNSIGNED_CHAR, bottom_block(RANK), 0, MPI_COMM_WORLD);
        }
        /* propagate coloumns */
        if((RANK%NPCOLS)%2 == 0){
            /* send right coloumn */
            for(int i=1; i<=BLOCKROWS; i++) temp[i-1] = block[i*BLOCKROWS - 1];
            MPI_Send(temp, BLOCKROWS, MPI_UNSIGNED_CHAR, right_block(RANK), 0, MPI_COMM_WORLD);
            MPI_Recv(left_clmn, BLOCKROWS, MPI_UNSIGNED_CHAR, left_block(RANK), 0, MPI_COMM_WORLD, &status);
            /* send left coloumn */
            for(int i=0; i<BLOCKROWS; i++) temp[i] = block[i*BLOCKROWS];
            MPI_Send(temp, BLOCKROWS, MPI_UNSIGNED_CHAR, left_block(RANK), 0, MPI_COMM_WORLD);
            MPI_Recv(right_clmn, BLOCKROWS, MPI_UNSIGNED_CHAR, left_block(RANK), 0, MPI_COMM_WORLD, &status);
        }else{
            MPI_Recv(left_clmn, BLOCKROWS, MPI_UNSIGNED_CHAR, left_block(RANK), 0, MPI_COMM_WORLD, &status);
            for(int i=1; i<=BLOCKROWS; i++) temp[i-1] = block[i*BLOCKROWS - 1];
            MPI_Send(temp, BLOCKROWS, MPI_UNSIGNED_CHAR, right_block(RANK), 0, MPI_COMM_WORLD);

            MPI_Recv(right_clmn, BLOCKROWS, MPI_UNSIGNED_CHAR, left_block(RANK), 0, MPI_COMM_WORLD, &status);
            for(int i=0; i<BLOCKROWS; i++) temp[i] = block[i*BLOCKROWS];
            MPI_Send(temp, BLOCKROWS, MPI_UNSIGNED_CHAR, left_block(RANK), 0, MPI_COMM_WORLD);
        }

        /* propagate corners in diagonal */
        if( (RANK/NPCOLS)%2 == 0 ){   /* blocks in even rows send first */
            MPI_Send(block, 1, MPI_UNSIGNED_CHAR, top_left_blk(RANK), 0, MPI_COMM_WORLD);
            MPI_Send(block + NPCOLS-1, 1, MPI_UNSIGNED_CHAR, top_right_blk(RANK), 0, MPI_COMM_WORLD);
            MPI_Send(block + NPCOLS*(NPROWS-1), 1, MPI_UNSIGNED_CHAR, btm_left_blk(RANK), 0, MPI_COMM_WORLD);
            MPI_Send(block + NPCOLS*NPROWS - 1, 1, MPI_UNSIGNED_CHAR, btm_right_blk(RANK), 0, MPI_COMM_WORLD);

            MPI_Recv(&btm_right, 1, MPI_UNSIGNED_CHAR, btm_right_blk(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&btm_left, 1, MPI_UNSIGNED_CHAR, btm_left_blk(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_right, 1, MPI_UNSIGNED_CHAR, top_right_blk(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_left, 1, MPI_UNSIGNED_CHAR, top_left_blk(RANK), 0, MPI_COMM_WORLD, &status);
        }else{
            MPI_Recv(&btm_right, 1, MPI_UNSIGNED_CHAR, btm_right_blk(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&btm_left, 1, MPI_UNSIGNED_CHAR, btm_left_blk(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_right, 1, MPI_UNSIGNED_CHAR, top_right_blk(RANK), 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&top_left, 1, MPI_UNSIGNED_CHAR, top_left_blk(RANK), 0, MPI_COMM_WORLD, &status);

            MPI_Send(block, 1, MPI_UNSIGNED_CHAR, top_left_blk(RANK), 0, MPI_COMM_WORLD);
            MPI_Send(block + NPCOLS-1, 1, MPI_UNSIGNED_CHAR, top_right_blk(RANK), 0, MPI_COMM_WORLD);
            MPI_Send(block + NPCOLS*(NPROWS-1), 1, MPI_UNSIGNED_CHAR, btm_left_blk(RANK), 0, MPI_COMM_WORLD);
            MPI_Send(block + NPCOLS*NPROWS - 1, 1, MPI_UNSIGNED_CHAR, btm_right_blk(RANK), 0, MPI_COMM_WORLD);
        }

        // printf("propagation done\n");

        // for(int p=0; p<NUM_PROC; p++){
            // if(RANK == 0){
            //     printf("RANK %d\n", RANK);
            //     printf(" top_1: %d   top_2: %d    top_3: %d    top_4: %d\n", top_left, top_row[0], top_row[1], top_right);
            //     printf("left_1: %d                           right_1: %d\n", left_clmn[0], right_clmn[0]);
            //     printf("left_2: %d                           right_2: %d\n", left_clmn[1], right_clmn[1]);
            //     printf(" btm_1: %d   btm_2: %d    btm_3: %d    btm_4: %d\n", btm_left, btm_row[0], btm_row[1], btm_right);
            // }
            // MPI_Barrier(MPI_COMM_WORLD);
        // }

        if(RANK == 0){
            tpend = MPI_Wtime();
            tprop += get_time(tpstart, tpend);
        }

        #pragma omp parallel shared(BLOCKROWS, BLOCKCOLS, block, top_left, top_row, top_right, left_clmn, right_clmn, btm_left, btm_row, btm_right)
        {
            // if(RANK == 0) printf("ciao sono il thread %d\n", omp_get_thread_num());
            // if(RANK == 0){
            //     for(int i=0; i<BLOCKSIZE; i++){
            //         for(int j=0; j<BLOCKSIZE; j++){
            //             printf("%d ", block[i*BLOCKSIZE+j]);
            //         }
            //         printf("\n");
            //     }
            // }

            #pragma omp for schedule(static) collapse(2)
            for(int i=0; i<BLOCKROWS; i++){
                for(int j=0; j<BLOCKCOLS; j++){
                    if(check_neighbours(block, BLOCKROWS, BLOCKCOLS, i, j, top_left, top_row, top_right, left_clmn, right_clmn, btm_left, btm_row, btm_right, RANK) == 1){ /* cell will be or remain alive */
                        if(*(((unsigned char*)block) + i*BLOCKCOLS + j) <= 127) *(((unsigned char*)block) + i*BLOCKCOLS + j) = 127;
                    }else{
                        if(*(((unsigned char*)block) + i*BLOCKCOLS + j) >= 128) *(((unsigned char*)block) + i*BLOCKCOLS + j) = 128;
                    }

                }
            }

            // if(RANK == 0){
            //     for(int i=0; i<BLOCKSIZE; i++){
            //         for(int j=0; j<BLOCKSIZE; j++){
            //             printf("%d ", block[i*BLOCKSIZE+j]);
            //         }
            //         printf("\n");
            //     }
            // }
            
            #pragma omp for schedule(static) collapse(2)
            for(int i=0; i<BLOCKROWS; i++){
                for(int j=0; j<BLOCKCOLS; j++){
                    if(*(((unsigned char*)block) + i*BLOCKCOLS + j) == 127) *(((unsigned char*)block) + i*BLOCKCOLS + j) = 255;
                    else if(*(((unsigned char*)block) + i*BLOCKCOLS + j) == 128) *(((unsigned char*)block) + i*BLOCKCOLS + j) = 0;
                }
            }

            // if(RANK == 0){
            //     for(int i=0; i<BLOCKSIZE; i++){
            //         for(int j=0; j<BLOCKSIZE; j++){
            //             printf("%d ", block[i*BLOCKSIZE+j]);
            //         }
            //         printf("\n");
            //     }
            // }
        }
        // printf("evolution done\n");
        

        if(s % SAVE == 0){
            MPI_Gatherv(block, BLOCKROWS*BLOCKCOLS, MPI_UNSIGNED_CHAR, board, counts, disps, blocktype, 0, MPI_COMM_WORLD); 
            
            if(RANK == 0) save_snap(board, DIM, MAXVAL, s);

            MPI_Scatterv(board, counts, disps, blocktype, block, BLOCKROWS*BLOCKCOLS, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        }
    }
    if(RANK == 0){
        tend = MPI_Wtime();
        *evo_T = get_time(tbegin, tend);
        *avg_propT = (int)tprop/STEPS;
    }
}
