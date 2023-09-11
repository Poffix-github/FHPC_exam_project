#include <gol.h>
#include <stdio.h>
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

#define MAX_NUMBER_OF_CELLS 8


/* Count the number of alive neighbours 
 */
char check_neighbours(const void* board, const int dim, const int i, const int j){
    const int off_sets[MAX_NUMBER_OF_CELLS][2] = {  {-1, -1}, {-1, 0}, {-1, 1},
                                                     {0, -1},           {0, 1},
                                                     {1, -1},  {1, 0},  {1, 1}  };
    int count=0;    /* number of alive cells around the ij-th */
    
    for(int z=0; z<MAX_NUMBER_OF_CELLS; z++){
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
void evolution_ordered(void* board, const int dim){
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(check_neighbours(board, dim, i, j) == 1){
                *(((unsigned short int*)board) + i*dim + j) = 255; /* board[i][j] = true, the cell is alive */
            }else{
                *(((unsigned short int*)board) + i*dim + j) = 0; /* board[i][j] = false, the cell is dead */
            }
        }
    }
}

/* Evolves the whole board once. 
 * The evolution is static, meaning the evaluation of the board is disentangled from the update. 
 */
void evolution_static(void* board, const int dim){
    /* Options:
     * - save list of cells to modify;
     * - mark cells to be modified;
     * */
    
    #pragma omp parallel for collapse(2)
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(check_neighbours(board, dim, i, j) == 1){ /* cell will be or remain alive */
                if(*(((unsigned char*)board) + i*dim + j) <= 127){  /* cell is currently dead */
                    *(((unsigned char*)board) + i*dim + j) = 127;
                }
            }else{ /* cell will be or remain dead */
                if(*(((unsigned char*)board) + i*dim + j) >= 128){ /* cell is currently alive */
                    *(((unsigned char*)board) + i*dim + j) = 128;
                }
            }
        }
    }
    
    #pragma omp parallel collapse(2)
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(*(((unsigned char*)board) + i*dim + j) == 127){
                *(((unsigned char*)board) + i*dim + j) = 255;
            }else{
                if(*(((unsigned char*)board) + i*dim + j) == 128){
                    *(((unsigned char*)board) + i*dim + j) = 0;
                }
            }
        }
    }
}
