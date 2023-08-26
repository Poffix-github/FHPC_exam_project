#include <stdio.h>

int check_neighbours(const int* board, const int dim, const int i, const int j){
    /* count the number of alive neighbours */
    const int off_sets[9] = {-1,0,1};
    int count=0;    /* number of alive cells around the ij-th */
    for(int z=0; z<3; z++){
        /* compute couple of coordinates */
        const int k = i + off_sets[z] < 0 ? dim-1 : (i + off_sets[z])%dim;
        for(int x=0; x<3; x++){
            if(k != 0 || x != 1){   /* if both conditions are false we are in the current cell, 
                                     * it is not a neighbour of itself so we skip the iteration 
                                     */
                const int l = j + off_sets[x] < 0 ? dim-1 : (j + off_sets[x])%dim;
                if(*(board + k*dim + l) >= 1){  /* essential not to treat integer as boolean */
                    count++;
                }
            }
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
void evolution_ordered(int* board, const int dim){
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(check_neighbours(board, dim, i, j) == 1){
                *(board + i*dim + j) = 1; /* board[i][j] = true, the cell is alive */
            }else{
                *(board + i*dim + j) = 0; /* board[i][j] = false, the cell is dead */
            }
        }
    }
}

/* Evolves the whole board once. 
 * The evolution is static, meaning the evaluation of the board is disentangled from the update. 
 */
void evolution_static(int* board, const int dim){
    /* Options:
     * - save list of cells to modify;
     * - mark cells to be modified;
     * */
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(check_neighbours(board, dim, i, j) == 1){ /* cell will be or remain alive */
                if(*(board + i*dim + j) <= 0){  /* cell is currently dead */
                    *(board + i*dim + j) = -1;
                }
            }else{ /* cell will be or remain dead */
                if(*(board + i*dim + j) >= 1){ /* cell is currently alive */
                    *(board + i*dim + j) = 2;
                }
            }
        }
    }
    
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(*(board + i*dim + j) == -1){
                *(board + i*dim + j) = 1;
            }else{
                if(*(board + i*dim + j) == 2){
                    *(board + i*dim + j) = 0;
                }
            }
        }
    }
    /* CELL NOMENCLATURE
     * 0:   dead
     * 1:   alive
     *-1:   is dead and will live next update
     * 2:   is alive and will die next update
     * 
     * NOTE: to check if a cell is alive or dead at the current state, 
     * it is necessary to check if the cell is >=1 or <=0 respectively.
     */
}

/* Performs a pretty printing of the given board. Specifies cell coordinates and 
 * only uses Ascii characters. 
 */
void print_board(const int* board, const int dim){
    for(int i=-1; i<dim; i++){
        for(int j=-1; j<dim; j++){
            /* coordinates */
            if(i == -1 && j == -1){ /* top left most cell */
                printf("  │");
            }else{
                if(i == -1){    /* coordinates of columns */
                    printf(" %d│", j);
                }else{
                    if (j == -1){   /* coordinates of rows */
                        printf("%d │", i);
                    }else{
                        /* cells status */
                        if(*(board + i*dim + j)){
                            printf("██│");
                        }else{
                            printf("  │");
                        }
                    }
                }
            }
        }
        printf("\n");
        /* grid, horizontal lines */
        for(int j=-1; j<dim; j++){
            printf("──┼");
        }
        printf("\n");
    }
}

/* Prints every cell of the given board in the most light way possible. */
void print_board_minimal(const int* board, const int dim){
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            printf("%d", *(board + i*dim + j)); /* just prints the content of the board */
        }
        printf("\n");
    }
}

int main(){
    int board[25] = {0,0,0,0,0,
                     0,0,0,0,0,
                     0,1,1,0,0,
                     0,0,1,0,0,
                     0,0,0,0,0};
    const int dim = 5; 
    printf("initial state\n");
    print_board(board, dim);
    for(int i=0; i<10; i++){
        printf("step %d\n", i+1);
        evolution_static(board, dim);
        print_board(board, dim);
    }
    return 0;
} 
