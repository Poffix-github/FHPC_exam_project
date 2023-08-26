#include "include/gol.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



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
    get_args(argc, argv, &action, &k, &e, &n, &s, &fname);
    
    if( action == INIT ){
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
        
        printf("reading image\n");
        read_pgm_image( &board, &maxval, &size, &size, fname);
        
        switch( maxval ){
            case -1:
                printf("I/O error in header\n"); break;
            case -2:
                printf("Memory not sufficient\n"); break;
            case -3:
                printf("I/O error in body\n"); break;
            default:
                /* swap endianism */
                if ( LITTLE_ENDIAN ){
                    swap_image( board, size, size, maxval);
                }
                
                printf("start evolution\n");
                for(int i=0; i<n; i++){
                    if (i % 100 == 0){
                    }
                    if( e == ORDERED ){
                        evolution_ordered(board, size);
                    }else{
                        printf("evolution step: %d\n", i);
                        evolution_static(board, size);
                    }
                    
                    if(i % s == 0){
                        /* swap endianism */
                        if ( LITTLE_ENDIAN ){
                            swap_image( board, size, size, maxval);
                        }
                        
                        char *snapshot_name;
                        size_t name_size = 29*sizeof(char);
                        snapshot_name = (char *)malloc(name_size);
                        snprintf(snapshot_name, name_size, "snapshots/snapshot_%05d.pgm", i);
                        
                        write_pgm_image( board, maxval, size, size, snapshot_name);
                    }
                }
        }
        
        free(board);
    }
    
    if ( fname != NULL ){
        free ( fname );
    }
}
