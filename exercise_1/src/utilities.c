#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <gol.h>
#include <omp.h>


#define INIT 1
#define RUN  2

char fname_deflt[] = "game_of_life.pgm";

void get_args(int argc, char **argv, int *action, int *k, int *e, int *n, int *s, char **fname){
/* parser for program parameters
 */
    *action = 0;
    char *optstring = "irk:e:f:n:s:";

    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
      
        case 'i':
            *action = INIT; break;
      
        case 'r':
            *action = RUN; break;
      
        case 'k':
            *k = atoi(optarg) < 100 ? 100 : atoi(optarg); break;
        case 'e':
            *e = atoi(optarg); break;

        case 'f':
            *fname = (char*)malloc( sizeof(optarg)+1 );
            sprintf(*fname, "%s", optarg );
            break;

        case 'n':
            *n = atoi(optarg); break;

        case 's':
            *s = atoi(optarg); break;

        default :
            printf("argument -%c not known\n", c ); break;
        }
    }
    
    if( *fname == NULL ){
        *fname = (char*)malloc( sizeof(fname_deflt)+1 );
        sprintf(*fname, "%s", fname_deflt );
    }
}

unsigned int easy_seed(char *string){
/* definse a sedd using a string,
 * it reads each character as an integer and sums them.
 */
    unsigned int seed=0;
    
    for(int i=(strlen(string)-1 > 10 ? 10 : strlen(string)-1); i>=0; i--){
        seed += string[i];
    }
    return seed;
}

void *random_board(const int size, char *string){
/* builds a board and fills it with random values.
 * There is a 0.1% chance that a cell is alive.
 */
    const short int   p = 200;
    char              *board;
    
    board = (char*)calloc( size*size, sizeof(char) );
    
    srand(easy_seed(string));
    #pragma omp parallel shared(board, size)
    {
        const int THREAD_NUM = omp_get_thread_num();
        srand(easy_seed(string)*THREAD_NUM);
        #pragma for schedule(static) collapse(2)
        for(int i=0; i<size; i++){
            for(int j=0; j<size; j++){
                board[i*size +j] = (unsigned char) (rand() % 1001 > p ? 0 : 255); 
                /* board[i][j] has p/10 % chance to be alive */
            }
        }
    }
    return (void*) board;
}

void save_snap(void* board, const int size, const int maxval, const int num){
    /* swap endianism */
    if ( LITTLE_ENDIAN ) swap_image( board, size, size, maxval);
                        
    char *snapshot_name;
    size_t name_size = 29*sizeof(char);
    snapshot_name = (char *)malloc(name_size);
    snprintf(snapshot_name, name_size, "snapshots/snapshot_%05d.pgm", num);
                        
    write_pgm_image( board, maxval, size, size, snapshot_name);
                        
    /* swap endianism */
    if ( LITTLE_ENDIAN ) swap_image( board, size, size, maxval);
}

int get_time(double tbegin, double tend){
    return (int)((tend - tbegin)*1000000.0);
}
