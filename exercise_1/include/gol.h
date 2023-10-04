#pragma once

/* evolution_parallel.c */

void evolution_ordered(void* board, const int DIM, const int STEPS, const int MAXVAL, const int SAVE, const int NUM_PROC, const int RANK, int* evo_T, int* avg_propT);

void evolution_static(void* board, const int DIM, const int STEPS, const int MAXVAL, const int SAVE, const int NUM_PROC, const int RANK, int* evo_T, int* avg_propT);


/* read_write_pgm_image_parallel.c */

void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);

void swap_image( void *image, int xsize, int ysize, int maxval );


/* utilities.c */

void get_args(int argc, char **argv, int *action, int *k, int *e, int *n, int *s, char **fname);

void *random_board(const int size, char *string);

void save_snap(void* board, const int size, const int maxval, const int num);

int get_time(double tbegin, double tend);
