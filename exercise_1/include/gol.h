#pragma once

/* evolution_parallel.c */

void evolution_ordered(void* board, const int dim, const int steps, const int maxval, const int save);

void evolution_static(void* board, const int DIM, const int STEPS, const int maxval, const int SAVE, const int num_proc, const int rank);


/* read_write_pgm_image_parallel.c */

void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);

void swap_image( void *image, int xsize, int ysize, int maxval );


/* utilities.c */

void get_args(int argc, char **argv, int *action, int *k, int *e, int *n, int *s, char **fname);

void *random_board(const int size, char *string);

void save_snap(void* board, const int size, const int maxval, const int num);
