#pragma once

/* evolution.c */

void evolution_ordered(void* board, const int dim);

void evolution_static(void* board, const int dim);


/* read_write_pgm_image.c */

void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);

void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);

void swap_image( void *image, int xsize, int ysize, int maxval );


/* utilities.c */

void get_args(int argc, char **argv, int *action, int *k, int *e, int *n, int *s, char **fname);

void *random_board(const int size, char *string);
