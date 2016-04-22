#ifndef _TREATMENT_CHECKBOARD_H_
#define _TREATMENT_CHECKBOARD_H_
/*
void compute_cell(int x, int y, int div4);
float *compute_seq(unsigned iterations);
float *compute_parallel(unsigned iterations);
float *compute_task(unsigned iterations);

void treatment_task(int x, int y);
*/


void cb_sand_init_homogeneous();
void cb_sand_init_center();

int cb_seq (int argc, char **argv);
int cb_parallel (int argc, char **argv);
int cb_parallel_task (int argc, char **argv);

#endif