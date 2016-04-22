#ifndef _TREATMENTS_H_
#define _TREATMENTS_H_

void compute_cell(int x, int y, int div4);
float *compute_seq(unsigned iterations);
float *compute_parallel(unsigned iterations);
float *compute_task(unsigned iterations);

void treatment_task(int x, int y);

int seq (int argc, char **argv);
int parallel (int argc, char **argv);
int parallel_task (int argc, char **argv);

#endif