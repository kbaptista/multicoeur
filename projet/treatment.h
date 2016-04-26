#ifndef _TREATMENT_H_
#define _TREATMENT_H_

void cb_sand_init_homogeneous();
void cb_sand_init_center();

int cb_seq (int argc, char **argv);
int cb_parallel (int argc, char **argv);
int cb_parallel_task (int argc, char **argv);

int dl_seq (int argc, char **argv);
int dl_parallel (int argc, char **argv);
int dl_parallel_task (int argc, char **argv);

#endif