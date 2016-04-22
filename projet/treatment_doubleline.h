#ifndef _TREATMENT_DOUBLELINE_H_
#define _TREATMENT_DOUBLELINE_H_

void dl_sand_init_homogeneous();
void dl_sand_init_center();

int dl_seq (int argc, char **argv);
int dl_parallel (int argc, char **argv);
int dl_parallel_task (int argc, char **argv);

#endif