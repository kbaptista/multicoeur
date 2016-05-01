#ifndef _SAND_H_
#define _SAND_H_

unsigned get (unsigned x, unsigned y);
static void print(int table);
static void sand_init_homogeneous();
static void sand_init_center();
void sand_init (int init);
static void copy(int table);

void treatment_compare(int argc, char **argv, int init_compare_res, int size_compare_res);
void size_compare(int argc, char **argv, int init_compare_res);
void init_compare(int argc, char **argv);

#endif