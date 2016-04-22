#include "display.h"
#include "treatment_checkboard.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


#define DIM 128
#define MAX_HEIGHT  4

#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned ocean[DIM][DIM][2];

// indice indiquant laquelle des deux tables est actuellement en lecture.
unsigned table; //init 0 implicite

//booléen qui change après le dernier traitement, indique donc si compute a encore du traitement à faire ou non.
static unsigned is_end = 1;

// ------------------------------------------------------------------------------
// -------------               Fonctions Utilitaires                -------------
// ------------------------------------------------------------------------------

/** 
  * Print matrices time after time on the shell.
  * Is adapted to print values until 999 without deformations.
*/
static void print(int table)
{
    for (int y = 0; y < DIM; y++)
    {
      for (int x = 0; x < DIM; x++) 
      {
        if(ocean[y][x][table] < 10)
        {
          printf("  %d ",ocean[y][x][table]);
        }
        else if (ocean[y][x][table] < 100)
        {
          printf(" %d ",ocean[y][x][table]);
        }
        else
        {
          printf("%d ",ocean[y][x][table]);
        }
      }
      printf("\n");
    }
    printf("\n");
}

/** 
  * initialise matrices with 5 grain in every case.
*/
void cb_sand_init_homogeneous()
{
//#pragma omp parallel for schedule(static) collapse(2)

  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][0] = 5;
      ocean[y][x][1] = 5;
    }
  }
}

/** 
  * initialise matrices with no one grain except in the middle case which have 100000 grain.
*/
void cb_sand_init_center()
{
  int center_value = 100000;
//#pragma omp parallel for schedule(static) collapse(2)
  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][0] = 0;
      ocean[y][x][1] = 0;
    }
  }
  ocean[DIM/2][DIM/2][0] = center_value;
  ocean[DIM/2][DIM/2][1] = center_value;
}


// ------------------------------------------------------------------------------
// -------------              Fonctions de traitement               -------------
// ------------------------------------------------------------------------------

// callback
unsigned cb_get (unsigned x, unsigned y)
{
  return ocean[y][x][table];
}

/** 
  * Copy all values from table to 1 - table
*/
static void copy(int table){
  #pragma omp parallel for collapse(2)
  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][1-table] = ocean[y][x][table];
    }
  }
}

/** 
  * Applied modifications into 1-table for ocean[x][y] sandpile.
*/
static void compute_cell(int x, int y, int div4){
  ocean[x][y][1-table] -= div4*4;
  ocean[x-1][y][1-table] += div4;
  ocean[x+1][y][1-table] += div4;
  ocean[x][y-1][1-table] += div4;
  ocean[x][y+1][1-table] += div4;
  is_end = 1;
}

/** 
  * Compute fonction for sequentiel treatment
*/
static float *cb_compute_seq(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x++){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
          int div4 = ocean[x][y][table]/4;
          compute_cell(x,y,div4);
        }
      } 
    }
    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}


/** 
  * Compute fonction for sequentiel treatment which alternate at every line
*/
static float *cb_compute_seq_alternance(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x=x+2){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
          int div4 = ocean[x][y][table]/4;
          compute_cell(x,y,div4);
        }
      } 
    }
    for (int x = 2; x < DIM-1; x=x+2){
      for (int y = DIM-2; y > 0; y--){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
          int div4 = ocean[x][y][table]/4;
          compute_cell(x,y,div4);
        }
      }
    }
    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for collapsed treatment
*/
float *cb_compute_parallel(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

/*
  for (unsigned i = 0; i < iterations; i++){
    #pragma omp parallel for collapse(2)
    for (int x = 1; x < DIM-1; x++){
      for (int y = 1; y < DIM-1; y++){
        //printf("%d et %d\n",omp_get_thread_num(),omp_get_num_threads());
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
            int div4 = ocean[y][x][table]/4;
            #pragma omp critical
            compute_cell(y, x, div4);
        }
      } 
    }
    */

  for (unsigned i = 0; i < iterations; i++){
    #pragma omp parallel for collapse(2)
    for (int x = 1; x < DIM-1; x=x+2){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
          int div4 = ocean[x][y][table]/4;
            #pragma omp critical
          compute_cell(x,y,div4);
        }
      } 
    }
    #pragma omp parallel for collapse(2)
    for (int x = 2; x < DIM-1; x=x+2){
      for (int y = DIM-2; y > 0; y--){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
          int div4 = ocean[x][y][table]/4;
            #pragma omp critical
          compute_cell(x,y,div4);
        }
      }
    }

    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Adapted treatment to one case, with task specificity
*/
static void treatment_task(int x, int y){
  int div4 = ocean[x][y][table]/4;

  //#pragma omp critical
  {
    compute_cell(x, y, div4);
  }
}

/** 
  * Compute fonction for task treatment
*/
static float *cb_compute_task(unsigned iterations)
{
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  // 1 passage 
/*
  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x ++){
      for (int y = 1; y < DIM-1; y++){
        #pragma omp task depend(in:ocean[x][y][table], ocean[x][y-1][table], ocean[x][y+1][table], ocean[x-1][y][table], ocean[x+1][y][table]) depend(out:ocean[x][y][table], ocean[x][y-1][table], ocean[x][y+1][table], ocean[x-1][y][table], ocean[x+1][y][table])
        if(ocean[x][y][table] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }
*/

    // 2 passage (pour optimiser le cache ?)
  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x =x+2){
      for (int y = 1; y < DIM-1; y++){
        #pragma omp task depend(in:ocean[x][y][table], ocean[x][y-1][table], ocean[x][y+1][table], ocean[x-1][y][table], ocean[x+1][y][table]) depend(out:ocean[x][y][table], ocean[x][y-1][table], ocean[x][y+1][table], ocean[x-1][y][table], ocean[x+1][y][table])
        if(ocean[x][y][table] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }

    for (int x = 2; x < DIM-1; x = x+2){
      for (int y = DIM-2; y > 0 ; y--){
        #pragma omp task depend(in:ocean[x][y][table], ocean[x][y-1][table], ocean[x][y+1][table], ocean[x-1][y][table], ocean[x+1][y][table]) depend(out:ocean[x][y][table], ocean[x][y-1][table], ocean[x][y+1][table], ocean[x-1][y][table], ocean[x+1][y][table])
        if(ocean[x][y][table] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }

    //#pragma omp taskwait
    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}

// ------------------------------------------------------------------------------
// ------------- Fonction initiale dépendante de la méthode appelée -------------
// ------------------------------------------------------------------------------

int cb_seq (int argc, char **argv)
{
  display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                cb_get,                // callback func
                cb_compute_seq);       // callback func
}

int cb_parallel (int argc, char **argv)
{
  //omp_set_nested(1);
  display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                cb_get,                // callback func
                cb_compute_parallel);  // callback func

}

int cb_parallel_task (int argc, char **argv)
{
  #pragma omp parallel

  #pragma omp single
  {
    display_init (argc, argv,
                  DIM,              // dimension ( = x = y) du tas
                  MAX_HEIGHT,       // hauteur maximale du tas
                  cb_get,              // callback func
                  cb_compute_task);    // callback func
  }
}