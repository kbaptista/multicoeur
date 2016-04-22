#include "treatments.h"

// ------------------------------------------------------------------------------
// -------------              Fonctions de traitement               -------------
// ------------------------------------------------------------------------------

// callback
unsigned get (unsigned x, unsigned y)
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
void compute_cell(int x, int y, int div4){
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
float *compute_seq(unsigned iterations){
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
float *compute_parallel(unsigned iterations){
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
void treatment_task(int x, int y){
  int div4 = ocean[x][y][table]/4;

  //#pragma omp critical
  {
    compute_cell(x, y, div4);
  }
}

/** 
  * Compute fonction for task treatment
*/
float *compute_task(unsigned iterations)
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

void compute_cl(){

}