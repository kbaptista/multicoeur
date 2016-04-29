#include "display.h"
#include "treatment.h"

#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define MAX_HEIGHT  4

#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned DIM;
unsigned *ocean;

static bool is_end = false;

// ------------------------------------------------------------------------------
// -------------               Utility Functions                    -------------
// ------------------------------------------------------------------------------

/** 
  * Print matrices line after line on the shell.
  * Is adapted to print values until 999 without deformations.
*/
static void print()
{
  for (int x = 0; x < DIM; x++)
  {
    for (int y = 0; y < DIM; y++) 
    {
      if(ocean[x*DIM+y]< 10)
      {
        printf("  %d ",ocean[x*DIM+y]);
      }
      else if (ocean[x*DIM+y] < 100)
      {
        printf(" %d ",ocean[x*DIM+y]);
      }
      else
      {
        printf("%d ",ocean[x*DIM+y]);
      }
    }
    printf("\n");
  }
}

/** 
  * initialize matrices with 5 squares in every square.
*/
void sand_init_homogeneous()
{
  for (int x = 0; x < DIM; x++)
  {
    for (int y = 0; y < DIM; y++)
    {
      if(x > 0 && x < DIM-1 && y > 0 && y < DIM-1)
        ocean[x*DIM+y] = 5;
      else
        
        ocean[x*DIM+y] = 0;
    }
  }
}

/** 
  * initialize matrices with no one grain except in the middle case which have 100000 grain.
*/
void sand_init_center()
{
  int center_value = 100000;
  for (int x = 0; x < DIM; x++)
  {
    for (int y = 0; y < DIM; y++)
    {
      ocean[x*DIM+y] = 0;
    }
  }
  ocean[DIM*DIM/2+DIM/2] = center_value;
}

// callback
unsigned get (unsigned x, unsigned y)
{
  return ocean[x*DIM+y];
}

// ------------------------------------------------------------------------------
// -------------              Treatments Functions               ----------------
// ------------------------------------------------------------------------------

/** 
  * Applied modifications into 1-table for ocean[x][y] sandpile.
*/
static inline void compute_cell(int x, int y, int div4)
{
  ocean[x*DIM+y]      -= div4*4;
  ocean[x*DIM+y+1]    += div4;
  ocean[x*DIM+y-1]    += div4;
  ocean[(x-1)*DIM+y]  += div4;
  ocean[(x+1)*DIM+y]  += div4;
  is_end = false;
}



/** 
  * Compute fonction for sequential treatment
*/
static inline float *compute_seq(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;
 
  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}


/*
TODO : ne peut pas fonctionner tel quel. revoir la conception de la ligne
*/
/*
static inline float *compute_seq_find_stabilize(unsigned iterations)
{
  if(init == false){
    init = true;
  }

  if(is_end == true)
  {
    return DYNAMIC_COLORING;
    free(tmp);
  }
  is_end = true;
 
  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      if(tmp[i] == true)
      {
        for (int y = 1; y < DIM-1; y++)
        {
          if(ocean[x*DIM+y] >= MAX_HEIGHT)
          {
            int div4 = ocean[x*DIM+y]/4;
            compute_cell(x,y,div4);
            if(ocean[(x+1)*DIM+y] >= MAX_HEIGHT)
              tmp[i+1] = true;
            if(ocean[(x-1)*DIM+y] >= MAX_HEIGHT)
              tmp[i-1] = true;
          }
        }
      }

    }
  }
  return DYNAMIC_COLORING;
}*/

// TODO : ne fonctionne pas ! plus le temps
//methode nécessitant deux matrices : une écriture une lecture
  /*
static inline float *compute_seq_reverse(unsigned iterations)
{
  if(init == false){
    tmp = malloc(sizeof(unsigned)*DIM*DIM);
    init = true;
  }

  if(is_end == true)
  {
    print();
    return DYNAMIC_COLORING;
  }
  is_end = false;
 
  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        int val = ocean[x*DIM+y]%4;
        val += ocean[(x-1)*DIM+y]/4;
        val += ocean[(x+1)*DIM+y]/4;
        val += ocean[x*DIM+y-1]/4;
        val += ocean[x*DIM+y+1]/4;

        if(val != ocean[x*DIM+y]){
          is_end = true;
        }
        tmp[x*DIM+y] = val;
      }
    }
    memcpy(ocean,tmp,DIM*DIM*sizeof(unsigned));
  }
  return DYNAMIC_COLORING;
}*/

static inline float *compute_seq_multipleline_expensif(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-5; y+=4)
      {
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
        if(ocean[x*DIM+y+1] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+1]/4;
          compute_cell(x,y+1,div4);
        }
        if(ocean[x*DIM+y+2] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+2]/4;
          compute_cell(x,y+2,div4);
        }
        if(ocean[x*DIM+y+3] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+3]/4;
          compute_cell(x,y+3,div4);
        }
      } 
    }
  }
  return DYNAMIC_COLORING;
}

static inline float *compute_seq_multipleline_introspectif(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  int table = 0;
  int table_2 = DIM*DIM;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-5; y+=4)
      {
          int val = ocean[table+x*DIM+y]%4 + 
                    ocean[table+(x-1)*DIM+y]/4 +
                    ocean[table+(x+1)*DIM+y]/4 +
                    ocean[table+x*DIM+y-1]/4 +
                    ocean[table+x*DIM+y+1]/4;
          ocean[table_2+x*DIM+y] = val;

          val = ocean[table+x*DIM+y+1]%4 + 
                ocean[table+(x-1)*DIM+y+1]/4 +
                ocean[table+(x+1)*DIM+y+1]/4 +
                ocean[table+x*DIM+y]/4 +
                ocean[table+x*DIM+y+2]/4;
          ocean[table_2+x*DIM+y+1] = val;

          val = ocean[table+x*DIM+y+2]%4 + 
                ocean[table+(x-1)*DIM+y+2]/4 +
                ocean[table+(x+1)*DIM+y+2]/4 +
                ocean[table+x*DIM+y+1]/4 +
                ocean[table+x*DIM+y+3]/4;
          ocean[table_2+x*DIM+y+2] = val;

          val = ocean[table+x*DIM+y+3]%4 + 
                ocean[table+(x-1)*DIM+y+3]/4 +
                ocean[table+(x+1)*DIM+y+3]/4 +
                ocean[table+x*DIM+y+2]/4 +
                ocean[table+x*DIM+y+4]/4;
          ocean[table_2+x*DIM+y+3] = val;

          if(ocean[table+x*DIM+y] >= MAX_HEIGHT || ocean[table+(x-1)*DIM+y] >= MAX_HEIGHT || 
            ocean[table+(x+1)*DIM+y] >= MAX_HEIGHT || ocean[table+x*DIM+y-1] >= MAX_HEIGHT || ocean[table+x*DIM+y+1] >= MAX_HEIGHT)
            is_end = false;
        /*
        
        if(ocean[x*DIM+y+1] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+1]/4;
          compute_cell(x,y+1,div4);
        }
        if(ocean[x*DIM+y+2] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+2]/4;
          compute_cell(x,y+2,div4);
        }
        if(ocean[x*DIM+y+3] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+3]/4;
          compute_cell(x,y+3,div4);
        }*/
      }
      //boucle nécessaire pour finir le traitement de la matrice si DIM-1 n'est pas multiple de 5
      if((DIM-1)%5 != 0){
        for (int y = (DIM-1)%5; y < DIM-1; y++){
          int val = ocean[table+x*DIM+y]%4 + 
                      ocean[table+(x-1)*DIM+y]/4 +
                      ocean[table+(x+1)*DIM+y]/4 +
                      ocean[table+x*DIM+y-1]/4 +
                      ocean[table+x*DIM+y+1]/4;
          ocean[table+x*DIM+y] = val;
          if(ocean[table+x*DIM+y] >= MAX_HEIGHT)
            is_end = false;
        }
      }
    }
    table = table == 0 ? DIM*DIM:0;
    table_2 = table_2 == 0 ? DIM*DIM:0;
  }

//si on termine les itérations avec écriture dans la seconde table, on doit recopier son contenu dans la première
  if(table != 0)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        ocean[table+x*DIM+y] = ocean[table_2+x*DIM+y];
      }
    }
  }

  print();
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for collapsed treatment 
*/
static inline float *compute_parallel_for(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
#pragma omp parallel shared(is_end,ocean)
    {
    unsigned ocean_private[DIM*DIM];

#pragma omp for schedule(static,4)
        for (int x = 1; x < DIM-1; x++)
        {
          for (int y = 1; y < DIM-1; y++)
          {
          int val = ocean_private[x*DIM+y]%4 + 
                    ocean_private[(x-1)*DIM+y]/4 +
                    ocean_private[(x+1)*DIM+y]/4 +
                    ocean_private[x*DIM+y-1]/4 +
                    ocean_private[x*DIM+y+1]/4;

            ocean_private[x*DIM+y] = val;

            if(ocean[x*DIM+y] >= MAX_HEIGHT)
              is_end = false;
          }
        }

#pragma omp for schedule(static,4)
      for (int x = 1; x < DIM-1; x++){
        for (int y = 1; y < DIM-1; y++){
          ocean[x*DIM+y] = ocean_private[x*DIM+y];
        }
      }
    }
  }

  return DYNAMIC_COLORING;
}

static inline float *compute_parallel_p_iteration(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  unsigned nb_lines = round((DIM-2) / omp_get_max_threads());

  //nb de lignes spécifiques au dernier thread : il prend la différence.
  unsigned last_thread_lines = nb_lines*omp_get_max_threads() != DIM-2 ? nb_lines+abs(DIM-2 - nb_lines*omp_get_max_threads()) : nb_lines;

#pragma omp parallel shared(is_end, iterations) firstprivate(nb_lines)
  {
    unsigned my_num = omp_get_thread_num();
    int table = 0;

    int begin = ((my_num)*nb_lines)-iterations; 

    unsigned ocean_private[DIM*((my_num==omp_get_max_threads()-1?last_thread_lines:nb_lines)+iterations*2)][2];

    for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines+iterations*2:nb_lines+iterations*2 ; x < final ; x++)
    {
      for(int y = 0 ; y < DIM ; y++)
      {
        if(begin+x > 0 && begin+x < DIM-1)
        {
          ocean_private[x*DIM+y][table] = ocean[(begin+x/*+1*/)*DIM+y];
          // if(my_num == 1 && y == 0)
          //   printf("  (%d)\n  ",x);
        }
        else{
          ocean_private[x*DIM+y][table] = 0;
        }
        ocean_private[x*DIM+y][1-table] = 0;
      }
    }

    if(my_num == 1){
      printf(" --- pre-iteration --- je traite %d\n",begin+iterations);
      for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines+iterations*2:nb_lines+iterations*2; x < final ; x++)
      {
        for(int y = 0 ; y < DIM ; y++)
        {
          printf("%d ",ocean_private[x*DIM+y][table]);
          printf("%d - ",ocean_private[x*DIM+y][1-table]);
        }
        printf("[%d]\n",x);
      }
    }

    //fin de l'initialisation



    for (unsigned i = 1; i <= iterations; i++)
    {
      int taille = (iterations-i)*2+(my_num==omp_get_max_threads()-1?last_thread_lines:nb_lines);
      for (int x = i; x < i+taille; x++)
      {
        for (int y = 1; y < DIM-1; y++)
        {
          int val = ocean_private[x*DIM+y][table]%4 + 
                    ocean_private[(x-1)*DIM+y][table]/4 +
                    ocean_private[(x+1)*DIM+y][table]/4 +
                    ocean_private[x*DIM+y-1][table]/4 +
                    ocean_private[x*DIM+y+1][table]/4;

          ocean_private[x*DIM+y][1-table] = val;

          if(ocean_private[x*DIM+y][table] >= MAX_HEIGHT)
            is_end = false;
        }
        
      }

      
      table = 1-table;
    }
//tmp
    // if(my_num == 10){
    //   printf(" --- post_iterations ---\n");
    //   for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines+iterations*2:nb_lines+iterations*2; x < final ; x++)
    //   {
    //     for(int y = 0 ; y < DIM ; y++)
    //     {
    //       printf("%d ",ocean_private[x*DIM+y][table]);
    //       printf("%d - ",ocean_private[x*DIM+y][1-table]);
    //     }
    //     printf("[%d]\n",x);
    //   }
    // }


    begin = my_num*nb_lines;
    for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines:nb_lines; x < final ; x++)
    {
      for(int y = 1 ; y < DIM-1 ; y++)
      {
        if(begin+x+1 > 0 && begin+x < DIM-1)
        {
#pragma omp critical
          ocean[(begin+x+1)*DIM+y] = ocean_private[(x+iterations)*DIM+y][table];
          // printf("\n- %d - %d -\n",(begin+x+1), ocean_private[(x+iterations)*DIM+y][table]);
        }
      }
    }
  }
  //tmp
    print();
  // exit(0);

  return DYNAMIC_COLORING;
}

/*
  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
  0   1   2   5   1   3   3   3   3   3   3   1   5   2   1   0 
  0   2   3   1   5   1   3   3   3   3   1   5   1   3   2   0 
  0   5   1   7   3   7   3   5   5   3   7   3   7   1   5   0 
  0   1   5   3   7   3   7   3   3   7   3   7   3   5   1   0 
  0   3   1   7   3   7   3   5   5   3   7   3   7   1   3   0 
  0   3   3   3   7   3   5   5   5   5   3   7   3   3   3   0 
  0   3   3   5   3   5   5   5   5   5   5   3   5   3   3   0 
  0   3   3   5   3   5   5   5   5   5   5   3   5   3   3   0 
  0   3   3   3   7   3   5   5   5   5   3   7   3   3   3   0 
  0   3   1   7   3   7   3   5   5   3   7   3   7   1   3   0 
  0   1   5   3   7   3   7   3   3   7   3   7   3   5   1   0 
  0   5   1   7   3   7   3   5   5   3   7   3   7   1   5   0 
  0   2   3   1   5   1   3   3   3   3   1   5   1   3   2   0 
  0   1   2   5   1   3   3   3   3   3   3   1   5   2   1   0 
  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
*/

/** 
  * Compute fonction for task treatment
*/
static inline float *compute_parallel_task(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
#pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}

// --------------------------------------------------------------------------------------------------------
// ----------------------------- Initials functions to adapt at arguments  --------------------------------
// --------------------------------------------------------------------------------------------------------

int display(int argc, char ** argv)
{
  switch((int)argv[4][0])
  {
    case 115 : // ascii s
      display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                get,                // callback func
                compute_seq);       // callback func
      break;
    case 102 : //ascii f
      omp_set_nested(1);
      display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                get,                // callback func
                compute_parallel_for);  // callback func
      break;
    case 116 : //ascii t
#pragma omp parallel

#pragma omp single
      {
        display_init (argc, argv,
                      DIM,              // dimension ( = x = y) du tas
                      MAX_HEIGHT,       // hauteur maximale du tas
                      get,              // callback func
                      compute_parallel_task);    // callback func
      }
      break;
    case 105 : //ascii i
      omp_set_nested(1);
      display_init (argc, argv,
              DIM,              // dimension ( = x = y) du tas
              MAX_HEIGHT,       // hauteur maximale du tas
              get,              // callback func
              compute_parallel_p_iteration);    // callback func
      break;
    default :
      printf("Unrecognize Algorithm. Please, read our manual by using ./sand\n");
      break;
  }
}

void without_display(float * (*compute_func_t) (unsigned iterations))
{
  // TODO : prendre une décision pour la valeur de interation par défault.
  do{
    compute_func_t(10);
  }
  while(is_end == false);
}

int performance(int argc, char ** argv)
{
  unsigned long temps;
  struct timeval t1, t2;

  gettimeofday(&t1,NULL);

  switch((int)argv[4][0])
  {
    case 115 : // ascii s
      printf("Sequential Algorithm\n");
      without_display(compute_seq);
      break;
    case 108 : //ascii l
      // on a besoin de deux matrices dans ce cas particulier
      ocean = realloc(ocean,sizeof(unsigned)*DIM*DIM*2);

      //initialisation de la seconde matrice
      for (int x = 0; x < DIM; x++)
      {
        for (int y = 0; y < DIM; y++)
        {
          ocean[DIM*DIM+x*DIM+y] = 0;
        }
      }

      without_display(compute_seq_multipleline);
      printf("Sequential Algorithm | 4 lines each loop turn\n");
      break;
    case 102 : //ascii f
      without_display(compute_parallel_for);
      printf("Parallel Algorithm\n");
      break;
    case 116 : //ascii t
      without_display(compute_parallel_task);
      printf("Parallel Task Algorithm\n");
      break;
    case 105 : //ascii i
      without_display(compute_parallel_p_iteration);
      printf("Parallel Task Algorithm\n");
      break;
    default :
      printf("Unrecognize Algorithm. Please, read our manual by using ./sand\n");
      break;
  }
  gettimeofday(&t2,NULL);
   
  temps = TIME_DIFF(t1,t2);
  printf("Algorithm time = %ld.%03ldms \n", temps/1000, temps%1000);
  //print();
  free(ocean);
}


void treatment(int argc, char ** argv)
{
  DIM = strtol(argv[3],NULL,10);
  if(!DIM)
  {
    printf("Unrecognize size. Please, read our manual by using ./sand\n");
    return;
  }
  
  // pas d'autres allocations durant le traitement,
  // ne sera pas free par nous à la fin si on entre dans display car aurait nécessité de modifier display 
  // MAIS sera nettoyé sans perte à la fermture du processus en tout les cas.
  ocean = malloc(sizeof(unsigned)*DIM*DIM); 

   //initialisation de la matrice
  switch((int)argv[2][0])
  {
    case 99 : //ascii c
      sand_init_center();
      printf("Centered Case\n");
      break;
    case 104: //ascii h
      sand_init_homogeneous();
      printf("Homogeneous Case\n");
      break;
    default :
      printf("Unrecognize configuration. Please, read our manual by using ./sand\n");
      break;
  }

  // objectif
  if(!atoi(argv[1]))
  {
    performance(argc,argv);
  }
  else
  {
    display(argc,argv);  
  }
}
