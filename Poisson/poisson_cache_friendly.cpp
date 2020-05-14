#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sys/time.h>
#include <omp.h>

#define TIME_RESOLUTION 1000000

const int block_length = 64;

long long unsigned initial_time;
long long unsigned init_update;
long long unsigned measure;
timeval t;

void start (void) {
    gettimeofday(&t, NULL);
    initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
    measure = 0;
}

void start_update(void){
  gettimeofday(&t, NULL);
  init_update = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}

void end_update(void){
  gettimeofday(&t, NULL);
  long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;

  measure += final_time - init_update;
}

long long unsigned stop (void) {
	gettimeofday(&t, NULL);
	long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;

	return final_time - initial_time;
}


void initBoundaries(double * arr, int size){

  for(int i = 0; i < size; i++){
    arr[i] = 0.0;
    arr[size*(size-1) + i] = 100.0;
  }

  for(int i = 0; i < size; i++){
    arr[i*size] = 100.0;
    arr[i*size + (size-1)] = 100.0;
    arr[size*(size-1) + i] = 100.0;
  }

  for(int i = 1; i < size-1; i++)
    for(int j = 1; j < size-1; j++){
      arr[i*size + j] = 50.0;
    }
}

double compute_diff_parallel(double * w, double * u, int size,int nthreads){
  double max = 0.0;

  #pragma omp parallel num_threads(nthreads)
  {
    double localMax = 0.0;
    double max_actual;
    #pragma omp for
    for(int i = 0; i < size; i++){
      for(int j = 0; j < size; j++){
        max_actual = fabs(w[i*size + j] - u[i*size + j]);
        if(max_actual > localMax)
          localMax = max_actual;
      }
    }

    #pragma omp critical
    {
      if (localMax > max){
        max = localMax;
      }
    }
  }
  //std::cout << max << std::endl;
  return max;
}

double compute_diff(double * w, double * u, int size){
  double max = 0.0;
  double max_actual;

  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      max_actual = fabs(w[i*size + j] - u[i*size + j]);
      if(max_actual > max)
        max = max_actual;
    }
  }

  return max;
}

/**

Implementação 3 com otimização. Não faz copia.

*/

int compute_block_row_parallel_optimized(double * w, int size,double p, double tolerance, int nthreads, int block_X, int block_Y){

  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;
  int I_BLOCK_SIZE;
  double localMax = 0.0,temp,max_actual;

  while(diff > tolerance){
    localMax = 0.0;
    diff = 0.0;
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp single
        start_update();

      #pragma omp for private(localMax,temp,max_actual)
      for(int i = 0; i < size; i += block_X){
        //primeiros 2 blocos Red à mão
        I_BLOCK_SIZE = (i + block_X == size) ? block_X - 1 : block_X;
        int ii = i == 0 ? 1 : i;
        int j = 0;
        int maxi = i + I_BLOCK_SIZE; // max i para os blocos Red
        int maxib = (i + block_X == size) ? i + I_BLOCK_SIZE : i + I_BLOCK_SIZE - 1;// max i para os blocos Black
        // BLOCO 1 RED
        for(; ii < maxi; ii++){
          int maxj = j + block_Y;
          int jj = 1 + ((ii + 1) & 1);
          for(; jj < maxj; jj += 2){
            temp = w[ii*size + jj];

            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

            max_actual = fabs(w[ii*size + jj] - temp);
            if (max_actual > localMax)
              localMax = max_actual;
          }
        }
        j += block_Y;
        ii = i == 0 ? 1 : i;
        //BLOCO 2 RED
        for(; ii < maxi; ii++){
          int maxj = j + block_Y;
          int jj = j + (ii & 1) ;
          for(; jj < maxj; jj += 2){
            temp = w[ii*size + jj];

            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

            max_actual = fabs(w[ii*size + jj] - temp);
            if (max_actual > localMax)
              localMax = max_actual;
          }
        }
        j += block_Y;
        int offset = 2 * block_Y;

        //loop genérico
        for(; j < size; j += block_Y){
          ii = i == 0 ? 1 : i + 1;
          // -1 por causa das bordas ! acertar no fim ! Block Black - 2
          for(; ii < maxib; ii++){
            int j_offset = j - offset;
            int maxj = j_offset + block_Y;
            int jj = j_offset == 0 ? 1 + (ii & 1) : j_offset + ((ii+1) & 1);
            for(; jj < maxj; jj += 2){
              temp = w[ii*size + jj];

              w[ii*size + jj] =
                aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

              max_actual = fabs(w[ii*size + jj] - temp);
              if (max_actual > localMax)
                localMax = max_actual;
            }
          }
          // Block Red
          ii = i == 0 ? 1 : i;
          for(; ii < maxi; ii++){
            int maxj = j + block_Y == size ? j + block_Y - 1 : j + block_Y ;
            int jj = j + (ii & 1);
            for(; jj < maxj; jj += 2){
              temp = w[ii*size + jj];

              w[ii*size + jj] =
                aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

              max_actual = fabs(w[ii*size + jj] - temp);
              if (max_actual > localMax)
                localMax = max_actual;
            }
          }
        }
        j -= offset;
        //ultimos 2 blocos BLACK
        ii = i == 0 ? 1 : i + 1;
        for(; ii < maxib; ii++){
          int maxj = j + block_Y;
          int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1);
          for(; jj < maxj; jj += 2){
            temp = w[ii*size + jj];

            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

            max_actual = fabs(w[ii*size + jj] - temp);
            if (max_actual > localMax)
              localMax = max_actual;
          }
        }
        j += block_Y;
        ii = i == 0 ? 1 : i + 1;
        for(; ii < maxib; ii++){
          int maxj = j + block_Y - 1;
          int jj = j + ((ii+1) & 1);
          for(; jj < maxj; jj += 2){
            temp = w[ii*size + jj];

            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

            max_actual = fabs(w[ii*size + jj] - temp);
            if (max_actual > localMax)
              localMax = max_actual;
          }
        }

        #pragma omp critical
        {
          if (localMax > diff){
            diff = localMax;
          }
        }
        localMax = 0.0;
      }
      //bordas
      #pragma omp for private(localMax,temp,max_actual)
      for(int i = 0; i < size - block_X; i+= block_X){
        for(int j = 0; j < size; j+= block_Y){
          int J_BLOCK_SIZE = (j + block_Y == size) ? block_Y - 1 : block_Y;
          int ii = i + block_X - 1;
          int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;
          int maxj = j + J_BLOCK_SIZE;
          for(; jj < maxj; jj += 2){
            temp = w[ii*size + jj];

            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

            max_actual = fabs(w[ii*size + jj] - temp);
            if (max_actual > localMax)
              localMax = max_actual;
          }

          ii = i + block_X;
          jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;
          for(; jj < maxj; jj += 2){
            temp = w[ii*size + jj];

            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;

            max_actual = fabs(w[ii*size + jj] - temp);
            if (max_actual > localMax)
              localMax = max_actual;
          }
        }

        #pragma omp critical
        {
          if (localMax > diff){
            diff = localMax;
          }
        }
        localMax = 0.0;
      }
    }

    end_update();

    iter++;
  }

  std::cout << "Computational time: " << measure << std::endl;
  return iter;
}

/**

Implementação 3: Versão paralela com o uso de blocos. Red e Black calculados na minha linha.

*/
int compute_block_row_parallel(double * w, int size,double p, double tolerance, int nthreads, int block_X, int block_Y){

  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;
  int I_BLOCK_SIZE;

  while(diff > tolerance){
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for
      for(int i = 0; i < size*size; i++)
        u[i] = w[i];

      #pragma omp single
        start_update();

      #pragma omp for
      for(int i = 0; i < size; i += block_X){
        //primeiros 2 blocos Red à mão
        I_BLOCK_SIZE = (i + block_X == size) ? block_X - 1 : block_X;
        int ii = i == 0 ? 1 : i;
        int j = 0;
        int maxi = i + I_BLOCK_SIZE; // max i para os blocos Red
        int maxib = (i + block_X == size) ? i + I_BLOCK_SIZE : i + I_BLOCK_SIZE - 1;// max i para os blocos Black
        // BLOCO 1 RED
        for(; ii < maxi; ii++){
          int maxj = j + block_Y;
          int jj = 1 + ((ii + 1) & 1);
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
        j += block_Y;
        ii = i == 0 ? 1 : i;
        //BLOCO 2 RED
        for(; ii < maxi; ii++){
          int maxj = j + block_Y;
          int jj = j + (ii & 1) ;
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
        j += block_Y;
        int offset = 2 * block_Y;

        //loop genérico
        for(; j < size; j += block_Y){
          ii = i == 0 ? 1 : i + 1;
          // -1 por causa das bordas ! acertar no fim ! Block Black - 2
          for(; ii < maxib; ii++){
            int j_offset = j - offset;
            int maxj = j_offset + block_Y;
            int jj = j_offset == 0 ? 1 + (ii & 1) : j_offset + ((ii+1) & 1);
            for(; jj < maxj; jj += 2){
              w[ii*size + jj] =
                aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
            }
          }
          // Block Red
          ii = i == 0 ? 1 : i;
          for(; ii < maxi; ii++){
            int maxj = j + block_Y == size ? j + block_Y - 1 : j + block_Y ;
            int jj = j + (ii & 1);
            for(; jj < maxj; jj += 2){
              w[ii*size + jj] =
                aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
            }
          }
        }
        j -= offset;
        //ultimos 2 blocos BLACK
        ii = i == 0 ? 1 : i + 1;
        for(; ii < maxib; ii++){
          int maxj = j + block_Y;
          int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1);
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
        j += block_Y;
        ii = i == 0 ? 1 : i + 1;
        for(; ii < maxib; ii++){
          int maxj = j + block_Y - 1;
          int jj = j + ((ii+1) & 1);
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
      }
      //bordas
      #pragma omp for
      for(int i = 0; i < size - block_X; i+= block_X){
        for(int j = 0; j < size; j+= block_Y){
          int J_BLOCK_SIZE = (j + block_Y == size) ? block_Y - 1 : block_Y;
          int ii = i + block_X - 1;
          int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;
          int maxj = j + J_BLOCK_SIZE;
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }

          ii = i + block_X;
          jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
      }
    }

    end_update();

    iter++;
    diff = compute_diff_parallel(w,u,size,nthreads);
  }

  std::cout << "Computational time: " << measure << std::endl;
  return iter;
}

/**

Implementação 3: Versão sequencial com o uso de blocos. Red e Black calculados na minha linha.

*/
int compute_block_row_sequencial(double * w, int size,double p, double tolerance){

  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;
  int I_BLOCK_SIZE;

  while(diff > tolerance){
    std::copy(w, w + size*size, u);

    for(int i = 0; i < size; i += block_length){
      //primeiros 2 blocos Red à mão
      I_BLOCK_SIZE = (i + block_length == size) ? block_length - 1 : block_length;
      int ii = i == 0 ? 1 : i;
      int j = 0;
      int maxi = i + I_BLOCK_SIZE; // max i para os blocos Red
      int maxib = (i + block_length == size) ? i + I_BLOCK_SIZE : i + I_BLOCK_SIZE - 1;// max i para os blocos Black
      // BLOCO 1 RED
      for(; ii < maxi; ii++){
        int maxj = j + block_length;
        int jj = 1 + ((ii + 1) & 1);
        for(; jj < maxj; jj += 2){
          w[ii*size + jj] =
            aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
        }
      }
      j += block_length;
      ii = i == 0 ? 1 : i;
      //BLOCO 2 RED
      for(; ii < maxi; ii++){
        int maxj = j + block_length;
        int jj = j + (ii & 1) ;
        for(; jj < maxj; jj += 2){
          w[ii*size + jj] =
            aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
        }
      }
      j += block_length;
      int offset = 2 * block_length;

      //loop genérico
      for(; j < size; j += block_length){
        ii = i == 0 ? 1 : i + 1;
        // -1 por causa das bordas ! acertar no fim ! Block Black - 2
        for(; ii < maxib; ii++){
          int j_offset = j - offset;
          int maxj = j_offset + block_length;
          int jj = j_offset == 0 ? 1 + (ii & 1) : j_offset + ((ii+1) & 1);
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
        // Block Red
        ii = i == 0 ? 1 : i;
        for(; ii < maxi; ii++){
          int maxj = j + block_length == size ? j + block_length - 1 : j + block_length ;
          int jj = j + (ii & 1);
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
      }
      j -= offset;
      //ultimos 2 blocos BLACK
      ii = i == 0 ? 1 : i + 1;
      for(; ii < maxib; ii++){
        int maxj = j + block_length;
        int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1);
        for(; jj < maxj; jj += 2){
          w[ii*size + jj] =
            aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
        }
      }
      j += block_length;
      ii = i == 0 ? 1 : i + 1;
      for(; ii < maxib; ii++){
        int maxj = j + block_length - 1;
        int jj = j + ((ii+1) & 1);
        for(; jj < maxj; jj += 2){
          w[ii*size + jj] =
            aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
        }
      }
    }

    for(int i = 0; i < size - block_length; i+= block_length){
      for(int j = 0; j < size; j+= block_length){
        int J_BLOCK_SIZE = (j + block_length == size) ? block_length - 1 : block_length;
        int ii = i + block_length - 1;
        int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;
        int maxj = j + J_BLOCK_SIZE;
        for(; jj < maxj; jj += 2){
          w[ii*size + jj] =
            aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
        }

        ii = i + block_length;
        jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;
        for(; jj < maxj; jj += 2){
          w[ii*size + jj] =
            aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
        }
      }
    }

    iter++;
    diff = compute_diff(w,u,size);
  }

  return iter;
}

/**

Nenhuma das implementações do relatório. Versão sequencial para uma versão naive do uso de blocos ainda com as iterações Red e Black separadas.

*/
int compute_block_row_sequencial_naive(double * w, int size,double p, double tolerance){
  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;
  int counter = 0;

  while(diff > tolerance){
    std::copy(w, w + size*size, u);

    for(int i = 0; i < size; i += block_length){
      for(int j = 0; j <  size; j += block_length){
        int I_BLOCK_SIZE = (i + block_length == size) ? block_length - 1 : block_length;
        int J_BLOCK_SIZE = (j + block_length == size) ? block_length - 1 : block_length;

        int maxi = i + I_BLOCK_SIZE;
        int ii = i == 0 ? 1 : i;
        for(; ii < maxi; ii++){
          int maxj = j + J_BLOCK_SIZE;
          int jj = j == 0 ? 1 + ((ii + 1) & 1) : j + (ii & 1) ;
          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
              counter++;
          }
        }
      }
    }
    for(int i = 0; i < size; i += block_length){
      for(int j = 0; j <  size; j += block_length){
        int I_BLOCK_SIZE = (i + block_length == size) ? block_length - 1 : block_length;
        int J_BLOCK_SIZE = (j + block_length == size) ? block_length - 1 : block_length;

        int maxi = i + I_BLOCK_SIZE;
        int ii = i == 0 ? 1 : i;
        for(; ii < maxi; ii++){
          int maxj = j + J_BLOCK_SIZE;
          int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;

          for(; jj < maxj; jj += 2){
            w[ii*size + jj] =
              aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
          }
        }
      }
    }

    iter++;
    diff = compute_diff(w,u,size);
  }

  return iter;
}

/**

Nenhuma das implementações do relatório. Versão paralela para uma versão naive do uso de blocos ainda com as iterações Red e Black separadas.

*/
int compute_block_row_parallel_naive(double * w, int size,double p, double tolerance, int nthreads){
  //Red points
  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;

  while(diff > tolerance){
    //std::copy(w, w + size*size, u);
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for
      for(int i = 0; i < size*size; i++)
        u[i] = w[i];

      #pragma omp single
        start_update();


      #pragma omp for
      for(int i = 0; i < size; i += block_length){
        for(int j = 0; j <  size; j += block_length){
          int I_BLOCK_SIZE = (i + block_length == size) ? block_length - 1 : block_length;
          int J_BLOCK_SIZE = (j + block_length == size) ? block_length - 1 : block_length;

          int maxi = i + I_BLOCK_SIZE;
          int ii = i == 0 ? 1 : i;
          for(; ii < maxi; ii++){
            int maxj = j + J_BLOCK_SIZE;
            int jj = j == 0 ? 1 + ((ii + 1) & 1) : j + (ii & 1) ;
            for(; jj < maxj; jj += 2){
              w[ii*size + jj] =
                aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
            }
          }
        }
      }

      #pragma omp for
      for(int i = 0; i < size; i += block_length){
        for(int j = 0; j <  size; j += block_length){
          int I_BLOCK_SIZE = (i + block_length == size) ? block_length - 1 : block_length;
          int J_BLOCK_SIZE = (j + block_length == size) ? block_length - 1 : block_length;

          int maxi = i + I_BLOCK_SIZE;
          int ii = i == 0 ? 1 : i;
          for(; ii < maxi; ii++){
            int maxj = j + J_BLOCK_SIZE;
            int jj = j == 0 ? 1 + (ii & 1) : j + ((ii+1) & 1) ;

            for(; jj < maxj; jj += 2){
              w[ii*size + jj] =
                aux*w[ii*size + jj] + p * ( w[(ii-1)*size + jj] + w[ii*size + jj - 1] + w[ii*size + jj + 1] + w[(ii+1)*size + jj] ) / 4.0;
            }
          }
        }
      }
    }

    end_update();

    iter++;
    diff = compute_diff_parallel(w,u,size,nthreads);
  }

  std::cout << "Computational time: " << measure << std::endl;
  return iter;
}

/**

Implementação 2 do relatório: Versão totalmente paralela e sem blocos com a implementação dada pelo professor

*/
int compute_noblock_row_parallel(double * w, int size,double p, double tolerance, int nthreads){

  //Red points
  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;

  while(diff > tolerance){
    //std::copy(w, w + size*size, u);

    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for
      for(int i = 0; i < size*size; i++)
        u[i] = w[i];

      #pragma omp single
        start_update();


      #pragma omp for
      for(int i = 1; i < size - 1; i++){
        for(int j = 1 + ((i+1) & 1); j < size - 1; j += 2){
          w[i*size + j] = aux*w[i*size + j] + p * ( w[(i-1)*size + j] + w[i*size + j - 1] + w[i*size + j + 1] + w[(i+1)*size + j] ) / 4.0;
        }
      }


      #pragma omp for
      for(int i = 1; i < size - 1; i++){
        for(int j = 1 + (i & 1); j < size - 1; j += 2){
          w[i*size + j] = aux*w[i*size + j] + p * ( w[(i-1)*size + j] + w[i*size + j - 1] + w[i*size + j + 1] + w[(i+1)*size + j] ) / 4.0;
        }
      }
    }
    end_update();

    iter++;
    diff = compute_diff_parallel(w,u,size,nthreads);
    //std::cout << diff << std::endl;
  }

  std::cout << "Computational time: " << measure << std::endl;

  return iter;
}


/**

Implementação 1 do relatório: Versão paralela e sem blocos com a implementação dada pelo professor.

*/

int compute_noblock_row_semi_parallel(double * w, int size,double p, double tolerance, int nthreads){

  //Red points
  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;

  while(diff > tolerance){
    std::copy(w, w + size*size, u);

    start_update();
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for
      for(int i = 1; i < size - 1; i++){
        for(int j = 1 + ((i+1) & 1); j < size - 1; j += 2){
          w[i*size + j] = aux*w[i*size + j] + p * ( w[(i-1)*size + j] + w[i*size + j - 1] + w[i*size + j + 1] + w[(i+1)*size + j] ) / 4.0;
        }
      }


      #pragma omp for
      for(int i = 1; i < size - 1; i++){
        for(int j = 1 + (i & 1); j < size - 1; j += 2){
          w[i*size + j] = aux*w[i*size + j] + p * ( w[(i-1)*size + j] + w[i*size + j - 1] + w[i*size + j + 1] + w[(i+1)*size + j] ) / 4.0;
        }
      }
    }
    end_update();

    iter++;
    diff = compute_diff(w,u,size);
    //std::cout << diff << std::endl;
  }

  std::cout << "Computational time: " << measure << std::endl;

  return iter;
}

/**

Implementação 1 do relatório: Versão sequencial e sem blocos com a implementação dada pelo professor.

*/

int  compute_noblock_sequencial(double * w, int size,double p, double tolerance){

  //Red points
  double * u = new double[size*size];
  double aux = 1.0 - p;
  double diff = 1.0 + tolerance;
  int iter = 0;

  while(diff > tolerance){
    std::copy(w, w + size*size, u);

    start_update();

    for(int i = 1; i < size - 1; i++){
      for(int j = 1 + ((i+1) & 1); j < size - 1; j += 2){
        w[i*size + j] = aux*w[i*size + j] + p * ( w[(i-1)*size + j] + w[i*size + j - 1] + w[i*size + j + 1] + w[(i+1)*size + j] ) / 4.0;
      }
    }

    for(int i = 1; i < size - 1; i++){
      for(int j = 1 + (i & 1); j < size - 1; j += 2){
        w[i*size + j] = aux*w[i*size + j] + p * ( w[(i-1)*size + j] + w[i*size + j - 1] + w[i*size + j + 1] + w[(i+1)*size + j] ) / 4.0;
      }
    }

    end_update();

    iter++;
    diff = compute_diff(w,u,size);
    //std::cout << diff << std::endl;
  }

  std::cout << "Computational time: " << measure << std::endl;

  return iter;
}

int main(int argc, char const *argv[]) {

  if (argc < 5){
    std::cout << "Erro nos argumentos" << std::endl;
    exit(1);
  }
  long long unsigned tt;

  int SIZE = atoi(argv[1]);
  double TOL = std::stod(argv[2]);
  int N_THREADS = atoi(argv[3]);
  int X_BLOCKS = atoi(argv[4]);
  int Y_BLOCKS = atoi(argv[5]);

  std::cout << "Config:" << std::endl;
  std::cout << "N: " << SIZE << std::endl;
  std::cout << "TOL: " << TOL << std::endl;
  std::cout << "THREADS: " << N_THREADS << std::endl;
  std::cout << "<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>" << std::endl << std::endl;


  double * w = new double[SIZE*SIZE];
  double * u = new double[SIZE*SIZE];
  double * tmp = new double[SIZE*SIZE];

  //Parâmetro de relaxação
  double p = 2/(1+sin(M_PI/(SIZE-1)));

  int n_block_per_line = SIZE / block_length;
  int n_blocks = (SIZE * SIZE) / (block_length * block_length);

  initBoundaries(w,SIZE);
  initBoundaries(tmp,SIZE);

  std::cout << p << std::endl;


  start();
  int iter1 = compute_noblock_row_parallel(tmp,SIZE,p,TOL,N_THREADS);
  tt = stop();
  std::cout << "Parallel Compute no Block: " << tt << " usecs" << std::endl;

  start();
  int iter = compute_block_row_parallel(w,SIZE,p,TOL,N_THREADS,X_BLOCKS,Y_BLOCKS);
  //int iter = compute_noblock_sequencial(w,SIZE,p,TOL);
  tt = stop();
  std::cout << "Sequencial: " << tt << " usecs" << std::endl;

  std::cout << "ITER: " << iter << std::endl;
  std::cout << "ITER1: " << iter1 << std::endl;


  bool flag = false;
  for(int i = 0; i < SIZE && !flag; i++){
    for(int j = 0; j < SIZE; j++){
      if ((w[i*SIZE + j] - tmp[i*SIZE + j]) != 0.0){
        std::cout <<  "ERRO: " << i << " - " << j << " : "  << w[i*SIZE + j] - tmp[i*SIZE + j] << std::endl;
        flag = true;
        break;
      }
    }
  }

  /**
  for(int i = 0; i < SIZE; i++){
    for(int j = 0; j < SIZE; j++){
      std::cout << w[i*SIZE + j] << " : ";
    }
    std::cout << std::endl;
  }

    std::cout << "          "  << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

  for(int i = 0; i < SIZE; i++){
      for(int j = 0; j < SIZE; j++){
        std::cout << tmp[i*SIZE + j] << " : ";
      }
      std::cout << std::endl;
    }*/

  return 0;
}
