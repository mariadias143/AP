
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <limits>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <functional>
#include <stdlib.h>
#include <sys/time.h>

#define TIME_RESOLUTION 1000000

long long unsigned initial_time;
timeval t;
//Parâmetro de relaxação
double p = 2/(1+sin(pi/(N-1)));

void start (void) {
	gettimeofday(&t, NULL);
	initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}

void printResults (long long unsigned tt, int n_execs, int type) {
  switch (type) {
    case 0:{
        std::cout << "Sequencial execution time: " << (tt / n_execs) << " usecs" << std::endl;
        break;
    }
    case 1:{
      std::cout << "OpenMP execution time: " << (tt / n_execs) << " usecs" << std::endl;
      break;
    }
    default:{
      std::cout << "error" << std::endl;
    }
  }
}

long long unsigned stop (void) {
	gettimeofday(&t, NULL);
	long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;

	return final_time - initial_time;
}


int computeDifference(int * u, int * w, int n){
  int * diff = new int[n*n]; 

  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      diff[i*n+j] = abs(w[i*n+j] - u[i*n+j]);
    
  int max = std::max_element(diff, diff + n*n);
  return max;
}


void poisson(int n, int tol){
	long long unsigned tt;

  int * w = new int[n*n]; // values of current iteration
  int * u = new int[n*n]; // last iteration values

  // Temperature is zero at the top and 100 on the other boundaries
  for(int i=0; i<n; i++){
    grid[i] = 0; //above border
    grid[i + n*(n-1)] = 100; //lower border
  }

  for(int j=0; j<n; j++){
    grid[j*n] = 100; //left border
    grid[n-1 + j*n] = 100; //right border
  }

  for(int i=1; i<n-1; i++)
    for(int j=1; j<n-1; i++)
      grid[i*n + j] = 50; //initial aproximation for interior points
  

  int diff = tol + 1;
  int n_iter = 0;
  while(diff > tol){
    std::copy(std::begin(w), std::end(w), std::begin(u)); // u = w;

    // Update red
    for(int i=1; i<n-1; i++)
      for(int j=1+remainder(i,2); j<n-1; j+=2)
        w[i*n + j] = (1-p)*w[i*n + j] + p * ( w[(i-1)*n+j] + w[i*n+(j-1)] + w[(i+1)*n+j] + w[i*n+(j+1)] ) / 4;

    // Update black
    for(int i=1; i<n-1; i++)
      for(int j=1+remainder(i+1,2); j<n-1; j+=2)
        w[i*n + j] = (1-p)*w[i*n + j] + p * ( w[(i-1)*n+j] + w[i*n+(j-1)] + w[(i+1)*n+j] + w[i*n+(j+1)] ) / 4;

    n_iter++;
    
  }

  diff = computeDifference(u,w,n);

  return n_iter;
}


int main(int argc, char const *argv[]) {

  int n = atoi(argv[1]);
  int tol = atoi(argv[2]);

  poisson(n,tol);

  return 0;
}
