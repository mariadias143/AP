
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <string.h>
#include <limits>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <functional>
#include <stdlib.h>
#include <sys/time.h>

using namespace std;

#define TIME_RESOLUTION 1000000

long long unsigned initial_time;
timeval t;

void start (void) {
    gettimeofday(&t, NULL);
    initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
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
    
  int max = *max_element(diff, diff + n*n);
  return max;
}

void compare(int* w, int* u, int n){
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++){
      if(w[i*n + j] != u[i*n + j]){
        std::cout << "i:" << i << " j:" << j << " w:" << w[i*n + j] << " u:" << u[i*n + j] << std::endl;
        exit(0);
      }
    }
}

int poissonSeq(int n, int tol, double p, int * w, int * u){
	long long unsigned tt;

  int diff = tol + 1;
  int n_iter = 0;
  while(diff > tol){
    std::copy(w, w + n*n, u); // u = w;
    //compare(w,u,n);
    
    // Update red
    for(int i=1; i<n-1; i++)
      for(int j=1+remainder(i,2); j<n-1; j+=2)
        w[i*n + j] = (1-p)*w[i*n + j] + p * ( w[(i-1)*n+j] + w[i*n+(j-1)] + w[(i+1)*n+j] + w[i*n+(j+1)] ) / 4;

    // Update black
    for(int i=1; i<n-1; i++)
      for(int j=1+remainder(i+1,2); j<n-1; j+=2)
        w[i*n + j] = (1-p)*w[i*n + j] + p * ( w[(i-1)*n+j] + w[i*n+(j-1)] + w[(i+1)*n+j] + w[i*n+(j+1)] ) / 4;

    n_iter++;
    diff = computeDifference(u,w,n);
  }

  return n_iter;
}


int poissonOMP(int n, int tol, double p, int * w, int * u, int n_threads){
	long long unsigned tt;

  int diff = tol + 1;
  int n_iter = 0;
  while(diff > tol){
    std::copy(w, w + n*n, u); // u = w;
    //compare(w,u,n);
    
    // Update red
    #pragma omp parallel num_threads(n_threads)
    {
    #pragma omp for
    for(int i=1; i<n-1; i++)
      for(int j=1+remainder(i,2); j<n-1; j+=2)
        w[i*n + j] = (1-p)*w[i*n + j] + p * ( w[(i-1)*n+j] + w[i*n+(j-1)] + w[(i+1)*n+j] + w[i*n+(j+1)] ) / 4;

    // Update black
    #pragma omp for
    for(int i=1; i<n-1; i++)
      for(int j=1+remainder(i+1,2); j<n-1; j+=2)
        w[i*n + j] = (1-p)*w[i*n + j] + p * ( w[(i-1)*n+j] + w[i*n+(j-1)] + w[(i+1)*n+j] + w[i*n+(j+1)] ) / 4;

    }
    n_iter++;
    diff = computeDifference(u,w,n);
  }

  return n_iter;
}


int main(int argc, char const *argv[]) {
  long long unsigned tt;

  int n = atoi(argv[1]);
  int tol = 5; 
  int nthreads = atoi(argv[2]);

  //Parâmetro de relaxação
  double p = 2/(1+sin(M_PI/(n-1)));

  int * w = new int[n*n]; // values of current iteration
  int * u = new int[n*n]; // last iteration values

  // Temperature is zero at the top and 100 on the other boundaries
  for(int i=0; i<n; i++){
    w[i] = 0; //above border
    w[i + n*(n-1)] = 100; //lower border
  }

  for(int j=0; j<n; j++){
    w[j*n] = 100; //left border
    w[n-1 + j*n] = 100; //right border
  }

  for(int i=1; i<(n-1); i++){
    for(int j=1; j<(n-1); j++){
      w[i*n + j] = 50; //initial aproximation for interior points
    }
  }

  start();
  int res1 = poissonSeq(n,tol,p,w,u);
  tt = stop();
  cout << "Sequencial: " << tt << " usecs" << endl;

  start();
  int res2 = poissonOMP(n,tol,p,w,u,nthreads);
  tt = stop();
  cout << "OMP: " << tt << " usecs" << endl;


  return 0;
}
