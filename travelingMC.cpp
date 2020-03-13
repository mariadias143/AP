
#include <iostream>
#include <random>
#include <cstdlib>
#include <math.h>

float dist(float * dists,int origin,int destination,int n_cities){
  return  dists[origin*n_cities + destination];
}

void randperm(int * res, int n){
  int hash_table [n] = {0};
  int n_missing = 0;

  while(n_missing < n){
    int rand_index;
    bool hasRand = false;

    while (!hasRand) {
      rand_index = rand() % n;
      if (hash_table[rand_index] == 0)
        hasRand = true;
    }
    res[n_missing++] = rand_index;
    hash_table[rand_index] = 1;
  }
}

void generateDists(float * dists, int n){
  float * x_coords = new float[n];
  float * y_coords = new float[n];

  for (int i = 0; i < n; i++) {
    x_coords[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/10.0));
    y_coords[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/10.0));
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j){
        dists[i*n + j] = -1.0;
      }
      else{
        dists[i*n + j] = sqrt( pow(x_coords[j] - x_coords[i], 2.0) + pow(y_coords[j] - y_coords[i], 2.0));
      }
    }
  }
}

void swap(int * res, int i, int j){
  int temp = res[i];
  res[i] = res[j];
  res[j] = temp;
}

float travelingMC(float * dists, int * res,int n){
  float Tdist = 0.0;

  randperm(res,n);
  Tdist = dist(dists,res[n-1],res[0],n);
  for (size_t i = 0; i < n-1; i++) {
    Tdist += dist(dists,res[i],res[i+1],n);
  }

  int i = 0;
  int rand_index,previous,next1,next2;
  while (i < 100) {
    rand_index = rand() % n;
    if (rand_index == 0){
      previous = n-1; next1 = 1; next2 = 2;
    }
    else if (rand_index == n-2){
      previous = n-3; next1 = n-1; next2 = 0;
    }
    else if (rand_index == n-1){
      previous = n-2; next1 = 0; next2 = 1;
    }
    else{
      previous = rand_index-1; next1 = rand_index+1; next2 = rand_index+2;
    }

    float delta =
      dist(dists,res[previous],res[next1],n) +
      dist(dists,res[rand_index],res[next2],n) -
      dist(dists,res[previous],res[rand_index],n) -
      dist(dists,res[next1],res[next2],n);

    if (delta < 0){
      swap(res,rand_index,next1);
      Tdist += delta;
      i = 0;
      std::cout << "entrei" << std::endl;
    }
    else{
      i++;
    }
  }

  return Tdist;
}


int main(int argc, char const *argv[]) {

  srand(time(NULL));
  int n = 100;
  int * res = new int[n];

  /*
  randperm(res,n);

  for (int i = 0; i < n; i++) {
    std::cout << res[i] << std::endl;
  }*/
  float * dists = new float [n*n];
  generateDists(dists, n);
  travelingMC(dists,res,n);

  return 0;
}
