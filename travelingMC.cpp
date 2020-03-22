#include <iostream>
#include <random>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <limits>
#include <algorithm>
#include <iterator>
#include <stdio.h>

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

float findNearestRoute(float * dists, int x, int * towns, int n, int visited, int * nextTown){
  float distance = 0.0;
  int selectedTown = towns[0];
  float d = dist(dists,x,towns[0],n);
  float d2;
  int ntowns = n-visited;

  for(int i=1; i<ntowns; i++){
    d2 = dist(dists,x,towns[i],n);
    if(d2 < d){
      d = d2;
      selectedTown = towns[i];
    }
  }
  distance = d;
  *nextTown = selectedTown;
  return distance;
}

int * setDiff(int x, int * vec, int size){
  int index = 0;
  int * res = new int[size];
  for(int i=0; i<size+1; i++){
    if(vec[i] != x){
      res[index] = vec[i];
      index++;
    }
  }
  return res;
}

float travellingGreedy(float * dists, int * res, int n){
  float Tdist = 0.0;
  int * town = new int[n];

  town[0] = rand() % n;
  res[0] = town[0];

  int visited = 1;
  int * to_visit = new int[n];
  
  for(int i=0; i<n; i++){
    to_visit[i] = i;
  }

  to_visit = setDiff(town[0],to_visit,n-1);

  float distance;
  int * nextTown = new int;

  for(int i=1; i<n; i++){
    distance = findNearestRoute(dists,town[i-1],to_visit,n,visited,nextTown);
    visited++;
    town[i] = *nextTown;
    res[i] = town[i];
    Tdist += distance;
    to_visit = setDiff(town[i],to_visit,n-visited);
  }
  
  Tdist += dist(dists,town[0],town[n],n);
  /*
  printf("Rota greedy\n");
  for(int i=0; i<n; i++){
    printf("%d -",res[i]);
  }
  printf("\n");
  */
  return Tdist;
}

float tourDistance(float * dists, int * route, int size){
  float Tdist = 0.0;
  for(int i=0; i<size-1; i++){
    Tdist += dist(dists, route[i],route[i+1],size);
  }
  return Tdist;
}


int * TwoOPTSwap(int * route, int i, int k, int size){
  int * new_route = new int[size];
  for ( int c = 0; c <= i - 1; ++c )    {
      new_route[c] = route[c];
  }

  int dec = 0;
  for ( int c = i; c <= k; ++c )    {
    new_route[c] = route[k-dec];
    dec++;
  }

  for ( int c = k + 1; c < size; ++c )    {
    new_route[c] = route[c];
  }

  return new_route;
}

float TwoOPT(float * dists, int n){
  int improve = 0;
  int * new_route = new int[n];
  int * route = new int[n];
  float best_distance = 0.0;
  float startingDistance = travellingGreedy(dists, route, n);

  while (improve < 20){

    best_distance = tourDistance(dists, route, n);
    for ( int i = 0; i < n - 1; i++ ){
      for ( int k = i + 1; k < n; k++){
        new_route = TwoOPTSwap(route, i, k, n);
        float new_distance = tourDistance(dists, new_route, n);

        if ( new_distance < best_distance ){
          improve = 0;
          route = new_route;
          best_distance = new_distance;
        }
      }
    }
    improve ++;
  }
  /*
  for(int i=0; i<n; i++){
    printf("%d -",route[i]);
  }
  printf("\n");
  */
  return best_distance;
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
    }
    else{
      i++;
    }
  }

  return Tdist;
}

float travelingSA(float * dists, int * res,int n){
  float Tdist = 0.0;

  randperm(res,n);
  Tdist = dist(dists,res[n-1],res[0],n);
  for (size_t i = 0; i < n-1; i++) {
    Tdist += dist(dists,res[i],res[i+1],n);
  }

  float randV;
  float temp = 1.0;
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

    randV = rand();
    float delta =
      dist(dists,res[previous],res[next1],n) +
      dist(dists,res[rand_index],res[next2],n) -
      dist(dists,res[previous],res[rand_index],n) -
      dist(dists,res[next1],res[next2],n);

    if (delta < 0 | exp(-delta/temp) >= randV) {
      swap(res,rand_index,next1);
      Tdist += delta;
      if (delta != 0.0)
        i = 0;
    }
    else{
      i++;
    }
    temp *= 0.999;
  }

  return Tdist;
}


int main(int argc, char const *argv[]) {

  srand(time(NULL));
  int n = atoi(argv[1]);
  int procs = atoi(argv[2]);
  int * res = new int[n];

  int * res_greedy = new int[n];
  int * res_twoOPT = new int[n];
  int * res_mc = new int[n];
  int * res_sa = new int[n];

  float minDist_twoOPT = std::numeric_limits<float>::max();
  float minDist_greedy = std::numeric_limits<float>::max();
  float minDist_mc = std::numeric_limits<float>::max();
  float minDist_sa = std::numeric_limits<float>::max();
  float tdist_ac = 0.0;
  float tdist_ac2 = 0.0;

  float * dists = new float [n*n];
  generateDists(dists, n);

  for (int i = 0; i < procs; i++) {
    tdist_ac = travellingGreedy(dists,res,n);

    if (tdist_ac < minDist_greedy){
      minDist_greedy = tdist_ac;
      memcpy(res,res_greedy,sizeof(int)*n);
    }
  }

  for (int i = 0; i < procs; i++) {
    tdist_ac2 = TwoOPT(dists,n);
    if (tdist_ac < minDist_twoOPT){
      minDist_twoOPT = tdist_ac2;
      memcpy(res,res_twoOPT,sizeof(int)*n);
    }
  }

  for (int i = 0; i < procs; i++) {
    tdist_ac = travelingMC(dists,res,n);

    if (tdist_ac < minDist_mc){
      minDist_mc = tdist_ac;
      memcpy(res,res_mc,sizeof(int)*n);
    }
  }

  for (int i = 0; i < procs; i++) {
    tdist_ac = travelingSA(dists,res,n);

    if (tdist_ac < minDist_sa){
      minDist_sa = tdist_ac;
      memcpy(res,res_sa,sizeof(int)*n);
    }
  }

  std::cout << "2-opt: " << minDist_twoOPT << std::endl;
  std::cout << "Greedy: " << minDist_greedy << std::endl;
  std::cout << "MC: " << minDist_mc << std::endl;
  std::cout << "SA: " << minDist_sa << std::endl;

  return 0;
}
