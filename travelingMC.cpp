
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
  int town[n] = { 0 };
  
  town[0] = rand() % n; 
  printf("Primeira cidade: %d\n",town[0]);

  int visited = 1;
  int * to_visit = new int[n];
  for(int i=0; i<n; i++){
    to_visit[i] = i;
  }

  to_visit = setDiff(town[0],to_visit,n-1);

  for(int i=0; i<n-1; i++){
    printf("%d -",to_visit[i]);
  }
  printf("\n");

  float distance;
  int * nextTown = new int;

  for(int i=1; i<n-1; i++){
    distance = findNearestRoute(dists,town[i-1],to_visit,n,visited,nextTown);
    printf("Next town is %d at a distance of %lf\n",*nextTown, distance);
    visited++;
    town[i] = *nextTown;
    Tdist += distance;
    to_visit = setDiff(town[i],to_visit,n-visited);
  }
  Tdist += dist(dists,town[n],town[0],n);
  return Tdist;
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


  srand(1);
  int n = 10;
  int * res = new int[n];



  /*
  randperm(res,n);

  for (int i = 0; i < n; i++) {
    std::cout << res[i] << std::endl;
  }*/

  
  float * dists = new float [n*n];
  generateDists(dists, n);

  travellingGreedy(dists,res,n);
  

  return 0;
}
