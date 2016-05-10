#include <sstream>
#include <fstream>
#include <iostream>
#include "cycleTimer.h"
#include "../src/eba.h"

#define MAX_N 500

int main(int argc, char** argv){
  if(argc != 2){
    std::cerr << "gib meh 1 file pls" << std::endl;
    return 1;
  }
  std::ifstream file(argv[1]);
  std::string line;

  getline(file,line);
  std::istringstream is(line);

  int s;
  is >> s;
  // bool **graph = new bool*[s];


  solution_t* solution = new solution_t[s];
  for(int i = 0; i < s; i++){
    // graph[i] = new bool[s];
    getline(file,line);
    is.clear();
    is.str(line);
    for(int j = 0; j < s; j++){
      int x;
      is >> x;
      if(x == 0)
        solution->graph[i][j] = 0;
      else
        solution->graph[i][j] = 1;
      std::cout << solution->graph[i][j] << " ";
    }
    std::cout << std::endl;
  }

  memset(solution, 0, sizeof(solution));
  memset(solution->matching, -1, sizeof(solution->matching));
  //dont think this is necesary
  // memset(solution->forest, 0, sizeof(solution->forest));
  // memset(solution->removed, 0, sizeof(solution->removed));
  // memset(solution->parity, 0, sizeof(solution->parity));
  // memset(solution->owner, 0, sizeof(solution->owner));
  for(int i = 0; i < s; i++){
    omp_init_lock(&solution->lock[i]);
  }
  omp_init_lock(&solution->blist_lock);


  double startTime = CycleTimer::currentSeconds();
  solve_eba(solution);
  double endTime = CycleTimer::currentSeconds();

  int* M = new int[s]();


  printf("Time taken:  %.4f ms\n", 1000.f * (endTime - startTime));

  return 0;
}
