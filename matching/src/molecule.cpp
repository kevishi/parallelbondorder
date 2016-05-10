// Simple Molecule class

#include "molecule.h"
#include <algorithm>
#include <iostream>
#include <set>
#include <map>
#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>
#include "cycleTimer.h"
#include <omp.h>
using namespace std;
extern int n;
typedef pair<int,int> edge;
typedef vector<int>* graph;

bool modified = 1;
int totalcount = 0;
int contentioncount = 0;
state *currentState;
struct pair_hash {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return v.first*31+v.second;
    }
};

bool sortAtomZ(const pair<Atom*,double> &a, const pair<Atom*,double> &b)
{   
  return (a.second < b.second); 
}

void Molecule::perceiveBonds() {
  // find the initial bonded atoms
  unsigned int max;
  double maxrad = 0.0;
  bool unset = false;
  Atom *atom,*nbr;
  vector<Atom*>::iterator i;
  vector<pair<Atom*,double> > zsortedAtoms;
  vector<double> rad;

  rad.resize(_atoms.size());

  for (i = _atoms.begin(); i < _atoms.end(); ++i) {
    pair<Atom*,double> entry(*i, (*i)->z());
    zsortedAtoms.push_back(entry);
  }
  sort(zsortedAtoms.begin(), zsortedAtoms.end(), sortAtomZ);
  max = zsortedAtoms.size();

  // find the maximum radius for the atoms
  for (unsigned int j = 0; j < max; ++j) {
    atom   = zsortedAtoms[j].first;
    rad[j] = atom->radius();
    maxrad = std::max(rad[j],maxrad);
  }

  int idx1, idx2;
  double d2,cutoff,zd;
  for (unsigned int j = 0; j < max ; ++j) {
    double maxcutoff = SQUARE(rad[j]+maxrad+0.45);
    atom = zsortedAtoms[j].first;
    for (unsigned int k = j + 1 ; k < max ; k++ ) {
      nbr = zsortedAtoms[k].first;

      // bonded if closer than elemental Rcov + tolerance
      cutoff = SQUARE(rad[j] + rad[k] + 0.45);
      // current z distance
      zd  = SQUARE(atom->z() - nbr->z());
      // bigger than max cutoff, which is determined using largest radius,
      // not the radius of k (which might be small, ie H, and cause an early  termination)
      // since we sort by z, anything beyond k will also fail
      if (zd > maxcutoff )
        break;

      d2  = SQUARE(atom->x() - nbr->x());
      if (d2 > cutoff)
        continue; // x's bigger than cutoff
      d2 += SQUARE(atom->y() - nbr->y());
      if (d2 > cutoff)
        continue; // x^2 + y^2 bigger than cutoff
      d2 += zd;
      if (d2 > cutoff)
        continue; // too far
      if (d2 < 0.16) // 0.4 * 0.4 = 0.16
        continue; // too close

      if (atom->isBonded(nbr))
        continue; // already handled

      addBond(atom,nbr,1);

    } // end inner loop
  } // end outer loop
} // end perceiveBonds

  // remove N-H bonds and x-N bonds where the N is connected to an H bond
void Molecule::preprocess(vector<Bond*>& relBonds, vector<Atom*>& relAtoms){
  vector<Bond*> unknowns;
  std::set<Atom*> no;
  set<Atom*> atoms;
  vector<Bond*>::iterator it;
  for(it = _bonds.begin(); it != _bonds.end(); it++){
    Bond* b = *it;
    if(b->start()->atomicNum() == 1 && b->end()->atomicNum() == 7){ //H-N
      no.insert(b->end());
    }else if(b->start()->atomicNum() == 7 && b->end()->atomicNum() == 1){ //N-H
      no.insert(b->start());
    } else if(b->start()->atomicNum() == 8){ //x-O
      no.insert(b->end());
    }else if(b->start()->atomicNum() == 8){ //O-x
      no.insert(b->end());
    } else{
      unknowns.push_back(b);
    }
  }

  for(it = unknowns.begin(); it != unknowns.end(); it++){
    Bond* b = *it;
    if(!(no.find(b->start()) != no.end() || no.find(b->end()) != no.end()) &&
       !(b->start()->atomicNum() == 1 || b->end()->atomicNum() == 1)){ //x-H, H-x
      relBonds.push_back(b);
      atoms.insert(b->start());
      atoms.insert(b->end());
    }
  }

  set<Atom*>::iterator it1;
  for(it1 = atoms.begin(); it1 != atoms.end(); it1++)
    relAtoms.push_back(*it1);

}




void unlockPath(std::vector<int> nodes, int end) {
  for(int i = 0; i<nodes.size(); i++){
    currentState->owner[nodes[i]] = -1;
    omp_unset_lock(&(currentState->lock[nodes[i]]));
  }
}

bool lockPath(std::vector<int> nodes) {
  int procId = omp_get_thread_num();
  for(int i = 0; i<nodes.size(); i++) {
    int u = nodes[i];
    totalcount++;
    if (omp_test_lock(&(currentState->lock[u]))) {
      currentState->owner[u] = procId;
      continue;
    } else {
      contentioncount++;
      while (procId > currentState->owner[u]) {
        if (omp_test_lock(&(currentState->lock[u]))) {
          currentState->owner[u] = procId;
          continue;
        }
      }
      if (procId < currentState->owner[u]) {
        unlockPath(nodes, u);
        return false;
      }
    }
  }
  return true;
}

bool findHead(int even, std::vector<int>& head_rev, int avoid) {
  if(even == avoid){
    return false;
  }

  if(currentState->matching[even] == -1){
    return true;
  }

  int a = currentState->matching[even];
  head_rev.push_back(a);

  for(int i = 0; i < n; i++){
    if (currentState->forest[a][i] && currentState->parity[i] == 2 && even != i) {
      head_rev.push_back(i);
      if(findHead(i,head_rev,avoid)){
        return true;
      }
      head_rev.pop_back();
    }
  }

  return false;
}

bool findTail(int prev, int even, std::vector<int>& path){
  if (currentState->matching[even] == -1) {

    path.push_back(even);
    path.push_back(prev);
    return true;
  }
  int a = currentState->matching[even];
  for(int i = 0; i < n; ++i){
    if (currentState->forest[a][i] && currentState->parity[i] == 2 && even != i) {
      if(findTail(a,i,path)){

        path.push_back(even);
        path.push_back(prev);
        return true;
      }
    }
  }

  return false;
}


bool checkPath(std::vector<int> path){
  if(currentState->matching[path[0]] != -1) {
    return false;
  }
  int old = -1;
  for(int i = 1; i<path.size()-1; i++) {
    old = path[i];
    i++;
    if(currentState->matching[old] != path[i] || currentState->matching[path[i]] != old)
      return false;
    if(i >= path.size()-1)
      return false;
  }
  if(currentState->matching[path[path.size()-1]] != -1)
    return false;
  return true;
}


bool findPath(int u, int v, std::vector<int>& path){
  // u - v -> ... -> root(v)
  std::vector<int> tail;
  if(!findTail(u,v,tail))
    return false;
  reverse(tail.begin(),tail.end());

  // u' -> ... -> root(u)
  int last = tail[tail.size() - 1];
  std::vector<int> head_rev;
  if(!findHead(u,head_rev,last))
    return false;

  //root(u) -> ... -> u -> v -> ... -> root(v)
  path = head_rev;
  reverse(path.begin(), path.end());
  path.insert(path.end(), tail.begin(), tail.end());

  if(!lockPath(path)){
    path.clear();
    return false;
  }

  if(checkPath(path)){
    return true;
  }else{
    unlockPath(path,-1);
    path.clear();
    return false;
  }
}

bool augmentPath(std::vector<int> nodes) {
  int old = nodes[0];
  currentState->matching[old] = nodes[1];
  currentState->matching[nodes[1]] = old;
  currentState->parity[old] = 0;
  currentState->parity[nodes[1]] = 0;
  currentState->owner[old] = -1;
  currentState->owner[nodes[1]] = -1;
  omp_unset_lock(&(currentState->lock[old]));
  omp_unset_lock(&(currentState->lock[nodes[1]]));
  int temp;
  int i = 1;
  while(i<nodes.size()-1){
    temp = nodes[i+1];
    i++;
    old = nodes[i];
    temp = nodes[i+1];
    i++;
    currentState->matching[old] = nodes[i];
    currentState->parity[old] = 0;
    currentState->matching[nodes[i]] = old;

    currentState->parity[nodes[i]] = 0;
    currentState->owner[old] = -1;
    currentState->owner[nodes[i]] = -1;
    omp_unset_lock(&(currentState->lock[old]));
    omp_unset_lock(&(currentState->lock[nodes[i]]));
  }
  return true;
}
int min(int x, int y) {
  return x<y?x:y;
}

void solve_justpre(state *s) {
  currentState = s;
  modified = 1;
  totalcount = 0;
  contentioncount = 0;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && s->graph[i][j]
            && s->matching[j] == -1
            && s->matching[i] == -1) {
          s->matching[i] = j;
          s->matching[j] = i;
      }
    }
  }
  return;
  
}

void solve_no_preprocess(state *s) {
  currentState = s;
  modified = 1;
  totalcount = 0;
  contentioncount = 0;
/*
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && s->graph[i][j]
            && s->matching[j] == -1
            && s->matching[i] == -1) {
          s->matching[i] = j;
          s->matching[j] = i;
      }
    }
    }*/

  while (modified) {
    for (int i = 0; i < n; i++) {
      s->parity[i] = 2*(s->matching[i] == -1);
    }
    memset(s->forest, 0, sizeof(s->forest));
    modified = 0;
    bool flag = 1;
    while (flag && !modified) {
      flag = 0;
#pragma omp parallel for schedule(dynamic)
      for (int v = 0; v < n; v++) {

        if (!modified && s->parity[v] ==2) {
          for (int w = 0; w < n; w++) {
            if (!modified && s->graph[v][w] && s->parity[w] != 1) {
              if(currentState->parity[w] < 1) {
                int x = currentState->matching[w];
                int low = min(w,x);

                if (omp_test_lock(&(currentState->lock[low]))) {
                  if(currentState->parity[w] < 1) {
                    int x = currentState->matching[w];
                    currentState->forest[v][w] = 1;
                    currentState->forest[w][v] = 1;
                    currentState->parity[w] = 1;
                    currentState->parity[x] = 2;
                    currentState->forest[w][x] = 1;
                    currentState->forest[x][w] = 1;
                  }
                  omp_unset_lock(&(currentState->lock[low]));
                }
              } else if (currentState->parity[w] == 2) {
                std::vector<int> path;
                if (findPath(v, w,path)) {
                  modified = augmentPath(path);
                }
              }
              flag = 1;
            }
          }
        }
      }
    }
  }
  
  return;
  
}

void solve(state *s) {
  currentState = s;
  modified = 1;

  totalcount = 0;
  contentioncount = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && s->graph[i][j]
            && s->matching[j] == -1
            && s->matching[i] == -1) {
          s->matching[i] = j;
          s->matching[j] = i;
      }
    }
    }

  while (modified) {
    for (int i = 0; i < n; i++) {
      s->parity[i] = 2*(s->matching[i] == -1);
    }
    memset(s->forest, 0, sizeof(s->forest));
    modified = 0;
    bool flag = 1;
    while (flag && !modified) {
      flag = 0;
#pragma omp parallel for schedule(static)
      for (int v = 0; v < n; v++) {

        if (!modified && s->parity[v] ==2) {
          for (int w = 0; w < n; w++) {
            if (!modified && s->graph[v][w] && s->parity[w] != 1) {
              if(currentState->parity[w] < 1) {
                int x = currentState->matching[w];
                int low = min(w,x);

                if (omp_test_lock(&(currentState->lock[low]))) {
                  if(currentState->parity[w] < 1) {
                    int x = currentState->matching[w];
                    currentState->forest[v][w] = 1;
                    currentState->forest[w][v] = 1;
                    currentState->parity[w] = 1;
                    currentState->parity[x] = 2;
                    currentState->forest[w][x] = 1;
                    currentState->forest[x][w] = 1;
                  }
                  omp_unset_lock(&(currentState->lock[low]));
                }
              } else if (currentState->parity[w] == 2) {
                std::vector<int> path;
                if (findPath(v, w,path)) {
                  modified = augmentPath(path);
                }
              }
              flag = 1;
            }
          }
        }
      }
    }
  }
  return;
}

void Molecule::doMatching() {
  // preprocessing

  vector<Bond*> relBonds;
  vector<Atom*> relAtoms;
  preprocess(relBonds, relAtoms);
  //matching
  // loop through the bonds
  // srand(time(NULL));
  vector<Bond*>::iterator it;
  std::map<Atom*, int> indexMap;

  for(int i = 0; i < relAtoms.size(); i++){
    indexMap[relAtoms[i]] = i;
  }
  unordered_map<edge,Bond*,pair_hash> bondMap;
  for(auto b: relBonds){
    edge e = make_pair(0,0);
    bondMap[e] = b;
    bondMap[make_pair(indexMap[b->start()],indexMap[b->end()])] = b;
    bondMap[make_pair(indexMap[b->end()],indexMap[b->start()])] = b;
  }

  graph G = new vector<int>[relAtoms.size()];
  for(auto e: relBonds){
    G[indexMap[e->start()]].push_back(indexMap[e->end()]);
    G[indexMap[e->end()]].push_back(indexMap[e->start()]);
  }

  //Assign weights w



  n = relAtoms.size();
  state *s = new state;
  memset(s, 0, sizeof(s));

  // s.matching = new short[n]();
  // s.graph = new bool*[n];
  // s.forest = new bool*[n];
  // for(int i = 0; i < n; i++){
    // s.graph[i] = new bool[n]();
    // s.forest[i] = new bool[n]();
  // }

  // s.removed = new bool[n]();
  // s.parent = new int[n]();
  // s.parity = new int[n]();

  // s.graph = G;
  // std::cout<<std::endl;

  for(int i = 0; i< n; i++) {
    omp_init_lock(&s->lock[i]);
  }
  for(it = relBonds.begin(); it < relBonds.end(); ++it) {
    int a = indexMap[(*it)->start()];
    int b = indexMap[(*it)->end()];
    int i = a < b ? a : b;
    int j = a < b ? b : a;
    s->graph[a][b] = 1;
    s->graph[b][a] = 1;
  }
  // for(int i = 0; i<relAtoms.size(); i++) {
       // std::cout << "["<<s.parity[i] <<"]";

  // }
  // std::cout<<std::endl;
  //memcpy(s.graph, adj, sizeof(bool) * sizeof(adj)*sizeof(adj));
  memset(s->matching, -1, sizeof(s->matching));
  memset(s->parity, 0, sizeof(s->parity));



  solve(s);
  //for (int i =  0; i < n; i++) {
  //  if (i % 5 == 0) {
  //    printf("\n");
  //  }
  //  printf("%d: %d\t", i, s.matching[i]);
  //}
  //solve_eba(s);


  for(it = relBonds.begin(); it != relBonds.end(); ++it){
    int a = indexMap[(*it)->start()];
    int b = indexMap[(*it)->end()];
    int i = a < b ? a : b;
    int j = a < b ? b : a;
    if(s->matching[i] == j ) {
      //std::cout << "Matching Found!: Bond:" << i << ","<<j << " Matching: "<<i<<","<<s.matching[i]<<std::endl;
      (*it)->setType(2);
    }
  }

  // for (it = _bonds.begin(); it < _bonds.end(); ++it)
    // std::cout << (*it)->start()->elementSymbol() << "(" << indexMap[(*it)->start()] << ":" << (*it)->start()<< ") - " << (*it)->end()->elementSymbol()  << "(" << indexMap[(*it)->end()] << ":"<< (*it)->end() << ") : " << (*it)->type() << std::endl ;
  // std::cout << "-----" << std::endl
// ;

  // SparseLU<SparseMatrix<scalar, ColMajor>, COLAMDOrdering<Index> >;
}
