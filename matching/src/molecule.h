// Simple Molecule class

#ifndef MOLECULE_H
#define MOLECULE_H

#define MAX_N 20000
#include <vector>
#include <omp.h>
//#include <unordered_map>

#include "atom.h"
#include "bond.h"
struct state {
  omp_lock_t lock[MAX_N];
  short owner[MAX_N];
  
  short matching[MAX_N];
  bool graph[MAX_N][MAX_N];
  bool forest[MAX_N][MAX_N];
  int parent[MAX_N];
  int parity[MAX_N];
  
};

void solve(state* s);
void solve_no_preprocess(state* s);
void solve_justpre(state* s);

class Molecule {
 public:
 Molecule() {};
 ~Molecule() {}

 void addAtom(Atom *atom) { _atoms.push_back(atom); }
 void addAtom(unsigned int element, double x, double y, double z)
 { Atom *a = new Atom(element, x, y, z, _atoms.size()); _atoms.push_back(a); }
 void addBond(Bond *bond) { _bonds.push_back(bond); }
 void addBond(Atom *a, Atom *b, unsigned short type = 0)
{
  Bond *bond = new Bond(a, b); bond->setType(type); _bonds.push_back(bond);
  a->addBond(bond);
  b->addBond(bond);
}

 unsigned int numberOfAtoms() { return _atoms.size(); }
 unsigned int numberOfBonds() { return _bonds.size(); }

 std::vector<Atom*> atoms() { return _atoms; }
 std::vector<Bond*> bonds() { return _bonds; }

 void clear() { _atoms.clear(); _bonds.clear(); }

 void perceiveBonds();
 void doMatching();

 protected:
  std::vector<Atom*> _atoms;
  std::vector<Bond*> _bonds;

 private:
  void preprocess(std::vector<Bond*>& relBonds, std::vector<Atom*>& relAtoms);
  /* void preprocess(std::unordered_map<Atom*, std::vector<Bond*>>& adjList); */
};

#endif
