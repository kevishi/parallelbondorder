// Simple Atom class

#ifndef ATOM_H
#define ATOM_H

#include <vector>
#include <math.h>

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

class Bond;

class Atom {
 public:
 Atom(unsigned short element) : _atomicNum(element), _x(0.0), _y(0.0), _z(0.0), _id(0) {};
 Atom(unsigned short e, double x, double y, double z, unsigned int id) : _atomicNum(e), _x(x), _y(y), _z(z), _id(id) {};
  ~Atom() {}

  unsigned short atomicNum() { return _atomicNum; }
  void setAtomicNum(unsigned short element) { _atomicNum = element; }

  char* elementSymbol();

  void setXYZ(double x, double y, double z) { _x = x; _y = y; _z = z; }
  double x() { return _x; }
  double y() { return _y; }
  double z() { return _z; }

  double distanceFrom(Atom &a)
  { return sqrt( SQUARE(_x - a._x) + SQUARE(_y - a._y) + SQUARE(_z - a._z) ); }
  double distanceSqFrom(Atom &a)
  { return (SQUARE(_x - a._x) + SQUARE(_y - a._y) + SQUARE(_z - a._z)); }

  unsigned int numberOfBonds() { return _bonds.size(); }
  unsigned int numberOfDoubleBonds();
  void addBond(Bond *bond) { _bonds.push_back(bond); }

  bool isBonded(Atom *a);

  double radius();

  unsigned int id() { return _id; }
  void setId(unsigned int id) { _id = id; }

  /* std::vector<Bond*> getBonds() {return _bonds;} */

 protected:
  double _x;
  double _y;
  double _z;
  unsigned short _atomicNum;
  unsigned int _id;

  std::vector<Bond*> _bonds;
};

#endif
