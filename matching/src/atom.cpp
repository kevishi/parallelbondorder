#include "atom.h"
#include "bond.h"

using namespace std;

unsigned int Atom::numberOfDoubleBonds()
{
  unsigned int count = 0;

  // loop through the bonds
  vector<Bond*>::iterator it;
  for (it = _bonds.begin(); it < _bonds.end(); ++it) {
    if ((*it)->type() == 2)
      ++count;
  }

  return count;
}

bool Atom::isBonded(Atom *a)
{
  bool bonded = false;
  if (a == NULL)
    return false;

  // loop through the bonds
  vector<Bond*>::iterator it;
  for (it = _bonds.begin(); it < _bonds.end(); ++it) {
    if ((*it)->neighbor(this) == a) {
      bonded = true;
      break;
    }
  }

  return bonded;
}

double Atom::radius()
{ // covalent radii for the elements
  switch (_atomicNum) {
  case 1:
    return 0.31;
    break;
  case 6:
    return 0.76;
    break;
  case 7:
    return 0.71;
    break;
  case 8:
    return 0.66;
    break;
  default:
    return 1.0;
  }
}

char* Atom::elementSymbol()
{
  char *sym = new char[4];

  switch (_atomicNum) {
  case 1:
    sym = "H";
    break;
  case 6:
    sym = "C";
    break;
  case 7:
    sym = "N";
    break;
  case 8:
    sym = "O";
    break;
  default:
    sym = "C";
  }

  return sym;
}
