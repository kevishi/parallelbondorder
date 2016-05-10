// Write SDF molecule files

#include "molecule.h"
#include "atom.h"

#include <iostream>
#include <fstream>

using namespace std;

#define BUFF_SIZE 1024

bool writeSDF(Molecule &mol, const char* filename)
{
  char buff[BUFF_SIZE];
  ofstream sdFile;
  sdFile.open(filename);

  if (!sdFile.is_open()) {
    cout << "Can't open file " << filename << endl;
    return false;
  }

  sdFile << filename << "\n";
  sdFile << " OpenBabel04091613483D\n"; // time/date stamp
  sdFile << "\n"; // comment line
  snprintf(buff, BUFF_SIZE, "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n",
           mol.numberOfAtoms(), mol.numberOfBonds());
  sdFile << buff;

  std::vector<Atom*> atoms = mol.atoms();
  for (std::vector<Atom*>::iterator i = atoms.begin(); i < atoms.end(); ++i) {
    snprintf(buff, BUFF_SIZE, "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n",
             (*i)->x(), (*i)->y(), (*i)->z(), (*i)->elementSymbol(),
             0,0,0,0,0,0,0,0,0,0,0,0);
    sdFile << buff;
  } // end atom loop

  std::vector<Bond*> bonds = mol.bonds();
  for (std::vector<Bond*>::iterator i = bonds.begin(); i < bonds.end(); ++i) {
    Atom *a1 = (*i)->start();
    Atom *a2 = (*i)->end();
    snprintf(buff, BUFF_SIZE, "%3d%3d%3d  0  0  0  0\n",
             a1->id() + 1, a2->id() + 1, (*i)->type());
    sdFile << buff;
  } // end bond loop

  // footer
  sdFile << "M  END\n";
  sdFile << "$$$$" << endl; // endl will perform a flush

  return true;
}
