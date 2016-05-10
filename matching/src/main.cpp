// Simple main.cpp
// for perfect matching of bonds from XYZ files

#include <vector>
#include <string>

#include <iostream>
#include <cstring>

#include "tokenize.h"
#include "molecule.h"
#include "atom.h"
#include "bond.h"

bool readXYZ(Molecule &mol, const char* filename);
bool writeSDF(Molecule &mol, const char* filename);

using namespace std;

int n = 1;

int main (int argc, char *argv[])
{
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " {filename}" << endl;
    return 1;
  }

  Molecule mol;
  // loop through the filenames
  for (unsigned int a = 1; a < argc; ++a) {
    if (!readXYZ(mol, argv[a]))
      cout << "Cannot read the XYZ file" << endl;

    mol.perceiveBonds();

    cout << "Molecule has " << mol.numberOfAtoms() << " atoms and " << mol.numberOfBonds() << " bonds." << endl;
    
    mol.doMatching();

    std::vector<Atom*> atoms = mol.atoms();
    unsigned int j = 0;
    for (std::vector<Atom*>::iterator i = atoms.begin(); i < atoms.end(); ++i, ++j) {
      if ((*i)->atomicNum() == 6 && (*i)->numberOfDoubleBonds() != 1) {
        // cout << " failed matching on atom " << j << endl;
        cout << " failed matching on atom " << *i << endl;
        cout << " has " << (*i)->numberOfDoubleBonds() << " double bonds" << endl;
        break;
      }
    } // end check loop

    // write an SD output
    char *filename = argv[a];
    // change extension
    char *pExt = strrchr(filename, '.');
    if (pExt != NULL)
      strcpy(pExt, ".sdf");
    else
      strcat(filename, ".sdf");

    if (!writeSDF(mol, filename))
      cout << "Cannot write the SDF file" << endl;

  } // end (loop through command-line args)

  return 0;
}
