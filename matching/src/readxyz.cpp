// How to read XYZ molecule files

#include "molecule.h"
#include "atom.h"
#include "tokenize.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

bool readXYZ(Molecule &mol, const char* filename)
{
  ifstream xyzFile;
  xyzFile.open(filename);

  if (!xyzFile.is_open()) {
    cout << "Can't open file " << filename << endl;
    return false;
  }

  // file opens OK
  mol.clear();

  string buffer;
  vector<string> vs;
  getline(xyzFile, buffer);  // parse # of atoms
  getline(xyzFile, buffer);  // comment-line, ignore
  while( getline(xyzFile, buffer) ) {
    vs = tokenize(buffer);
    if (vs.size() < 4)
      break; // incomplete line for some reason

    int element = 0;
    double x, y, z;
    switch (vs[0].c_str()[0]) {
    case 'H':
      element = 1;
      break;
    case 'C':
      element = 6;
      break;
    case 'N':
      element = 7;
      break;
    case 'O':
      element = 8;
      break;
    case 'S':
      element = 16;
      break;
    default: // leave as element 0, something is strange
      element = 0;
    }

    x = atof( vs[1].c_str() );
    y = atof( vs[2].c_str() );
    z = atof( vs[3].c_str() );

    mol.addAtom(element, x, y, z);
  } // end while loop (reading atoms)
  // done

  // success if there was actually something in the file
  return (mol.numberOfAtoms() > 0);
}
