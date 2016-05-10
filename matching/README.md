# Parallel Bond Order Assignment by Graph Theory

In chemistry and molecular biology, one often of a molecule as a graph, with atoms as the vertices and bonds as the bidirectional edges. Usually, we prefer to visualize different bond types: single, double, triple, etc. as a way of understanding their reactivity.

![matching](/images/matching.png)

In many cases, however, we are not given the edges of the graphs or their types. We must *perceive* them. Determining the connectivity is usually easy. Bonded atoms are close to each other, so simple distance constraints allow us to quickly find vertices.

Assigning bond types requires finding a [perfect weighted match](https://en.wikipedia.org/wiki/Matching_%28graph_theory%29) on the graph, assigning single and double bonds such that:

* All carbon atoms have only one double bond
* Pyrrole NH groups have no double bonds
* Pyridine N groups have only one double bond
* Hydrogens are ignored
* Edges out of a ring to an oxygen are assumed to be double bonds

In this project, you will optimize bond order assignment (aka Kekule assignment), preparing a correct implementation of an efficient parallel perfect matching algorithm. The baseline implementations in existing open source cheminformatics libraries are correct, but are exponential in the number of fused rings (cycles). Initially, this wasn't a problem because most molecules are simple. Then someone tried a nanotube with hundreds of rings. Oops.

[Labute](papers/Labute-2005.pdf) suggests a maximum matching algorithm, but the method is serial and too complicated anyway. Perfect matching algorithms can, of course, be parallel (e.g., [Chawla](papers/Chawla-2004.pdf) or [Mulmuley](papers/Mulmuley-1987.pdf)).

Several target molecules are included through GitHub, ranging in size from 6-10 atoms (tiny) to 500+ atoms.

Please contact me (geoffh@pitt.edu) if you have any questions about this.

##Deliverables:

- You can implement the matching algorithm in any language, although C++ is strongly preferred.
- The implementation should be highly parallel. OpenMP is strongly preferred.
- You must be willing to release your final code under a GPLv2, LGPL, or BSD (preferred) open source license

## Testing

The src/ code uses CMake and a C++ compiler.
- mkdir build
- cd build
- cmake ..
- make -j3
- ./pmatch ../molecules/*.xyz
