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

Several target molecules are included through GitHub, ranging in size from 6-10 atoms (tiny) to 500+ atoms.

## Testing

The src/ code uses CMake and a C++ compiler.
- mkdir build
- cd build
- cmake ..
- make -j3
- ./pmatch ../molecules/*.xyz
