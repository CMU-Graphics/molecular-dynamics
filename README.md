# molecular-dynamics

![Molecule visualization](interface.png)

Skeleton code and visualization for a basic molecular dynamics simulator, provided as part of the 15-462/662 final exam at CMU.  Visualization via [Polyscope](http://polyscope.run/).

The code will build and run as-is, but to make the simulator work you will have to implement the routines marked `TODO` in `src/MolecularDynamics.cpp`.  The correct behavior for these routines is described in the exam found in `15462_Final_FA2021.pdf`.

Several test files can be found in the `data/` subdirectory.

Note that you do _not_ have to build and run this version of the code for the exam (which can be more difficult to compile, due to the visualization routines).  It is provided "just for fun."  See Piazza for the official exam code.

**Disclaimer**: _This simulator is an extremely simplified version of a real molecular dynamics simulator, intended purely for pedagogical purposes.  It should not be used for serious scientific work._

### Build the code

**Unix-like machines**: configure (with cmake) and compile
```
cd molecular-dynamics
mkdir build
cd build
cmake ..
make -j6
```

**Windows / Visual Studio**

Install CMake, and use either the CMake GUI or the command line interface (as on unix) to generate a Visual Studio solution.  Build the solution with Visual Studio.

### Run the code
```
./bin/md /path/to/a/molecule.mol
```

### Edit the code

Copy your solutions into `src/MolecularDynamics.cpp`.  Modify the main file `src/main.cpp` to change UI/visualization.  `CMakeLists.txt` contains a few comments for adding additional files.

