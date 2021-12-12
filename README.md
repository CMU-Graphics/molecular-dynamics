# molecular-dynamics

![Molecule visualization](interface.png)

Skeleton code and visualization for a basic molecular dynamics simulator, provided as part of the 15-462/662 final exam at CMU.  Visualization via [Polyscope](http://polyscope.run/).

**IMPORTANT**: You do _not_ have to build and run this version of the code for the exam (which can be more difficult to compile, due to the visualization routines).  It is provided "just for fun."  See Piazza for the official exam code.

**NOTE**: To use this code, you should copy your solution _routines_ into `src/MolecularDynamics.cpp`.  Copying in the whole file will not work, since this version of the code has additional routines (such as `Molecule::write()`) that are not present in the official exam code.

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

