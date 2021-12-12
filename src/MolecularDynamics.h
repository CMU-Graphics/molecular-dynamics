// MolecularDynamics.h
// CMU 15-452/662 (Fall 2021)

#pragma once

#include <string>
#include <vector>

#include "Vector3.h"

// The Atom class stores all the basic information for a single
// atom.  It does not store any information about how atoms are
// connected by bonds.  Note that mass and charge are constants
// that are read in from a file---you do not need to set them or
// change them (though you may need to use them).
class Atom {
   public:
      double mass;
      double charge;
      Vector3 position;
      Vector3 velocity;
      Vector3 force;
};

// The Bond class represents a single bond between a pair of
// atoms i and j, as well as the separation distance this
// bond ideally achieves at equilibrium.  The separation
// distance is a constant read in from a file---you do not
// need to set it or change it.
class Bond {
   public:
      // indices of atoms as 0-based indices
      // into the list of atom positions
      int i, j;

      // equilibrium separation distance
      double r;
};

// The Corner class represents a triple of atoms i, j, k
// that seek to achieve an equilibrium angle theta.  This
// angle is a constant read in from a file---you do not
// need to set it or change it.
class Corner {
   public:
      // indices of atoms as 0-based indices
      // into the list of atom positions
      int i, j, k;

      // equilibrium bond angle between ji and jk,
      // given as an exterior angle in radians
      double theta;
};

// A Molecule is a collection of Atoms, connected by Bonds
// of given lengths. It also describes a collection of Corners
// where the molecule seeks to achieve given angles.
class Molecule {

   public:
      // properties and current state of each atom
      std::vector<Atom> atoms;

      // equilibrium separation distances
      std::vector<Bond> bonds;

      // equilibrium bond angles
      std::vector<Corner> corners;

      // read molecule description from a file
      // (note: sets all velocities and forces to zero)
      void read( std::string filename );

      // write current configuration to a molecule
      // description file
      void write( std::string filename );

      // returns the total energy of the system
      double totalEnergy();

      // steps the configuration of the molecule forward in
      // time by a duration t, using symplectic Euler
      void integrate( double t );

      // for debugging purposes, compares the derivatives
      // computed in code to finite differences of the energy
      void checkDerivatives();

   protected:

      // returns the potential energy of the system
      double potentialEnergy();

      // returns the kinetic energy of the system
      double kineticEnergy();

      // returns the total energy due to deviation of
      // bonded atoms from their equilibrium separation
      // distance
      double stretchPotential();

      // returns the total energy due to deviation of
      // bonded atoms from their equilibrium bend angles
      double bendPotential();

      // returns the total energy due to nonlocal
      // interactions, including both the Coulomb and
      // Lennard-Jones potentials
      double nonlocalPotential();

      // computes the force on each atom due to the
      // total potential energy
      void updateForces();

      // adds the gradient of the stretch, bend, or
      // nonlocal potential (with respect to the atom
      // positions) to the overall force
      void  addStretchForce();
      void     addBendForce();
      void addNonlocalForce();

      // returns the Coulomb potential for two charges
      // c1, c2 separated by a distance r
      double Coulomb( double r, double c1, double c2 );

      // returns the derivative of Coulomb potential w.r.t. r
      double dCoulomb_dr( double r, double c1, double c2 );

      // returns the Lennard-Jones potential for a
      // separation distance r
      double LennardJones( double r );

      // returns the derivative of Lennard-Jones potential w.r.t. r
      double dLennardJones_dr( double r );
};

