// MolecularDynamics.cpp
// CMU 15-452/662 (Fall 2021)

#include "MolecularDynamics.h"

// General tips:
//
//    - A good approach to writing this kind of code is to
//         * FIRST implement a given potential E
//         * THEN make a copy of the potential code, and modify
//           it to implement the derivative dE
//      Apart from just saving you time, doing it this way makes sure
//      that the potential and gradient code agree on all the details
//      (e.g., which elements get accessed, constants, etc.).  It
//      also provides an opportunity to double-check your derivatives:
//      as you modify your potential code to get the gradient code, you
//      can ask yourself, "is this right?"
//
//    - Your code will become *much* easier to read if, inside the loop,
//      you make temporary variables for things like the atom positions,
//      rather than always accessing them directly from the array.
//
//    - Please also use short variable names that reflect the
//      notation in the exam, e.g., "xi", "xj" instead of "atomPosition1",
//      "atomPosition2".  This will help the TAs verify that you implemented
//      everything correctly---and/or give you partial credit!
//
//    - Don't forget to incorporate the appropriate constant(s) from the top
//      of this file, as well as any relevant physical constants stored
//      with each atom.
//
//    - The Vector3 class provides an implementation of a 3D vector, and
//      is extremely similar to Vec3 from Scotty3D (or any other C++ graphics
//      library); see Vector3.h.
//
//    - Finally, search for the string "TODO" before submitting,
//      to make sure that everything is TODONE!
//

#include <fstream>
#include <sstream>

// Physical constants
const double kStretch = 10.; // stretch strength
const double kBend    = 10.; // bend strength
const double kCoulomb = 10.; // Coulomb strength
const double sigma    = 1.; // Lennard-Jones diameter
const double epsilon  = 1.; // Lennard-Jones well depth

// returns the total energy of the system
double Molecule::totalEnergy()
{
   return potentialEnergy() + kineticEnergy();
}

// returns the potential energy of the system
double Molecule::potentialEnergy()
{
   double U = 0.;

   U += stretchPotential();
   U += bendPotential();
   U += nonlocalPotential();

   return U;
}

// returns the kinetic energy of the system
double Molecule::kineticEnergy()
{
   double K = 0.;

   // add up kinetic energies 1/2 mv^2 of each atom
   for( Atom a : atoms )
   {
      K += a.mass * a.velocity.norm2() / 2.;
   }

   return K;
}

// steps the configuration of the molecule forward in
// time by a duration t, using symplectic Euler
void Molecule::integrate( double t )
{
   // TODO Implement a single time step of symplectic Euler.
   //
   // Tips:
   //    - Think carefully about the order in which quantities get updated.
   //    - Don't forget to evaluate the forces!
}

// returns the total energy due to deviation of bonded atoms
// from their equilibrium separation distance
double Molecule::stretchPotential()
{
   // TODO Sum up and return the stretch potential due to all bonds.
   //
   // Tips:
   //    - Look at Molecule::KineticEnergy for an example of looping over an
   //      array and adding up a quantity.  (Though this time you don't want
   //      to loop over the `atoms` array!)

   return 0.;
}

// returns the total energy due to deviation of bonded
// atoms from their equilibrium bend angles
double Molecule::bendPotential()
{
   // TODO Sum up and return the bend potential due to all corners.
   //
   // Tips:
   //    - Think about how big the bond angle can possibly be, and
   //      choose a formula for computing this angle accordingly.
   //    - If you observe any numerical difficulty (e.g., INF or NaN),
   //      you might try clamping the inputs to any math functions to
   //      a valid input range.

   return 0.;
}

// returns the total energy due to nonlocal interactions,
// including both the Coulomb and Lennard-Jones potentials
double Molecule::nonlocalPotential()
{
   // TODO Sum up and return the nonlocal potentials due to all distinct pairs of atoms.
   //
   // Tips:
   //    - Make sure to call the subroutines Coulomb() and LennardJones(),
   //      rather than implementing these functions inline.
   //    - If you get bogus values, you may be forgetting to skip something.
   
   return 0.;
}

// computes the force on each atom due to the
// total potential energy
void Molecule::updateForces()
{
   // initialize all forces to zero
   for( Atom& a : atoms )
   {
      a.force = Vector3{ 0., 0., 0. };
   }

   addStretchForce();
   addBendForce();
   addNonlocalForce();
}

// adds the gradient of the stretch potential
// potential to the overall force
void Molecule::addStretchForce()
{
   // TODO Compute the stretch forces due to each bond, and
   // add these forces to each atom's total force (Atom::force).
   //
   // Tips:
   //    - Make sure not to overwrite the force!  I.e., make
   //      sure to use += not just =.
   //    - Remember that the potential term for _one_ bond
   //      will contribute a force to _two_ different atoms.
}

// adds the gradient of the bend potential
// potential to the overall force
void Molecule::addBendForce()
{
   // TODO Compute the bend forces due to each corner, and
   // add these forces to each atom's total force (Atom::force).
   //
   // Tips:
   //    - Make sure not to overwrite the force!  I.e., make
   //      sure to use += not just =.
   //    - Remember that the potential term for _one_ corner
   //      will contribute a force to _three_ different atoms.
   //    - Don't forget to take advantage of the expression
   //      we gave you for the gradient of an angle with!
}

// adds the gradient of the nonlocal potentials
// (both Coulomb and Lennard-Jones) to the overall force
void Molecule::addNonlocalForce()
{
   // TODO Compute the nonlocal forces due to each pair of
   // distinct atoms, and add these forces to each atom's
   // total force (Atom::force).
   //
   // Tips:
   //    - Make sure not to overwrite the force!  I.e., make
   //      sure to use += not just =.
   //    - Don't forget to add both the Coulomb and
   //      Lennard-Jones terms, making use of the subroutines
   //      dCoulomb_dr() and dLennardJones_dr().
   //    - Remember that each pairwise potential term
   //      will contribute a force to _two_ different atoms.
}

// returns the Coulomb potential for two charges
// c1, c2 separated by a distance r
double Molecule::Coulomb( double r, double c1, double c2 )
{
   // TODO Return the Coulomb potential for two
   // charges c1, c2 separated by a distance r.
   //
   // Tips:
   //    - This should be a one-liner.
   //    - Don't forget the constant.

   return 0.;
}

// returns the derivative of Coulomb potential w.r.t. r
double Molecule::dCoulomb_dr( double r, double c1, double c2 )
{
   // TODO Return the derivative of the Coulomb
   // potential with respect to the distance r
   //
   // Tips:
   //    - This should also be a one-liner.

   return 0.;
}

// returns the Lennard-Jones potential for a
// separation distance r
double Molecule::LennardJones( double r )
{
   // TODO Return the Lennard-Jones potential for
   // a separation distance r.
   //
   // Tips:
   //    - This should be a one- or two-liner.

   return 0.;
}

// returns the derivative of Lennard-Jones potential w.r.t. r
double Molecule::dLennardJones_dr( double r )
{
   // TODO Return the derivative of the Lennard-Jones
   // potential with respect to the distance r.
   //
   // Tips:
   //    - This should also be a one- or two- liner.

   return 0.;
}

// read molecule description from a file
// (note: sets all velocities and forces to zero)
void Molecule::read( std::string filename )
{
   std::ifstream in( filename.c_str() );

   if( !in.is_open() )
   {
      std::cerr << "Warning: could not read from file " << filename << std::endl;
      return;
   }

   atoms.clear();
   bonds.clear();
   corners.clear();

   std::string s, token;
   while( getline( in, s ) )
   {
      std::stringstream line( s );

      // skip blank and comment lines
      if( s.length() == 0 || s[0] == '#' )
      {
         continue;
      }

      // vertex line --- read: x y z mass charge
      if( s[0] == 'v' )
      {
         Atom a;
         double x, y, z;

         line >> token >> x >> y >> z >> a.mass >> a.charge;
         a.position = Vector3{ x, y, z };
         a.velocity = Vector3{ 0., 0., 0. };
         a.force    = Vector3{ 0., 0., 0. };

         atoms.push_back( a );
      }

      // edge line --- read: index1 index2 restLength
      if( s[0] == 'e' )
      {
         Bond b;

         line >> token >> b.i >> b.j >> b.r;

         if( b.i==b.j )
         {
            std::cout << "WARNING: bond has repeated atom!" << std::endl;
         }

         bonds.push_back( b );
      }
      
      // edge line --- read: index1 index2 restLength
      if( s[0] == 'a' )
      {
         Corner c;

         line >> token >> c.i >> c.j >> c.k >> c.theta;
         
         // convert from degrees to radians
         c.theta *= M_PI/180.;
         
         // convert from interior to exterior
         c.theta = M_PI - c.theta;

         if( c.i==c.j || c.j==c.k || c.k==c.i )
         {
            std::cout << "WARNING: corner has repeated atom!" << std::endl;
         }

         corners.push_back( c );
      }
   }
}

// write current configuration to a molecule
// description file
void Molecule::write( std::string filename )
{
   std::ofstream out( filename.c_str() );

   if( !out.is_open() )
   {
      std::cerr << "Warning: could not write to file " << filename << std::endl;
      return;
   }

   out << "# ATOMS (x y z mass charge)" << std::endl;
   for( Atom a : atoms )
   {
      Vector3 p = a.position;
      out << "v ";
      out << p.x << " " << p.y << " " << p.z << " ";
      out << a.mass << " " << a.charge << std::endl;
   }
   out << std::endl;

   out << "# BONDS (index1 index2 length)" << std::endl;
   for( Bond b : bonds )
   {
      out << "e " << b.i << " " << b.j << " " << b.r << std::endl;
   }
   out << std::endl;

   out << "# CORNERS (index1 index2 index3 angle)" << std::endl;
   for( Corner c : corners )
   {
      out << "a " << c.i << " " << c.j << " " << c.k << " ";
      out << (180.*(M_PI - c.theta)/M_PI) << std::endl;
   }
   out << std::endl;
}

// for debugging purposes, compares the derivatives
// computed in code to finite differences of the energy
void Molecule::checkDerivatives()
{
   // To validate our code for the potential energy derivatives, we're going to
   // take finite differences of the value returned by the energy function,
   // with respect to perturbations of the atoms' coordinates.  This method is
   // reliable in the sense that it always gives us an approximation of the
   // derivatives of the expression computed in the energy code---no matter
   // what that code looks like.  On the other hand, it's inaccurate and very
   // slow to evaluate (since we have to evaluate the energy function twice for
   // every atom).  So, we prefer closed-form expressions---and merely want to
   // check here that those closed-form expressions are correct!

   const double eps = 1e-7; // perturbation size
   double maxError;

   std::cout << "Checking stretch forces..." << std::endl;
   double E0stretch = stretchPotential(); // current energy
   for( Atom& a : atoms ) a.force = Vector3{ 0., 0., 0. };
   addStretchForce();
   maxError = 0.;
   for( Atom& a : atoms )
   {
      Vector3 dEdx;
      for( int i = 0; i < 3; i++ )
      {
         a.position[i] += eps;
         double Estretch = stretchPotential();
         dEdx[i] = ( Estretch - E0stretch ) / eps;
         a.position[i] -= eps;
      }
      double error = ( dEdx - (-a.force)).norm();
      maxError = std::max( maxError, error );
   }
   std::cout << "   worst error: " << maxError << std::endl;

   std::cout << "Checking bend forces..." << std::endl;
   double E0bend = bendPotential(); // current energy
   for( Atom& a : atoms ) a.force = Vector3{ 0., 0., 0. };
   addBendForce();
   maxError = 0.;
   for( Atom& a : atoms )
   {
      Vector3 dEdx;
      for( int i = 0; i < 3; i++ )
      {
         a.position[i] += eps;
         double Ebend = bendPotential();
         dEdx[i] = ( Ebend - E0bend ) / eps;
         a.position[i] -= eps;
      }
      double error = ( dEdx - (-a.force)).norm();
      maxError = std::max( maxError, error );
   }
   std::cout << "   worst error: " << maxError << std::endl;

   std::cout << "Checking nonlocal forces..." << std::endl;
   double E0nonlocal = nonlocalPotential(); // current energy
   for( Atom& a : atoms ) a.force = Vector3{ 0., 0., 0. };
   addNonlocalForce();
   maxError = 0.;
   for( Atom& a : atoms )
   {
      Vector3 dEdx;
      for( int i = 0; i < 3; i++ )
      {
         a.position[i] += eps;
         double Enonlocal = nonlocalPotential();
         dEdx[i] = ( Enonlocal - E0nonlocal ) / eps;
         a.position[i] -= eps;
      }
      double error = ( dEdx - (-a.force)).norm();
      maxError = std::max( maxError, error );
   }
   std::cout << "   worst error: " << maxError << std::endl;
}

