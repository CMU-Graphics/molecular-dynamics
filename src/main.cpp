#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "args/args.hxx"
#include "imgui.h"

#include "MolecularDynamics.h"
Molecule molecule;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;
polyscope::PointCloud *psCloud;
polyscope::CurveNetwork *psNetwork;

// Some algorithm parameters
float timestep = .01;
bool runIntegration = false;

std::vector<Vector3> getPositions( Molecule& m )
{
   std::vector<Vector3> p;

   for( Atom a : m.atoms )
   {
      p.push_back( a.position );
   }

   return p;
}

std::vector<std::vector<size_t>> getEdges( Molecule& m )
{
   std::vector<std::vector<size_t>> e;

   for( Bond b : m.bonds )
   {
      std::vector<size_t> ij;
      ij.push_back( b.i );
      ij.push_back( b.j );
      e.push_back( ij );
   }

   return e;
}

void updateMolecule()
{
  psCloud->updatePointPositions( getPositions(molecule) );
  psNetwork->updateNodePositions( getPositions(molecule) );
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  ImGui::SliderFloat("time step", &timestep, 0.001, 1.);

  ImGui::Checkbox("integrate", &runIntegration);
  if( runIntegration || ImGui::Button("take step") )
  {
     molecule.integrate( timestep );
     updateMolecule();
     std::cout << "total energy: " << molecule.totalEnergy() << std::endl;
  }

  if( ImGui::Button("heat up") )
  {
     for( Atom& a : molecule.atoms )
     {
        a.velocity += 0.1 * Vector3{
           polyscope::randomReal( -1., 1. ),
           polyscope::randomReal( -1., 1. ),
           polyscope::randomReal( -1., 1. )
        };
     }
     std::cout << "total energy: " << molecule.totalEnergy() << std::endl;
  }


  if( ImGui::Button("write to molecule.mol") )
  {
     molecule.write( "molecule.mol" );
     std::cout << "Wrote to file molecule.mol" << std::endl;
  }

  if (ImGui::Button("check derivatives")) {
     molecule.checkDerivatives();
  }
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("Molecular Dynamics Simulator");
  args::Positional<std::string> inputFilename(parser, "molecule", "A molecule file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a molecule file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  molecule.read( args::get(inputFilename) );

  // Register the molecule with polyscope
  psCloud = polyscope::registerPointCloud( "atoms", getPositions(molecule) );
  psNetwork = polyscope::registerCurveNetwork( "bonds", getPositions(molecule), getEdges(molecule) );

  // Set visualization options
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  polyscope::options::groundPlaneHeightFactor = 0.3;

  // Set size of atoms and bonds
  psCloud->setPointRadius( .025 );
  psNetwork->setRadius( .008 );

  // Color atoms according to CPK coloring
  std::vector<Vector3> color;
  for( Atom a : molecule.atoms )
  {
     int m = (int) round(a.mass);
     Vector3 c;
     switch( m )
     {
        case 1: // hydrogen
           c = Vector3{ .9, .9, .9 };
           break;
        case 12: // carbon
           c = Vector3{ .1, .1, .1 };
           break;
        case 14: // nitrogen
           c = Vector3{ 0., .0, .9 };
           break;
        case 16: // oxygen
           c = Vector3{ .9, 0., 0. };
           break;
         default:
           c = Vector3{ 0., 0., 0. };
     }
     color.push_back( c );
  }
  psCloud->addColorQuantity( "mass", color );
  psCloud->getQuantity( "mass" )->setEnabled( true );
  psCloud->setMaterial( "wax" );
  psNetwork->setColor( glm::vec3( .5, .5, .5 ) );
  psNetwork->setMaterial( "wax" );

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
