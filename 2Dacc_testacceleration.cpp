//******************** FILE: 2DACC_TESTACCELERATION.CPP ********************
//
// Description: Source file used to the the particle acceleration functionalities
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: May 2012
// Last update: May 2012
//
//******************************************************************************

// Standard header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <ctime>
#include <algorithm>

// Project specific header files
#include "./beams.hpp"
#include "./particles.hpp"
#include "./constants.hpp"
#include "./general.hpp"

// Boost header files
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>



using namespace std;


// Main function
int main() {

  /*// FEATURE TO BE CHARACTERIZED
  // Off-axis acceleration of a single electron

  // Define the driver beams
  double lambda0 = 800e-9;
  double w0 = 4e-6;
  //double ka = 8.0;
  double zR = PI*pow(w0,2.0)/lambda0;
  //CNonparaxialTM01 npTM01(lambda0,1e9,ka,125,PI);
  CParaxialTM01 pTM01(lambda0,1e15,w0,10e-15,0);
  //CBeam * NonParaxialBeam = &npTM01;
  CBeam * ParaxialBeam = &pTM01;

  // Define the particle to be accelerated
  vector<double> InitialCoords(4,0.0);
  //InitialCoords[0] = lambda0/10.0;
  //InitialCoords[1] = 1.7*lambda0;
  CParticle theElectron("Electron",-QE,ME,InitialCoords);

  // Open output files
  ofstream out1, out2;
  out1.open("./dat/test_acceleration/acceleration2D_nonparaxial.dat", ios::out);
  out1 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out1 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out1 << "#" << endl;
  out1 << "# DESCRIPTION" << endl;
  out1 << "# On-axis electron acceleration test" << endl;
  out1 << "#" << endl;
  theElectron.printparameters(out1);
  out1 << "#" << endl;
  NonParaxialBeam->printparameters(out1);
  out1 << "#" << endl;

  out2.open("./dat/test_acceleration/acceleration2D_paraxial.dat", ios::out);
  out2 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out2 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out2 << "#" << endl;
  out2 << "# DESCRIPTION" << endl;
  out2 << "# On-axis electron acceleration test" << endl;
  out2 << "#" << endl;
  theElectron.printparameters(out2);
  out2 << "#" << endl;
  ParaxialBeam->printparameters(out2);
  out2 << "#" << endl;

  // Particle acceleration
  //theElectron.accelerate2D(NonParaxialBeam,-75e-15,1e-16,100e-15,out1);
  //theElectron.reset_coordinates(InitialCoords);
  theElectron.accelerate2D(ParaxialBeam,-50e-15,1e-15,15e-12,out2);

  // Close output files
  //out1 << endl;
  //out1 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  //out1.close();
  out2 << endl;
  out2 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out2.close();*/


  // FEATURE TO BE CHARACTERIZED
  // Acceleration of a cloud of electron

  // Definition of the driver beams
  double lambda0 = 800e-9;
  double w0 = 3e-6;
  double ka = 277.6;
  double phi0 = 0.0;
  CNonparaxialTM01 npTM01(lambda0,4.0e14,ka,277.0,phi0);
  CParaxialTM01 pTM01(lambda0,4.0e14,w0,10e-15,phi0);
  CBeam * NonParaxialBeam = &npTM01;
  CBeam * ParaxialBeam = &pTM01;

  // Open output files
  ofstream out1, out2;
  out1.open("./dat/test_acceleration/acceleration2D_nonparaxial.dat", ios::out);
  out1 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out1 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out1 << "#" << endl;
  out1 << "# DESCRIPTION" << endl;
  out1 << "# Acceleration of a cloud of electrons" << endl;
  out1 << "#" << endl;
  NonParaxialBeam->printparameters(out1);
  out1 << "#" << endl;
  out1 << "# ri [m]  zi [m]  Ki [eV]  rf [m]  zf [m]  Kf [eV]" << endl;

  out2.open("./dat/test_acceleration/acceleration2D_paraxial.dat", ios::out);
  out2 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out2 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out2 << "#" << endl;
  out2 << "# DESCRIPTION" << endl;
  out2 << "# Acceleration of a cloud of electrons" << endl;
  out2 << "#" << endl;
  out2 << "#" << endl;
  ParaxialBeam->printparameters(out2);
  out2 << "#" << endl;
  out2 << "# ri [m]  zi [m]  Ki [eV]  rf [m]  zf [m]  Kf [eV]" << endl;


  // Parametrization of the random number generator
  typedef boost::mt19937 engine_type; // Mersenne-Twister rule
  engine_type eng;
  eng.seed(42);
  typedef boost::normal_distribution<> normal_dist_type;
  typedef boost::variate_generator< engine_type&, normal_dist_type > normal_dist_rng;
  normal_dist_rng coord_r_rng(eng,normal_dist_type(0,80e-9));
  normal_dist_rng coord_z_rng(eng,normal_dist_type(0,80e-9));

  // Particle acceleration
  int N = 1000;

  // Define the particles to be accelerated
  vector<double> InitialCoords(4,0.0);
  CParticle theElectron1("Electron",-QE,ME,InitialCoords);
  CParticle theElectron2("Electron",-QE,ME,InitialCoords);

  // Generate the initial coordinates vector
  vector<double> InitialR (N,0); // Radial coordinate vector
  vector<double> InitialZ (N,0); // Axial coordinate vector
  for (int n = 0; n<N; n++) {

    InitialR[n] = coord_r_rng();
    InitialZ[n] = coord_z_rng();

  } // end for
  sort(InitialR.begin(),InitialR.end()); // Sort the radial coordinates vector

  // Particle acceleration loop
  for (int n = 0; n<N; n++) {

    cout << n << endl;

    // Set initial coordinates
    InitialCoords[0] = InitialR[n];
    InitialCoords[1] = InitialZ[n];
    theElectron1.reset_coordinates(InitialCoords);
    theElectron2.reset_coordinates(InitialCoords);

    // Write initial coordinates
    out1 << theElectron1.coordinates[0]/1e-9 << "  " << theElectron1.coordinates[1]/1e-9 << "  " << theElectron1.get_kineticenergy() << "  ";
    out2 << theElectron2.coordinates[0]/1e-9 << "  " << theElectron2.coordinates[1]/1e-9 << "  " << theElectron2.get_kineticenergy() << "  ";

    // Accelerate the particles
    //theElectron1.accelerate2D(NonParaxialBeam,-50e-15,15e-12);
    theElectron2.accelerate2D(ParaxialBeam,-50e-15,15e-12);

    // Write final coordinates
    out1 << theElectron1.coordinates[0]/1e-9 << "  " << (theElectron1.coordinates[1]-0.0045)/1e-9 << "  " << theElectron1.get_kineticenergy();
    out2 << theElectron2.coordinates[0]/1e-9 << "  " << (theElectron2.coordinates[1]-0.0045)/1e-9 << "  " << theElectron2.get_kineticenergy();
    out1 << endl;
    out2 << endl;

    // Separate index in output file
    //if ((n+1) % 200 == 0) {
    //  out1 << endl << endl;
    //  out2 << endl << endl;
    //} // end if

  } // end for

    // Close output files
  out1 << endl;
  out1 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out1.close();
  out2 << endl;
  out2 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out2.close();


  /*// FEATURE TO BE CHARACTERIZED
  // Reproduction of the 1D results with the 2D code

  // Definition of the driver beam
  CNonparaxialTM01 npTM01(800e-9,2.0e15,1.0,1.0,0.0);
  CBeam * NonParaxialBeam = &npTM01;

  // Define the particle to be accelerated
  vector<double> InitialCoords(4,0.0);
  CParticle theElectron("Electron",-QE,ME,InitialCoords);

  // Open output files
  ofstream out1;
  out1.open("./dat/test_acceleration/acceleration2D_energy_z0phi.dat", ios::out);
  out1 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out1 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out1 << "#" << endl;
  out1 << "# DESCRIPTION" << endl;
  out1 << "# Energy gain dependence on the z0 and phi0 parameters" << endl;
  out1 << "#" << endl;
  NonParaxialBeam->printparameters(out1);
  out1 << "#" << endl;
  theElectron.printparameters(out1);
  out1 << "#" << endl;
  out1 << "# z0/a  phi0 [rad]  K [eV]" << endl;

  // Parameter sweep
  for (double z0a = -35.0; z0a<= 35.0; z0a+=0.5) {

    cout << z0a << endl;

    for (double phi0 = 0.0; phi0 <= 2.0*PI; phi0 += PI/50.0) {

      // Set the beam's and particle's parameters
      npTM01.resetarrays_phase(phi0);
      InitialCoords[1] = z0a*npTM01.charlength_z;
      theElectron.reset_coordinates(InitialCoords);

      // Particle acceleration
      double tini = -abs(z0a*npTM01.charlength_z/C0) - 10e-15;
      theElectron.accelerate2D(NonParaxialBeam,tini,10e-12);

      // Write results to file
      out1 << z0a << "  " << phi0 << "  " << theElectron.get_kineticenergy() << endl;

    } // end for

    out1 << endl;

  } // end for


  // Close output files
  out1 << endl;
  out1 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out1.close();*/


  /*// FEATURE TO BE CHARACTERIZED
  // Acceleration of a cloud of electron: generation of synchronized counter-propagating electron beams

  // Definition of the driver beam
  CNonparaxialTM01 npTM01(800e-9,2.0e15,1.0,1.0,PI);
  CBeam * NonParaxialBeam = &npTM01;

  // Open output files
  ofstream out1;
  out1.open("./dat/test_acceleration/acceleration2D_counterpropagating.dat", ios::out);
  out1 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out1 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out1 << "#" << endl;
  out1 << "# DESCRIPTION" << endl;
  out1 << "# Acceleration of a cloud of electrons: generation of synchronized counterpropagating electron beams" << endl;
  out1 << "#" << endl;
  NonParaxialBeam->printparameters(out1);
  out1 << "#" << endl;
  out1 << "# ri [m]  zi [m]  Ki [eV]  rf [m]  zf [m]  Kf [eV]" << endl;

  // Parametrization of the random number generator
  typedef boost::mt19937 engine_type; // Mersenne-Twister rule
  engine_type eng;
  eng.seed(42);
  typedef boost::normal_distribution<> normal_dist_type;
  typedef boost::variate_generator< engine_type&, normal_dist_type > normal_dist_rng;
  normal_dist_rng coord_r_rng(eng,normal_dist_type(0,1e-6));
  normal_dist_rng coord_z_rng(eng,normal_dist_type(-30.0*npTM01.charlength_z,1e-6));

  // Particle acceleration
  int N = 500;

  // Define the particles to be accelerated
  vector<double> InitialCoords(4,0.0);
  CParticle theElectron("Electron",-QE,ME,InitialCoords);

  // Generate the initial coordinates vector
  vector<double> InitialR (N,0); // Radial coordinate vector
  vector<double> InitialZ (N,0); // Axial coordinate vector
  for (int n = 0; n<N; n++) {

    InitialR[n] = coord_r_rng();
    InitialZ[n] = coord_z_rng();

  } // end for
  sort(InitialZ.begin(),InitialZ.end()); // Sort the radial coordinates vector

  // Particle acceleration loop
  for (int n = 0; n<N; n++) {

    cout << n << endl;

    // Set initial coordinates
    InitialCoords[0] = InitialR[n];
    InitialCoords[1] = InitialZ[n];
    theElectron.reset_coordinates(InitialCoords);

    // Write initial coordinates
    out1 << theElectron.coordinates[0] << "  " << theElectron.coordinates[1] << "  " << theElectron.get_kineticenergy() << "  ";

    // Accelerate the particles
    double tini = -abs(theElectron.coordinates[1]/C0) - 10e-15;
    theElectron.accelerate2D(NonParaxialBeam,tini,15e-12);

    // Write final coordinates
    out1 << theElectron.coordinates[0] << "  " << theElectron.coordinates[1] << "  " << theElectron.get_kineticenergy();
    out1 << endl;

    // Separate index in output file
    if ((n+1) % (N/5) == 0) {
      out1 << endl << endl;
    } // end if

  } // end for

    // Close output files
  out1 << endl;
  out1 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out1.close();*/

  /*// FEATURE TO BE CHARACTERIZED
  // Off-axis acceleration of a single electron

  // Define the driver beams
  CNonparaxialTM01 npTM01(800e-9,2.0e14,124.366,277,0.0933*PI);
  CParaxialTM01 pTM01(800e-9,2.0e14,2e-6,10e-15,0.0933*PI);
  CBeam * NonParaxialBeam = &npTM01;
  CBeam * ParaxialBeam = &pTM01;

  // Define the particle to be accelerated
  vector<double> InitialCoords(4,0.0);
  InitialCoords[0] = 0.2e-6;
  CParticle theElectron("Electron",-QE,ME,InitialCoords);

  // Open output files
  ofstream out1, out2;
  out1.open("./dat/test_acceleration/acceleration2D_nonparaxial.dat", ios::out);
  out1 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out1 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out1 << "#" << endl;
  out1 << "# DESCRIPTION" << endl;
  out1 << "# Off-axis electron acceleration test" << endl;
  out1 << "#" << endl;
  theElectron.printparameters(out1);
  out1 << "#" << endl;
  NonParaxialBeam->printparameters(out1);
  out1 << "#" << endl;

  out2.open("./dat/test_acceleration/acceleration2D_paraxial.dat", ios::out);
  out2 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out2 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out2 << "#" << endl;
  out2 << "# DESCRIPTION" << endl;
  out2 << "# Off-axis electron acceleration test" << endl;
  out2 << "#" << endl;
  theElectron.printparameters(out2);
  out2 << "#" << endl;
  ParaxialBeam->printparameters(out2);
  out2 << "#" << endl;

  // Particle acceleration
  theElectron.accelerate2D(NonParaxialBeam,-100e-15,5e-16,10e-12,out1);
  theElectron.reset_coordinates(InitialCoords);
  theElectron.accelerate2D(ParaxialBeam,-100e-15,5e-16,10e-12,out2);

  // Close output files
  out1 << endl;
  out1 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out1.close();
  out2 << endl;
  out2 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out2.close();*/

  return 0;
}
