//******************** FILE: PARTICLES.HPP ********************
//
// Description: This header file contains the particle class definition as well as the definition of its derived class.
//              It also contains the definition of the ODE systems in 1 and 2 dimensions.
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: May 2012
// Last update: May 2012
//
//*********************************************************

#ifndef PARTICLES_HPP_INCLUDED
#define PARTICLES_HPP_INCLUDED

// Standard header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>

// GSL header files
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

// Project specific header files
#include "./constants.hpp"

using namespace std;


// ******************* PROTOTYPES ***********************************
int odesys1D(double,const double *, double *, void *); // Function odesys1D prototype
int odesys2D(double,const double *, double *, void *); // Function odesys2D prototype



//******************** PARTICLE CLASS DEFINITION ********************
// CLASS CParticle define fully general particle objects. See derived classes CElectron and CProton for specific instances
class CParticle {
  protected:
    const string name; // Name of the particle
    const double charge; // Charge of the particle
    const double mass; // Mass of the particle
  public:
    CParticle (string,double,double,vector<double>); // Constructor
    vector<double> coordinates; // Coordinate vector of the particle: [r,z,vr,vz]
    void reset_coordinates (vector<double>); // Reset the particle's coordinates
    double get_kineticenergy (); // Compute kinetic energy
    double get_totalenergy (); // Compute total energy
    void accelerate1D (CBeam *,double,double); // On-axis acceleration, no output stream
    void accelerate1D (CBeam *,double,double,double,ostream&); // On-axis acceleration, write to given output stream
    void accelerate2D (CBeam *,double,double); // Two-dimensional acceleration, no output stream
    void accelerate2D (CBeam *,double,double,double,ostream&); // Two-dimensional acceleration, write to given output stream
    void printparameters (ostream&); // Print the particle parameters to output stream
    friend int odesys1D(double,const double *, double *, void *); // Friendship granted to the odesys1D function
    friend int odesys2D(double,const double *, double *, void *); // Friendship granted to the odesys2D function
};


// STRUCT BeamParticlePair definition
// Used as a parameter container in the odesys1D and odesys2D functions
struct BeamParticlePair {
  CBeam * beam; // Pointer to the beam object
  CParticle * particle; // Pointer to the particle object
} ;


// CONSTRUCTOR CParticle definition
// Set the member attributes to the prescribed values
//
// Input
//   - initialname: name of the particle
//   - initialcharge: charge of the particle
//   - initialmass: mass of the particle
//   - initialcoords: initial coordinates [r,z,vr,vz]
//
CParticle::CParticle (string initialname, double initialcharge, double initialmass, vector<double> initialcoords)
  :  name(initialname),
     charge(initialcharge),
     mass(initialmass)
{

  coordinates = initialcoords;

} // End constructor CParticle definition


// MEMBER FUNCTION reset_coordinates(vector<double>) definition
// Reset the coordinates of a particle to given values
//
// Input
//   - newcoords: new coordinates [r,z,vr,vz]
// Output
//     none
//
void CParticle::reset_coordinates(vector<double> newcoords) {

  coordinates = newcoords;

} // End member function reset_coordinates definition


// MEMBER FUNCTION get_kineticenergy() definition
// Compute the kinetic energy of the particle
//
// Input
//     none
// Output
//   - kineticenergy: kinetic energy of the particle
//
double CParticle::get_kineticenergy () {

  double vr = coordinates[2];
  double vz = coordinates[3];
  double gamma = pow(1.0-(pow(vr,2.0)+pow(vz,2.0))/pow(C0,2.0), -0.5);
  double kineticenergy = (gamma-1.0)*mass*pow(C0,2.0)/QE; // Total kinetic energy in eV

  return kineticenergy;

} // End member function get_kineticenergy definition


// MEMBER FUNCTION get_totalenergy() definition
// Compute the total energy of the particle
//
// Input
//     none
// Output
//   - totalenergy: total energy of the particle
//
double CParticle::get_totalenergy () {

  double vr = coordinates[2];
  double vz = coordinates[3];
  double gamma = pow(1.0-(pow(vr,2.0)+pow(vz,2.0))/pow(C0,2.0), -0.5);
  double totalenergy =gamma*mass*pow(C0,2.0)/QE; // Total energy in eV

  return totalenergy;

} // End member function get_totalenergy definition


// MEMBER FUNCTION accelerate1D (Cbeam*,double,double) definition (NO OUTPUT)
// Function that simulate the 1D laser-driven acceleration of the particle by a given driver beam. No trajectory output.
//
// Input
//   - beam: pointer to the corresponding beam object
//   - tstart: start time (the beam passes throught its waist at z=0, t=0)
//   - tstop: stop time (the beam passes throught its waist at z=0, t=0)
// Output
//     none
//
void CParticle::accelerate1D (CBeam * beam, double tstart,double tstop) {

  // Set integrator parameters
  double t = tstart;
  double dt = 1e-17;
  double t_step = 5e-15;
  double t_max = tstop;
  const double eps_abs = 1e-20;
  const double eps_rel = 1e-12;
  BeamParticlePair params;
  params.beam = beam;
  params.particle = this;

  // Set GSL odeiv parameters
  const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45; // Runge-Kutta-Fehlberg 4-5 stepper
  gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type,4);
  gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
  gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (4);
  gsl_odeiv_system sys = {odesys1D, NULL, 4, &params};

  // Numerical integration of the system of ODEs
  int status = GSL_SUCCESS;
  double t_target = t;
  for (; t_target <= t_max; t_target += t_step ) {
    while (t < t_target) {
      status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,coordinates.data());
      if (status != GSL_SUCCESS)
        break;
    } // end while
    if (status != GSL_SUCCESS)
      break;
  } // end for

} // end member function accelerate1D definition


// MEMBER FUNCTION accelerate1D (Cbeam*,double,double,double,ofsteam&) definition (WITH OUTPUT)
// Function that simulate the 1D laser-driven acceleration of the particle by a given driver beam. With trajectory output.
//
// Input
//   - beam: pointer to the corresponding beam object
//   - tstart: start time (the beam passes throught its waist at z=0, t=0)
//   - tincr: time step between each trajectory point
//   - tstop: stop time (the beam passes throught its waist at z=0, t=0)
//   - out: valide ofstream object for output
// Output
//     none
//
void CParticle::accelerate1D (CBeam * beam, double tstart, double tincr, double tstop, ostream& out) {

  // Set integrator parameters
  double t = tstart;
  double dt = 1e-17;
  double t_step = tincr;
  double t_max = tstop;
  const double eps_abs = 1e-20;
  const double eps_rel = 1e-12;
  BeamParticlePair params;
  params.beam = beam;
  params.particle = this;

  // Set GSL odeiv parameters
  const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45; // Runge-Kutta-Fehlberg 4-5 stepper
  gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type,4);
  gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
  gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (4);
  gsl_odeiv_system sys = {odesys1D, NULL, 4, &params};

  // Numerical integration of the system of ODEs
  out << "# t   r [m]   z [m]   vr/c   vz/c   K [MeV]" << endl;
  int status = GSL_SUCCESS;
  double t_target = t;
  for (; t_target <= t_max; t_target += t_step ) {
    while (t < t_target) {
      status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,coordinates.data());
      if (status != GSL_SUCCESS)
        break;
    } // end while
    if (status != GSL_SUCCESS)
      break;
    // Write trajectory to file
    out << t << "  " << coordinates[0] << "  " << coordinates[1];
    out << "  " << coordinates[2]/C0 << "  " << coordinates[3]/C0 << "  " << this->get_kineticenergy() << endl;
  } // end for

} // end member function accelerate1D definition


// MEMBER FUNCTION accelerate2D (Cbeam*,double,double) definition (NO OUTPUT)
// Function that simulate the 2D laser-driven acceleration of the particle by a given driver beam. No trajectory output.
//
// Input
//   - beam: pointer to the corresponding beam object
//   - tstart: start time (the beam passes throught its waist at z=0, t=0)
//   - tstop: stop time (the beam passes throught its waist at z=0, t=0)
// Output
//     none
//
void CParticle::accelerate2D (CBeam * beam, double tstart,double tstop) {

  // Set integrator parameters
  double t = tstart;
  double dt = 1e-17;
  double t_step = 5e-15;
  double t_max = tstop;
  const double eps_abs = 1e-20;
  const double eps_rel = 1e-12;
  BeamParticlePair params;
  params.beam = beam;
  params.particle = this;

  // Set GSL odeiv parameters
  const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45; // Runge-Kutta-Fehlberg 4-5 stepper
  gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type,4);
  gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
  gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (4);
  gsl_odeiv_system sys = {odesys2D, NULL, 4, &params};

  // Numerical integration of the system of ODEs
  int status = GSL_SUCCESS;
  double t_target = t;
  for (; t_target <= t_max; t_target += t_step ) {
    while (t < t_target) {
      status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,coordinates.data());
      if (status != GSL_SUCCESS)
        break;
    } // end while
    if (status != GSL_SUCCESS)
      break;
  } // end for

} // end member function accelerate2D definition


// MEMBER FUNCTION accelerate2D (Cbeam*,double,double,double,ofstream&) definition (WITH OUTPUT)
// Function that simulate the 2D laser-driven acceleration of the particle by a given driver beam. With trajectory output.
//
// Input
//   - beam: pointer to the corresponding beam object
//   - tstart: start time (the beam passes throught its waist at z=0, t=0)
//   - tincr: time step between each trajectory point
//   - tstop: stop time (the beam passes throught its waist at z=0, t=0)
//   - out: valide ofstream object for output
// Output
//     none
//
void CParticle::accelerate2D (CBeam * beam, double tstart, double tincr, double tstop, ostream& out) {

  // Set integrator parameters
  double t = tstart;
  double dt = 1e-17;
  double t_step = tincr;
  double t_max = tstop;
  const double eps_abs = 1e-20;
  const double eps_rel = 1e-12;
  BeamParticlePair params;
  params.beam = beam;
  params.particle = this;

  // Set GSL odeiv parameters
  const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45; // Runge-Kutta-Fehlberg 4-5 stepper
  gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type,4);
  gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
  gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (4);
  gsl_odeiv_system sys = {odesys2D, NULL, 4, &params};

  // Numerical integration of the system of ODEs
  out << "# t   r [m]   z [m]  vr/c   vz/c   K [MeV]" << endl;
  int status = GSL_SUCCESS;
  double t_target = t;
  for (; t_target <= t_max; t_target += t_step ) {
    while (t < t_target) {
      status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,coordinates.data());
      if (status != GSL_SUCCESS)
        break;
    } // end while
    if (status != GSL_SUCCESS)
      break;
    // Write trajectory to file
    vector<complex<double> > fields = beam->fields(coordinates[0],coordinates[1],t);
    out << t << "  " << coordinates[0] << "  " << coordinates[1];
    out << "  " << coordinates[2]/C0 << "  " << coordinates[3]/C0 << "  " << this->get_kineticenergy();
    out << "  " << real(fields[0]) << "  " << real(fields[1]) << "  " << real(fields[2]) << endl;
  } // end for

} // end member function accelerate2D definition


// MEMBER FUNCTION printparameters(ofstream) definition
// Print the particle parameters to a given output stream
//
// Input
//   - out: valid ofstream object (output stream)
// Output
//     none
//
void CParticle::printparameters(ostream& out) {

  out << "# PARTICLE PARAMETERS" << endl;
  out << "# Name = " << name << endl;
  out << "# Mass = " << mass << " kg" << endl;
  out << "# Charge = " << charge << " C" << endl;

} // end member function printparameters(ofstream) definition



//******************** SYSTEM OF ODES DEFINITIONS ********************

// FUNCTION odesys1D (double, const double *, double *, void *) definition
// System of ODEs describing the dynamics of a particle subject to the Lorentz force in 1 dimension
//
// Input
//   - t: current time
//   - y[]: pointer to the coordinates container
//   - f[]: pointer to the derivatives container
//    - params: parameters (to be cast in a BeamParticlePair structure)
// Output
//   - returns GSL_SUCCESS if successful
//
int odesys1D(double t,const double y[], double f[], void * params){

  // Cast parameters into the appropriate structure
  BeamParticlePair& p = *static_cast<BeamParticlePair* >(params);
  CBeam * theBeam = p.beam;
  CParticle * theParticle = p.particle;

  // Get the current axial coordinates (recall that r=0 and vr=0 in 1D)
  double z = y[1];
  double vz = y[3];

  // Compute longitudinal electric field at the wanted position
  complex<double> Ez = theBeam->Ezfield(0.0,z,t);

  // Update the derivatives container
  f[0] = 0.0; // dr/dt
  f[1] = vz; // dz/dt
  f[2] = 0.0; // dvr/dt
  f[3] = (theParticle->charge/theParticle->mass)*real(Ez)*pow(1.0-pow(vz,2.0)/pow(C0,2.0),1.5);// dvz/dt

  return GSL_SUCCESS;

} // end function odesys1D definition


// FUNCTION odesys2D (double, const double *, double *, void *) definition
// System of ODEs describing the dynamics of a particle subject to the Lorentz force in 2 dimensions
//
// Input
//   - t: current time
//   - y[]: pointer to the coordinates container
//   - f[]: pointer to the derivatives container
//    - params: parameters (to be cast in a BeamParticlePair structure)
// Output
//   - returns GSL_SUCCESS if successful
//
int odesys2D(double t,const double y[], double f[], void * params){

  // Cast parameters into the appropriate structure
  BeamParticlePair& p = *static_cast<BeamParticlePair* >(params);
  CBeam * theBeam = p.beam;
  CParticle * theParticle = p.particle;

  // Get the current coordinates
  double r = y[0];
  double z = y[1];
  double vr = y[2];
  double vz = y[3];

  // Compute the electromagnetic fields
  vector<complex<double> > theFields = theBeam->fields(r,z,t);
  double Er = real(theFields[0]);
  double Ez = real(theFields[1]);
  double Bphi = MU0*real(theFields[2]);

  // Update the derivatives container
  double gamma = pow(1.0-(pow(vr,2.0)+pow(vz,2.0))/pow(C0,2.0),-0.5);
  f[0] = vr; // dr/dt
  f[1] = vz; // dz/dt
  f[2] = (theParticle->charge/(gamma*theParticle->mass))*(Er - vz*Bphi - vr*(vr*Er+vz*Ez)/pow(C0,2.0)); // dvr/dt
  f[3] = (theParticle->charge/(gamma*theParticle->mass))*(Ez + vr*Bphi - vz*(vr*Er+vz*Ez)/pow(C0,2.0)); // dvz/dt

  return GSL_SUCCESS;

} // end function odesys2D definition


#endif // PARTICLES_HPP_INCLUDED
