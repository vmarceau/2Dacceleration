//******************** FILE: BEAMS.HPP ********************
//
// Description: This header file contains the beam class definition as well as the definition of its derived class.
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: April 2012
// Last update: July 2012
//
//*********************************************************

#ifndef BEAMS_HPP_INCLUDED
#define BEAMS_HPP_INCLUDED

// Standard header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

// GSL header files
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

// Project specific header files
#include "./constants.hpp"
#include "./general.hpp"

using namespace std;


//******************** BEAM CLASS DEFINITION ********************
// CLASS CBeam is an abstract base class. No instances of it can be created
class CBeam {
  protected:
    double lambda0; // Wavelength of maximum spectral amplitude [m]
    double k0; // Wavevector of maximum spectral amplitude [1/m]
    double omega0; // Angular frequency of maximum spectral amplitude [rad/s]
  public:
    double charlength_r; // Characteristic length in the radial direction
    double charlength_z; // Characteristic length in the axial direction
    virtual vector<complex<double> > fields (double,double,double) =0; // Electromagnetic fields
    virtual complex<double> Ezfield (double,double,double) =0; // Longitudinal electric field
    virtual void printparameters (ostream&) =0; // Print beam parameters
    double errormap_fixedtime (double,double,double,double,double,double,double,double); // Compute Maxwell's equations error
                                                                                         // parameter map at fixed time (no output)
    double errormap_fixedtime (double,double,double,double,double,double,double,double,ostream&); // Compute Maxwell's equations
                                                                                                  // error parameter map at fixed
                                                                                                  // time (full output)
};


// MEMBER FUNCTION errormap_fixedtime (double,double,double,double,double,double,double,double) definition
// Compute the Maxwell's equation error parameter map at a given time (no output)
//
// Input
//   - t: time
//   - maxr: upper boundary for the r coordinate
//   - stepr: grid step for the r coordinate
//   - maxz: upper and lower boundaries for the z coordinate
//   - stepz: grid step for the z coordinate
//   - hr: step size for the computation of the r derivative
//   - hz: step size for the computation of the z derivative
//   - ht: step size for the computation of the t derivative
//   - out: valid ofstream object (output stream)
// Output
//   - epsilonMmax: maximum value of the error parameter
//
double CBeam::errormap_fixedtime (double t, double maxr, double stepr, double maxz, double stepz, double hhr, double hhz, double hht) {

  // Fields and derivatives
  double Er, Ez, Bphi;
  double dErdz, dErdt, dEzdr, dEzdt, dBphidr, dBphidz, dBphidt;

  // Make sure that ht is a representable number
  volatile double tempt = 1e-15*(t/1e-15 + hht/1e-15);
  double ht = 1e-15*(tempt/1e-15 - t/1e-15);

  // Define result container
  double Ermax = 0.0;
  double Ezmax = 0.0;
  double epsilonMmax = 0.0;
  int sizer = static_cast<int>(maxr/stepr);
  int sizez = static_cast<int>(2*maxz/stepz + 1);
  double epsilonM[sizer][sizez];

  // Iterate over map and calculate the unscaled error parameter
  for (int nr = 0; nr<sizer; nr++) {

    double r = (nr+1)*stepr; // Current r coordinate
    volatile double tempr = r + hhr; // Make sure that hr is a representable number
    double hr = tempr - r;

    for (int nz = 0; nz<sizez; nz++) {

      double z = nz*stepz - maxz; // Current z coordinate
      volatile double tempz = z + hhz; // Make sur that hz is a representable number
      double hz = tempz - z;


      // Compute the fields around the given point
      vector<complex<double> > theFields = this->fields(r,z,t);
      vector<complex<double> > theFieldsRback = this->fields(r-hr,z,t);
      vector<complex<double> > theFieldsRfront = this->fields(r+hr,z,t);
      vector<complex<double> > theFieldsZback = this->fields(r,z-hz,t);
      vector<complex<double> > theFieldsZfront = this->fields(r,z+hz,t);
      vector<complex<double> > theFieldsTback = this->fields(r,z,(t/1e-15-ht/1e-15)*1e-15);
      vector<complex<double> > theFieldsTfront = this->fields(r,z,(t/1e-15+ht/1e-15)*1e-15);

      // Compute the derivatives with three point rule
      dErdz = (real(theFieldsZfront[0]) - real(theFieldsZback[0]))/(2*hz);
      dErdt = (real(theFieldsTfront[0]) - real(theFieldsTback[0]))/(2*ht);
      dEzdr = (real(theFieldsRfront[1]) - real(theFieldsRback[1]))/(2*hr);
      dEzdt = (real(theFieldsTfront[1]) - real(theFieldsTback[1]))/(2*ht);
      dBphidr = MU0*(real(theFieldsRfront[2]) - real(theFieldsRback[2]))/(2*hr);
      dBphidz = MU0*(real(theFieldsZfront[2]) - real(theFieldsZback[2]))/(2*hz);
      dBphidt = MU0*(real(theFieldsTfront[2]) - real(theFieldsTback[2]))/(2*ht);

      // Compute field amplitudes
      Er = real(theFields[0]);
      Ez = real(theFields[1]);
      Bphi = MU0*real(theFields[2]);

      // Compute unscaled error parameter (SI units)
      epsilonM[nr][nz] = abs(dErdz-dEzdr+dBphidt) + C0*sqrt(pow(dBphidz+dErdt/pow(C0,2.0),2.0) + pow(Bphi/r+dBphidr-dEzdt/pow(C0,2.0),2.0));

      // Check for maximum field amplitude
      if (abs(Er) > Ermax)
        Ermax = Er;
      if (abs(Ez) > Ezmax)
        Ezmax = Ez;

      // Check for maximum error parameter
      if (epsilonM[nr][nz] > epsilonMmax)
        epsilonMmax = epsilonM[nr][nz];

    } // end for
  } // end for


  // Compute scaling factor
  double scale;
  if (Ermax > Ezmax)
    scale = lambda0/Ermax;
  else
    scale = lambda0/Ezmax;

  return epsilonMmax*scale;

} // end member function errormap_fixedtime (double,double,double,double,double,double,double,double) definition


// MEMBER FUNCTION errormap_fixedtime (double,double,double,double,double,double,double,double,ostream&) definition
// Compute the Maxwell's equation error parameter map at a given time (full output)
//
// Input
//   - t: time
//   - maxr: upper boundary for the r coordinate
//   - stepr: grid step for the r coordinate
//   - maxz: upper and lower boundaries for the z coordinate
//   - stepz: grid step for the z coordinate
//   - hr: step size for the computation of the r derivative
//   - hz: step size for the computation of the z derivative
//   - ht: step size for the computation of the t derivative
//   - out: valid ofstream object (output stream)
// Output
//   - epsilonMmax: maximum value of the error parameter
//
double CBeam::errormap_fixedtime (double t, double maxr, double stepr, double maxz, double stepz, double hhr, double hhz, double hht, ostream& out) {

  // Fields and derivatives
  double Er, Ez, Bphi;
  double dErdz, dErdt, dEzdr, dEzdt, dBphidr, dBphidz, dBphidt;

  // Make sure that ht is a representable number
  volatile double tempt = 1e-15*(t/1e-15 + hht/1e-15);
  double ht = 1e-15*(tempt/1e-15 - t/1e-15);


  // Define result container
  double Ermax = 0.0;
  double Ezmax = 0.0;
  double epsilonMmax = 0.0;
  int sizer = static_cast<int>(maxr/stepr);
  int sizez = static_cast<int>(2*maxz/stepz + 1);
  double epsilonM[sizer][sizez];

  // Iterate over map and calculate the unscaled error parameter
  for (int nr = 0; nr<sizer; nr++) {

    double r = (nr+1)*stepr; // Current r coordinate
    volatile double tempr = r + hhr; // Make sure that hr is a representable number
    double hr = tempr - r;

    for (int nz = 0; nz<sizez; nz++) {

      double z = nz*stepz - maxz; // Current z coordinate
      volatile double tempz = z + hhz; // Make sur that hz is a representable number
      double hz = tempz - z;

      // Compute the fields around the given point
      vector<complex<double> > theFields = this->fields(r,z,t);
      vector<complex<double> > theFieldsRback = this->fields(r-hr,z,t);
      vector<complex<double> > theFieldsRfront = this->fields(r+hr,z,t);
      vector<complex<double> > theFieldsZback = this->fields(r,z-hz,t);
      vector<complex<double> > theFieldsZfront = this->fields(r,z+hz,t);
      vector<complex<double> > theFieldsTback = this->fields(r,z,(t/1e-15-ht/1e-15)*1e-15);
      vector<complex<double> > theFieldsTfront = this->fields(r,z,(t/1e-15+ht/1e-15)*1e-15);

      // Compute the derivatives with three point rule
      dErdz = (real(theFieldsZfront[0]) - real(theFieldsZback[0]))/(2*hz);
      dErdt = (real(theFieldsTfront[0]) - real(theFieldsTback[0]))/(2*ht);
      dEzdr = (real(theFieldsRfront[1]) - real(theFieldsRback[1]))/(2*hr);
      dEzdt = (real(theFieldsTfront[1]) - real(theFieldsTback[1]))/(2*ht);
      dBphidr = MU0*(real(theFieldsRfront[2]) - real(theFieldsRback[2]))/(2*hr);
      dBphidz = MU0*(real(theFieldsZfront[2]) - real(theFieldsZback[2]))/(2*hz);
      dBphidt = MU0*(real(theFieldsTfront[2]) - real(theFieldsTback[2]))/(2*ht);

      // Compute field amplitudes
      Er = real(theFields[0]);
      Ez = real(theFields[1]);
      Bphi = MU0*real(theFields[2]);

      // Compute unscaled error parameter (SI units)
      epsilonM[nr][nz] = abs(dErdz-dEzdr+dBphidt) + C0*sqrt(pow(dBphidz+dErdt/pow(C0,2.0),2.0) + pow(Bphi/r+dBphidr-dEzdt/pow(C0,2.0),2.0));

      // Check for maximum field amplitude
      if (abs(Er) > Ermax)
        Ermax = Er;
      if (abs(Ez) > Ezmax)
        Ezmax = Ez;

      // Check for maximum error parameter
      if (epsilonM[nr][nz] > epsilonMmax)
        epsilonMmax = epsilonM[nr][nz];

    } // end for
  } // end for


  // Compute scaling factor
  double scale;
  if (Ermax > Ezmax)
    scale = lambda0/Ermax;
  else
    scale = lambda0/Ezmax;

  // Output results
  out << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out << "#" << endl;
  out << "# DESCRIPTION" << endl;
  out << "# Maxwell's equations error parameter map at t = " << t << endl;
  out << "#" << endl;
  this->printparameters(out);
  out << "#" << endl;
  out << "# FIELD PARAMETERS" << endl;
  out << "# Ermax = " << Ermax << " , Ezmax = " << Ezmax;
  out << " , scaling factor = " << scale << " , epsilonMmax = " << epsilonMmax*scale << endl;
  out << "#" << endl;
  out << "# r/lambda0  z/lambda0  epsilonM" << endl;


  for (int nr = 0; nr<sizer; nr++) {
    double r = (nr+1)*stepr; // Current r coordinate
    for (int nz = 0; nz<sizez; nz++) {
      double z = nz*stepz - maxz; // Current z coordinate
      out << r/lambda0 << "  " << z/lambda0 << "  " << epsilonM[nr][nz]*scale << endl;
    } // end for
    out << endl;
  } // end for

  return epsilonMmax*scale;

} // end member function errormap_fixedtime (double,double,double,double,double,double,double,double,ostream&) definition



//******************** NONPARAXIALTM01 CLASS DEFINITION ********************
// CLASS CNonparaxialTM01 is derived from class Cbeam
// Class whose instances represent ultrashort and nonparaxial TM01 pulsed beams
class CNonparaxialTM01: public CBeam {
  protected:
    double E0; // Amplitude parameter [V/m]
    double a; // Confocal parameter [m]
    double s; // Spectral width parameter []
    double phi0; // Pulse phase [rad]
    vector<complex<double> > G0; // G0 factors container []
    double w0; // Beam waist size at wavelength lambda0 [m]
    double zR; // Rayleigh range at wavelength lambda0 [m]
    double Sz2PIr (double); // Calculate the z-component of the Poynting vector times 2*PI*r at beam waist
    static double static_Sz2PIr(double,void *); // Static function used to call the Sz2PIr function for a given object
  public:
    CNonparaxialTM01 (double,double,double); // Constructor for DUMMY PULSE
    CNonparaxialTM01 (double,double,double,double,double); // Constructor for FULLY PARAMETRIZED PULSE
    vector<complex<double> > fields (double,double,double); // Compute the electromagnetic fields
    complex<double> Ezfield (double,double,double); // Compute the longitudinal electric field Ez
    double power (); // Calculate the peak power of the pulse
    void reset_phase(double); // Reset the pulse phase phi0 to a given value
    void printparameters (ostream&); // Print the beam parameters to output stream
};


// CONSTRUCTOR CNonparaxialTM01 (dummy case) definition
// The lambda0, a, and s member attributes are set to the prescribed values, while E0=1 and phi0=0 are automatically set.
//
// Input
//   - wavelength: value for the lambda0 member attribute
//   - ka: value for the a member attribute times k0
//   - spectral: value for the s member attribute
//
CNonparaxialTM01::CNonparaxialTM01 (double wavelength, double ka, double spectral){

  // Set the values of the pulse parameters
  lambda0 = wavelength;
  k0 = 2.0*PI/lambda0;
  omega0 = k0*C0;
  a = ka/k0;
  s = spectral;
  E0 = 1.0;
  phi0 = 0.0;
  charlength_r = a;
  charlength_z = a;

  // Set the G0 factors
  // The G0 factors are defined as (LaTeX notation):
  // G_0^{(n)} = e^{-j\phi_0} \frac{\Gamma(s+n+1)}{\Gamma(s+1)} (\frac{j\omega_0}{s})^n
  G0.resize(3);
  G0[0] = exp(-I*phi0);
  G0[1] = G0[0]*(s+1.0)*(I*omega0/s);
  G0[2] = G0[1]*(s+2.0)*(I*omega0/s);

  // Set beam waist size w0 and Rayleigh range zR
  w0 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0*a,2.0))-1)/k0;
	zR = k0*pow(w0,2.0)/2;

} // End constructor CNonparaxialTM01 (dummy case) definition


// CONSTRUCTOR CNonparaxialTM01 (fully parametrized case) definition
// The lambda0, E0, a, s, and phi0 member attributes are set to the prescribed values.
//
// Input
//   - wavelength: value for the lambda0 member attribute
//   - peakpower: peak power of the pulse, used to set the E0 member attribute
//   - ka: value for the a member attribute times k0
//   - spectral: value for the s member attribute
//   - phase: value for the phi0 member attribute
//
CNonparaxialTM01::CNonparaxialTM01 (double wavelength, double peakpower, double ka, double spectral, double phase){

  // Create dummy pulse equivalent
  CNonparaxialTM01 npTM01dummy (wavelength,ka,spectral);

  // Set the values of the pulse parameters
  lambda0 = wavelength;
  k0 = 2.0*PI/lambda0;
  omega0 = k0*C0;
  E0 = sqrt(peakpower/npTM01dummy.power());
  a = ka/k0;
  s = spectral;
  phi0 = phase;
  charlength_r = a;
  charlength_z = a;

  // Set the G0 factors
  // The G0 factors are defined as (LaTeX notation):
  // G_0^{(n)} = e^{-j\phi_0} \frac{\Gamma(s+n+1)}{\Gamma(s+1)} (\frac{j\omega_0}{s})^n
  G0.resize(3);
  G0[0] = exp(-I*phi0);
  G0[1] = G0[0]*(s+1.0)*(I*omega0/s);
  G0[2] = G0[1]*(s+2.0)*(I*omega0/s);


  // Set beam waist size w0 and Rayleigh range zR

  w0 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0*a,2.0))-1)/k0;
	zR = k0*pow(w0,2.0)/2;

} // End constructor CNonparaxialTM01 (fully parametrized case) definition


// MEMBER FUNCTION fields(double,double,double) definition
// The electromagnetic fields of the nonparaxial TM01 pulse are calculated according to the closed-form expressions presented
// in the PhD thesis of Alexandre April.
//
// Input
//   - r: radial coordinate
//   - z: longitudinal coordinate
//   - t: time
// Output
//   - npTM01fields[0]: analytic (complex) value of the Er field component
//   - npTM01fields[1]: analytic (complex) value of the Ez field component
//   - npTM01fields[2]: analytic (complex) value of the Hphi field component
//
vector<complex<double> > CNonparaxialTM01::fields(double r,double z,double t)
{

  // Check if the fields are to be computed near a possibly singular point
  if ( abs(z) < 1e-10 && r/a > 0.95 && r/a < 1.1 ) { // If yes, proceed to an interpolation of the field at the desired position

    // Precomputation of the psi0 amplitude
    double psi0 = -exp(0.5)*pow(a*C0/omega0,1.5)*E0;
    cout << psi0 << endl;

    // Compute the list of field points necessary for the interpolation
    double rn[10] = {0.8*a,0.85*a,0.9*a,0.925*a,0.95*a,1.1*a,1.15*a,1.2*a,1.25*a,1.3*a};
    double REr[10], IEr[10], REz[10], IEz[10], RHphi[10], IHphi[10];
    for (int n=0; n<=9; n++) {
      complex<double> ztilde = z + I*a;
      complex<double> Rtilde = sqrt(pow(rn[n],2.0) + pow(ztilde,2.0));
      complex<double> sintheta = rn[n]/Rtilde;
      complex<double> costheta = ztilde/Rtilde;
      complex<double> tplus = t + Rtilde/C0 + I*a/C0;
      complex<double> tminus = t - Rtilde/C0 + I*a/C0;
      complex<double> fplusbase = 1.0-I*omega0*tplus/s;
      complex<double> fplus0 = pow(fplusbase,-(s+1.0));
      complex<double> fplus1 = fplus0/fplusbase;
      complex<double> fplus2 = fplus1/fplusbase;
      complex<double> fminusbase = 1.0-I*omega0*tminus/s;
      complex<double> fminus0 = pow(fminusbase,-(s+1.0));
      complex<double> fminus1 = fminus0/fminusbase;
      complex<double> fminus2 = fminus1/fminusbase;
      complex<double> Gminus0 = G0[0]*(fplus0 - fminus0);
      complex<double> Gplus1 = G0[1]*(fplus1 + fminus1);
      complex<double> Gminus1 = G0[1]*(fplus1 - fminus1);
      complex<double> Gplus2 = G0[2]*(fplus2 + fminus2);
      complex<double> Gminus2 = G0[2]*(fplus2 - fminus2);
      complex<double> Gfactor = Gminus0/pow(Rtilde,2.0) - Gplus1/(C0*Rtilde) + Gminus2/(3.0*pow(C0,2.0));
      // Radial electric field component Er
      complex<double> Er = 3.0*psi0*costheta*sintheta*Gfactor/Rtilde;
      REr[n] = real(Er);

      IEr[n] = imag(Er);
      // Longitudinal electric field component Ez
      complex<double> Ez = (2.0*psi0/Rtilde)*( 0.5*(3.0*pow(costheta,2.0)-1.0)*Gfactor - Gminus2/(3*pow(C0,2.0))  );
      REz[n] = real(Ez);
      IEz[n] = imag(Ez);
      // Azimutal magnetic field component Hphi
      complex<double> Hphi = psi0*sintheta*EPS0*( Gminus1/Rtilde - Gplus2/C0 )/Rtilde;
      RHphi[n] = real(Hphi);
      IHphi[n] = imag(Hphi);
    } // end for

    // Interpolate using cubic spline from GSL scientific library
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc ();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
    gsl_interp_accel *acc3 = gsl_interp_accel_alloc ();
    gsl_interp_accel *acc4 = gsl_interp_accel_alloc ();
    gsl_interp_accel *acc5 = gsl_interp_accel_alloc ();
    gsl_interp_accel *acc6 = gsl_interp_accel_alloc ();
    gsl_spline *spline1 = gsl_spline_alloc (gsl_interp_cspline, 10);
    gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, 10);
    gsl_spline *spline3 = gsl_spline_alloc (gsl_interp_cspline, 10);
    gsl_spline *spline4 = gsl_spline_alloc (gsl_interp_cspline, 10);
    gsl_spline *spline5 = gsl_spline_alloc (gsl_interp_cspline, 10);
    gsl_spline *spline6 = gsl_spline_alloc (gsl_interp_cspline, 10);
    gsl_spline_init (spline1, rn, REr, 10);
    gsl_spline_init (spline2, rn, IEr, 10);
    gsl_spline_init (spline3, rn, REz, 10);
    gsl_spline_init (spline4, rn, IEz, 10);
    gsl_spline_init (spline5, rn, RHphi, 10);
    gsl_spline_init (spline6, rn, IHphi, 10);
    double RErint = gsl_spline_eval (spline1, r, acc1);
    double IErint = gsl_spline_eval (spline2, r, acc2);
    double REzint = gsl_spline_eval (spline3, r, acc3);
    double IEzint = gsl_spline_eval (spline4, r, acc4);
    double RHphiint = gsl_spline_eval (spline5, r, acc5);
    double IHphiint = gsl_spline_eval (spline6, r, acc6);
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_spline_free (spline3);
    gsl_spline_free (spline4);
    gsl_spline_free (spline5);
    gsl_spline_free (spline6);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_interp_accel_free (acc3);
    gsl_interp_accel_free (acc4);
    gsl_interp_accel_free (acc5);
    gsl_interp_accel_free (acc6);

    // Return the components of the electromagnetic field
    vector<complex<double> > npTM01fields (3,complex<double>(0,0));
    npTM01fields[0] = complex<double>(RErint,IErint);
    npTM01fields[1] = complex<double>(REzint,IEzint);
    npTM01fields[2] = complex<double>(RHphiint,IHphiint);
    return npTM01fields;

  } // end if
  else { // If not, compute the fields directly according to the analytical expressions

    // Compute various quantities that will be necessary in the computation of the field components
    double psi0 = -exp(0.5)*pow(a*C0/omega0,1.5)*E0;
    complex<double> ztilde = z + I*a;
    complex<double> Rtilde = sqrt(pow(r,2.0) + pow(ztilde,2.0));
    complex<double> sintheta = r/Rtilde;
    complex<double> costheta = ztilde/Rtilde;
    complex<double> tplus = t + Rtilde/C0 + I*a/C0;
    complex<double> tminus = t - Rtilde/C0 + I*a/C0;
    complex<double> fplusbase = 1.0-I*omega0*tplus/s;
    complex<double> fplus0 = pow(fplusbase,-(s+1.0));
    complex<double> fplus1 = fplus0/fplusbase;
    complex<double> fplus2 = fplus1/fplusbase;
    complex<double> fminusbase = 1.0-I*omega0*tminus/s;
    complex<double> fminus0 = pow(fminusbase,-(s+1.0));
    complex<double> fminus1 = fminus0/fminusbase;
    complex<double> fminus2 = fminus1/fminusbase;
    complex<double> Gminus0 = G0[0]*(fplus0 - fminus0);
    complex<double> Gplus1 = G0[1]*(fplus1 + fminus1);
    complex<double> Gminus1 = G0[1]*(fplus1 - fminus1);
    complex<double> Gplus2 = G0[2]*(fplus2 + fminus2);
    complex<double> Gminus2 = G0[2]*(fplus2 - fminus2);
    complex<double> Gfactor = Gminus0/pow(Rtilde,2.0) - Gplus1/(C0*Rtilde) + Gminus2/(3.0*pow(C0,2.0));

    // Radial electric field component Er
    complex<double> Er = 3.0*psi0*costheta*sintheta*Gfactor/Rtilde;

    // Longitudinal electric field component Ez
    complex<double> Ez = (2.0*psi0/Rtilde)*( 0.5*(3.0*pow(costheta,2.0)-1.0)*Gfactor - Gminus2/(3*pow(C0,2.0))  );

    // Azimutal magnetic field component Hphi
    complex<double> Hphi = psi0*sintheta*EPS0*( Gminus1/Rtilde - Gplus2/C0 )/Rtilde;


    // Return the components of the electromagnetic field
    vector<complex<double> > npTM01fields (3,complex<double>(0,0));
    npTM01fields[0] = Er;
    npTM01fields[1] = Ez;
    npTM01fields[2] = Hphi;

    return npTM01fields;

  } // end else

} // End member function fields(double,double,double) definition


// MEMBER FUNCTION Ezfield(double,double,double) definition
// The longitudinal electric field of the nonparaxial TM01 pulse is calculated according to the closed-form expression presented
// in the PhD thesis of Alexandre April.
//
// Input
//   - r: radial coordinate
//   - z: longitudinal coordinate
//   - t: time
// Output
//   - Ez: analytic (complex) value of the Ez field component
//
complex<double> CNonparaxialTM01::Ezfield(double r,double z,double t)
{

  // Compute various quantities that will be necessary in the computation of the field component
  double psi0 = -exp(0.5)*pow(a*C0/omega0,1.5)*E0;
  complex<double> ztilde = z + I*a;
  complex<double> Rtilde = sqrt(pow(r,2.0) + pow(ztilde,2.0));
  complex<double> costheta = ztilde/Rtilde;
  complex<double> tplus = t + Rtilde/C0 + I*a/C0;
  complex<double> tminus = t - Rtilde/C0 + I*a/C0;
  complex<double> fplusbase = 1.0-I*omega0*tplus/s;
  complex<double> fplus0 = pow(fplusbase,-(s+1.0));
  complex<double> fplus1 = fplus0/fplusbase;
  complex<double> fplus2 = fplus1/fplusbase;
  complex<double> fminusbase = 1.0-I*omega0*tminus/s;
  complex<double> fminus0 = pow(fminusbase,-(s+1.0));
  complex<double> fminus1 = fminus0/fminusbase;
  complex<double> fminus2 = fminus1/fminusbase;
  complex<double> Gminus0 = G0[0]*(fplus0 - fminus0);
  complex<double> Gplus1 = G0[1]*(fplus1 + fminus1);
  complex<double> Gminus2 = G0[2]*(fplus2 - fminus2);
  complex<double> Gfactor = Gminus0/pow(Rtilde,2.0) - Gplus1/(C0*Rtilde) + Gminus2/(3.0*pow(C0,2.0));

  // Longitudinal electric field component Ez
  complex<double> Ez = (2.0*psi0/Rtilde)*( 0.5*(3.0*pow(costheta,2.0)-1.0)*Gfactor - Gminus2/(3*pow(C0,2.0))  );

  return Ez;

} // End member function Ezfield(double,double,double) definition


// MEMBER FUNCTION power(void) definition
// Returns the peak power of the pulse by integrating the z-component of the Poynting vector at z=0,t=0 in the transverse plane
//
// Input
//     none
// Output:
//   - Ppeak: peak power of the pulse
//
double CNonparaxialTM01::power() {

  // Create pulse equivalent with unitary amplitude and null phase
  CNonparaxialTM01 npTM01equiv (lambda0,k0*a,s);

  // Create integration workspace
  gsl_integration_workspace * intworkspace = gsl_integration_workspace_alloc(10000);

  // Define the function to be integrated and the result containers
  double Ppeak, error;
  double upperbound = 20.0*npTM01equiv.w0;
  gsl_function F;
  F.function = &CNonparaxialTM01::static_Sz2PIr;
  F.params = &npTM01equiv;

  // Perform numerical integration using GSL QAGS adaptive integration with known singular points
  gsl_integration_qags(&F,0,upperbound,0,1e-6,10000,intworkspace,&Ppeak,&error);

  // Scale power to the pulse under investigation
  Ppeak = pow(E0,2.0)*Ppeak;

  // Free allocated workspace
  gsl_integration_workspace_free(intworkspace);

  // Return peak power of the pulse
  return Ppeak;

} // end member function power(void) definition


// MEMBER FUNCTION static_Sz2PIr(double,void*) definition
// Static member function used to call the member function Sz2PIr(double,void*)
//
// Input
//   - r: radial coordinate
//   - params: Pointer to a CNonparaxial object
// Output
//   - Sz: z-component of the Poynting vector times 2*PI*r at z=0, t=0.
//
double CNonparaxialTM01::static_Sz2PIr(double r, void * params) {

  CNonparaxialTM01 * beamPtr = static_cast<CNonparaxialTM01*>(params);
  return beamPtr->Sz2PIr(r);

} // end member function static_Sz2PIr(double,void*) definition


// MEMBER FUNCTION Sz2PIr(double) definition
// Integrand function used in the power(void) member function
//
// Input
//   - r: radial coordinate
// Output
//   - Sz: z-component of the Poynting vector times 2*PI*r at z=0, t=0.
//
double CNonparaxialTM01::Sz2PIr(double r) {

  // Get electromagnetic fields
  vector<complex<double> > thefields = this->fields(r,0.0,0.0);

  // Compute Sz times 2*PI*r
  double Sz = 2.0*PI*r*real(thefields[0])*real(thefields[2]);

  return Sz;

} // end member function Sz2PIr(double) definition


// MEMBER FUNCTION reset_phase(double) definition
// Reset the pulse phase to a given value
//
// Input
//   - newphi0: new value of the pulse phase
// Output
//     none
//
void CNonparaxialTM01::reset_phase(double newphi0) {

  phi0 = newphi0;

  // Reset the G0 factors
  G0[0] = exp(-I*phi0);
  G0[1] = G0[0]*(s+1.0)*(I*omega0/s);
  G0[2] = G0[1]*(s+2.0)*(I*omega0/s);

} // end member function reset_phase(double) definition


// MEMBER FUNCTION printparameters(ofstream) definition
// Print the pulsed beam parameters to a given output stream
//
// Input
//   - out: valid ofstream object (output stream)
// Output
//     none
//
void CNonparaxialTM01::printparameters(ostream& out) {

  out << "# BEAM PARAMETERS" << endl;
  out << "# Type: Ultrashort and nonparaxial TM01 pulsed beam" << endl;
  out << "# Wavelength of maximum spectral amplitude: lambda0 = " << lambda0 << endl;
  out << "# Peak power: Ppeak = " << this->power() << endl;
  out << "# Amplitude parameter: E0 = " << E0 << endl;
  out << "# Confocal parameter: k0a = " << 2.0*PI*a/lambda0 << endl;
  out << "# Spectral width parameter: s = " << s << endl;
  out << "# Pulse phase: phi0 = " << phi0 << endl;
  out << "# Beam waist size at wavelength lambda0: w0 = " << w0 << endl;
  out << "# Rayleigh range at wavelength lambda0: zR = " << zR << endl;

} // end member function printparameters(ofstream) definition



//******************** PARAXIALTM01 CLASS DEFINITION ********************
// CLASS CParaxialTM01 is derived from class Cbeam
// Class whose instances represent paraxial few-cycle TM01 pulsed beams
class CParaxialTM01: public CBeam {
  protected:
    double E0; // Amplitude parameter [V/m]
    double w0; // Beam waist size [m]
    double zR; // Rayleigh range [m]
    double T; // Pulse duration [s]
    double phi0; // Pulse phase [rad]
  public:
    CParaxialTM01 (double,double,double,double,double); // Constructor for FULLY PARAMETRIZED PULSE
    vector<complex<double> > fields (double,double,double); // Compute the electromagnetic fields
    complex<double> Ezfield (double,double,double); // Compute the longitudinal electric field Ez
    double power (void); // Calculate the peak power of the pulse
    void reset_phase (double); // Reset the pulse phase phi0 to a given value
    void printparameters (ostream&); // Print the beam parameters to output stream
};


// CONSTRUCTOR CParaxialTM01 (fully parametrized case) definition
// The lambda0, E0, a, s, and phi0 member attributes are set to the prescribed values.
//
// Input
//   - wavelength: value for the lambda0 member attribute
//   - peakpower: peak power of the pulse, used to set the E0 member attribute
//   - waistsize: beam waist size w0 of the pulse
//   - duration: duration T of the pulse
//   - phase: value for the phi0 member attribute
//
CParaxialTM01::CParaxialTM01 (double wavelength, double peakpower, double waistsize, double duration, double phase){

  // Set the values of the pulse parameters
  lambda0 = wavelength;
  k0 = 2.0*PI/lambda0;
  omega0 = k0*C0;
  w0 = waistsize;
  zR = k0*pow(w0,2.0)/2.0;
  E0 = sqrt(k0*ETA0*exp(-1.0)*peakpower/(PI*zR));
  T = duration;
  phi0 = phase;
  charlength_r = zR;
  charlength_z = w0;

} // End constructor CParaxialTM01 (fully parametrized case) definition


// MEMBER FUNCTION fields(double,double,double) definition
// The electromagnetic fields of the paraxial TM01 pulse are calculated according to the closed-form expressions presented
// in the MSc thesis of Pierre-Louis Fortin.
//
// Input
//   - r: radial coordinate
//   - z: longitudinal coordinate
//   - t: time
// Output
//   - pTM01fields[0]: analytic (complex) value of the Er field component
//   - pTM01fields[1]: analytic (complex) value of the Ez field component
//   - pTM01fields[2]: analytic (complex) value of the Hphi field component
//
vector<complex<double> > CParaxialTM01::fields(double r,double z,double t)
{

  // Compute various quantities that will be necessary in the computation of the electric field components
  complex<double> qtilde = z+I*zR;
  complex<double> amplitude = E0*exp(0.5)*pow(I*zR/qtilde,2.0)*sqrt(k0/zR)*exp(-I*k0*pow(r,2.0)/(2.0*qtilde));
  complex<double> pulse_shape = exp(-pow(t-z/C0,2.0)/pow(T,2.0));
  complex<double> propagator = exp(I*(omega0*t - k0*z - phi0));
  complex<double> common_factor = amplitude*pulse_shape*propagator;

  // Radial electric field component
  complex<double> Er = common_factor*r;

  // Longitudinal electric field component
  complex<double> Ez = common_factor*(-2.0*I/k0)*(1.0-I*k0*pow(r,2.0)/(2.0*qtilde));

  // Azimutal magnetic field component
  complex<double> Hphi = Er/ETA0;

  // Return the components of the electromagnetic field
  vector<complex<double> > pTM01fields (3,complex<double>(0,0));
  pTM01fields[0] = Er;
  pTM01fields[1] = Ez;
  pTM01fields[2] = Hphi;
  return pTM01fields;

} // End member function fields(double,double,double) definition


// MEMBER FUNCTION Ezfield(double,double,double) definition
// The longitudinal electric field of the paraxial TM01 pulse is calculated according to the closed-form expression presented
// in the MSc thesis of Pierre-Louis Fortin.
//
// Input
//   - r: radial coordinate
//   - z: longitudinal coordinate
//   - t: time
// Output
//   - Ez: analytic (complex) value of the Ez field component
//
complex<double> CParaxialTM01::Ezfield(double r,double z,double t)
{

  // Compute various quantities that will be necessary in the computation of the electric field component
  complex<double> qtilde = z+I*zR;
  complex<double> amplitude = E0*exp(0.5)*pow(I*zR/qtilde,2.0)*sqrt(k0/zR)*exp(-I*k0*pow(r,2.0)/(2.0*qtilde));
  complex<double> pulse_shape = exp(-pow(t-z/C0,2.0)/pow(T,2.0));
  complex<double> propagator = exp(I*(omega0*t - k0*z - phi0));
  complex<double> common_factor = amplitude*pulse_shape*propagator;

  // Longitudinal electric field component
  complex<double> Ez = common_factor*(-2.0*I/k0)*(1.0-I*k0*pow(r,2.0)/(2.0*qtilde));

  return Ez;

} // End member function Ezfield(double,double,double) definition


// MEMBER FUNCTION power(void) definition
// Returns the peak power of the paraxial pulse
//
// Input
//     none
// Output:
//   - Ppeak: peak power of the pulse
//
double CParaxialTM01::power(void) {

  // Compute the analytical expression for the peak power of the pulse
  double Ppeak = (PI/ETA0)*pow(E0,2.0)*exp(1.0)*(zR/k0);

  return Ppeak;

} // end member function power(void) definition


// MEMBER FUNCTION reset_phase(double) definition
// Reset the pulse phase to a given value
//
// Input
//   - newphi0: new value of the pulse phase
// Output
//     none
//
void CParaxialTM01::reset_phase(double newphi0) {

  phi0 = newphi0;

} // end member function reset_phase(double) definition


// MEMBER FUNCTION printparameters(ofstream) definition
// Print the pulsed beam parameters to a given output stream
//
// Input
//   - out: valid ofstream object (output stream)
// Output
//     none
//
void CParaxialTM01::printparameters(ostream& out) {

  out << "# BEAM PARAMETERS" << endl;
  out << "# Type: Few-cycle paraxial TM01 pulsed beam" << endl;
  out << "# Central wavelength: lambda0 = " << lambda0 << endl;
  out << "# Peak power: Ppeak = " << this->power() << endl;
  out << "# Amplitude parameter: E0 = " << E0 << endl;
  out << "# Beam waist size: w0 = " << w0 << endl;
  out << "# Rayleigh range: zR = " << zR << endl;
  out << "# Duration: T = " << T << endl;
  out << "# Pulse phase: phi0 = " << phi0 << endl;

} // end member function printparameters(ofstream) definition


#endif // BEAMS_HPP_INCLUDED
