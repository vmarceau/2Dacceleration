//******************** FILE: 2DACC_BEAMCHARACTERIZATION.CPP ********************
//
// Description: Source file that contains algorithms to characterize the beam objects
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: April 2012
// Last update: April 2012
//
//******************************************************************************

// Standard header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <ctime>

// Project specific header files
#include "./beams.hpp"
#include "./constants.hpp"
#include "./general.hpp"


using namespace std;


// Main function
int main() {

  // FEATURE TO BE CHARACTERIZED
  // Density of electric energy versus radial coordinate r and time t at beam waist

  // Definition of the pulsed beams to be used
  /*double lambda0 = 800e-9;
  CNonparaxialTM01 npTM01(lambda0,1.0e12,124.366,277.0,0.0);
  CParaxialTM01 pTM01(lambda0,1.0e12,2.0e-6,10e-15,0.0);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out1;
  out1.open("./dat/test_beamcharacterization/pnpTM01_we_vs_rt.dat", ios::out);
  out1 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out1 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out1 << "#" << endl;
  out1 << "# DESCRIPTION" << endl;
  out1 << "# Electric field energy density at beam waist versus radial coordinate and time" << endl;
  out1 << "#" << endl;
  beam1->printparameters(out1);
  out1 << "#" << endl;
  beam2->printparameters(out1);
  out1 << "#" << endl;
  out1 << "# r/lambda0  t*1e-15  w_e(npTM01)  w_e(pTM01)" << endl;

  for (double r=0.0*lambda0; r<=20.0*lambda0; r+=lambda0/25.0) {
    for (double t=-30.0e-15; t<=30.0e-15; t+=1.0e-15/25.0) {
      vector<complex<double> > fields1 = beam1->fields(r,0.0*beam1->charlength_z,t);
      vector<complex<double> > fields2 = beam2->fields(r,0.0*beam1->charlength_z,t);
      double w_e1 = EPS0*(pow(real(fields1[0]),2.0) + pow(real(fields1[1]),2.0))/2.0;
      double w_e2 = EPS0*(pow(real(fields2[0]),2.0) + pow(real(fields2[1]),2.0))/2.0;
      out1 << r/lambda0 << "  " << t/1e-15 << "  " << w_e1 << "  " << w_e2 << endl;
    } // end for
    out1 << endl;
  } // end for

  // Close output file
  out1 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out1.close();*/


  // FEATURE TO BE CHARACTERIZED
  // Electromagnetic field components at beam waist versus radial coordinate

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  //double w0 = 2e-6;
  //double ka = 124.366;
  //double zR = PI*pow(w0,2.0)/lambda0;
  CNonparaxialTM01 npTM01(lambda0,2e9,1.2,20,PI/7.0);
  //CParaxialTM01 pTM01(lambda0,2.0e14,w0,10e-15,0.0933*PI);
  CBeam * beam1 = &npTM01;
  //CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out2;
  out2.open("./dat/test_beamcharacterization/pnpTM01_fields_vs_r.dat", ios::out);
  out2 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out2 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out2 << "#" << endl;
  out2 << "# DESCRIPTION" << endl;
  out2 << "# Electromagnetic field components at beam waist at t=0" << endl;
  out2 << "#" << endl;
  beam1->printparameters(out2);
  out2 << "#" << endl;
  //beam2->printparameters(out2);
  out2 << "#" << endl;
  out2 << "# r/lambda0  npTM01[Er,Ez,Hphi]  pTM01[Er,Ez,Hphi]" << endl;

  for (double r=-3.0*lambda0; r<=3.0*lambda0; r+=lambda0/10000) {
      vector<complex<double> > fields1 = beam1->fields(r,4.428*lambda0,1.437e-14);
      //vector<complex<double> > fields2 = beam2->fields(r,0.0,0.0);
      out2 << r/lambda0 << "  " << real(fields1[0]) << "  " << real(fields1[1]) << "  " << MU0*real(fields1[2]) << endl;
      //out2 << "  " << real(fields2[0]) << "  " << real(fields2[1]) << "  " << real(fields2[2]) << endl;
  } // end for

  // Close output file
  out2 << endl;
  out2 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out2.close();

  /*// FEATURE TO BE CHARACTERIZED
  // Electromagnetic field components at pulse peak versus r and z coordinates

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double w0 = 2e-6;
  double ka = 124.366;
  double zR = PI*pow(w0,2.0)/lambda0;
  CNonparaxialTM01 npTM01(lambda0,1.0e12,ka,277,0.0);
  CParaxialTM01 pTM01(lambda0,1.0e12,w0,10e-15,0.0);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out3;
  out3.open("./dat/test_beamcharacterization/pnpTM01_fields_rz.dat", ios::out);
  out3 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out3 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out3 << "#" << endl;
  out3 << "# DESCRIPTION" << endl;
  out3 << "# Electromagnetic field components at pulse peak versus r and z coordinates" << endl;
  out3 << "#" << endl;
  beam1->printparameters(out3);
  out3 << "#" << endl;
  beam2->printparameters(out3);
  out3 << "#" << endl;
  out3 << "# z/zR  r/w0  npTM01[Er,Ez,Hphi]  pTM01[Er,Ez,Hphi]" << endl;

  for (double z=-1.0*zR; z<=1.0*zR; z+=zR/500.0) {
    for (double r=0; r<=5.0*w0; r+=w0/50.0) {
      vector<complex<double> > fields1 = beam1->fields(r,z,0*zR/C0);
      vector<complex<double> > fields2 = beam2->fields(r,z,0*zR/C0);
      out3 << z/zR << "  " << r/w0 << "  " << abs(fields1[0]) << "  " << abs(fields1[1]) << "  " << abs(fields1[2]);
      out3 << "  " << abs(fields2[0]) << "  " << abs(fields2[1]) << "  " << abs(fields2[2]) << endl;
    } // end for
    out3 << endl;
  } // end for

  // Close output file
  out3 << endl;
  out3 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out3.close();*/


  /*// FEATURE TO BE CHARACTERIZED
  // Electromagnetic field components versus r and z coordinates averaged over time

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double w0 = 2e-6;
  double ka = 124.366;
  double zR = PI*pow(w0,2.0)/lambda0;
  CNonparaxialTM01 npTM01(lambda0,1.0e12,ka,277,0.0);
  CParaxialTM01 pTM01(lambda0,1.0e12,w0,10e-15,0.0);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out4;
  out4.open("./dat/test_beamcharacterization/pnpTM01_fields_rz.dat", ios::out);
  out4 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out4 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out4 << "#" << endl;
  out4 << "# DESCRIPTION" << endl;
  out4 << "# Electromagnetic field components versus r and z coordinates averaged over time" << endl;
  out4 << "#" << endl;
  beam1->printparameters(out4);
  out4 << "#" << endl;
  beam2->printparameters(out4);
  out4 << "#" << endl;
  out4 << "# r/w0  z/zR  npTM01[<E^2>]  pTM01[<E^2>]" << endl;

  // Results container
  int N = 201;
  vector<double> sample (N,0.0);
  vector<vector<double> > EsquaredP, EsquaredNP;
  EsquaredP.resize(N);
  EsquaredNP.resize(N);
  for (int n=0; n<N; n++) {
    EsquaredP[n] = sample;
    EsquaredNP[n] = sample;
  }

  for (double t=-100e-14; t<=100e-14; t+=1e-14) {
    cout << t << endl;
    for (int indexR=0; indexR<N; indexR+=1) {
      for (int indexZ=0; indexZ<N; indexZ+=1) {
        double r = indexR*w0/40.0;
        double z = -5.0*zR + indexZ*zR/20.0;
        vector<complex<double> > fields1 = beam1->fields(r,z,t);
        vector<complex<double> > fields2 = beam2->fields(r,z,t);
        EsquaredNP[indexR][indexZ] += pow(abs(fields1[0]),2.0) + pow(abs(fields1[1]),2.0);
        EsquaredP[indexR][indexZ] += pow(abs(fields2[0]),2.0) + pow(abs(fields2[1]),2.0);
      } // end for
    } // end for
  } // end for

  // Write results
  for (int indexR=0; indexR<N; indexR+=1) {
    for (int indexZ=0; indexZ<N; indexZ+=1) {
      double r = indexR*w0/40.0;
      double z = -5.0*zR + indexZ*zR/20.0;
      out4 << r/w0 << "  " << z/zR << "  " << EsquaredNP[indexR][indexZ] << "  " << EsquaredP[indexR][indexZ] << endl;
    } // end for
    out4 << endl;
  } // end for

  // Close output file
  out4 << endl;
  out4 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out4.close();*/


  /*// FEATURE TO BE CHARACTERIZED
  // Electromagnetic force components at pulse peak versus r and z coordinates

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double w0 = 2e-6;
  double ka = 124.366;
  double zR = PI*pow(w0,2.0)/lambda0;
  CNonparaxialTM01 npTM01(lambda0,2.0e14,ka,277,PI/2.0);
  CParaxialTM01 pTM01(lambda0,2.0e14,w0,10e-15,PI/2.0);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out5;
  out5.open("./dat/test_beamcharacterization/pnpTM01_force_rz.dat", ios::out);
  out5 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out5 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out5 << "#" << endl;
  out5 << "# DESCRIPTION" << endl;
  out5 << "# Electromagnetic force components at pulse peak versus r and z coordinates" << endl;
  out5 << "# Force computed for an electron with velocity 0.995c along the z axis" << endl;
  out5 << "#" << endl;
  beam1->printparameters(out5);
  out5 << "#" << endl;
  beam2->printparameters(out5);
  out5 << "#" << endl;
  out5 << "# z/zR  r/w0  npTM01[Fr,Fz]  pTM01[Fr,Fz,Hphi]" << endl;

  double vz = 0.999*C0;
  double gamma = pow(1.0-pow(vz/C0,2.0),-0.5);
  for (double z=-5.0*zR; z<=5.0*zR; z+=zR/20.0) {
    for (double r=0*w0; r<=5.0*w0; r+=w0/40.0) {
      vector<complex<double> > fields1 = beam1->fields(r,z,z/C0);
      vector<complex<double> > fields2 = beam2->fields(r,z,z/C0);
      double npFr = -QE*(real(fields1[0]) - vz*MU0*real(fields1[2]));
      double npFz = -QE*real(fields1[1]);
      double pFr = -QE*(real(fields2[0]) - vz*MU0*real(fields2[2]));
      double pFz = -QE*real(fields2[1]);
      out5 << z/zR << "  " << r/w0 << "  " << npFr << "  " << npFz << "  " << pFr << "  " << pFz << endl;
    } // end for
    out5 << endl;
  } // end for

  // Close output file
  out5 << endl;
  out5 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out5.close();*/


  /*// FEATURE TO BE CHARACTERIZED
  // Maxwell's equation error parameter map

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double k0 = 2.0*PI/lambda0;
  double w0 = 2.5*lambda0;
  double ka = k0*w0*sqrt(1+pow(0.5*k0*w0,2.0));
  CNonparaxialTM01 npTM01(lambda0,2.0e14,ka,277,0.0);
  CParaxialTM01 pTM01(lambda0,2.0e14,w0,10e-15,0.0);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Declaration of the output files
  ofstream outnp,outp;
  outnp.open("./dat/test_beamcharacterization/npTM01_errormap_ka124.dat", ios::out);
  outp.open("./dat/test_beamcharacterization/pTM01_errormap_ka124.dat", ios::out);
  double eps1 = beam1->errormap_fixedtime(0.0,10*lambda0,0.025*lambda0,10*lambda0,0.025*lambda0,1e-11,1e-11,1e-20,outnp);
  double eps2 = beam2->errormap_fixedtime(0.0,10*lambda0,0.025*lambda0,10*lambda0,0.025*lambda0,1e-11,1e-11,1e-20,outp);*/

  /*// FEATURE TO BE CHARACTERIZED
  // Pulse power check

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  CNonparaxialTM01 npTM01(lambda0,2.0e14,10,100,0.0);
  CBeam * beam1 = &npTM01;

  double sum1 = 0.0;
  double sum2 = 0.0;
  double deltar = lambda0/10000.0;
  for (double r=0; r<=1.0*lambda0; r+=deltar) {
      vector<complex<double> > fields = beam1->fields(r,0.0,0.0);
      sum1 += 2.0*PI*r*real(fields[0])*real(fields[2])*deltar;
      sum2 += PI*r*real(fields[0]*conj(fields[2]))*deltar;
  } // end for
  cout << "Max instantaneous power = " << sum1 << " W " << endl;
  cout << "Half instantaneous power = " << sum1/2.0 << " W" << endl;
  cout << "Average power = " << sum2 << " W" << endl;*/

  return 0;
}
