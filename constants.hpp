//******************** FILE: CONSTANTS.HPP ********************
//
// Description: This header file contains the definition of constants frequently used in the project.
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: April 2012
// Last update: April 2012
//
//*************************************************************

#ifndef CONSTANTS_HPP_INCLUDED
#define CONSTANTS_HPP_INCLUDED

// Standard header files
#include <complex>

// GSL header files
#include <gsl/gsl_const_mksa.h>


// Imaginary unit []
#ifndef I
#define I (complex<double>(0,1.0))
#endif

// PI constant []
#ifndef PI
#define PI M_PI
#endif

// Velocity of light in free space [m/s]
#ifndef C0
#define C0 GSL_CONST_MKSA_SPEED_OF_LIGHT
#endif

// Permittivity of free space [F/m]
#ifndef EPS0
#define EPS0 GSL_CONST_MKSA_VACUUM_PERMITTIVITY
#endif

// Permeability of free space [s^2/F*m]
#ifndef MU0
#define MU0 GSL_CONST_MKSA_VACUUM_PERMEABILITY
#endif

// Impedance of free space [Ohms]
#ifndef ETA0
#define ETA0 (MU0*C0)
#endif

// Electron mass [kg]
#ifndef ME
#define ME GSL_CONST_MKSA_MASS_ELECTRON
#endif

// Proton mass [kg]
#ifndef MP
#define MP GSL_CONST_MKSA_MASS_PROTON
#endif

// Unit charge [C]
#ifndef QE
#define QE GSL_CONST_MKSA_ELECTRON_CHARGE
#endif


#endif // CONSTANTS_HPP_INCLUDED
