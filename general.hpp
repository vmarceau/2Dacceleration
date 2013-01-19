//******************** FILE: GENERAL.HPP ********************
//
// Description: This header file contains the definition of functions of general purpose
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: May 2012
// Last update: May 2012
//
//*********************************************************

#ifndef GENERAL_HPP_INCLUDED
#define GENERAL_HPP_INCLUDED

// Standard header files
#include <vector>
#include <ctime>

using namespace std;

// Function GetDate() definition
// Returns a string with the current date in YYYY-MM-DD format
//
// Input
//     none
// Output
//   - strTime: Current date in YYYY-MM-DD format
//
string GetDate() {
  char cptime[50];
  time_t now = time(NULL);
  strftime(cptime, 50, "%F", localtime(&now)); //uses short month name
  string strTime = cptime;
  return strTime;
} // End function GetDate() definition


// Function GetTime() definition
// Returns a string with the current time in HH:MM:SS format
//
// Input
//     none
// Output
//   - strTime: Current time in HH:MM:SS format
//
string GetTime() {
  char cptime[50];
  time_t now = time(NULL);
  strftime(cptime, 50, "%H:%M:%S", localtime(&now));
  string strTime = cptime;
  return strTime;
} // End function GetTime() definition


#endif // GENERAL_HPP_INCLUDED

