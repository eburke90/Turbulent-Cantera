/**
 *  @file stringUtils.h
 *       Contains declarations for string manipulation functions
 *       within Cantera.
 */

/*
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_STRINGUTILS_H
#define CT_STRINGUTILS_H

#include "ct_defs.h"

#include <string>

namespace Cantera {

  class Phase;
  class ThermoPhase;

  //! Convert a double into a c++ string
  /*!
   *  This routine doesn't assume a formatting. You
   *  must supply the formatting
   *
   * @param x double to be converted
   * @param fmt   Format to be used (printf style)
   */
  std::string fp2str(const double x, const std::string &fmt);

  //! Convert a double into a c++ string
  /*!
   * The default format to use is equivalent to the default
   * format used by printf's %g formatting.
   *
   * @param x double to be converted
   */
  std::string fp2str(const double x);

  //!  Convert an int to a string using a format converter
  /*!
   *  @param n          int to be converted
   *  @param fmt        format converter for an int int the printf command
   */
  std::string int2str(const int n, const std::string &fmt);

  //!  Convert an int to a string 
  /*!
   *  @param n          int to be converted
   */
  std::string int2str(const int n);

  //! Strip the leading and trailing white space
  //! from a string
  /*!
   *  The command isprint() is used to determine printable
   *  characters.
   *
   *    @param   s       Input string
   *    @return  Returns a copy of the string, stripped
   *             of leading and trailing white space
   */
  std::string stripws(const std::string &s);

  //! Strip non-printing characters wherever they are
  /*!
   *   @param s        Input string
   *   @return         Returns a copy of the string,
   *                   stripped of all non-printing characters.
   */
  std::string stripnonprint(const std::string &s);

  //! Cast a copy of a string to lower case
  /*!
   *   @param s        Input string
   *   @return         Returns a copy of the string,
   *                   with all characters lowercase.
   */
  std::string lowercase(const std::string &s);

  //! Parse a composition string into a map consisting of individual key:composition
  //! pairs.
  /*!
   *  The composition is a double.
   * Example
   *
   *  Input is
   *
   *    "fire:0   ice:1   snow:2"
   *
   *  Output is
   *             x["fire"] = 0
   *             x["ice"]  = 1
   *             x["snow"] = 2
   *
   *     @param ss   original string consisting of multiple key:composition
   *                 pairs on multiple lines
   *     @param x    Output map consisting of a composition
   *                 map, which is a string to double map
   */
  void parseCompString(const std::string &ss, Cantera::compositionMap& x);


  //! Parse a composition string into individual key:composition
  //! pairs
  /*!
   *
   *     @param ss   original string consisting of multiple key:composition
   *                 pairs on multiple lines
   *     @param w    Output vector consisting of single key:composition
   *                 items in each index.
   */
  void split(const std::string &ss, std::vector<std::string>& w);

  //! Interpret a string as a list of floats, and convert it to a vector
  //! of floats
  /*!
   *   @param str     String input vector
   *   @param a       Output pointer to a vector of floats
   *   @param delim   character delimiter. Defaults to a space 
   *   @return        Returns the number of floats found and converted
   */
  int fillArrayFromString(const std::string& str, doublereal* const a,
			  const char delim = ' ');

 
  //!  Generate a logfile name based on an input file name
  /*!
   *   It tries to find the basename. Then, it appends a .log
   *   to it.
   *
   *   @param infile      Input file name
   *
   *  @return Returns a logfile name
   */
  std::string logfileName(const std::string& infile);

  
  //! Get the file name without the path or extension
  /*!
   *   @param fullPath   Input file name consisting
   *                     of the full file name
   *
   *  @return Returns the basename
   */
  std::string getBaseName(const std::string& fullPath);

  //! Translate a string into one integer value
  /*!
   *  No error checking is done on the conversion. The c stdlib function
   *  atoi() is used.
   *
   *  @param val   String value of the integer
   *
   *  @return      Returns an integer
   */
  int intValue(std::string val);

  //! Translate a string into one doublereal value
  /*!
   *  No error checking is done on the conversion. The c stdlib function
   *  atof() is used.
   *
   *  @param val   String value of the double
   *
   *  @return      Returns a doublereal value
   */
  doublereal fpValue(std::string val);

  //! Translate a string into one doublereal value
  /*!
   *  Error checking is carried on the conversion. 
   *
   *  @param val   String value of the double
   *
   *  @return      Returns a doublereal value
   */
  doublereal fpValueCheck(std::string val);

  //! Line wrap a string via a copy operation
  /*!
   *   @param  s    Input string to be line wrapped
   *   @param  len  Length at which to wrap. The 
   *                default is 70.
   */
  std::string wrapString(const std::string &s, 
			 const int len=70);

  //! Routine strips off white space from a c character string
  /*!
   *     This routine strips off blanks and tabs (only leading and trailing
   *     characters) in 'str'.  On return, it returns the number of
   *     characters still included in the string (excluding the null character).
   *
   *      Comments are excluded -> All instances of the comment character, '!',
   *                               are replaced by NULL character thereby terminating
   *                               the string
   *
   *     Parameter list:
   *
   * @param  str   On output 'str' contains the same characters as on
   *               input except the leading and trailing white space and
   *               comments have been removed.
   */
  int stripLTWScstring(char str[]);

  //! Translate a char string into a single double
  /*!
   * atofCheck is a wrapper around the C stdlib routine atof().
   * It does quite a bit more error checking than atof() or
   * strtod(), and is quite a bit more restrictive.
   *
   *   First it interprets both E, e, d, and D as exponents.
   *   atof() only interprets e or E as an exponent character.
   *
   *   It only accepts a string as well formed if it consists as a 
   *   single token. Multiple words will produce an error message
   *
   *   It will produce an error for NAN and inf entries as well,
   *   in contrast to atof() or strtod().
   *   The user needs to know that a serious numerical issue
   *   has occurred.
   *
   *   It does not accept hexadecimal numbers.
   *
   *  @param dptr  pointer to the input c string
   *  @return      Returns the double
   *
   * On any error, it will throw a CanteraError signal.
   */
  doublereal atofCheck(const char * const dptr);


  //! Interpret one or two token string as a single double
  /*!
   *   This is similar to atof(). However, the second token
   *   is interpreted as an MKS units string and a conversion
   *   factor to MKS is applied.
   *
   *   Example
   *  " 1.0 atm"
   *
   *   results in the number 1.01325e5 
   *
   *   @param strSI string to be converted. One or two tokens
   *
   *   @return returns a converted double  
   */
  doublereal strSItoDbl(const std::string& strSI); 

  //! This function separates a string up into tokens
  //! according to the location of white space.
  /*!
   *  White space includes the new line character. tokens
   *  are stripped of leading and trailing white space.
   *
   *  The separate tokens are returned in a string vector, v.
   *
   *  @param oval   String to be broken up
   *  @param v     Output vector of tokens.
   */
  void tokenizeString(const std::string& oval,
                      std::vector<std::string>& v);

}

#endif
