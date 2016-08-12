/*
 * Potential.cc
 *
 *  Created on: Jul 8, 2016
 *      Author: dbpc
 */

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

#include "Potential.h"

/**
 * Header Interface for Class Point
 *
class Potential
{
  public:
	int p;
    int DEFAULT_P = 12;

    Potential() { this->p = DEFAULT_P; };
	Potential(int p)
	 :
     p(p)
	{};
	int getP() { return p;};
	void setP(int p) { this->p = p; };
	std::vector<std::complex<double> > getSR(std::complex<double> from,
			                                 std::complex<double> to,
			                                 const std::vector<std::complex<double> > &sCoeff);
	std::vector<std::complex<double> > getSS(std::complex<double> from,
			                                 std::complex<double> to,
			                                 const std::vector<std::complex<double> > &sCoeff);
	std::vector<std::complex<double> > getRR(std::complex<double> from,
			                                 std::complex<double> to,
			                                 const std::vector<std::complex<double> > &rCoeff);

	std::vector<std::complex<double> > getRCoeff(std::complex<double> xi, std::complex<double> xstar);
	std::vector<std::complex<double> > getSCoeff(std::complex<double> xi, std::complex<double> xstar);

    std::vector<std::complex<double> > getRVector(std::complex<double> y, std::complex<double> xstar);
    std::vector<std::complex<double> > getSVector(std::complex<double> y, std::complex<double> xstar);

    std::complex<double>               direct(std::complex<double> yj, std::complex<double> xi);

};
*/

std::vector<std::complex<double> > Potential::getSCoeff(std::complex<double> xi,
		                                                std::complex<double> xstar)
{
  std::vector<std::complex<double> > ans(p);
  ans[0].real(1.0); ans[0].imag(0.0);                    // first coefficient of S-expansions is 1
  for (int i=1; i<p; ++i)
  {
    ans[i] = -1.0 * pow((xi - xstar), i)/ ((double) i);  // type cast on i
  }
  return ans;
}



/**
 * Note that SS transformation/translation matrix we work with here is the
 * transpose of the matrix from the notes (also shown in notes - Main.cc).
 * We follow the work of Yang Wang's Master Thesis.  See the notes in Main.cc
 * for the S|S translation and its development
 */
std::vector<std::complex<double> > Potential::getSR(std::complex<double> from, std::complex<double> to,
		                                            const std::vector<std::complex<double> > &sCoeff)
{
  std::complex<double> t = to - from;
  // SR transformation matrix with size p x p
  std::vector<std::vector<std::complex<double> > > sr(p, std::vector<std::complex<double> >(p));

  // std::log(z) computes the complex natural logarithm of a complex value z
  // according to cppreference - the returned value is in the range of [-i \pi, i \pi]
  //
  // Recall:
  //   the notes in Main.cc file about complex logarithmic function
  //   Let u + iv = w = ln(z) where u,v are real numbers and we is complex
  //   By the definition of the complex natural logarithm e^{w} = z
  //   Since z is complex, let z = re^{i \theta} with r = ||z||
  //   Now, z = e^{w} implies re^{i \theta} = e^{u+iv} = e^u e^{iv}
  //   Therefore, r = e^u and i \theta = iv
  //   Hence, u + iv = w = ln(z) implies ln(r) + i \theta = ln(z)
  //   That is, ln(z) = ln(||z||) + i \theta
  //   Note that \theta (argument of z) is determined only up multiples of 2 \pi
  //   and is the angle made from the positive x-axis to the ray that extends
  //   from the pole (origin) to the point z = x + iy in the complex plane.
  //   The largest this rotation can be is 360 degrees
  //   \theta is defined to be somewhere in the interval [-\pi, \pi].
  //
  // Ex: z = 0 + 1i = x + iy
  //     r = ||z|| = 0^2 + 1^2 = 1
  //     tan( \theta ) = y/x = 1/0 (undefined - on positive y-axis) implies \theta = \pi / 2
  //
  //     ln(z) = ln(||z||) + i \theta
  //           = ln(1) + i \pi / 2
  //           = 0 + i \pi / 2
  //           = i \pi / 2
  //
  // Ex: z = -1.0 - 0.0i = x + iy
  //     r = ||z|| = (-1)^2 + (-0.0)^2 = 1
  //     tan( \theta ) = y/x = (-0.0)/(-1.0) implies \theta = - \pi
  //
  //     ln(z) = ln(||z||) + i \theta
  //           = ln(1) - i \pi
  //           = 0 - i \pi
  //           = - i \pi
  //
  // Ex: z = -1.0 + 0.0i = x + iy
  //     r = ||z|| = (-1.0)^2 + (0.0)^2 = 1
  //     tan( \theta ) = y/x = 0.0/(-1.0) implies \theta = \pi
  //
  //     ln(z) = ln(||z||) + i \theta
  //           = ln(1) + i \pi
  //           = 0 + i \pi
  //           = i \pi
  sr[0][0] = std::log(t);

  // first row and column of SR translation matrix
  // see Main.cc notes for matrix (transpose)
  sr[1][0] = 1.0/t;
  for (int i=2; i<p; ++i)  // first column (following Wang)
    sr[i][0] = -1.0 * sr[i-1][0] * double(i-1) / ( double(i) * t );
  for (int j=1; j<p; ++j)  // first row
    sr[0][j] = 1.0 / std::pow(t,j);

  // Working down from row to row following Yang Wang's thesis
  // Rows of Wang's thesis corresponds to columns of tranformation
  // matrix in Main.cc doxygen html output
  // Therefore, we are working across from column to column w.r.t. Main.cc illustration
  // We work this by row w.r.t. Wang's thesis (as in code below)
  //
  // Ex: Going from row 1 to 2 for count starting at 0
  //               (row 2 to row 3 for count starting at 1)
  //  i = 2 (row 2)
  //
  //    j = 1 (column 1 - count starts at 0)
  //        - multiply by a negative
  //        - divide by a power of t
  //
  //    j = 2 (column 2 - count starts at 0)
  //        - multiply by a negative
  //        - divide by a power of t
  //        - divide by 2 = i
  //        - multiply by 3 = i + j - 1
  //
  //    j = 3 (column 3)
  //        - multiply by a negative
  //        - divide by a power of t
  //        - divide by 2 = i
  //        - multiply by 4 = i + j - 1
  //
  //    j = 4 (column 4)
  //        - multiply by a negative
  //        - divide by a power of t
  //        - divide by 2 = i
  //        - multiply by 5 = i + j - 1
  //
  // Ex: Going from row 2 to 3 for count starting at 0
  //               (row 3 to row 4 for count starting at 1)
  //  i = 3 (row 3)
  //
  //    j = 1 (column 1 - count starts at 0)
  //        - multiply by a negative
  //        - divide by a power of t
  //
  //    j = 2 (column 2 - count starts at 0)
  //        - multiply by a negative
  //        - divide by a power of t
  //        - divide by 3 = i
  //        - multiply by 4 = i + j - 1
  //
  //    j = 3 (column 3)
  //        - multiply by a negative
  //        - divide by a power of t
  //        - divide by 3 = i
  //        - multiply by 5 = i + j - 1
  //
  //    j = 4 (column 4)
  //        - multiply by a negative
  //        - divide by a power of t
  //        - divide by 3 = i              // affectively multiplying
  //        - multiply by 6 = i + j - 1    // by 2 here
  //
  for (int i=1; i<p; ++i)
    for (int j=1; j<p; ++j)
      sr[i][j] = (-1.0) * sr[i-1][j] * double(i+j-1) / double(i);

  std::vector<std::complex<double> > ans(p);


  // Performing the translation:
  // Following the notes in Main.cc, the translation of the series to the new location (new series)
  // results in a change of the coefficients of the original series (and a change of center)
  // The new coefficients are obtained through a linear combination of the coefficients of
  // the series centered at the old location (old series).  In Main.cc, we can see that the
  // linear combination for each coefficient is a row/vector multiply where the row comes from
  // the matrix (transpose matrix in notes) that has just been created above.  The vector is
  // made up of the coefficients of the old series
  for (int i=0; i<p; ++i)
  {
    ans[i] = 0.0;
    for (int j=0; j<p; ++j)
      ans[i] += sr[i][j] * sCoeff[j];
  }

  return ans;

}


/**
 * Note that SS transformation/translation matrix we work with here is the
 * transpose of the matrix from the notes (also shown in notes - Main.cc).
 * We follow the work of Yang Wang's Master Thesis.  See the notes in Main.cc
 * for the S|S translation and its development
 */
std::vector<std::complex<double> > Potential::getSS(std::complex<double> from, std::complex<double> to,
		                                            const std::vector<std::complex<double> > &sCoeff)
{
  std::complex<double> t = to - from;
  // SS transformation matrix with size p x p
  std::vector<std::vector<std::complex<double> > > ss(p, std::vector<std::complex<double> >(p));

  // main diagonal of ss (SS transformation/translation matrix)
  for (int i=0; i<p; ++i)
      ss[i][i] = 1.0;

  // first colum of ss
  // The first column is defined recursively where multiplication by (double)(i-1)
  // cancels out the previous division
  ss[1][0] = t;                     // first column (0) second row (1)
  for (int i=2; i<p; ++i)           // first column rest of the rows
      ss[i][0] = ss[i-1][0] * (double)(i-1) * (-1.0) / (double)(i);

  // The lower triangular part of ss matrix working per row from left to right
  // note that the coefficient progression along a row are just the progression
  // of combinations for that rows number minus 2 (say for row p we are looking at
  // \f$ \begin{pmatrix} p-2 \\ k \end{pmatrix} \f$ where \f$k\f$ progresses from
  // 1 to p-2).
  // Example: Row 6
  //          i = 5
  //            j = i-1 = 4 (first term to left of diagonal of ones
  //              ss[i][j] = (-1.0)*ss[i][j+1] * (double)j / (double)(i-j)
  //              ss[5][4] = (-1.0)*1.0*4.0/(5.0-4.0) = -4.0
  //            j = i-2 = 3 (second term from left)
  //              ss[i][j] = (-1.0)*ss[i][j+1] * (double)j / (double)(i-j)
  //              ss[5][3] = (-1.0)*(-4.0)*(3.0)/(5.0-3.0) = 12.0/2.0 = 6.0
  //            j = i-3 = 5-3 = 2
  //              ss[i][j] = (-1.0)*ss[i][j+1] * (double)j / (double)(i-j)
  //              ss[5][2] = (-1.0)*(6.0)*(2.0)/(5.0-2.0) = -12.0/3.0 = -4.0
  //            j = i-4 = 5-4 = 1
  //              ss[i][j] = (-1.0)*ss[i][j+1] * (double)j / (double)(i-j)
  //              ss[5][1] = (-1.0)*(-4.0)*(1.0)/(5.0-1.0) = 4.0/4.0 = 1.0
  // More explanation can be found in Main.cc at S|S translation
  for (int i=1; i<p; ++i)           // all rows after first row (row 0)
      for (int j=i-1; j>=1; --j)    // lower triangular elements
          ss[i][j]= (-1.0) * ss[i][j+1] * ((double)j/(double)(i-j));

  // The entries in the upper triangular part of the ss matrix are all zero
  // (all entries to the right of the main diagonal of ones have value zero)
  for (int i=0; i<p; ++i)           // all rows
      for (int j=i+1; j<p; ++j)     // upper triangular part of p x p matrix
          ss[i][j]=0.0;

  // We now have the ss translation matrix

  // Performing the translation:
  // Following the notes in Main.cc, the translation of the series to the new location (new series)
  // results in a change of the coefficients of the original series (and a change of center)
  // The new coefficients are obtained through a linear combination of the coefficients of
  // the series centered at the old location (old series).  In Main.cc, we can see that the
  // linear combination for each coefficient is a row/vector multiply where the row comes from
  // the matrix (transpose matrix in notes) that has just been created above.  The vector is
  // made up of the coefficients of the old series
  std::vector<std::complex<double> > ans(p);
  for (int i=0; i<p; ++i)
  {
    ans[i] = 0.0;                                   // initializing ans
    for (int j=0; j<p; ++j)
    {
      ans[i] = ans[i] + (sCoeff[j] * ss[i][j]);     // row/vector multiply
    }
  }

  return ans;

}

std::vector<std::complex<double> > Potential::getRR(std::complex<double> from, std::complex<double> to,
		                                            const std::vector<std::complex<double> > &rCoeff)
{
  std::complex<double>  t = to - from;
  std::vector<std::vector<std::complex<double> > > rr(p, std::vector<std::complex<double> >(p));
  for (int i=0; i<p; ++i)
    rr[i][i] = 1.0;            // main diagonal
  for (int j=1; j<p; ++j)
    rr[0][j] = rr[0][j-1] * t; // first row

  /**
   *  In the for loop
   *    rr[i][j] = rr[i-1][j] * double(j-i+1) / (t * double(i))
   *  calls the previous row value above entry (rows correspond to columns in Main.cc write-up)
   *  The count starts on one since the first row is done (above)
   *
   *  Example:
   *    i = 1:
   *     j = 2:
   *       rr[i][j] = rr[i-1][j] * double(j-i+1) / (t * double(i))
   *                = t^2 * 2 / (t * 1) = 2 t
   *       Since
   *       - rr[i-1][j] = rr[1-1][2] = rr[0][2] = t^2
   *       - double(j-i+1) = double(2-1+1) = double(2) = 2
   *       - double(i) = 1
   *
   *    i = 1:
   *     j = 3:
   *       rr[i][j] = rr[i-1][j] * double(j-i+1) / (t * double(i))
   *                = t^3 * 3 / (t * 1) = 3 t^2
   *       Since
   *       - rr[i-1][j] = rr[1-1][3] = rr[0][3] = t^3
   *       - double(j-i+1) = double(3-1+1) = double(3) = 3
   *       - double(i) = 1
   *
   *    i = 2:
   *     j = 3:
   *       rr[i][j] = rr[i-1][j] * double(j-i+1) / (t * double(i))
   *                = 3t^2 * 2 / (t * 2) = 3 t
   *       Since
   *       - rr[i-1][j] = rr[2-1][3] = rr[1][3] = 3t^2
   *       - double(j-i+1) = double(3-2+1) = double(2) = 2
   *       - double(i) = 2
   *
   *    i = 2:
   *     j = 4:
   *       rr[i][j] = rr[i-1][j] * double(j-i+1) / (t * double(i))
   *                = 4t^3 * 3 / (t * 2) = 6 t^2
   *       Since
   *       - rr[i-1][j] = rr[2-1][4] = rr[1][4] = 4t^3
   *       - double(j-i+1) = double(4-2+1) = double(3) = 3
   *       - double(i) = 2
   *
   */
  for (int i=1; i<p; ++i)      // upper triangular part of matrix
    for (int j=i+1; j<p; j++)  // (columns after main diagonal)
      rr[i][j]
           = rr[i-1][j] * double(j-i+1) / (t * double(i)) ;

  for (int i=1; i<p; i++)      // lower triangular part of matrix
    for (int j=0; j<i; j++)    // (columns before main diagonal)
      rr[i][j] = 0.0;

  // Performing the translation:
  // Following the notes in Main.cc, the translation of the series to the new location (new series)
  // results in a change of the coefficients of the original series (and a change of center)
  // The new coefficients are obtained through a linear combination of the coefficients of
  // the series centered at the old location (old series).  In Main.cc, we can see that the
  // linear combination for each coefficient is a row/vector multiply where the row comes from
  // the matrix (transpose matrix in Main.cc notes) that has just been created above.  The vector
  // of the matrix-vector product is made up of the coefficients of the old series.  The
  // For each row i of the matrix, the row-vector products are given below.
  std::vector<std::complex<double> > ans(p);
  for (int i=0; i<p; i++)
  {
    ans[i] = 0.0;
    for (int j=0; j<p; j++)
    {
      ans[i] = ans[i] + ( rr[i][j] * rCoeff[j] );
    }
  }

  return ans;

}

// powers of the R-expansion power series (see Main.cc discussion for details).
std::vector<std::complex<double> > Potential::getRVector(std::complex<double> y, std::complex<double> xstar)
{
  std::vector<std::complex<double> > rVec(p);
  for (int i=0; i<p; ++i)
    rVec[i] = std::pow(y - xstar, i);
  return rVec;
}

// direct (exact) calculation of the potential acting on yj by xi
std::complex<double> Potential::direct(std::complex<double> yj, std::complex<double> xi)
{
  std::complex<double> ans = std::log(yj-xi);
  return ans;
}
