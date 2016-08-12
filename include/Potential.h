/*
 * Potential.h
 *
 *  Created on: Jul 13, 2016
 *      Author: dbpc
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <vector>
#include <complex>

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




#endif /* POTENTIAL_H_ */
