/*
 * MatVec.h
 *
 *  Created on: Jan 12, 2015
 *      Author: maud
 */

#ifndef SRC_LIB_MATVEC_H_
#define SRC_LIB_MATVEC_H_

#include <array>
#include <ostream>
#include <initializer_list>
#include <cmath>

class MatVec {
private:
	std::array<double, 3> Vec;

public:
	MatVec();                      // standard Initialisierung
	MatVec(const MatVec& other);   // Copy Constructor
	MatVec(MatVec&& other);        // Move Constructor

	auto begin() -> decltype(Vec.begin()); // Iterators
	auto end() -> decltype(Vec.end());

	double& operator [](int i);    // Elementweiser Zugriff
	const double& operator [](int i) const;

	MatVec& operator =(const MatVec& other);  // copy assignment
	MatVec& operator =(MatVec&& other);       // move assignment

	bool operator ==(const MatVec& other) const;
	bool operator !=(const MatVec& other) const;

	MatVec& operator +=(const MatVec& other);
	MatVec& operator -=(const MatVec& other);
	MatVec& operator *=(const double& scale);
	MatVec& operator /=(const double& scale);

	MatVec operator +(const MatVec& other) const;   // Elementare Operationen
	MatVec operator -(const MatVec& other) const;
	MatVec operator *(const double& scale) const;
	MatVec operator /(const double& scale) const;
	MatVec operator -();                      // Negativer Vektor



	double operator *(const MatVec& other) const;   // Skalarprodukt

	double norm2();
	double norm();

	std::ostream& print(std::ostream& os) const;


};

std::ostream& operator << (std::ostream& os, const MatVec& some) {
	return some.print(os);
}


#endif /* SRC_LIB_MATVEC_H_ */
