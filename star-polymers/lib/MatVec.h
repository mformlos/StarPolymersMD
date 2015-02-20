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
#include <limits>
#include <initializer_list>
#include <cmath>

class MatVec {
private:
	std::array<double, 3> Vec;

public:
	MatVec() = default;                      // standard Initialisierung
	MatVec(std::initializer_list<double> elements); //elementweise Initialisierung
	MatVec(const MatVec& other) = default;   // Copy Constructor
	MatVec(MatVec&& other) = default;        // Move Constructor
	~MatVec() = default;

	auto begin() -> decltype(Vec.begin()); // Iterators
	auto end() -> decltype(Vec.end());

	double& operator [](int i);    // Elementweiser Zugriff
	const double& operator [](int i) const;

	MatVec& operator =(const MatVec& other) = default;  // copy assignment
	MatVec& operator =(MatVec&& other) = default;       // move assignment

	friend void swap(MatVec& meca, MatVec& mecb) {
		std::swap(meca.Vec, mecb.Vec);
	}

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

	MatVec operator /(const std::array<double,3>& other) const; //elementweise division
	MatVec operator %(const std::array<double,3>& other) const; //elementweise multiplikation


	double operator *(const MatVec& other) const;   // Skalarprodukt

	double norm2();
	double norm();

	std::ostream& print(std::ostream& os) const;

	template<class UnaryFunction>
	UnaryFunction operator()(UnaryFunction func) {
		for (auto& el : Vec) func( el );
		return func;
	}

};

std::ostream& operator <<(std::ostream& os, const MatVec& some);

MatVec floor(MatVec mec);

MatVec round(MatVec mec);

double min(MatVec mec);

#endif /* SRC_LIB_MATVEC_H_ */
