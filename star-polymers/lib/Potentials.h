/*
 * Potentials.h
 *
 *  Created on: Feb 19, 2015
 *      Author: maud
 */

#ifndef LIB_POTENTIALS_H_
#define LIB_POTENTIALS_H_

#include <cmath>
#include <iostream>

inline double Fene_Potential(double r2) {
	double potential{};
	try {
		if (r2 == 0.0) throw 1;
		potential = -15.0*2.25*log(1.0 - r2/2.25);
	}
	catch(...) { std::cout << "r = 0 !" << std::endl; }
	return potential;
}

inline double Fene_Force(double r2){
	double force{};
	try {
		if (r2 == 0.0) throw 1;
		force = -30.0/(1.0 - (r2/2.25));
	}
	catch(...) { std::cout << "r = 0 !" << std::endl; }
	return force;
}

inline double Fene_Anchor_Potential(double r2) {
	double potential{};
	try {
		if (r2 == 0.0) throw 1;
		potential = -15.0*20.25*log(1.0 - r2/20.25);
	}
	catch(...) { std::cout << "r = 0 !" << std::endl; }
	return potential;
}

inline double Fene_Anchor_Force(double r2) {
	double force{};
	try {
		if (r2 == 0.0) throw 1;
		force = -30.0/(1.0 - (r2/20.25));
	}
	catch(...) { std::cout << "r = 0 !" << std::endl; }
	return force;
}

inline double TypeAA_Potential(double r2) {
	double potential{};
	if (r2 < 1.0594631) {
		double rm24 { 1.0 / r2 };
		rm24 *= rm24*rm24;
		rm24 *= rm24;
		rm24 *= rm24;
		potential = 4.0*rm24*(rm24 - 1.0) + 1.0;
	}
	return potential;
}

inline double TypeAA_Force(double r2) {
	double force{};
	if (r2 < 1.0594631) {
		double rm2 { 1.0 / r2 };
		double rm24 { rm2*rm2*rm2 };
		rm24 *= rm24;
		rm24 *= rm24;
		force = 192.0*rm2*rm24*(rm24 - 0.5);
	}
	return force;
}

inline double TypeBB_Potential(double r2, double lambda) {
	double potential{};
	if (r2 < 1.0594631) {
		double rm24 { 1.0 / r2 };
		rm24 *= rm24*rm24;
		rm24 *= rm24;
		rm24 *= rm24;
		potential = 4.0*rm24*(rm24 - 1.0) + 1.0 - lambda;
	}
	else if (r2 < 2.25) {
		double rm24 { 1.0 / r2 };
		rm24 *= rm24*rm24;
		rm24 *= rm24;
		rm24 *= rm24;
		potential = 4.0*lambda*(rm24*(rm24 - 1.0) - 0.0000593);
	}
	return potential;
}

inline double TypeBB_Force(double r2, double lambda) {
	double force{};
	if (r2 < 1.0594631) {
		double rm2 { 1.0 / r2 };
		double rm24 { rm2*rm2*rm2 };
		rm24 *= rm24;
		rm24 *= rm24;
		force = 192.0*rm2*rm24*(rm24 - 0.5);
	}
	else if (r2 < 2.25) {
		double rm2 { 1.0 / r2 };
		double rm24 { rm2*rm2*rm2 };
		rm24 *= rm24;
		rm24 *= rm24;
		force = 192.0*lambda*rm2*rm24*(rm24 - 0.5);
	}
	return force;
}

#endif /* LIB_POTENTIALS_H_ */
