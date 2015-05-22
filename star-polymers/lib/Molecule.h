/*
 * Molecule.h
 *
 *  Created on: Dec 22, 2014
 *      Author: maud
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <vector>

#include "Particle.h"
#include "Rand.h"
#include <../eigen/Eigen/Dense>



class Molecule {

protected:
	double Ekin;

public:
	double Epot;
	unsigned NumberOfMonomers;
	unsigned AType;
	unsigned BType;
	unsigned Arms;
	std::vector<MDParticle> Monomers;

	Molecule(unsigned);
	Molecule(unsigned, double);

	MDParticle& operator [](int i);    // Elementweiser Zugriff
	const MDParticle& operator [](int i) const;

	void initialize_straight_chain(unsigned A, unsigned B, double Temperature, double bond = 1.01);
	void initialize_open_star(unsigned A, unsigned B, unsigned Arms, double Temperature, double Bond = 1.01, double AnchorBond = 2.0);
	double calculate_Ekin();

	std::ostream& print(std::ostream& os) const;
	template<class UnitaryFunc>
	UnitaryFunc unitary(UnitaryFunc&& func) const;

	template<class UnaryFunc, class BinaryFunc>
	void operator() (UnaryFunc& ufunc, BinaryFunc& bfunc) const;


};

std::ostream& operator <<(std::ostream& os, const Molecule& some);



#endif /* MOLECULE_H_ */
