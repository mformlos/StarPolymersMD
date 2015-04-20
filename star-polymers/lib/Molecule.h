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



class Molecule {

protected:
	double Ekin;

public:
	double Epot;
	unsigned NumberOfMonomers;
	std::vector<Particle> Monomers;

	Molecule(unsigned);
	Molecule(unsigned, double);

	Particle& operator [](int i);    // Elementweiser Zugriff
	const Particle& operator [](int i) const;

	void initialize_straight_chain(unsigned A, unsigned B, double bond, double temperature);
	double calculate_Ekin();

	std::ostream& print(std::ostream& os) const;


};

std::ostream& operator <<(std::ostream& os, const Molecule& some);



#endif /* MOLECULE_H_ */
