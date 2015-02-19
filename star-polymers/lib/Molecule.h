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



class Molecule {

protected:
	double Ekin;
	double Epot;

public:
	unsigned NumberOfMonomers;
	std::vector<Particle> Monomers;

	Molecule(unsigned);
	Molecule(unsigned, double);

	void initialize_straight_chain(double);
	double calculate_Ekin();


	std::ostream& print(std::ostream& os) const;


};

std::ostream& operator <<(std::ostream& os, const Molecule& some);



#endif /* MOLECULE_H_ */
