/*
 * Molecule.h
 *
 *  Created on: Dec 22, 2014
 *      Author: maud
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <vector>

#include "../src/lib/Particle.h"



class Molecule {

protected:
	double Ekin;

public:
	unsigned NumberOfMonomers;
	std::vector<Particle> Monomers;

	Molecule(unsigned);
	Molecule(unsigned, double);

	void initialize_straight_chain();
	double calculate_Ekin();

};




#endif /* MOLECULE_H_ */
