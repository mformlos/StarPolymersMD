/*
 * Box.h
 *
 *  Created on: Feb 19, 2015
 *      Author: maud
 */

#ifndef LIB_BOX_H_
#define LIB_BOX_H_

#include <array>
#include <vector>
#include <forward_list>
#include "Molecule.h"
#include "MatVec.h"

class Box {
protected:
	double SystemTime;
	double Temperature;
	std::array<double,3> Size;
	std::vector<Molecule> Molecules;
	std::vector<Particle> Fluid;
	std::vector<std::vector<std::vector<std::forward_list<Particle*>>>> CellList;
	std::forward_list<Particle*> VerletList;
	//Thermostat Thermo;  //to be declared



public:
	Box(double Lx, double Ly, double Lz, double temperature);

	void add_chain(unsigned N, double mass, double bondLength);


	void propagate(double dt);

	void calculate_forces();

	std::ostream& print_molecules(std::ostream& os) const;

	// outputs
};



#endif /* LIB_BOX_H_ */
