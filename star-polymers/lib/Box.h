/*
 * Box.h
 *
 *  Created on: Feb 19, 2015
 *      Author: maud
 */

#ifndef LIB_BOX_H_
#define LIB_BOX_H_

#include <iostream>
#include <array>
#include <vector>
#include <forward_list>
#include "Molecule.h"
#include "MatVec.h"
#include "Potentials.h"

//class Thermostat;

class Box {
	friend class Thermostat_None;
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

	MatVec& wrap (MatVec& pos);
	MatVec wrap (MatVec&& pos);
	void wrap(Particle& part);
	void wrap(Molecule& mol);
	void wrap();

	MatVec relative_position(Particle& one, Particle& two);

	void propagate(double dt);

	void calculate_forces();

	std::ostream& print_molecules(std::ostream& os) const;
	std::ostream& print_Epot(std::ostream& os) const;
	// outputs
};



#endif /* LIB_BOX_H_ */
