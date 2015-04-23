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
#include "Molecule.h"
#include "MatVec.h"
#include "Potentials.h"

//class Thermostat;

class Box {
	friend class Thermostat_None;
	friend class Lowe_Andersen;
	friend class Nose_Hoover;
	friend class Andersen;
protected:
	double SystemTime;
	double Temperature;
	double Lambda;
	double Cutoff;
	double VerletRadius;
	double VerletRadius2;
	unsigned NumberOfMonomers;



	std::array<double,3> Size;
	std::array<double,3> CellSize;
	std::vector<Molecule> Molecules;
	std::vector<Particle> Fluid;
	std::vector<std::vector<std::vector<std::forward_list<Particle*>>>> CellList;
	//Thermostat Thermo;  //to be declared



public:
	Box(double Lx, double Ly, double Lz, double temperature, double Lambda);

	void add_chain(unsigned N, double mass, double bondLength);
	void add_chain(unsigned A, unsigned B, double mass, double bondLength);

	MatVec& wrap (MatVec& pos);
	MatVec wrap (MatVec&& pos);
	void wrap(Particle& part);
	void wrap(Molecule& mol);
	void wrap();

	MatVec relative_position(Particle& one, Particle& two);

	void propagate(double dt);

	void calculate_forces(bool calc_epot = false);
	void calculate_forces_verlet();
	double calculate_ekin();
	double calculate_radius_of_gyration();
	unsigned numberOfMonomers();

	void update_VerletLists();
	void check_VerletLists();

	std::ostream& print_molecules(std::ostream& os) const;
	std::ostream& print_Epot(std::ostream& os) const;
	std::ostream& print_Ekin(std::ostream& os);
	std::ostream& print_Temperature(std::ostream& os);
	std::ostream& print_radius_of_gyration(std::ostream& os);
	// outputs
};



#endif /* LIB_BOX_H_ */
