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
#include "Functions.h"

//class Thermostat;

class Box {
	friend class Thermostat_None;
	friend class Lowe_Andersen;
	friend class Nose_Hoover;
	friend class Andersen;
	friend class MPC;
protected:
	double SystemTime;
	double Temperature;
	double Lambda;
	double Cutoff;
	double VerletRadius;
	double VerletRadius2;
	unsigned NumberOfMonomers;


	std::array<double,3> Size;
	std::array<int,3> CellSize;
	std::array<double,3> CellSideLength;
	std::vector<Molecule> Molecules;
	std::vector<MPCParticle> Fluid;
	std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>> CellList;
	//Thermostat Thermo;  //to be declared



public:
	Box(double Lx, double Ly, double Lz, double temperature, double Lambda);

	void add_chain(unsigned N, double mass, double bondLength);
	void add_chain(unsigned A, unsigned B, double mass, double bondLength);

	Vector3d& wrap (Vector3d& pos);
	Vector3d wrap (Vector3d&& pos);
	void wrap(Particle& part);
	void wrap(Molecule& mol);
	void wrap();

	Vector3d relative_position(Particle& one, Particle& two);

	void propagate(double dt);

	void calculate_forces(bool calc_epot = false);
	void calculate_forces_verlet(bool calc_epot = false);
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


	template<class UnitaryFunc>
	UnitaryFunc unitary(UnitaryFunc&& func) const;
};



#endif /* LIB_BOX_H_ */
