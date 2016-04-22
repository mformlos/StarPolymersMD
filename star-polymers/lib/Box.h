/*
 * Box.h
 *
 *  Created on: Feb 19, 2015
 *      Author: maud
 */

#ifndef LIB_BOX_H_
#define LIB_BOX_H_

#include <cstdio>
#include <iostream>
#include <array>
#include <vector>
#include <list>
#include <tuple>
#include "Molecule.h"
#include "MatVec.h"
#include "Potentials.h"
#include "Functions.h"
#include <csignal>
#include "Analysis.h"
#include "Vel_Autocorr.h"


//class Thermostat;

class MPC;

class Box {
	friend class Thermostat_None;
	friend class Lowe_Andersen;
	friend class Nose_Hoover;
	friend class Andersen;
	friend class MPC;
	friend class VelocityX;
//protected:
public:
	double SystemTime;
	double Temperature;
	double Lambda;
	double Cutoff;
	double VerletRadius;
	double VerletRadius2;
	unsigned NumberOfMonomers;

	std::array<int,3> BoxSize;
	std::array<int,3> CellSize;
	std::array<double,3> CellSideLength;
	std::vector<Molecule> Molecules;
	std::vector<Molecule> Molecules_com_reference_frame;
	Vector3d COM_Pos;
	Vector3d COM_Vel;
	std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>> CellList;
	//Thermostat Thermo;  //to be declared



public:
	Box(int Lx, int Ly, int Lz, double temperature, double Lambda);


	void add_chain(unsigned A, unsigned B, double Mass, double Bond = 1.01);
	void add_star(unsigned A, unsigned B, unsigned Arms, double Mass, double Bond = 1.01, double AnchorBond = 2.0);
	void add_star(string filename, unsigned A, unsigned B, unsigned Arms, double Mass, bool set_zero = false);
	void add_gaussian(unsigned N, double Mass);
	void add_gaussian(string filename, unsigned N, double Mass, bool set_zero = false);

	void resize(double Lx, double Ly, double Lz);
	void set_center_of_mass_to_zero(Molecule&);
	void center_of_mass_reference_frame();


	void update_VerletLists();
	void check_VerletLists();

	void calculate_forces(bool calc_epot = false);
	void calculate_forces_verlet(bool calc_epot = false);
	void calculate_forces_gaussian(bool calc_epot = false);


	double calculate_ekin();
	double calculate_epot(MDParticle&, MDParticle&);
	double calculate_radius_of_gyration();
	std::tuple<double, Matrix3d> calculate_gyration_tensor();
	Vector3d calculate_rotation_frequency();
	std::list<unsigned> calculate_patches();
	std::vector<double> calculate_patches_new();


	unsigned numberOfMonomers();

	std::ostream& print_molecules(std::ostream& os) const;
	void print_PDB(FILE*, int step);
	void print_PDB_with_velocity(FILE*, int step);
	std::ostream& print_Epot(std::ostream& os) const;
	std::ostream& print_Ekin(std::ostream& os);
	std::ostream& print_Temperature(std::ostream& os);
	std::ostream& print_Ekin_and_Temperature(std::ostream& os);
	std::ostream& print_radius_of_gyration(std::ostream& os);
	std::ostream& print_center_of_mass(std::ostream& os);

	// outputs


	template<class UnitaryFunc>
	UnitaryFunc unitary(UnitaryFunc&& func) const;

	template<class UnaryFunc, class BinaryFunc>
	void operator() (UnaryFunc& ufunc, BinaryFunc& bfunc) const;


	template<class UnaryFunc>
	void operator() (UnaryFunc& ufunc) const {
		auto first = Molecules.cbegin(), last = Molecules.cend();
		for(; first != last; ++first) {
			ufunc(*first);
		}
	}

};



#endif /* LIB_BOX_H_ */
