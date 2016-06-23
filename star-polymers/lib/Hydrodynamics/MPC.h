/*
 * MPC.h
 *
 *  Created on: Apr 24, 2015
 *      Author: maud
 */

#ifndef LIB_MPC_H_
#define LIB_MPC_H_

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include "Functions.h"
#include "Particle.h"
#include "Molecule.h"
#include "Rand.h"
#include "Analysis.h"
#include "Hydrodynamics.h"
#include <csignal>

class Box;
class MPC : public Hydrodynamics{
	friend class Box;
	friend class VelocityX;
protected:
	double Temperature;
	double Shear;
	unsigned step_update;
	unsigned NumberOfCells;
	unsigned NumberOfMPCParticles;
	double c;
	double s;
	double delrx;
	bool shear_on;
	bool angular_momentum;
	std::array<int,3> BoxSize;
	std::vector<std::vector<MPCParticle*>> MPCCellList;
	std::vector<unsigned> MPCCellListFluidParticles;
public:
	std::vector<MPCParticle> Fluid;
	MPC(Box&, double, int N_c, double aShear = 0., unsigned step = 1, bool angular_mom = false);

	void initialize();
	void initialize(string filename, std::array<int,3> Box_old);

	//MPC routine:
	void step(const long int& t, const double& dt);
	void streaming(const double& dt);
	void sort();
	void sort(std::vector<std::vector<MPCParticle*>>&);
	void collide(unsigned, const Vector3d&);
	void collide(const std::vector<MPCParticle*>& List, const Vector3d& CMV);
	void thermostat(unsigned i, const Vector3d& CMV);
	void thermostat(const std::vector<MPCParticle*>&, const Vector3d&);
	void periodic_image_box(unsigned);
	void undo_periodic_image_box(unsigned);
	void check_bounds();
	inline void shiftParticles(const Vector3d& Shift);
	void LEBC(Particle&);
	void LEBC(Vector3d&, Vector3d&);
	Vector3d& wrap(Vector3d&);
	Vector3d wrap(Vector3d&&);
	void wrap(Particle&);
	void wrap(Vector3d&, Vector3d&);

	void LEBC_to_zero(Vector3d&, Vector3d&);
	Vector3d& wrap_to_zero(Vector3d&);
	void wrap_to_zero(Vector3d&, Vector3d&);
	void wrap_to_zero(Particle&);

    Vector3d relative_position(Particle& one, Particle& two);
    Vector3d relative_position(const Vector3d, const Vector3d);




	//getter methods:
	inline void calculateCMV(unsigned Index, Vector3d&);
	inline void calculateCMV(std::vector<MPCParticle*>& List, Vector3d& CMV );
	inline void calculateFluidVelocity(unsigned Index, Vector3d&);
	inline void calculateCMP(unsigned Index, Vector3d&);
	inline void calculateAngular(unsigned Index, Vector3d& Angular, const Vector3d& CMV, const Vector3d& CMP, const Matrix3d& Rotation);
	inline void calculateAngularMomentum(unsigned Index, Vector3d&, const Vector3d& CMV, const Vector3d& CMP);
	double calculateEkinInCell(unsigned Index);
	double calculateEkinTotal();
	double calculateCurrentTemperature();
	unsigned filledCells();
	unsigned NumberOfParticles() const {return NumberOfMPCParticles;};

	template<class UnitaryFunc>
	UnitaryFunc unitary(UnitaryFunc& func) const {
		return for_each(Fluid.cbegin(), Fluid.cend(), func);
	}

	template<class UnaryFunc, class BinaryFunc>
	void operator() (UnaryFunc& ufunc, BinaryFunc& bfunc) const {
		auto first = Fluid.cbegin(), last = Fluid.cend();
		auto second = first;

		for(; first != last; ++first) {
			ufunc( *first );
			for( second = first + 1; second != last; ++second )
				bfunc( *first, *second );
		}
	}

	template<class UnaryFunc>
	void operator() (UnaryFunc& ufunc) const {
		auto first = Fluid.cbegin(), last = Fluid.cend();
		for(; first != last; ++first) {
			ufunc( *first );
		}
	}

	void print_fluid(FILE*, int, int, int);
	void print_fluid_with_coordinates(FILE*, int, int, int);
	void print_fluid_complete(FILE*);

};



#endif /* LIB_MPC_H_ */
