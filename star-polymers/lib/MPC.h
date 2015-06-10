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
#include "Functions.h"
#include "Particle.h"
#include "Molecule.h"
#include "Rand.h"
#include "Box.h"
#include "Analysis.h"
#include "Hydrodynamics.h"

class MPC : public Hydrodynamics{
protected:
	double Temperature;
	double Shear;
	unsigned NumberOfCells;
	unsigned NumberOfMPCParticles;
	double c;
	double s;
	double delrx;
	bool shear_on;
	std::array<int,3> BoxSize;
	std::vector<std::vector<MPCParticle*>> MPCCellList;
	std::vector<unsigned> MPCCellListFluidParticles;
public:
	std::vector<MPCParticle> Fluid;
	MPC(Box&, double, double aShear = 0.);

	void initialize();

	//MPC routine:
	void step(const double& dt);
	void streaming(const double& dt);
	void sort();
	void sort(std::vector<std::vector<MPCParticle*>>&);
	void collide(unsigned, const Vector3d&);
	void collide(const std::vector<MPCParticle*>& List, const Vector3d& CMV);
	void thermostat(unsigned i, const Vector3d& CMV);
	void thermostat(const std::vector<MPCParticle*>&, const Vector3d&);
	inline void shiftParticles(const Vector3d& Shift);
	inline void LEBC(Particle&);
	inline Vector3d& wrap(Vector3d&);
	inline Vector3d wrap(Vector3d&&);
	inline void wrap(Particle&);





	//getter methods:
	inline void calculateCMV(unsigned Index, Vector3d&);
	inline void calculateCMV(std::vector<MPCParticle*>& List, Vector3d& CMV );
	inline void calculateCMP(unsigned Index, Vector3d&);
	double calculateEkinInCell(unsigned Index);
	double calculateEkinTotal();
	double calculateCurrentTemperature();
	unsigned filledCells();

	/*template<class UnitaryFunc>
	UnitaryFunc unitary(UnitaryFunc&& func) const;

	template<class UnaryFunc, class BinaryFunc>
	void operator() (UnaryFunc& ufunc, BinaryFunc& bfunc) const;*/

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

};



#endif /* LIB_MPC_H_ */
