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
#include "Rand.h"
#include "Box.h"

class MPC {
protected:
	Box& SimBox;
	double Temperature;
	unsigned NumberOfCells;
	unsigned NumberOfMPCParticles;
	double c;
	double s;
	std::array<unsigned,3> BoxSize;
	std::vector<MPCParticle> Fluid;

public:
	MPC(Box&, double, int, int, int);

	void initializeMPC();

	//MPC routine:
	void MPCstep(double dt);
	void streaming(double dt);
	void sort();
	void collide(unsigned, const Vector3d&);
	void thermostat(unsigned i, const Vector3d& CMV);
	void shiftParticles(Vector3d& Shift);

	//getter methods:
	inline void calculateCMV(unsigned Index, Vector3d&);
	inline void calculateCMP(unsigned Index, Vector3d&);
	double calculateEkinInCell(unsigned Index);
	double calculateEkinTotal();
	double calculateCurrentTemperature();
	unsigned filledCells();
};



#endif /* LIB_MPC_H_ */
