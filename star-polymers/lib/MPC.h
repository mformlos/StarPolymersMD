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
#include "Box.h"
#include "Functions.h"
#include "Particle.h"
#include "Rand.h"

class MPC {
public:
	Box& SimBox;
	double Temperature;
	unsigned NumberOfCells;
	unsigned NumberOfMPCParticles;
	double c;
	double s;

	std::array<unsigned,3> BoxSize;
	std::vector<MPCParticle> Fluid;


	MPC(Box&, double, int, int, int);

	void initializeMPC();

	//MPC routine:
	void streaming(double dt);
	void sort();
	void collide();
	void thermostat();
	void MPCstep(double dt);

	//getter methods:
	inline void calculateCMV(unsigned Index, Vector3d&);
	inline void calculateCMP(unsigned Index, Vector3d&);
	double calculateEkinInCell(unsigned Index);
	double calculateEkinTotal();
	double calculateCurrentTemperature();
	int filledCells();
};



#endif /* LIB_MPC_H_ */
