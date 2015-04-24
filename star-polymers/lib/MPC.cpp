/*
 * MPC.cpp
 *
 *  Created on: Apr 24, 2015
 *      Author: maud
 */

#include "MPC.h"

MPC::MPC(Box& box, double aTemperature, int Lx, int Ly, int Lz) :
SimBox(box),
Temperature(aTemperature) {
	NumberOfCells = Lx*Ly*Lz - 1;
	NumberOfMPCParticles = 10*(NumberOfCells+1);
	BoxSize[0] = Lx;
	BoxSize[1] = Ly;
	BoxSize[2] = Lz;
	c = cos(130.0*M_PI/180.0);
	s = sin(130.0*M_PI/180.0);
	Fluid.reserve(NumberOfMPCParticles);
	for ( unsigned i = 0 ; i < NumberOfMPCParticles ; i++ ) {
		Fluid.push_back(MPCParticle(1.0));
	}
}

void MPC::initializeMPC() {
	Vector3d CMV(0., 0., 0.);
	for (auto& part : Fluid) {
		for (unsigned i = 0 ; i < 3; i++) {
			part.Position(i) = BoxSize[i]*Rand::real_uniform();
			part.Velocity(i) = Rand::real_uniform() - 0.5;
		}
		CMV += part.Velocity;
	}
	CMV /= NumberOfMPCParticles;
	for (auto& part : Fluid) {
		part.Velocity -= CMV;
	}
	//TODO: scale Velocities according to temperature
}

//MPC routine:

void MPC::streaming(double dt) {
	for (auto& part : Fluid) {
		part.Position += part.Velocity*dt;
	}
}

void MPC::sort() {
	for (auto& part : Fluid) {
		part.CellIndex = (int)part.Position(0) + BoxSize[0]*(int)part.Position(1)+BoxSize[0]*BoxSize[1]*(int)part.Position(2);
	}
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.CellIndex = (int)mono.Position(0) + BoxSize[0]*(int)mono.Position(1)+BoxSize[0]*BoxSize[1]*(int)mono.Position(2);
		}
	}
}

void MPC::collide() {
	Vector3d CMV;
	double phi { };
	double theta { };
	Vector3d RotationAxis;
	Matrix3d RotationMatrix;
	for (unsigned i = 0; i <= NumberOfCells; i++) {
		calculateCMV(i, CMV);
		phi = 2.*M_PI*(Rand::real_uniform()-0.5);
		theta = 2.*(Rand::real_uniform()-0.5);
		RotationAxis(0) = sqrt(1-theta*theta)*cos(phi);
		RotationAxis(1) = sqrt(1-theta*theta)*sin(phi);
		RotationAxis(2) = theta;
		RotationMatrix(0,0) = RotationAxis(0)*RotationAxis(0) + (1 - RotationAxis(0)*RotationAxis(0))*c;
		RotationMatrix(0,1) = RotationAxis(0)*RotationAxis(1)*(1 - c) - RotationAxis(2)*s;
		RotationMatrix(0,2) = RotationAxis(0)*RotationAxis(2)*(1 - c) + RotationAxis(1)*s;
		RotationMatrix(1,0) = RotationAxis(0)*RotationAxis(1)*(1 - c) + RotationAxis(2)*s;
		RotationMatrix(1,1) = RotationAxis(1)*RotationAxis(1) + (1 - RotationAxis(1)*RotationAxis(1))*c;
		RotationMatrix(1,2) = RotationAxis(1)*RotationAxis(2)*(1 - c) - RotationAxis(0)*s;
		RotationMatrix(2,0) = RotationAxis(0)*RotationAxis(2)*(1 - c) - RotationAxis(1)*s;
		RotationMatrix(2,1) = RotationAxis(1)*RotationAxis(2)*(1 - c) + RotationAxis(0)*s;
		RotationMatrix(2,2) = RotationAxis(2)*RotationAxis(2) + (1 - RotationAxis(2)*RotationAxis(2))*c;

		for (auto& part : Fluid ) {
			if (part.CellIndex == i) part.Velocity = CMV + RotationMatrix*(part.Velocity - CMV);
		}
		for (auto& mol : SimBox.Molecules ) {
			for (auto& mono : mol.Monomers) {
				if (mono.CellIndex == i) mono.Velocity = CMV + RotationMatrix*(mono.Velocity - CMV);
			}
		}
	}
}
	void thermostat();
	void MPCstep(double dt);


inline void MPC::calculateCMV(unsigned Index, Vector3d& CMV) {
	CMV(0) = CMV(1) = CMV(2) = 0.;
	double totalMass { };
	for (auto& part : Fluid) {
		if (part.CellIndex == Index) {
			CMV += part.Velocity*part.Mass;
			totalMass += part.Mass;
		}
	}
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			if (mono.CellIndex == Index) {
				CMV += mono.Velocity*mono.Mass;
				totalMass += mono.Mass;
			}
		}
	}
	CMV /= totalMass;
}

inline void MPC::calculateCMP(unsigned Index, Vector3d& CMP) {
	CMP(0) = CMP(1) = CMP(2) = 0.;
	double totalMass { };
	for (auto& part : Fluid) {
		if (part.CellIndex == Index) {
			CMP += part.Position*part.Mass;
			totalMass += 1.;
		}
	}
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			if (mono.CellIndex == Index) {
				CMP += mono.Position*mono.Mass;
				totalMass += mono.Mass;
			}
		}
	}
	CMP /= totalMass;
}

double MPC::calculateEkinInCell(unsigned Index) { }
double MPC::calculateEkinTotal() {}
double MPC::calculateCurrentTemperature() {}
int MPC::filledCells() {}


