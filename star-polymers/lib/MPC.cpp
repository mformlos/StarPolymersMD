/*
 * MPC.cpp
 *
 *  Created on: Apr 24, 2015
 *      Author: maud
 */

#include "MPC.h"

MPC::MPC(Box& box, double aTemperature, double aShear) :
SimBox { box },
Temperature { aTemperature },
Shear { aShear } {
	NumberOfCells = SimBox.BoxSize[0]*SimBox.BoxSize[1]*SimBox.BoxSize[2] - 1;
	NumberOfMPCParticles = 10*(NumberOfCells+1);
	BoxSize[0] = SimBox.BoxSize[0];
	BoxSize[1] = SimBox.BoxSize[1];
	BoxSize[2] = SimBox.BoxSize[2];
	c = cos(130.0*M_PI/180.0);
	s = sin(130.0*M_PI/180.0);
	delrx = 0.;
	if (Shear != 0.) shear_on = true;
	else shear_on = false;
	Fluid.reserve(NumberOfMPCParticles);
	for ( unsigned i = 0 ; i < NumberOfMPCParticles ; i++ ) {
		Fluid.push_back(MPCParticle(1.0));
	}
}

void MPC::initializeMPC() {
	Vector3d CMV(0., 0., 0.);
	//double ekin { };
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
		//ekin += part.Velocity.squaredNorm();
	}
	//double vel_scale = sqrt(3.*(double)NumberOfMPCParticles*Temperature/(ekin*2.0));
	/*for (auto& part : Fluid) {
		part.Velocity *= vel_scale;
	}*/
}

//MPC routine:
void MPC::MPCstep(const double& dt) {
	delrx += Shear*BoxSize[1]*dt;
	delrx -= BoxSize[0]*floor(delrx/BoxSize[0]);
	streaming(dt);
	Vector3d Shift(Rand::real_uniform() - 0.5, Rand::real_uniform() - 0.5, Rand::real_uniform() - 0.5);
	shiftParticles(Shift);
	sort();
	//#pragma omp parallel num_threads(3)
	//{
    	//#pragma omp for
		for (unsigned Index = 0; Index <= NumberOfCells; ++Index) {
			Vector3d CMV { };
			calculateCMV(Index, CMV);
			thermostat(Index, CMV);
			collide(Index, CMV);
		}
	//}
	Shift = -Shift;
	shiftParticles(Shift);
}

void MPC::streaming(const double& dt) {
	for (auto& part : Fluid) {
		part.Position += part.Velocity*dt;
		if (shear_on) LEBC(part);
		else wrap(part);
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

void MPC::collide(unsigned Index, const Vector3d& CMV) {
	double phi { };
	double theta { };
	Vector3d RotationAxis { };
	Matrix3d RotationMatrix { };
	phi = 2.*M_PI*(Rand::real_uniform());
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
		if (part.CellIndex == Index) part.Velocity = CMV + RotationMatrix*(part.Velocity - CMV);
	}
	for (auto& mol : SimBox.Molecules ) {
		for (auto& mono : mol.Monomers) {
			if (mono.CellIndex == Index) mono.Velocity = CMV + RotationMatrix*(mono.Velocity - CMV);
		}
	}
}

void MPC::thermostat(unsigned index, const Vector3d& CMV) {
	double ekinOld { };
	double scaling { };
	unsigned count { };
	for (auto& part : Fluid) {
		if (part.CellIndex == index) {
			ekinOld += part.Mass*(part.Velocity - CMV).squaredNorm();
			count++;
		}
	}
	if (count < 2) return;
	ekinOld *= 0.5;
	scaling = pow((Rand::real_gamma(3.*(count-1)*Temperature/2.)/ekinOld), 0.5);
	for (auto& part : Fluid) {
		if (part.CellIndex == index) {
			part.Velocity = scaling*part.Velocity + (1-scaling)*CMV;
		}
	}
}

inline void MPC::shiftParticles(const Vector3d& Shift) {
	for (auto& part : Fluid) {
		part.Position += Shift;
		if(shear_on) LEBC(part);
		else wrap(part);
	}
	for (auto& mol : SimBox.Molecules ) {
		for (auto& mono : mol.Monomers) {
			mono.Position += Shift;
			if(shear_on) LEBC(mono);
			else wrap(mono);
		}
	}
}



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

double MPC::calculateEkinInCell(unsigned Index) {
	double ekin { };
	Vector3d CMV { };
	calculateCMV(Index, CMV);
	for (auto& part : Fluid) {
		if (part.CellIndex == Index) ekin += part.Mass*(part.Velocity - CMV).squaredNorm();
	}
	ekin *= 0.5;
	return ekin;
}

double MPC::calculateEkinTotal() {
	double ekin { };
	for (unsigned i = 0; i <= NumberOfCells; i++) {
		ekin += calculateEkinInCell(i);
	}
	return ekin;
}

double MPC::calculateCurrentTemperature() {
	double currenttemp { };
	currenttemp = 2.*calculateEkinTotal()/(3.*(NumberOfMPCParticles - filledCells()));
	return currenttemp;
}

unsigned MPC::filledCells() {
	unsigned count { };
	for (unsigned i = 0; i <= NumberOfCells; i++) {
		for (auto& part : Fluid) {
			if (part.CellIndex == i) {
				count++;
				break;
			}
		}
	}
	return count;
}

inline void MPC::LEBC(Particle &part) {
	double cy { floor(part.Position(1) / BoxSize[1]) };
	part.Position(0) -= BoxSize[0]*floor(part.Position(0)/BoxSize[0]);
	part.Position(0) -= cy*delrx;
	part.Position(0) -= BoxSize[0]*floor(part.Position(0)/BoxSize[0]);
	part.Position(1) -= BoxSize[1]*cy;
	part.Position(2) -= BoxSize[2]*floor(part.Position(2)/BoxSize[2]);
	part.Velocity(0) -= cy*Shear*BoxSize[1];
}

inline Vector3d& MPC::wrap(Vector3d& pos) {
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= floor(pos(i)/BoxSize[i]) * BoxSize[i];
	}
	return pos;
}

inline Vector3d MPC::wrap(Vector3d&& pos) {
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= floor(pos(i)/BoxSize[i]) * BoxSize[i];
	}
	return pos;
}

inline void MPC::wrap(Particle& part) {
	part.Position = wrap(part.Position);
}




