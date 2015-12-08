/*
 * MPC.cpp
 *
 *  Created on: Apr 24, 2015
 *      Author: maud
 */

#include "MPC.h"

MPC::MPC(Box& box, double aTemperature, int N_c, double aShear, unsigned step, bool angular_mom) :
Hydrodynamics { box },
Temperature { aTemperature },
Shear { aShear },
step_update {step},
angular_momentum { angular_mom }{
	NumberOfCells = SimBox.BoxSize[0]*SimBox.BoxSize[1]*SimBox.BoxSize[2] - 1;
	NumberOfMPCParticles = N_c*(NumberOfCells+1);
	BoxSize[0] = SimBox.BoxSize[0];
	BoxSize[1] = SimBox.BoxSize[1];
	BoxSize[2] = SimBox.BoxSize[2];
	MPCCellList = std::vector<std::vector<MPCParticle*>>(NumberOfCells+1, std::vector<MPCParticle*>());
	MPCCellListFluidParticles = std::vector<unsigned>(NumberOfCells+1, 0);
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

void MPC::initialize() {
	Vector3d CMV(0., 0., 0.);
	//double ekin { };
	for (auto& part : Fluid) {
		for (unsigned i = 0 ; i < 3; i++) {
			part.Position(i) = BoxSize[i]*(Rand::real_uniform() - 0.5);
			part.Velocity(i) = Rand::real_uniform() - 0.5;
		}
		wrap(part);
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

void MPC::initialize(string filename) {
	ifstream file {filename};
	string line { };
	Vector3d pos {Vector3d::Zero()};
	Vector3d vel {Vector3d::Zero()};
	Vector3d CMV(0., 0., 0.);
	int count {};
	for (auto& part : Fluid) {
		if (file >> pos(0) >> pos(1) >> pos(2) >> vel(0) >> vel(1) >> vel(2)) {
			if (fabs(round(pos(0)/BoxSize[0])) >= 1.0 || (fabs(round(pos(1)/BoxSize[1])) >= 1.0) || fabs(round(pos(2)/BoxSize[2])) >= 1.0) {
				for (unsigned i = 0 ; i < 3; i++) {
					part.Position(i) = BoxSize[i]*(Rand::real_uniform() - 0.5);
					part.Velocity(i) = Rand::real_uniform() - 0.5;
				}
				count ++;
			}
			else {
				part.Position = pos;
				part.Velocity = vel;
			}
		}
		else {
			for (unsigned i = 0 ; i < 3; i++) {
				part.Position(i) = BoxSize[i]*(Rand::real_uniform() - 0.5);
				part.Velocity(i) = Rand::real_uniform() - 0.5;
			}
		}
		wrap(part);
		CMV += part.Velocity;
	}
	std::cout << count << std::endl;
	CMV /= NumberOfMPCParticles;
	for (auto& part : Fluid) {
		part.Velocity -= CMV;
	}
}

//MPC routine:
void MPC::step(const long int& t, const double& dt) {
	delrx += Shear*BoxSize[1]*dt;
	delrx -= BoxSize[0]*round(delrx/BoxSize[0]);
	if (t%step_update) return;
	else {
		streaming(dt);
		Vector3d Shift(Rand::real_uniform() - 0.5, Rand::real_uniform() - 0.5, Rand::real_uniform() - 0.5);
		shiftParticles(Shift);
		sort();

		#pragma omp parallel
		{
			#pragma omp for
			for (unsigned Index = 0; Index <= NumberOfCells; ++Index) {
				if (MPCCellListFluidParticles[Index] < MPCCellList[Index].size()) periodic_image_box(Index);
				Vector3d CMV { };
				calculateCMV(Index, CMV);
				thermostat(Index, CMV);
				collide(Index, CMV);
				if (MPCCellListFluidParticles[Index] < MPCCellList[Index].size()) undo_periodic_image_box(Index);
			}
		}
		Shift = -Shift;
		shiftParticles(Shift);
	}
}

void MPC::periodic_image_box(unsigned Index) {
	std::vector<MPCParticle*>::iterator iter = MPCCellList[Index].begin()+MPCCellListFluidParticles[Index];
	int periodic_box = (int)round(((*iter) -> Position(1))/BoxSize[1]);
	for (iter = MPCCellList[Index].begin() ; iter < MPCCellList[Index].begin() + MPCCellListFluidParticles[Index]; iter++) {
		(*iter) -> Velocity(1) += Shear*BoxSize[1]*periodic_box;
	}
}

void MPC::undo_periodic_image_box(unsigned Index) {
	std::vector<MPCParticle*>::iterator iter = MPCCellList[Index].begin()+MPCCellListFluidParticles[Index];
	int periodic_box = (int)round(((*iter) -> Position(1))/BoxSize[1]);
	for (iter = MPCCellList[Index].begin() ; iter < MPCCellList[Index].begin() + MPCCellListFluidParticles[Index]; iter++) {
		(*iter) -> Velocity(1) -= Shear*BoxSize[1]*periodic_box;
	}
}


void MPC::streaming(const double& dt) {
	for (auto& part : Fluid) {
		part.Position += part.Velocity*dt;
		wrap(part);
	}
}



void MPC::sort() {
	for (auto& element : MPCCellList) {
		element.clear();
	}
	for (auto& number : MPCCellListFluidParticles) {
		number = 0;
	}
	for (auto& part : Fluid) {
		part.CellIndex = (int)(part.Position(0)+BoxSize[0]*0.5) + BoxSize[0]*(int)(part.Position(1)+BoxSize[1]*0.5)+BoxSize[0]*BoxSize[1]*(int)(part.Position(2)+BoxSize[2]*0.5);
		//std::cout << part.CellIndex << std::endl;
		try {
			if (part.CellIndex < 0 || part.CellIndex > NumberOfCells) throw 20;
			MPCCellList[part.CellIndex].push_back(&part);
			MPCCellListFluidParticles[part.CellIndex]++;
		}
		catch(...) {std::cout << "CellIndex of fluid out of range! " << part.CellIndex << " " << part.Position.transpose() << std::endl;
			exit(0); }
	}
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			Particle WrappedMonomer{mono.Position, mono.Velocity, mono.Mass};
			wrap(WrappedMonomer);
			mono.CellIndex = (int)(WrappedMonomer.Position(0) + BoxSize[0]*0.5) + BoxSize[0]*(int)(WrappedMonomer.Position(1)+BoxSize[1]*0.5) + BoxSize[0]*BoxSize[1]*(int)(WrappedMonomer.Position(2) + BoxSize[2]*0.5);
			//std::cout << mono.CellIndex << std::endl;
			try {
				if (mono.CellIndex < 0 || mono.CellIndex > NumberOfCells) throw 20;
				MPCCellList[mono.CellIndex].push_back(&mono);
			}
			catch(...) {std::cout << "CellIndex of monomer out of range! " << mono.CellIndex << std::endl;
				exit(0);}
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
	Vector3d CMP(0., 0., 0.);
	Vector3d Angular(0., 0., 0.);
	bool amc { angular_momentum && MPCCellList[Index].size() > 2 };
	if (amc) {
		calculateCMP(Index, CMP);
		calculateAngular(Index, Angular, CMV, CMP, RotationMatrix);
	}

	for (auto& part : MPCCellList[Index] ) {
		part -> Velocity = CMV + RotationMatrix*(part -> Velocity - CMV);
		if (amc) part -> Velocity += Angular.cross(part -> Position - CMP);
	}
}


void MPC::thermostat(unsigned Index, const Vector3d& CMV){
	double ekinOld { };
	double scaling { };
	unsigned count { };
	for (unsigned i = 0; i < MPCCellListFluidParticles[Index]; i++) {
		MPCParticle * part {MPCCellList[Index][i]};
		ekinOld += part -> Mass * (part -> Velocity -CMV).squaredNorm();
		count++;
	}

	if (count < 2) return;
	ekinOld *= 0.5;
	scaling = sqrt(Rand::real_gamma(3.*(count-1)*Temperature/2.)/ekinOld);
	for (unsigned i = 0; i < MPCCellListFluidParticles[Index]; i++) {
		MPCParticle * part {MPCCellList[Index][i]};
		part -> Velocity = scaling * part -> Velocity + (1-scaling)*CMV;
	}

}

inline void MPC::shiftParticles(const Vector3d& Shift) {
	for (auto& part : Fluid) {
		part.Position += Shift;
		wrap(part);
	}
	for (auto& mol : SimBox.Molecules ) {
		for (auto& mono : mol.Monomers) {
			mono.Position += Shift;
		}
	}
}


inline void MPC::calculateCMV(unsigned Index, Vector3d& CMV ){
	CMV(0) = CMV(1) = CMV(2) = 0.;
	double totalMass { };
	for (auto& part : MPCCellList[Index]) {
		CMV += part -> Velocity * part -> Mass;
		totalMass += part -> Mass;
	}
	if (totalMass > 0) CMV /= totalMass;
}

inline void MPC::calculateCMP(unsigned Index, Vector3d& CMP) {
	CMP(0) = CMP(1) = CMP(2) = 0.;
	double totalMass { };
	for (auto& part : MPCCellList[Index]) {
		CMP += part -> Position * part -> Mass;
		totalMass += part -> Mass;
	}
	if (totalMass > 0) CMP /= totalMass;
}

inline void MPC::calculateFluidVelocity(unsigned Index, Vector3d& FluidVelocity) {
	FluidVelocity(0) = FluidVelocity(1) = FluidVelocity(2) = 0.;
	unsigned i {};
	for (i = 0; i < MPCCellListFluidParticles[Index]; i++) {
		MPCParticle * part {MPCCellList[Index][i]};
	    FluidVelocity += part -> Velocity;
	}
	if (i > 0) FluidVelocity /= (double)i;
}

inline void MPC::calculateAngular(unsigned Index, Vector3d& Angular, const Vector3d& CMV, const Vector3d& CMP, const Matrix3d& Rotation) {
	Vector3d Cvec(0., 0., 0.);
	Matrix3d InertiaTensor = Matrix3d::Zero();
	for (auto& part : MPCCellList[Index]) {
		Cvec += (part -> Position - CMP).cross(Rotation*(part -> Velocity - CMV) - part -> Velocity);
		InertiaTensor(0,0) += pow(part -> Position(1) - CMP(1), 2) + pow(part -> Position(2) - CMP(2), 2);
		InertiaTensor(1,1) += pow(part -> Position(0)-CMP(0), 2) + pow(part -> Position(2)-CMP(2), 2);
		InertiaTensor(2,2) += pow(part -> Position(0)-CMP(0), 2) + pow(part -> Position(1)-CMP(1), 2);
		InertiaTensor(0,1) -= (part -> Position(0) - CMP(0))*(part -> Position(1)-CMP(1));
		InertiaTensor(0,2) -= (part -> Position(0) - CMP(0))*(part -> Position(2)-CMP(2));
		InertiaTensor(1,2) -= (part -> Position(1) - CMP(1))*(part -> Position(2)-CMP(2));
	}
	InertiaTensor(1,0) = InertiaTensor(0,1);
	InertiaTensor(2,0) = InertiaTensor(0,2);
	InertiaTensor(2,1) -= InertiaTensor(1,2);
	Angular = InertiaTensor.llt().solve(-Cvec);
}

inline void MPC::calculateAngularMomentum(unsigned Index, Vector3d& AngularMomentum, const Vector3d& CMP) {
	for (auto& part : MPCCellList[Index]) {
		AngularMomentum += (part -> Position - CMP).cross(part -> Velocity);
	}
}


double MPC::calculateEkinInCell(unsigned Index) {
	double ekin { };
	Vector3d CMV { };
	calculateCMV(Index, CMV);
	for (unsigned i = 0; i < MPCCellListFluidParticles[Index]; i++) {
		MPCParticle * part {MPCCellList[Index][i]};
		ekin += part -> Mass*(part -> Velocity - CMV).squaredNorm();
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
		if (MPCCellListFluidParticles[i] > 0) count++;
	}
	return count;
}


void MPC::LEBC(Particle &part) {
	double cy { round(part.Position(1) / BoxSize[1]) };
	part.Position(0) -= BoxSize[0]*round(part.Position(0)/BoxSize[0]);
	part.Position(0) -= cy*delrx;
	part.Position(0) -= BoxSize[0]*round(part.Position(0)/BoxSize[0]);
	part.Position(1) -= BoxSize[1]*cy;
	part.Position(2) -= BoxSize[2]*round(part.Position(2)/BoxSize[2]);
	part.Velocity(0) -= cy*Shear*BoxSize[1];
}

void MPC::LEBC(Vector3d& pos, Vector3d& vel) {
	double cy { round(pos(1) / BoxSize[1]) };
	pos(0) -= BoxSize[0]*round(pos(0)/BoxSize[0]);
	pos(0) -= cy*delrx;
	pos(0) -= BoxSize[0]*round(pos(0)/BoxSize[0]);
	pos(1) -= BoxSize[1]*cy;
	pos(2) -= BoxSize[2]*round(pos(2)/BoxSize[2]);
	vel(0) -= cy*Shear*BoxSize[1];
}

void MPC::LEBC_to_zero(Vector3d& pos, Vector3d& vel) {
	double cy { round(pos(1) / BoxSize[1]) };
		pos(0) -= BoxSize[0]*round(pos(0)/BoxSize[0]);
		pos(0) -= cy*delrx;
		pos(0) -= BoxSize[0]*round(pos(0)/BoxSize[0]);
		pos(1) -= BoxSize[1]*cy;
		pos(2) -= BoxSize[2]*round(pos(2)/BoxSize[2]);
		vel(0) -= cy*Shear*BoxSize[1];
}

Vector3d MPC::relative_position(Particle& one, Particle& two) {
	Vector3d rel_pos {two.Position - one.Position};
	double cy {round(rel_pos(1)/BoxSize[1])};
	rel_pos(0) -= BoxSize[0]*round(rel_pos(0)/BoxSize[0]);
	rel_pos(0) -= cy*delrx;
	rel_pos(0) -= BoxSize[0]*round(rel_pos(0)/BoxSize[0]);
	rel_pos(1) -= BoxSize[1]*cy;
	rel_pos(2) -= BoxSize[2]*round(rel_pos(2)/BoxSize[2]);
	return rel_pos;
}

Vector3d MPC::relative_position(const Vector3d pos_one, const Vector3d pos_two) {
	Vector3d rel_pos {pos_two - pos_one};
	double cy {round(rel_pos(1)/BoxSize[1])};
	rel_pos(0) -= BoxSize[0]*round(rel_pos(0)/BoxSize[0]);
	rel_pos(0) -= cy*delrx;
	rel_pos(0) -= BoxSize[0]*round(rel_pos(0)/BoxSize[0]);
	rel_pos(1) -= BoxSize[1]*cy;
	rel_pos(2) -= BoxSize[2]*round(rel_pos(2)/BoxSize[2]);
	return rel_pos;
}

Vector3d& MPC::wrap(Vector3d& pos) {
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= round(pos(i)/BoxSize[i]) * BoxSize[i];
	}
	return pos;
}

Vector3d MPC::wrap(Vector3d&& pos) {
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= round(pos(i)/BoxSize[i]) * BoxSize[i];
	}
	return pos;
}

void MPC::wrap(Particle& part) {
	if (shear_on) LEBC(part);
	else part.Position = wrap(part.Position);
}

void MPC::wrap(Vector3d& pos, Vector3d& vel) {
	if (shear_on) LEBC(pos, vel);
	else wrap(pos);
}

Vector3d& MPC::wrap_to_zero(Vector3d& pos){
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= round(pos(i)/BoxSize[i]) * BoxSize[i];
	}
	return pos;
}

void MPC::wrap_to_zero(Particle& part){
	if (shear_on) LEBC_to_zero(part.Position, part.Velocity);
	else wrap_to_zero(part.Position);
}

void MPC::wrap_to_zero(Vector3d& pos, Vector3d& vel){
	if (shear_on) LEBC_to_zero(pos, vel);
	else wrap_to_zero(pos);
}


/*void MPC::print_fluid(FILE* fluid_file, int step, int z_start, int z_stop) {
	fprintf(fluid_file, "TIME     %d \n", step);
	for (int z = z_start; z < z_stop; z++) {
		fprintf(fluid_file, "Z %d \n", z);
		for (int x = 0; x < BoxSize[0]; x++) {
			for (int y = 0; y < BoxSize[1]; x++) {
				int Index = x + BoxSize[0]*y+BoxSize[0]*BoxSize[1]*z;
				Vector3d CMV { };
				calculateFluidVelocity(Index, CMV);
				fprintf(fluid_file, "%f %f %f ", CMV(0), CMV(1), CMV(2));
			}
			fprintf(fluid_file, "\n");
		}
		fprintf(fluid_file, "\n");
	}
	fprintf(fluid_file, "\n");
}*/

void MPC::print_fluid(FILE* fluid_file, int step, int z_start, int z_stop) {
	fprintf(fluid_file, "TIME %d \n", step);
	for (int z = z_start; z <= z_stop; z++) {
		fprintf(fluid_file, "Z %d \n", z );
		for(int x = 0; x < BoxSize[0]; x++) fprintf(fluid_file, "%f ", x+0.5);
		fprintf(fluid_file, "\n");
		for (int y = 0; y < BoxSize[1]; y++) fprintf(fluid_file, "%f ", y+0.5);
		fprintf(fluid_file, "\n");
		for(int x = 0; x < BoxSize[0]; x++) {
			for (int y = 0; y < BoxSize[1]; y++) {
				int Index = x + BoxSize[0]*y+BoxSize[0]*BoxSize[1]*z;
				Vector3d CMV { };
				calculateFluidVelocity(Index, CMV);
				fprintf(fluid_file, "%f %f ", CMV(0), CMV(1));
			}
			fprintf(fluid_file, "\n");
		}
		fprintf(fluid_file, "\n");
	}
	fprintf(fluid_file, "\n");
}

void MPC::print_fluid_with_coordinates(FILE* fluid_file, int step, int z_start, int z_stop) {
	fprintf(fluid_file, "TIME     %d \n", step);
	for (int z = z_start; z < z_stop; z++) {
		for (int x = 0; x < BoxSize[0]; x++) {
			for (int y = 0; y < BoxSize[1]; y++) {
				int Index = x + BoxSize[0]*y+BoxSize[0]*BoxSize[1]*z;
				Vector3d CMV { };
				calculateFluidVelocity(Index, CMV);
				fprintf(fluid_file, "%f %f %f %f %f %f \n", x + 0.5, y + 0.5, z + 0.5, CMV(0), CMV(1), CMV(2));
			}
		}
	}
	fprintf(fluid_file, "\n");
}

void MPC::print_fluid_complete(FILE* fluid_file) {
    Vector3d pos {Vector3d::Zero()};
    Vector3d vel {Vector3d::Zero()};
    Vector3d shift_pos {SimBox.COM_Pos};
    /*for (unsigned i = 0; i < 3; i++) {
    	shift_pos(i) -= BoxSize[i];
    }*/
	for (auto& part : Fluid) {

		/*pos = part.Position - shift_pos;
		vel = part.Velocity - SimBox.COM_Vel;
		wrap(pos, vel);*/
		fprintf(fluid_file, "%f %f %f %f %f %f \n", pos(0), pos(1), pos(2), vel(0), vel(1), vel(2));
	}
}




