/*
 * Andersen_stochastic.cpp
 *
 *  Created on: Mar 30, 2016
 *      Author: maud
 */

#include "Andersen_stochastic.h"

const std::string AndersenStochastic::Name = "Andersen stochastic";
AndersenStochastic::AndersenStochastic(Box& box, double dt, double T, double nu):
		Thermostat(box, dt),
		TargetTemperature(T),
		Nu(nu),
		Sigma(sqrt(T)) {
		update_temp();
		dtime(dt);
		SimBox.update_VerletLists();
		SimBox.calculate_forces_verlet();

	}

void AndersenStochastic::update_temp() {}

void AndersenStochastic::dtime(double& dt) {
	Thermostat::dtime(dt);
	DeltaTHalf = DeltaT * 0.50;
}

void AndersenStochastic::propagate(bool calc_epot) {
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			if (mono.Mass == 1.0) mono.Velocity += mono.Force*DeltaTHalf;
			else mono.Velocity += mono.Force*(DeltaTHalf/mono.Mass);
			mono.Position += mono.Velocity*DeltaT;
		}
	}
	//SimBox.wrap();
 	SimBox.check_VerletLists();
 	//SimBox.calculate_forces(calc_epot);
	SimBox.calculate_forces_verlet(calc_epot);

	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			if (mono.Mass == 1.0) mono.Velocity += mono.Force*DeltaTHalf;
			else mono.Velocity += mono.Force*(DeltaTHalf/mono.Mass);
		}
	}
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			if (Rand::real_uniform() < Nu*DeltaT){
				for (unsigned i = 0; i < 3; i++) {
					if (mono.Mass == 1.0) mono.Velocity(i) = Rand::real_normal(0.0, Sigma);
					else mono.Velocity(i) = Rand::real_normal(0.0, Sigma/sqrt(mono.Mass));
				}
			}
		}
	}
}


std::string AndersenStochastic::name() const {return Name;}

std::string AndersenStochastic::info() const {
	std::string str{"Thermostat "};
	str += Name;
	str += " dtime ";
	str += std::to_string ( DeltaT );
	return str;
}





