
#include "Thermostat_None.h"

const std::string Thermostat_None::Name = "None";
Thermostat_None::Thermostat_None(Box& box, double dt, bool gaussian):
		Thermostat(box, dt) {
		update_temp();
		dtime(dt);
		SimBox.update_VerletLists();
		if (gaussian) SimBox.calculate_forces_gaussian();
		else SimBox.calculate_forces_verlet();
	}

void Thermostat_None::update_temp() { }

void Thermostat_None::dtime(double dt) {
	Thermostat::dtime(dt);
	DeltaTHalf = DeltaT * 0.50;
}

void Thermostat_None::propagate(bool calc_epot) {
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
			mono.Position += mono.Velocity*DeltaT;
		}
	}

	//SimBox.wrap();
 	SimBox.check_VerletLists();
	SimBox.calculate_forces_verlet(calc_epot);

	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
		}
	}
}

void Thermostat_None::propagate_gaussian(bool calc_epot) {
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
			mono.Position += mono.Velocity*DeltaT;
		}
	}

	//SimBox.wrap();
 	SimBox.check_VerletLists();
	SimBox.calculate_forces_gaussian(calc_epot);

	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
		}
	}
}


std::string Thermostat_None::name() const {return Name;}

std::string Thermostat_None::info() const {
	std::string str{"Thermostat "};
	str += Name;
	str += " dtime ";
	str += std::to_string ( DeltaT );
	return str;
}


