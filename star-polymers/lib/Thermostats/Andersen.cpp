
#include "Andersen.h"

const std::string Andersen::Name = "Andersen";
Andersen::Andersen(Box& box, double dt, double T, unsigned step):
		Thermostat(box, dt),
		TargetTemperature(T),
		UpdateStep(step),
		Step { } {
		update_temp();
		dtime(dt);
		SimBox.update_VerletLists();
		SimBox.calculate_forces_verlet();

	}

void Andersen::update_temp() {}

void Andersen::dtime(double& dt) {
	Thermostat::dtime(dt);
	DeltaTHalf = DeltaT * 0.50;
}

void Andersen::propagate(bool calc_epot) {
	Step++;
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
	if (!(Step%UpdateStep)) {
		for (auto& mol : SimBox.Molecules) {
			for (auto& mono : mol.Monomers) {
				for (unsigned i = 0; i < 3; i++) {
					if (mono.Mass == 1.0) mono.Velocity(i) = Rand::real_normal(0.0, sqrt(TargetTemperature));
					else mono.Velocity(i) = Rand::real_normal(0.0, sqrt(TargetTemperature/mono.Mass));
				}
			}
		}
	}
}

void Andersen::propagate_gaussian(bool calc_epot) {
	Step++;
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
	SimBox.calculate_forces_gaussian(calc_epot);

	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			if (mono.Mass == 1.0) mono.Velocity += mono.Force*DeltaTHalf;
			else mono.Velocity += mono.Force*(DeltaTHalf/mono.Mass);
		}
	}
	if (!(Step%UpdateStep)) {
		for (auto& mol : SimBox.Molecules) {
			for (auto& mono : mol.Monomers) {
				for (unsigned i = 0; i < 3; i++) {
					if (mono.Mass == 1.0) mono.Velocity(i) = Rand::real_normal(0.0, sqrt(TargetTemperature));
					else mono.Velocity(i) = Rand::real_normal(0.0, sqrt(TargetTemperature/mono.Mass));
				}
			}
		}
	}
}

std::string Andersen::name() const {return Name;}

std::string Andersen::info() const {
	std::string str{"Thermostat "};
	str += Name;
	str += " dtime ";
	str += std::to_string ( DeltaT );
	return str;
}


