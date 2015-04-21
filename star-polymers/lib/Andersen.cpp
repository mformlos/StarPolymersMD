
#include "Andersen.h"

const std::string Andersen::Name = "Andersen";
Andersen::Andersen(Box& box, double dt, double T, unsigned step):
		Thermostat(box, dt),
		TargetTemperature(T),
		UpdateStep(step),
		Step { } {
		update_temp();
		dtime(dt);
	}

void Andersen::update_temp() { }

void Andersen::dtime(double dt) {
	Thermostat::dtime(dt);
	DeltaTHalf = DeltaT * 0.50;
}

void Andersen::propagate() {
	Step++;
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
			mono.Position += mono.Velocity*DeltaT;
		}
	}
	SimBox.wrap();
	SimBox.calculate_forces();
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
		}
	}
	if (!(Step%UpdateStep)) {
		for (auto& mol : SimBox.Molecules) {
			for (auto& mono : mol.Monomers) {
				for (unsigned i = 0; i < 3; i++) {
					mono.Velocity[i] = Rand::real_normal(0.0, sqrt(TargetTemperature/mono.Mass));
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


