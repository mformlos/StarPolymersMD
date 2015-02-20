#include "Lowe_Andersen.h"



using namespace std;

const string Lowe_Andersen::Name = "Lowe_Andersen";

Lowe_Andersen::Lowe_Andersen(Box &box, double dt, double temp, double nu) :
		Thermostat(box, dt),
		TargetTemperature(temp),
		Nu(nu) {
	update_temp();
	dtime(dt);
}

void Lowe_Andersen::update_temp() {  // for all masses equal
	Sigma = sqrt(2*TargetTemperature / SimBox.Molecules[0].Monomers[0].Mass);
}

void Lowe_Andersen::dtime(double dt) {
	Thermostat::dtime(dt);
	DeltaTHalf = DeltaT * 0.5;
	NuDt = Nu*DeltaT*exp(-DeltaT*Nu);
}

void Lowe_Andersen::collide(Particle& one, Particle& two) {
	double reduced_mass = one.Mass*two.Mass/(one.Mass + two.Mass);
	double sigma = sqrt(TargetTemperature/reduced_mass);
	double therm_v = Rand::real_normal(0, sigma);
	MatVec unit_sep = SimBox.relative_position(one, two);
	unit_sep /= unit_sep.norm();
	MatVec velocity_diff = two.Velocity - one.Velocity;
	MatVec dv = unit_sep*(therm_v - velocity_diff*unit_sep);
	one.Velocity += dv*(reduced_mass/one.Mass);
	two.Velocity -= dv*(reduced_mass/two.Mass);
}

void Lowe_Andersen::propagate() {
	// velocity verlet
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
			mono.Position += mono.Velocity*DeltaT;
		}
	}

	for (unsigned i = 0; i < SimBox.Molecules.size(); i++) {
		for (unsigned j = 0; j < SimBox.Molecules[i].Monomers.size(); j++) {
			if (NuDt < Rand::real_uniform()) continue;
			for (unsigned l = j+1; j < SimBox.Molecules[i].Monomers.size(); l++) {
				collide(SimBox.Molecules[i].Monomers[j], SimBox.Molecules[i].Monomers[l]);
			}
			for (unsigned k = i+1; k < SimBox.Molecules.size(); k++) {
				for (unsigned l = 0; l < SimBox.Molecules[k].Monomers.size(); l++) {
					collide(SimBox.Molecules[i].Monomers[j], SimBox.Molecules[k].Monomers[l]);
				}
			}
		}
	}

	SimBox.calculate_forces();

	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
		}
	}
}

string Lowe_Andersen::name() const {return Name;}

string Lowe_Andersen::info() const {
	string str{"Thermostat "};
	str += Name;
	str += " dtime ";
	str += to_string ( DeltaT );
	str += " nu ";
	str += to_string ( Nu );
return str;
}
