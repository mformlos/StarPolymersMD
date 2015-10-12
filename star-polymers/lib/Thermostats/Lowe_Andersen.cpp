#include "Lowe_Andersen.h"



using namespace std;

const string Lowe_Andersen::Name = "Lowe_Andersen";

Lowe_Andersen::Lowe_Andersen(Box &box, double dt, double temp, double nu, double r) :
		Thermostat(box, dt),
		TargetTemperature(temp),
		Nu(nu),
		InteractionRadius{r} {
	update_temp();
	dtime(dt);
	SimBox.calculate_forces();
}

void Lowe_Andersen::update_temp() {  // for all masses equal
}

void Lowe_Andersen::dtime(double dt) {
	Thermostat::dtime(dt);
	DeltaTHalf = DeltaT * 0.5;
	NuDt = Nu*DeltaT*exp(-DeltaT*Nu);
	std::cout << "NuDt = " << NuDt << std::endl;
}

void Lowe_Andersen::collide(Particle& one, Particle& two) {
	Vector3d unit_sep = SimBox.relative_position(one, two);
	double dist = unit_sep.norm();
	if (dist < InteractionRadius) {
		unit_sep /= dist;
		double reduced_mass { (one.Mass == two.Mass) ? 0.5*one.Mass : one.Mass*two.Mass/(one.Mass + two.Mass) };
		double sigma = sqrt(TargetTemperature/(reduced_mass));// TODO: why does this work?
		double therm_v { sigma*Rand::real_normal()};
		Vector3d velocity_diff = two.Velocity - one.Velocity;
		Vector3d dv = unit_sep*(therm_v + velocity_diff.dot(unit_sep));
		one.Velocity += dv*(reduced_mass/one.Mass);
		two.Velocity -= dv*(reduced_mass/two.Mass);
	}
}

void Lowe_Andersen::propagate(bool calc_epot) {
	// velocity verlet
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
			mono.Position += mono.Velocity*DeltaT;
		}
	}
	SimBox.wrap();
	SimBox.check_VerletLists();
    SimBox.calculate_forces_verlet(calc_epot);

    for (auto& mol: SimBox.Molecules) {
    	for (auto& first : mol.Monomers) {
    		for(auto& second: first.VerletList) {
    			if (NuDt*0.5 < Rand::real_uniform()) continue;
    			collide(first, *second);
    		}
    	}
    }


	/*for (unsigned i = 0; i < SimBox.Molecules.size(); i++) {
		for (unsigned j = 0; j < SimBox.Molecules[i].Monomers.size(); j++) {
			for (unsigned l = j+1; l < SimBox.Molecules[i].Monomers.size(); l++) {
				if (NuDt < Rand::real_uniform()) continue;
				//std::cout << j << " " << l << std::endl;
				collide(SimBox.Molecules[i].Monomers[j], SimBox.Molecules[i].Monomers[l]);
			}

			for (unsigned k = i+1; k < SimBox.Molecules.size(); k++) {
				for (unsigned l = 0; l < SimBox.Molecules[k].Monomers.size(); l++) {
					collide(SimBox.Molecules[i].Monomers[j], SimBox.Molecules[k].Monomers[l]);
				}
			}
		}
	}*/

	//SimBox.calculate_forces(calc_epot);

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
