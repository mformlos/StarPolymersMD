#include "Nose_Hoover.h"

const std::string Nose_Hoover::Name = "Nose_Hoover";

Nose_Hoover::Nose_Hoover(Box& box, double dt, double temp, double a_q1, double a_q2 )
	: Thermostat(box, dt),
	  TargetTemperature(temp)
	, q1(a_q1)
	, q2(a_q2)
	, xi1(0.0)
	, xi2(0.0)
	, nuxi1(0.0)
	, nuxi2(0.0)
	, g1(0.0)
	, g2(0.0) {
	DeltaTHalf = 0.5*dt;
	DeltaT4 = 0.5*DeltaTHalf;
	DeltaT8 = 0.5*DeltaT4;
	update_temp();
	dtime(dt);
	SimBox.update_VerletLists();
	SimBox.calculate_forces_verlet();
}

Nose_Hoover::~Nose_Hoover() {}

void Nose_Hoover::pos_vel() {
	for (auto& mol : SimBox.Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
			mono.Position += mono.Velocity*DeltaT;
		}
	}

	SimBox.wrap();
	SimBox.check_VerletLists();
	SimBox.calculate_forces_verlet();

	for (auto& mol : SimBox.Molecules) {
			for (auto& mono : mol.Monomers) {
				mono.Velocity += (mono.Force/mono.Mass)*DeltaTHalf;
		}
	}
}

void Nose_Hoover::chain() {
	double ekin = SimBox.calculate_ekin();
	g2 = (q1*nuxi1*nuxi1 - TargetTemperature) / q2;
	nuxi2 += g2*DeltaT4;				// L_G2
	nuxi1 *= exp(-nuxi2*DeltaT8);	// L_nuxi1
	g1 = (2.0*ekin - (3.*SimBox.NumberOfMonomers - 1)*TargetTemperature) / q1;
	nuxi1 += g1*DeltaT4;				// L_G1
	nuxi1 *= exp(-nuxi2*DeltaT8);	// L_nuxi1
	xi1 += nuxi1*DeltaTHalf;			// L_xi
	xi2 += nuxi2*DeltaTHalf;			// L_xi
	double s = exp(-nuxi1*DeltaTHalf);

	for (auto& mol : SimBox.Molecules) {
		for (auto& mono: mol.Monomers) {
			mono.Velocity *= s;// L_Cv
		}
	}

	nuxi1 *= exp(-nuxi2*DeltaT8);	// L_nuxi1
	ekin *= s*s;
	g1 = (2.0*ekin - (3.*SimBox.NumberOfMonomers - 1)*TargetTemperature) / q1;
	nuxi1 += g1*DeltaT4;				// L_G1
	nuxi1 *= exp(-nuxi2*DeltaT8);	// L_nuxi1
	g2 = (q1*nuxi1*nuxi1 - TargetTemperature) / q2;
	nuxi2 += g2*DeltaT4;				// L_G2
}

void Nose_Hoover::propagate() {
	chain();
	pos_vel();
	chain();
}

void Nose_Hoover::update_temp(){
}

std::string Nose_Hoover::name() const {return Name;}

std::string Nose_Hoover::info() const {
	std::string str{"Thermostat "};
	str += Name;
	str += " dtime ";
	str += std::to_string ( DeltaT );
	str += " q1 ";
	str += std::to_string ( q1 );
	str += " q2 ";
	str += std::to_string ( q2 );
return str;
}
