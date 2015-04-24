#include "Molecule.h"

Molecule::Molecule(unsigned N) :
	Ekin { },
	Epot { },
	NumberOfMonomers { N }
	{
		Monomers.reserve(N);
	}

Molecule::Molecule(unsigned N, double mass) :
		Ekin { },
		Epot { },
		NumberOfMonomers { N }
		{
			Monomers.reserve(N);
			for ( unsigned i = 0 ; i < N ; i++ ) {
			    Monomers.push_back(MDParticle(mass, 0));
			}
		}

MDParticle& Molecule::operator [](int i) {return Monomers[i];}

const MDParticle& Molecule::operator [](int i) const {return Monomers[i];}

double Molecule::calculate_Ekin() {
	Ekin = 0.0;
	for (auto& m : Monomers) {
		Ekin += m.Mass*(m.Velocity.dot(m.Velocity));
	}
	Ekin /= 2.0;
	return Ekin;
}

void Molecule::initialize_straight_chain(unsigned A, unsigned B, double bondLength, double temperature) {
	Vector3d tempPos{Vector3d::Zero()};
	Vector3d average_vel{Vector3d::Zero()};
	for (unsigned i = A; i < (A+B); i++) Monomers[i].AmphiType = 1;
	for (unsigned i = 0; i < NumberOfMonomers; i++) {
		Monomers[i].Position = tempPos;
		for (unsigned j = 0; j < 3; j++) {
			Monomers[i].Velocity(j) = Rand::real_uniform() -0.5;
		}
		average_vel += Monomers[i].Velocity;
		if (i != (NumberOfMonomers - 1)) {
 			Monomers[i].set_neighbor(Monomers[i+1]);
		}
		tempPos(0) += bondLength;
	}
	average_vel /= NumberOfMonomers;
	for (auto& mono : Monomers) {
		mono.Velocity -= average_vel;
	}
	calculate_Ekin();
	std::cout<< "ekin before: " << Ekin;
	double vel_scale = sqrt(3.*(double)NumberOfMonomers*temperature/(Ekin*2.0));
	std::cout << " vel scale: " << vel_scale;
	for (auto& mono : Monomers) {
		mono.Velocity *= vel_scale;
	}
	std::cout << " ekin after: " << calculate_Ekin() << std::endl;
}

std::ostream& Molecule::print(std::ostream& os) const {
	for (unsigned i = 0; i < NumberOfMonomers; i++) {
		os << Monomers[i].Position.transpose() <<" ; " << Monomers[i].Velocity.transpose() << " ; " << Monomers[i].Force.transpose() << '\n';
	}
	return os;
}

std::ostream& operator <<(std::ostream& os, const Molecule& some) {
	return some.print(os);
}
