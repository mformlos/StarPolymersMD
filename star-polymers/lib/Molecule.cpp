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
				//Particle newp(mass, 0, 0);
				//Monomers[i] = newp;
			    Monomers.push_back(Particle(mass, 0, 0));
			}
		}

Particle& Molecule::operator [](int i) {return Monomers[i];}

const Particle& Molecule::operator [](int i) const {return Monomers[i];}

double Molecule::calculate_Ekin() {
	Ekin = 0.0;
	for (auto& m : Monomers) {
		Ekin += m.Mass*(m.Velocity*m.Velocity);
	}
	Ekin /= 2.0;
	return Ekin;
}

void Molecule::initialize_straight_chain(unsigned A, unsigned B, double bondLength, double temperature) {
	MatVec tempPos{};
	MatVec average_vel{};
	for (unsigned i = A; i < (A+B); i++) Monomers[i].AmphiType = 1;
	for (unsigned i = 0; i < NumberOfMonomers; i++) {
		Monomers[i].Position = tempPos;
		Monomers[i].Velocity = MatVec{ Rand::real_uniform() -0.5, Rand::real_uniform() -0.5, Rand::real_uniform() -0.5};
		average_vel += Monomers[i].Velocity;
		if (i != 0) {
			Monomers[i].set_neighbor(Monomers[i-1]);
		}
		if (i != (NumberOfMonomers - 1)) {
 			Monomers[i].set_neighbor(Monomers[i+1]);
		}
		tempPos[0] += bondLength;
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
		os << Monomers[i].Position <<" ; " << Monomers[i].Velocity << " ; " << Monomers[i].Force << '\n';
	}
	return os;
}

std::ostream& operator <<(std::ostream& os, const Molecule& some) {
	return some.print(os);
}
