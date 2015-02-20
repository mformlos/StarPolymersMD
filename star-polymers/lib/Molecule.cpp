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

void Molecule::initialize_straight_chain(double bondLength) {
	MatVec tempPos{};
	for (unsigned i = 0; i < NumberOfMonomers; i++) {
		Monomers[i].Position = tempPos;
		if (i != 0) {
			Monomers[i].set_neighbor(Monomers[i-1]);
		}
		if (i != (NumberOfMonomers - 1)) {
 			Monomers[i].set_neighbor(Monomers[i+1]);
		}
		tempPos[0] += bondLength;
	}
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
