#include "../src/lib/Molecule.h"

Molecule::Molecule(unsigned N) :
	NumberOfMonomers { N },
	Ekin { }
	{
	for ( int i = 0 ; i < N ; i++ ) {
		Monomers.push_back(Particle());
	}
}

Molecule::Molecule(unsigned N, double m) :
		NumberOfMonomers { N },
		Ekin { }
		{
		for ( int i = 0 ; i < N ; i++ ) {
			Monomers.push_back(Particle(m, 0, 0));
		}
}

double Molecule::calculate_Ekin() {
	Ekin = 0.0;
	for (auto& m : Monomers) {
		Ekin += m.Mass*m.Velocity*m.Velocity;
	}
	Ekin /= 2.0;
	return Ekin;
}



