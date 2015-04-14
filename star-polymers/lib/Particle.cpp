#include "Particle.h"

Particle::Particle(double aMass, bool aAmphiType, bool aGhost) :
    Position { },
    Velocity { },
	Force { },
    Mass { aMass },
    AmphiType {aAmphiType },
    Ghost { aGhost },
	Neighbors { },
	VerletList { } {}

Particle::Particle(MatVec aPosition, MatVec aVelocity, double aMass, bool aAmphiType, bool aGhost) :
  Position { aPosition },
  Velocity { aVelocity },
  Force { },
  Mass { aMass },
  AmphiType {aAmphiType },
  Ghost { aGhost },
  Neighbors { },
  VerletList { } {}


void Particle::set_neighbor(Particle& neighbor) {
	Neighbors.push_front(&neighbor);
}

double Particle::mass() const {return(Mass);}

void Particle::clear_VerletList() {
	for (auto pointer : VerletList) {
		delete(pointer);
	}
	VerletList.clear();
}



