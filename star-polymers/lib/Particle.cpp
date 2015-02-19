#include "Particle.h"

Particle::Particle() :
  Position { },
  Velocity { },
  Force { },
  Mass { 0. },
  AmphiType { },
  Ghost { },
  Neighbors { } {}

Particle::Particle(double aMass, bool aAmphiType, bool aGhost) :
    Position { },
    Velocity { },
	Force { },
    Mass { aMass },
    AmphiType {aAmphiType },
    Ghost { aGhost },
	Neighbors { } {}

Particle::Particle(MatVec aPosition, MatVec aVelocity, double aMass, bool aAmphiType, bool aGhost) :
  Position { aPosition },
  Velocity { aVelocity },
  Force { },
  Mass { aMass },
  AmphiType {aAmphiType },
  Ghost { aGhost },
  Neighbors { } {}

void Particle::set_neighbor(Particle &neighbor) {
	Neighbors.push_front(&neighbor);
}

double Particle::mass() const {return(Mass);}



