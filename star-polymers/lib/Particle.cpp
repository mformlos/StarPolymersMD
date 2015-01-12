#include "../src/lib/Particle.h"

Particle::Particle() :
  Position { },
  Velocity { },
  Mass { 0. },
  AmphiType { },
  Ghost { } {}

Particle::Particle(double aMass, bool aAmphiType, bool aGhost) :
    Position { },
    Velocity { },
    Mass { aMass },
    AmphiType {aAmphiType },
    Ghost { aGhost } {}

Particle::Particle(MatVec aPosition, MatVec aVelocity, double aMass, bool aAmphiType, bool aGhost) :
  Position { aPosition },
  Velocity { aVelocity },
  Mass { aMass },
  AmphiType {aAmphiType },
  Ghost { aGhost } {}


double Particle::mass() const {return(Mass);}


