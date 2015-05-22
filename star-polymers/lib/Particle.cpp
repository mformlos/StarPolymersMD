#include "Particle.h"

Particle::Particle(double aMass) :
    Position { Vector3d::Zero() },
    Velocity { Vector3d::Zero() },
	Mass { aMass } { }

Particle::Particle(Vector3d aPosition, Vector3d aVelocity, double aMass) :
	Position { aPosition },
	Velocity { aVelocity },
	Mass { aMass } {}

MPCParticle::MPCParticle(double aMass) :
	Particle ( aMass ),
	CellIndex { } {}

MPCParticle::MPCParticle(Vector3d aPosition, Vector3d aVelocity, double aMass) :
	Particle(aPosition, aVelocity, aMass),
	CellIndex { } {}

MDParticle::MDParticle(double aMass, bool aAmphiType) :
	MPCParticle ( aMass ),
	Force { Vector3d::Zero() },
	VerletPosition { Vector3d::Zero() },
	AmphiType {aAmphiType },
	Anchor { },
	Neighbors { },
	VerletList { } {}

MDParticle::MDParticle(Vector3d aPosition, Vector3d aVelocity, double aMass, bool aAmphiType) :
	MPCParticle(aPosition, aVelocity, aMass),
	VerletPosition { Vector3d::Zero() },
    AmphiType {aAmphiType },
	Anchor { },
	Neighbors { },
	VerletList { } {}


void MDParticle::set_neighbor(MDParticle& neighbor) {
	Neighbors.push_front(&neighbor);
}


void MDParticle::clear_VerletList() {
	VerletList.clear();
}




