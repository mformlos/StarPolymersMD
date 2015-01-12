/*
 * Particle.h
 *
 *  Created on: Dec 22, 2014
 *      Author: maud
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>


using namespace std;

class Particle {

 public:
  //members
  MatVec Position;
  MatVec Velocity;
  double Mass;
  bool AmphiType;
  bool Ghost;

  //Constructor
  Particle();
  Particle(double, bool, bool);
  Particle(MatVec, MatVec, double, bool, bool);

  double mass() const;

};




#endif /* PARTICLE_H_ */
