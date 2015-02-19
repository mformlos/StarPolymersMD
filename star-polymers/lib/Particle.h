/*
 * Particle.h
 *
 *  Created on: Dec 22, 2014
 *      Author: maud
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <forward_list>
#include "MatVec.h"

class Particle {

 public:
  //members
  MatVec Position;
  MatVec Velocity;
  MatVec Force;
  double Mass;
  bool AmphiType;
  bool Ghost;
  std::forward_list<Particle*> Neighbors;


  //Constructor
  Particle();
  Particle(double, bool, bool);
  Particle(MatVec, MatVec, double, bool, bool);

  void set_neighbor(Particle &neighbor);

  void print_neighbor() {
	  for (auto& m : Neighbors) {
		  std::cout << m->Position << std::endl;
	  }
  }

  double mass() const;

};
#endif /* PARTICLE_H_ */
