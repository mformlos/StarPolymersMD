/*
 * Particle.h
 *
 *  Created on: Dec 22, 2014
 *      Author: maud
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <algorithm>
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
  Particle() = default;
  Particle(double, bool, bool);
  Particle(MatVec, MatVec, double, bool, bool);
  Particle(const Particle& other) = default;
  Particle(Particle&& other) = default;
  ~Particle() = default;

  Particle& operator =(const Particle& mec) = default;
  Particle& operator =(Particle && mec) = default;

  void set_neighbor(Particle& neighbor);

  void print_neighbor() {
	  for (auto& m : Neighbors) {
		  std::cout << m->Position << std::endl;
	  }
  }

  double mass() const;

};
#endif /* PARTICLE_H_ */
