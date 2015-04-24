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
#include <../eigen/Eigen/Dense>
using namespace Eigen;

class Particle {
public:
	Vector3d Position;
	Vector3d Velocity;
	double Mass;

	//Constructors
	//Particle() = default;
	Particle(double);
	Particle(Vector3d, Vector3d, double);
	//Particle(const Particle& other) = default;
	//Particle(Particle&& other) = default;
	~Particle() = default;

	//Particle& operator =(const Particle& mec) = default;
	//Particle& operator =(Particle && mec) = default;


};

class MPCParticle : public Particle {
public:
	unsigned CellIndex;

	//MPCParticle()=default;
	MPCParticle(double);
	MPCParticle(Vector3d, Vector3d, double);
	//MPCParticle& operator =(const MPCParticle& mec) = default;
	//MPCParticle& operator =(MPCParticle && mec) = default;
};

class MDParticle : public MPCParticle {

 public:
  //members

	Vector3d Force;
	Vector3d VerletPosition;
	bool AmphiType;
	std::forward_list<MDParticle*> Neighbors;
	std::forward_list<MDParticle*> VerletList;


	//Constructor

	MDParticle(double, bool);
	MDParticle(Vector3d, Vector3d, double, bool);
	//MDParticle(const MDParticle& other) = default;
	//MDParticle(MDParticle&& other) = default;
	//~MDParticle() = default;

	//MDParticle& operator =(const MDParticle& mec) = default;
	//MDParticle& operator =(MDParticle && mec) = default;

	void set_neighbor(MDParticle& neighbor);
	void clear_VerletList();

	void print_neighbor() {
		for (auto& m : Neighbors) {
			std::cout << m->Position.transpose() << std::endl;
		}
	}


};
#endif /* PARTICLE_H_ */
