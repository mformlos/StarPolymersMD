/*
 * Hydrodynamics_None.h
 *
 *  Created on: Jun 9, 2015
 *      Author: maud
 */

#ifndef LIB_HYDRODYNAMICS_NONE_H_
#define LIB_HYDRODYNAMICS_NONE_H_

#include "Hydrodynamics.h"

class Hydrodynamics_None : public Hydrodynamics {
public:
	Hydrodynamics_None(Box&);
	void initialize();
	void step(const double& dt);
	double calculateCurrentTemperature();
};

Hydrodynamics_None::Hydrodynamics_None(Box& box) : Hydrodynamics(box) {}
void Hydrodynamics_None::initialize() {}
void Hydrodynamics_None::step(const double& dt) {}
double Hydrodynamics_None::calculateCurrentTemperature() {return 0.0;}





#endif /* LIB_HYDRODYNAMICS_NONE_H_ */
