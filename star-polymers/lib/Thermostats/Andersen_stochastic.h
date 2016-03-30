/*
 * Andersen_stochastic.h
 *
 *  Created on: Mar 30, 2016
 *      Author: maud
 */

#ifndef STAR_POLYMERS_LIB_THERMOSTATS_ANDERSEN_STOCHASTIC_H_
#define STAR_POLYMERS_LIB_THERMOSTATS_ANDERSEN_STOCHASTIC_H_

#include "Thermostat.h"
#include "Box.h"
#include "Rand.h"
#include <cmath>
#include <string>

class AndersenStochastic : public Thermostat {
private:
	double TargetTemperature;
	double DeltaTHalf;
	double Nu;
	double Sigma;
	static const std::string Name;
public:
	AndersenStochastic(Box& box, double dt, double T, double Nu);
	void update_temp();
	void dtime(double& delta_time);
	void propagate(bool calc_epot = false);
	std::string name() const;
	std::string info() const;
};




#endif /* STAR_POLYMERS_LIB_THERMOSTATS_ANDERSEN_STOCHASTIC_H_ */
