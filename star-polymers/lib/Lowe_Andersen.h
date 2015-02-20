#ifndef LOWE_ANDERSEN_H
#define LOWE_ANDERSEN_H

#include "Thermostat.h"


#include "Thermostat.h"
#include "Box.h"
#include "Rand.h"
#include <cmath>
#include <string>

class Lowe_Andersen : public Thermostat {
private:
	double TargetTemperature;
	double DeltaTHalf;
	double Nu, NuDt, Sigma;
	static const std::string Name;
public:
	double time;

	Lowe_Andersen(Box& box, double dt, double temp, double nu);

	void dtime(double dt);
	void update_temp();
	void collide(Particle& one, Particle& two);
	void propagate();
	std::string name() const;
	std::string info() const;
};

#endif
