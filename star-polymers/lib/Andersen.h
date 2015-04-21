#ifndef LIB_ANDERSEN_H_
#define LIB_ANDERSEN_H_

#include "Thermostat.h"
#include "Box.h"
#include "Rand.h"
#include <cmath>
#include <string>

class Andersen : public Thermostat {
private:
	double TargetTemperature;
	double DeltaTHalf;
	unsigned UpdateStep;
	unsigned Step;
	static const std::string Name;
public:
	Andersen(Box& box, double dt, double T, unsigned step);
	void update_temp();
	void dtime(double delta_time);
	void propagate(bool calc_epot = false);
	std::string name() const;
	std::string info() const;
};

#endif
