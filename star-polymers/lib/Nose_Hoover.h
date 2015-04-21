#ifndef NOSE_HOOVER_H
#define NOSE_HOOVER_H


#include "Thermostat.h"
#include "Box.h"
#include <string>
#include <cmath>

class Nose_Hoover :public Thermostat {

private:
	double TargetTemperature;

	double q1, q2, xi1, xi2, nuxi1, nuxi2, g1, g2, DeltaTHalf, DeltaT4, DeltaT8;
	void chain();
	void pos_vel(bool);
	static const std::string Name;

public:
	Nose_Hoover(Box& box, double dt, double temp, double a_q1, double a_q2);
	~Nose_Hoover();
	void propagate(bool calc_epot = false);
	void update_temp();
	std::string name() const;
	std::string info() const;
};


#endif
