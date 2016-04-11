#ifndef LIB_THERMOSTAT_NONE_H_
#define LIB_THERMOSTAT_NONE_H_

#include "Thermostat.h"
#include "Box.h"
#include <cmath>
#include <string>

class Thermostat_None : public Thermostat {
private:
	double DeltaTHalf;
	static const std::string Name;
public:
	Thermostat_None(Box& box, double dt, bool gaussian = false);
	void update_temp();
	void dtime(double delta_time);
	void propagate(bool calc_epot = false);
	void propagate_gaussian(bool calc_epot = false);
	std::string name() const;
	std::string info() const;
};

#endif
