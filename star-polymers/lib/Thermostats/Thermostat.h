#ifndef LIB_THERMOSTAT_H_
#define LIB_THERMOSTAT_H_

#include <string>
#include "Box.h"

class Thermostat {
protected:
	Box& SimBox;
	double DeltaT;
	~Thermostat() = default;
public:
	Thermostat() = default;
	Thermostat(Box& box, double dt);
	double dtime() const;
	virtual void dtime(double new_dtime);
	virtual void update_temp() = 0; // für mögliche Erwärmung
	virtual void propagate(bool calc_epot = false) = 0;
	virtual void propagate_gaussian(bool calc_epot = false) = 0;
	virtual std::string name() const = 0;
	virtual std::string info() const = 0;
};

#endif
