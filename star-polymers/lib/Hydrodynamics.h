#ifndef LIB_HYDRODYNAMICS_H_
#define LIB_HYDRODYNAMICS_H_

#include <string>
#include "Box.h"

class Hydrodynamics {
protected:
	Box& SimBox;
	~Hydrodynamics() = default;
public:
	Hydrodynamics() = default;
	Hydrodynamics(Box& box);
	virtual void initialize() {};
	virtual void step(const double& dt) {};
	virtual double calculateCurrentTemperature() {return 0.0;};

};

#endif
