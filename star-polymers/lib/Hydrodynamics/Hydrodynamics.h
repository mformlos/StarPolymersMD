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
	virtual void step(const long int& t, const double& dt) {};
	virtual double calculateCurrentTemperature() {return 0.0;};
	virtual void print_fluid(FILE*, int, int, int) {};
	virtual void print_fluid_with_coordinates(FILE*, int, int, int) {};



};

#endif
