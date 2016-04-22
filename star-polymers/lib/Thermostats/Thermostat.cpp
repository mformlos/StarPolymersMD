#include "Thermostat.h"

Thermostat::Thermostat(Box& box, double dt) :
		SimBox ( box ),
		DeltaT { dt } { }

double Thermostat::dtime() const {
	return DeltaT;
}

void Thermostat::dtime(double& new_dtime) {
	DeltaT = new_dtime;
}

