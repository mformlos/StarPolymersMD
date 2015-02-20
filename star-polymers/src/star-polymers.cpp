//============================================================================
// Name        : star-polymers.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <iostream>
#include <forward_list>
#include "Box.h"
#include "Thermostat_None.h"

int main() {

	Box box(10., 10., 5., 1.);
	box.add_chain(5, 1., 1.01);

	box.print_molecules(std::cout);
	box.calculate_forces();
	box.print_Epot(std::cout);

	Thermostat *thermostat{};
	thermostat = new Thermostat_None{ box, 0.01 };
	thermostat->propagate();

	box.print_molecules(std::cout);
}
