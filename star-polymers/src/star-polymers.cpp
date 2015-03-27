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
#include "Lowe_Andersen.h"

int main() {

	Box box(10., 10., 5., 1.);
	box.add_chain(5, 1., 1.01, 1.);

	box.print_molecules(std::cout);
	box.calculate_forces();
	box.print_Epot(std::cout);
	box.print_Ekin(std::cout);

	Thermostat *thermostat{};
	thermostat = new Thermostat_None{ box, 0.001 };
	//thermostat = new Lowe_Andersen{ box, 0.01, 1., 1./(0.01*5), 5 };

	for (int n = 0; n < 50; n++) {
		thermostat->propagate();
		box.print_Epot(std::cout);
		box.print_Ekin(std::cout);
		std::cout << '\n';
	    box.print_molecules(std::cout);
	}

	box.print_molecules(std::cout);
}
