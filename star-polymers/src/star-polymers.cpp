//============================================================================
// Name        : star-polymers.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <forward_list>
#include "Box.h"
#include "Thermostat_None.h"
#include "Lowe_Andersen.h"

int main() {

	Box box(10., 10., 5., 1.);
	box.add_chain(5, 1., 1.01, 1.);

	ofstream temp_file;
	temp_file.open("temperature.dat", ios::out | ios::trunc);

	box.print_molecules(std::cout);
	box.calculate_forces();
	box.print_Epot(std::cout);
	box.print_Ekin(std::cout);

	Thermostat *thermostat{};
	//thermostat = new Thermostat_None{ box, 0.001 };
	thermostat = new Lowe_Andersen{ box, 0.001, 1., 500., 7.0 };

	for (int n = 0; n < 1000; n++) {
		thermostat->propagate();
		std::cout << n << " ";
		box.print_Epot(std::cout);
		box.print_Ekin(std::cout);
		box.print_Ekin(temp_file);
		temp_file << "\n";
		std::cout << '\n';
	    // box.print_molecules(std::cout);
	}

	box.print_molecules(std::cout);
}
