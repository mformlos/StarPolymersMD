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
#include "Nose_Hoover.h"

int main() {

	Box box(30., 30., 30., 0.5, 0.9);
	box.add_chain(20, 1., 1.01);

	ofstream temp_file;
	temp_file.open("temperature.dat", ios::out | ios::trunc);

	box.print_molecules(std::cout);
	box.calculate_forces();
	box.print_Epot(std::cout);
	box.print_Ekin(std::cout);

	Thermostat *thermostat{};
	//thermostat = new Thermostat_None{ box, 0.001 };
	//thermostat = new Lowe_Andersen{ box, 0.001, 1., 20.0, 7.0 };
	thermostat = new Nose_Hoover{box, 0.001, 0.5, 1., 1.};

	for (int n = 0; n < 1e7; n++) {
		thermostat->propagate();
		if ( n > 1e5 && !(n%1000)) {
			std::cout << n << " ";
			box.print_Epot(std::cout);
			box.print_Ekin(std::cout);
			temp_file << n << " ";
			box.print_Epot(temp_file);
			box.print_Ekin(temp_file);
			box.print_Temperature(temp_file);
			temp_file << "\n";
			std::cout << '\n';
		}
	    // box.print_molecules(std::cout);
	}

	box.print_molecules(std::cout);
}
