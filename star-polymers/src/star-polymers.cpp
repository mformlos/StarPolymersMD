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

int main(void) {

	Box box(10., 10., 5., 1.);
	box.add_chain(5, 1., 1.01);

	box.print_molecules(std::cout);
	std::cout << "git test";
	box.calculate_forces();
}
