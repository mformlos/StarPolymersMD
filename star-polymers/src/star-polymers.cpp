//============================================================================
// Name        : star-polymers.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <iostream>
#include "MatVec.h"

int main(void) {

	MatVec a{};
	MatVec b{};

	a[0] = 1;
	a[1] = 2;
	a[2] = 3;

	b[0] = b[1] = b[2] = 2;

	std::cout << "a : " << std::endl;
	std::cout << "b : " << std::endl;

	std::cout << "Skalarprodukt a*b = " << a*b << std::endl;


}
