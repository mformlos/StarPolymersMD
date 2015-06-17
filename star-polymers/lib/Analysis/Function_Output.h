/*
 * Histogramm.h
 *
 *  Created on: May 6, 2015
 *      Author: maud
 */

#ifndef LIB_FUNCTION_OUTPUT_H_
#define LIB_FUNCTION_OUTPUT_H_



#include <map>
#include <iterator>
#include <cmath>
#include <iostream>

//Histogramm symmetrisch um 0
//frei in Klasseneinteilung
class Function_Output {
protected:
	std::map<double, double> function;
	std::map<double, double>::iterator function_iter;

	double width, center_left, offset, scale;
	unsigned count;

	void add(double x, double func_x);
public:
	Function_Output();
	Function_Output(double width, double offset = 0.0, bool center = true);
	Function_Output(double width, bool center, double offset = 0.0);

	void operator ()(double x, double func_x);

	bool output(std::ostream& os);

	void output_reset();

	unsigned get_count() const {return count;}
};


#endif /* LIB_FUNCTION_OUTPUT_H_ */
