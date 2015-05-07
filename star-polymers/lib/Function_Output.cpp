/*
 * Function_Output.cpp
 *
 *  Created on: May 6, 2015
 *      Author: maud
 */
#include "Function_Output.h"

void Function_Output::add(double x, double func_x) {
	double klass { floor( (x - offset) / width + center_left) * width + offset };
	function[klass] += func_x;
	function_iter = function.begin();
	count++;
	scale = 1. / (count);
}

Function_Output::Function_Output() :
		Function_Output { 1, 0.0, true } { }

Function_Output::Function_Output(double width, double offset, bool center) :
		function { }, function_iter { },
		width { width }, center_left{ center ? 0.5 : 0.0}, offset{offset},
		scale { }, count { } {}

Function_Output::Function_Output(double width, bool center, double offset) :
		Function_Output{width, offset, center} {}

void Function_Output::operator ()(double x, double func_x) {
	add(x, func_x);
}

bool Function_Output::output(std::ostream& os) {
	os.precision(8);
	os << std::scientific;
	os << function_iter->first << '\t';
	os << (function_iter->second) * scale << '\t';
	os << std::flush;
	auto function_iter_buf = function_iter;
	if (++function_iter_buf == function.end())
		return false;

	++function_iter;
	return true;
}
void Function_Output::output_reset() {
	function_iter = function.begin();
}


