/*
 * Velocity_x.h
 *
 *  Created on: May 6, 2015
 *      Author: maud
 */

#ifndef LIB_VELOCITY_X_H_
#define LIB_VELOCITY_X_H_

#include "Analysis.h"
#include "Function_Output.h"
#include <../eigen/Eigen/Dense>

class VelocityX: public Analysis<Particle> {
protected:
	Function_Output vel_x_average;
	Vector3d x_hat;

public:
	VelocityX(double width, double offset = 0.0, bool center = true) :
		vel_x_average {width, offset, center},
		x_hat {1.0, 0.0, 0.0} {}

	VelocityX() :
		vel_x_average { },
		x_hat {1.0, 0.0, 0.0}{ }

	~VelocityX() = default;

	void operator() (const Particle& part) {
		vel_x_average(part.Position(1), part.Velocity.dot(x_hat));
	}
	double value() {return 1.0; }
	std::ostream& print_result(std::ostream& os) {
		bool out {};
	  	do {
		  	out = vel_x_average.output(os);
		  	os << '\n';
	  	} while(out);
	  	vel_x_average.output_reset();
		return os;
	}
};




#endif /* LIB_VELOCITY_X_H_ */
