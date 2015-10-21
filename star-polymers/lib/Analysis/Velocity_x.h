/*
 * Velocity_x.h
 *
 *  Created on: May 6, 2015
 *      Author: maud
 */

#ifndef LIB_VELOCITY_X_H_
#define LIB_VELOCITY_X_H_

#include "Analysis.h"
#include <map>
#include <iterator>
#include "Function_Output.h"
#include <../eigen/Eigen/Dense>

class VelocityX: public Analysis<Particle> {
protected:
	std::map<double, double> vel_x_average;
	std::map<double, double>::iterator vel_x_average_iter;
	std::map<double, double> vel_x_average_count;
	std::map<double, double>::iterator vel_x_average_count_iter;

	double width;
	//Function_Output vel_x_average;

public:
	VelocityX(double a_width = 0.5) :
		vel_x_average { },
		vel_x_average_iter { },
		vel_x_average_count { },
		vel_x_average_count_iter { },
		width {a_width} {}

	void initialize(string filename) {
		ifstream input {filename};
		double pos{ }, vel { };
		int count { };
	    while(input >> pos >> vel >> count) {
	    	double y {floor(pos/width)*width};
	    	vel_x_average[y] += vel*count;
	    	vel_x_average_count[y] += count;
	    }
    	vel_x_average_iter = vel_x_average.begin();
		vel_x_average_count_iter = vel_x_average_count.begin();
	}


	void operator() (const Particle& part) {
		double y {floor(part.Position(1)/ width) * width };
		vel_x_average[y] += part.Velocity(0);
		vel_x_average_iter = vel_x_average.begin();
		vel_x_average_count[y]++;
		vel_x_average_count_iter = vel_x_average_count.begin();
	}

	double value() {return 1.0; }
	std::ostream& print_result(std::ostream& os) {
		bool out {true};
		do {
			os.precision(8);
			os << std::scientific;
			os << vel_x_average_iter->first << ' ';
			os << (vel_x_average_iter->second) / (vel_x_average_count_iter -> second) <<  " " << vel_x_average_count_iter -> second << '\n';
			os << std::flush;
			auto vel_x_average_iter_buf = vel_x_average_iter;
			if (++vel_x_average_iter_buf == vel_x_average.end()) out = false;
			++vel_x_average_iter;
			++vel_x_average_count_iter;

		} while(out);
		vel_x_average_iter = vel_x_average.begin();
		vel_x_average_count_iter = vel_x_average_count.begin();
		return os;
	}



};




#endif /* LIB_VELOCITY_X_H_ */
