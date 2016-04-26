/*
 * Vel_distribution.h
 *
 *  Created on: Apr 26, 2016
 *      Author: maud
 */

#ifndef STAR_POLYMERS_LIB_ANALYSIS_VEL_DISTRIBUTION_H_
#define STAR_POLYMERS_LIB_ANALYSIS_VEL_DISTRIBUTION_H_

#include "Analysis.h"
#include <map>
#include <iterator>
#include <../eigen/Eigen/Dense>


class VelocityDistribution: public Analysis<Particle> {
protected:

	std::map<double, double> vel_dist;
	std::map<double, double>::iterator vel_dist_iter;
	unsigned count;
	double width;

public:
	VelocityDistribution(double a_width = 0.5) :
		vel_dist { },
		vel_dist_iter { },
		count { },
		width {a_width} {}



	void operator() (const Particle& part) {
		double x {floor(part.Velocity.norm()/ width) * width };
		vel_dist[x] += 1.;
		vel_dist_iter = vel_dist.begin();
		count++;
	}

	double value() {return 1.0; }
	std::ostream& print_result(std::ostream& os) {
		bool out {true};
		do {
			os.precision(8);
			os << std::scientific;
			os << vel_dist_iter->first << ' ';
			os << (vel_dist_iter->second) / ((double)count*width) << '\n';
			os << std::flush;
			auto vel_dist_iter_buf = vel_dist_iter;
			if (++vel_dist_iter_buf == vel_dist.end()) out = false;
			++vel_dist_iter;
		} while(out);
		vel_dist_iter = vel_dist.begin();
		os << '\n';
		return os;
	}
};




#endif /* STAR_POLYMERS_LIB_ANALYSIS_VEL_DISTRIBUTION_H_ */
