/*
 * Analysis.h
 *
 *  Created on: Apr 28, 2015
 *      Author: maud
 */

#ifndef STAR_POLYMERS_LIB_ANALYSIS_H_
#define STAR_POLYMERS_LIB_ANALYSIS_H_

#include <ostream>

template<class ... Args>
class Analysis {
public:
	virtual ~Analysis() {}
	virtual void operator() (const Args& ... args){}
	virtual double value(){return -1.0;}
	virtual std::ostream& print_result(std::ostream& os) {return os;}
};


#endif /* STAR_POLYMERS_LIB_ANALYSIS_H_ */
