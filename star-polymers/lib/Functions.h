#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>

template<typename T>
inline T my_modulus(T a, T b) {
	if (a < 0) a += b;
	else if (a >= b) a -= b;
	return a;
}
#endif
