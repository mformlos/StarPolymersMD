#include "MatVec.h"

	auto MatVec::begin() -> decltype(Vec.begin()) {
		return Vec.begin();
	}

	auto MatVec::end() -> decltype(Vec.end()) {
		return Vec.end();
	}

	MatVec::MatVec() : Vec{} {
	};

	MatVec::MatVec(std::initializer_list<double> elements) : MatVec { } {
		auto el_begin = elements.begin(), el_end = elements.end();
		auto vec_begin = Vec.begin(), vec_end = Vec.end();
		for (; el_begin != el_end && vec_begin != vec_end; ++el_begin, ++vec_begin) {
			*vec_begin = *el_begin;
		}
	}
	MatVec::MatVec(const MatVec& other) : Vec{} {
		for (unsigned i = 0; i < 3; i++) {
			Vec[i] = other.Vec[i];
		}
	}

	MatVec::MatVec(MatVec&& other) : Vec{} {
		Vec.swap(other.Vec);
	}

	double& MatVec::operator [](int i) {
		return Vec[i];
	}

	const double& MatVec::operator [](int i) const {
		return Vec[i];
	}


	MatVec& MatVec::operator =(const MatVec& other) {
		Vec = other.Vec;
		return *this;
	}

	MatVec& MatVec::operator =(MatVec&& other) {
		Vec.swap(other.Vec);
		return *this;
	}

	bool MatVec::operator ==(const MatVec& other) const {
		return Vec == other.Vec;
	}

	bool MatVec::operator !=(const MatVec& other) const {
		return Vec != other.Vec;
	}

	MatVec& MatVec::operator +=(const MatVec& other) {
		for (unsigned i = 0; i < 3; i++) {
			Vec[i] += other.Vec[i];
		}
		return *this;
	}

	MatVec& MatVec::operator -=(const MatVec& other) {
		for (unsigned i = 0; i < 3; i++) {
			Vec[i] -= other.Vec[i];
		}
		return *this;
	}

	MatVec& MatVec::operator *=(const double& scale) {
		for (unsigned i = 0; i < 3; i++) {
			Vec[i] *= scale;
		}
		return *this;
	}

	MatVec& MatVec::operator /=(const double& scale) {
		for (unsigned i = 0; i < 3; i++) {
			Vec[i] /= scale;
		}
		return *this;
	}

	MatVec MatVec::operator +(const MatVec& other) const {
		MatVec result{*this};
		result += other;
		return result;
	}

	MatVec MatVec::operator -(const MatVec& other) const {
		MatVec result{*this};
		result -= other;
		return result;
	}


	MatVec MatVec::operator *(const double& scale) const {
		MatVec result{*this};
		result *= scale;
		return result;
	}

	MatVec MatVec::operator /(const double& scale) const {
		MatVec result{*this};
		result /= scale;
		return result;
	}

	double MatVec::operator *(const MatVec& other) const {
		double result{};
		for (unsigned i = 0; i < 3; i++) {
			result += Vec[i]*other[i];
		}
		return result;
	}

	double MatVec::norm2() {
		double result{};
		for (unsigned i = 0; i < 3; i++) {
			result += Vec[i]*Vec[i];
		}
		return result;
	}

	double MatVec::norm() {
		return sqrt(norm2());
	}

	std::ostream& MatVec::print(std::ostream& os) const {
		for (auto& el : Vec) {
		os << el << " ";
		}
		return os;
}

std::ostream& operator <<(std::ostream& os, const MatVec& some) {
	return some.print(os);
}
