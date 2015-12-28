#include "Molecule.h"

Molecule::Molecule(unsigned N) :
	Ekin { },
	Epot { },
	NumberOfMonomers { N },
	AType { },
	BType { },
	Arms { }
	{
		Monomers.reserve(N);
	}

Molecule::Molecule(unsigned N, double mass) :
		Ekin { },
		Epot { },
		NumberOfMonomers { N },
		AType { },
		BType { },
		Arms { }
		{
			Monomers.reserve(N);
			for ( unsigned i = 0 ; i < N ; i++ ) {
			    Monomers.push_back(MDParticle(mass, 0));
			}
		}

Molecule::Molecule(const Molecule& other) :
		Ekin {other.Ekin},
		Epot {other.Epot},
		NumberOfMonomers { other.NumberOfMonomers },
		AType {other.AType},
		BType {other.BType},
		Arms {other.Arms} {
			Monomers.reserve(NumberOfMonomers);
			for(unsigned i = 0; i < NumberOfMonomers; i++ ) {
				Monomers.push_back(other.Monomers[i]);
			}
 		}

Molecule& Molecule::operator = (const Molecule& other) {
	Ekin = other.Ekin;
	Epot = other.Epot;
	NumberOfMonomers = other.NumberOfMonomers;
	AType = other.AType;
	BType = other.BType;
	Arms = other.Arms;
	Monomers.clear();
	Monomers.reserve(NumberOfMonomers);
	for(unsigned i = 0; i < NumberOfMonomers; i++ ) {
		Monomers.push_back(other.Monomers[i]);
	}
	return *this;
}

MDParticle& Molecule::operator [](int i) {return Monomers[i];}

const MDParticle& Molecule::operator [](int i) const {return Monomers[i];}

double Molecule::calculate_Ekin() {
	Ekin = 0.0;
	Vector3d com_vel {0., 0., 0.,};
	for (auto& m : Monomers) {
		com_vel += m.Velocity;
	}
	com_vel /= NumberOfMonomers;
	for (auto& m : Monomers) {
		Vector3d rel_vel {m.Velocity - com_vel};
		//Ekin += m.Mass*(m.Velocity.dot(m.Velocity));
		Ekin += m.Mass*(rel_vel.dot(rel_vel));
	}
	Ekin /= 2.0;
	return Ekin;
}

Vector3d Molecule::calculate_center_of_mass() {
	Vector3d com {Vector3d::Zero()};
	for (auto& mono : Monomers) {
		com += mono.Position;
	}
	com /= NumberOfMonomers;
	return com;
}

Vector3d Molecule::calculate_center_of_mass_velocity() {
	Vector3d com {Vector3d::Zero()};
	for (auto& mono : Monomers) {
		com += mono.Velocity;
	}
	com /= NumberOfMonomers;
	return com;
}


void Molecule::initialize_straight_chain(unsigned A, unsigned B, double temperature, double bondLength) {
	Vector3d tempPos{Vector3d::Zero()};
	Vector3d average_vel{Vector3d::Zero()};
	AType = A;
	BType = B;
	for (unsigned i = A; i < (A+B); i++) Monomers[i].AmphiType = 1;
	for (unsigned i = 0; i < NumberOfMonomers; i++) {
		Monomers[i].Position = tempPos;
		for (unsigned j = 0; j < 3; j++) {
			Monomers[i].Velocity(j) = Rand::real_uniform() -0.5;
		}
		average_vel += Monomers[i].Velocity;
		if (i != (NumberOfMonomers - 1)) {
 			Monomers[i].set_neighbor(Monomers[i+1]);
		}
		tempPos(0) += bondLength;
	}
	average_vel /= NumberOfMonomers;
	for (auto& mono : Monomers) {
		mono.Velocity -= average_vel;
	}
	calculate_Ekin();
	std::cout<< "ekin before: " << Ekin;
	double vel_scale = sqrt(3.*(double)NumberOfMonomers*temperature/(Ekin*2.0));
	std::cout << " vel scale: " << vel_scale;
	for (auto& mono : Monomers) {
		mono.Velocity *= vel_scale;
	}
	std::cout << " ekin after: " << calculate_Ekin() << std::endl;
}

void Molecule::initialize_open_star(unsigned A, unsigned B, unsigned arms, double Temperature, double Bond, double AnchorBond) {
	Vector3d direction{Vector3d::Zero()};
	Vector3d average_vel{Vector3d::Zero()};
	AType = A;
	BType = B;
	Arms = arms;
	double phi { }, theta { }, scale { };
	unsigned n {(unsigned)ceil(sqrt(Arms))};
	unsigned arm_count {0};
	unsigned index { };
	unsigned monomers_per_arm {A+B};
	Monomers[0].Anchor = 1; //set anchor monomer;
	for (unsigned i = 0; i < n; i++) {
		phi = 2.*M_PI*i/n;
		for (unsigned j = 1; j < n+1; j++) {
			if (arm_count >= Arms) break;
			theta = M_PI*j/(n+1);
			direction(0) = cos(phi)*sin(theta);
			direction(1) = sin(phi)*sin(theta);
			direction(2) = cos(theta);
			for (unsigned k = 1; k <= monomers_per_arm; k++) {
				index = arm_count*monomers_per_arm+k;
				if (k > A) Monomers[index].AmphiType = 1;
				if (k==1) {
					Monomers[index].Position = AnchorBond*direction;
					Monomers[0].set_neighbor(Monomers[index]);
					Monomers[index].set_neighbor(Monomers[index+1]);
				}
				else {
					Monomers[index].Position = (AnchorBond+(k-1)*Bond)*direction;
					if (k != monomers_per_arm) Monomers[index].set_neighbor(Monomers[index+1]);
				}
				for (unsigned l = 0; l < 3; l++) {
					Monomers[index].Velocity(l) = Rand::real_uniform() -0.5;
				}
				average_vel += Monomers[index].Velocity;

			}
			arm_count++;
		}
	}
	average_vel /= NumberOfMonomers;
	for (auto& mono : Monomers) {
		mono.Velocity -= average_vel;
	}
	calculate_Ekin();
	scale = sqrt(3.*(double)NumberOfMonomers*Temperature/(Ekin*2.0));
	for (auto& mono : Monomers) {
		mono.Velocity *= scale;
	}
}

void Molecule::star_from_file(string filename, unsigned A, unsigned B, unsigned arms, bool set_zero) {
	ifstream input {filename};
	AType = A;
	BType = B;
	Arms = arms;
	Monomers[0].Anchor = 1;
	string line { };
	std::getline(input, line);
	std::cout << line << endl;
	for (auto& mono : Monomers) {
		std::getline(input, line);

		mono.Position(0) = stod(line.substr(31, 7));
		mono.Position(1) = stod(line.substr(39, 7));
		mono.Position(2) = stod(line.substr(47, 7));
		mono.Velocity(0) = stod(line.substr(55, 7));
		mono.Velocity(1) = stod(line.substr(63, 7));
		mono.Velocity(2) = stod(line.substr(71, 7));
	}
	unsigned monomers_per_arm {AType+BType};
	for (unsigned i = 0; i < Arms; i++) {
		for (unsigned k = 1; k <= monomers_per_arm; k++) {
			unsigned index { i*monomers_per_arm+k };
			if (k > A) Monomers[index].AmphiType = 1;
			if (k==1) {
				Monomers[0].set_neighbor(Monomers[index]);
				Monomers[index].set_neighbor(Monomers[index+1]);
			}
			else {
				if (k != monomers_per_arm) Monomers[index].set_neighbor(Monomers[index+1]);
			}
		}
	}
	if (set_zero) {
		Vector3d com {calculate_center_of_mass()};
        Vector3d com_vel {calculate_center_of_mass_velocity()};
        std::cout << "center of mass before: " << com << std::endl;
        for (auto& mono : Monomers) {
        	mono.Position -= com;
        	mono.Velocity -= com_vel;
        }
	}
	/*Center = BoxCenter - Monomers[0].Position;
	for (auto& mono : Monomers) {
		mono.Position += Center;
	}*/
}

std::ostream& Molecule::print(std::ostream& os) const {
	for (unsigned i = 0; i < NumberOfMonomers; i++) {
		os << Monomers[i].Position.transpose() <<" ; " << Monomers[i].Velocity.transpose() << " ; " << Monomers[i].Force.transpose() << '\n';
	}
	return os;
}

std::ostream& operator <<(std::ostream& os, const Molecule& some) {
	return some.print(os);
}

template<class UnitaryFunc>
UnitaryFunc Molecule::unitary(UnitaryFunc&& func) const {
	return for_each(Monomers.cbegin(), Monomers.cend(), func);
}

template<class UnaryFunc, class BinaryFunc>
void Molecule::operator() (UnaryFunc& ufunc, BinaryFunc& bfunc) const {
	auto first = Monomers.cbegin(), last = Monomers.cend();
	auto second = first;

	for(; first != last; ++first) {
		ufunc( *first );
		for( second = first + 1; second != last; ++second )
			bfunc( *first, *second );
	}
}
