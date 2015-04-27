#include "Box.h"

Box::Box(double Lx, double Ly, double Lz, double temperature, double lambda) :
	SystemTime { },
	Temperature { temperature },
	Lambda { lambda },
	Cutoff { 1.5 },
	VerletRadius { 2.0 },
	VerletRadius2 { 4.0 },
	NumberOfMonomers { } {
		Size[0] = Lx;
		Size[1] = Ly;
		Size[2] = Lz;
		CellSize[0] = int(Size[0]/VerletRadius);
		CellSize[1] = int(Size[1]/VerletRadius);
		CellSize[2] = int(Size[2]/VerletRadius);
		CellSideLength[0] = Size[0]/(double)CellSize[0];
		CellSideLength[1] = Size[1]/(double)CellSize[1];
		CellSideLength[2] = Size[2]/(double)CellSize[2];
		CellList = std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>>(CellSize[0], std::vector<std::vector<std::forward_list<MDParticle*>>>(CellSize[1], std::vector<std::forward_list<MDParticle*>>(CellSize[2], std::forward_list<MDParticle*>())));
}



inline Vector3d& Box::wrap(Vector3d& pos) {
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= floor(pos(i)/Size[i]) * Size[i];
	}
	return pos;
}

inline Vector3d Box::wrap(Vector3d&& pos) {
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= floor(pos(i)/Size[i]) * Size[i];
	}
	return pos;
}

void Box::wrap(Particle& part) {
	part.Position = wrap(part.Position);
}

inline void Box::wrap(Molecule& mol) {
	for (auto& mono : mol.Monomers) wrap(mono);
}

void Box::wrap() {
	for (auto& mol : Molecules) wrap(mol);
}

Vector3d Box::relative_position(Particle& one, Particle& two) {
	Vector3d result;
	result = two.Position - one.Position;
	for (unsigned i = 0; i < 3; i++) {
		result(i) -= round(result(i)/Size[i]) * Size[i];
	}
	return result;
}

void Box::add_chain(unsigned N, double mass, double bondLength) {
	Molecules.push_back(Molecule{N, mass});
	Molecules.back().initialize_straight_chain(N, 0, bondLength, Temperature);
	wrap(Molecules.back());
	NumberOfMonomers += N;
}

void Box::add_chain(unsigned A, unsigned B, double mass, double bondLength) {
	Molecules.push_back(Molecule{A+B, mass});
	Molecules.back().initialize_straight_chain(A, B, bondLength, Temperature);
	wrap(Molecules.back());
	NumberOfMonomers += (A+B);
}


unsigned Box::numberOfMonomers() {
	unsigned N { };
	for(auto& mol: Molecules) {
		N += mol.Monomers.size();
	}
	return N;
}
std::ostream& Box::print_molecules(std::ostream& os) const {
	for (auto& mol : Molecules) {
		os << mol << std::endl;
	}
	return os;
}

std::ostream& Box::print_Epot(std::ostream& os) const {
	double PotentialEnergy { };
	for (auto& mol : Molecules) {
		PotentialEnergy += mol.Epot;
	}
	PotentialEnergy /= NumberOfMonomers;
	os << PotentialEnergy<< " ";
	return os;
}

double Box::calculate_ekin() {
	double KineticEnergy { };
	for (auto& mol : Molecules) {
		KineticEnergy += mol.calculate_Ekin();
	}

	return KineticEnergy;
}

double Box::calculate_radius_of_gyration() {
	double r_gyr { };
	Vector3d distance;
	for (auto& mol : Molecules) {
		double r_gyr_mol { };
		for (unsigned i = 0; i < mol.Monomers.size(); i++) {
			for (unsigned j = i+1; j < mol.Monomers.size(); j++) {
				distance = relative_position(mol[i], mol[j]);
				r_gyr_mol += distance.dot(distance);
			}
		}
		r_gyr_mol /= (mol.Monomers.size()*mol.Monomers.size());
		r_gyr_mol = sqrt(r_gyr_mol);
		r_gyr += r_gyr_mol;
	}
	r_gyr /= Molecules.size();
	return r_gyr;
}


std::ostream& Box::print_Ekin(std::ostream& os) {
	double KineticEnergy { };
	for (auto& mol : Molecules) {
		KineticEnergy += mol.calculate_Ekin();
	}
	KineticEnergy /= NumberOfMonomers;
	os << KineticEnergy << " ";
	return os;
}

std::ostream& Box::print_Temperature(std::ostream& os) {
	double Temperature { };
	for (auto& mol : Molecules) {
		Temperature += mol.calculate_Ekin()/mol.NumberOfMonomers;
	}
	Temperature *= 2./3.;
	os << Temperature << " ";
	return os;
}

std::ostream& Box::print_radius_of_gyration(std::ostream& os) {
	os << calculate_radius_of_gyration() << " ";
	return os;
}

void Box::calculate_forces(bool calc_epot) {
	double radius2 { };
	double force_abs { };
	Vector3d distance;
	Vector3d force;
	for (auto& mol : Molecules) {
		if (calc_epot) mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
		for (unsigned i = 0; i < mol.Monomers.size(); i++) {
			for (unsigned j = i + 1; j < mol.Monomers.size(); j++) {
				distance = relative_position(mol[i], mol[j]);
				radius2 = distance.dot(distance);
				if (mol[i].AmphiType == 1 && mol[j].AmphiType == 1) { // BB Type
					if (calc_epot) mol.Epot += TypeBB_Potential(radius2, Lambda);
					force_abs = TypeBB_Force(radius2, Lambda);
					force = distance*force_abs;
					mol[i].Force -= force;
					mol[j].Force += force;
				}
				else { // AA Type
					if (calc_epot) mol.Epot += TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					force = distance*force_abs;
					mol[i].Force -= force;
					mol[j].Force += force;
				}
			}
			for (auto& neighbor : mol[i].Neighbors) { //Fene bonds
				distance = relative_position(mol[i], *neighbor);
				radius2 = distance.dot(distance);
				if (calc_epot) mol.Epot += Fene_Potential(radius2);
				force_abs = Fene_Force(radius2);
				force = distance*force_abs;
				mol[i].Force -= force;
				neighbor -> Force += force;

			}
		}
	}
}

void Box::calculate_forces_verlet(bool calc_epot) {
	double radius2 { };
	double force_abs { };
	Vector3d distance { };
	Vector3d force { };
	for (auto& mol : Molecules) {
		mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
	}
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {

			for (auto& other : mono.VerletList) {
				distance = relative_position(mono, *other);
				radius2 = distance.dot(distance);
				if (mono.AmphiType == 1 && other -> AmphiType == 1) { // BB Type
					if (calc_epot) mol.Epot += 0.5*TypeBB_Potential(radius2, Lambda);
					force_abs = TypeBB_Force(radius2, Lambda);
					force = distance*force_abs;
					mono.Force -= force;
				}
				else { // AA Type
					if (calc_epot) mol.Epot += 0.5*TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					force = distance*force_abs;
					mono.Force -= force;
				}
			}
			for (auto& neighbor : mono.Neighbors) { //Fene bonds
				distance = relative_position(mono, *neighbor);
				radius2 = distance.dot(distance);
				if (calc_epot) mol.Epot += Fene_Potential(radius2);
				force_abs = Fene_Force(radius2);
				force = distance*force_abs;
				mono.Force -= force;
				neighbor -> Force += force;
			}
		}
	}
}

void Box::update_VerletLists() {
	std::array<int,3> CellNumber { };
	double radius2 { };
	Vector3d distance;
	//clear all Lists;
	for (auto& sheet : CellList) {
		for (auto& row : sheet) {
			for (auto& list : row) {
				list.clear();
			}
		}
	}
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.clear_VerletList();
		}
	}

	//sort into CellLists
	for (auto& mol: Molecules) {
		for (auto& mono : mol.Monomers) {
			for (int i = 0; i < 3; i++) {
				CellNumber[i] = (int)(mono.Position(i)/CellSideLength[i]);
			}
			mono.VerletPosition = mono.Position;
			CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono);
		}
	}

	//make VerletLists
	int p { }, q { }, r { };
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			for (int i = 0; i < 3; i++) {
				CellNumber[i] = (int)(mono.Position(i)/CellSideLength[i]);
			}
			for (int j = CellNumber[0]-1; j < CellNumber[0]+2; j++) {
				for (int k = CellNumber[1]-1; k < CellNumber[1]+2; k++) {
					for (int l = CellNumber[2]-1; l < CellNumber[2]+2; l++) {

						p = my_modulus(j, CellSize[0]);
						q = my_modulus(k, CellSize[1]);
						r = my_modulus(l, CellSize[2]);

						for (auto& other : CellList[p][q][r]) {
							if (other == &mono) continue;
							distance = relative_position(mono, *other);
							radius2 = distance.dot(distance);
							if (radius2 <= VerletRadius2) {
								mono.VerletList.push_front(other);
							}
						}
					}
				}
			}
		}
	}
}

void Box::check_VerletLists() {
	Vector3d displacement { };
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			displacement = mono.Position - mono.VerletPosition;
			for (unsigned i = 0; i < 3; i++) {
				displacement(i) -= round(displacement(i)/Size[i]) * Size[i];
			}
			if (displacement.norm() > (VerletRadius - Cutoff)*0.5) {
				update_VerletLists();
				return;
			}
		}
	}
}

