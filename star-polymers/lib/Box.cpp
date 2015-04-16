#include "Box.h"

Box::Box(double Lx, double Ly, double Lz, double temperature, double lambda) :
	SystemTime { },
	Temperature { temperature },
	Lambda { lambda },
	Cutoff { 1.5 },
	VerletRadius { 2.0 },
	VerletRadius2 { 4.0 },
	NumberOfMonomers { },
	count { } {
		Size[0] = Lx;
		Size[1] = Ly;
		Size[2] = Lz;
		CellSize[0] = int(Size[0]/Cutoff);
		CellSize[1] = int(Size[1]/Cutoff);
		CellSize[2] = int(Size[2]/Cutoff);
		CellList = std::vector<std::vector<std::vector<std::forward_list<Particle*>>>>(CellSize[0], std::vector<std::vector<std::forward_list<Particle*>>>(CellSize[1], std::vector<std::forward_list<Particle*>>(CellSize[2], std::forward_list<Particle*>())));
}



inline MatVec& Box::wrap(MatVec& pos) {
	return pos -=  floor(pos/Size) % Size;
}

inline MatVec Box::wrap(MatVec&& pos) {
	return pos -=  floor(pos/Size) % Size;
}

inline void Box::wrap(Particle& part) {
	part.Position = wrap(part.Position);
}

inline void Box::wrap(Molecule& mol) {
	for (auto& mono : mol.Monomers) wrap(mono);
}

void Box::wrap() {
	for (auto& mol : Molecules) wrap(mol);
}

MatVec Box::relative_position(Particle& one, Particle& two) {
	MatVec result { };
	result = two.Position - one.Position;
	result -= round(result/Size) % Size;
	return result;
}

void Box::add_chain(unsigned N, double mass, double bondLength) {
	Molecules.push_back(Molecule{N, mass});
	Molecules.back().initialize_straight_chain(bondLength, Temperature);
	wrap(Molecules.back());
	calculate_forces();
	NumberOfMonomers += N;
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

void Box::calculate_forces() {
	double radius2 { };
	double force_abs { };
	MatVec distance { };
	MatVec force { };
	count = 0.0;
	for (auto& mol : Molecules) {
		mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
		for (unsigned i = 0; i < mol.Monomers.size(); i++) {
			for (unsigned j = i + 1; j < mol.Monomers.size(); j++) {
				distance = relative_position(mol[i], mol[j]);
				radius2 = distance*distance;
				if (mol[i].AmphiType == 1 && mol[j].AmphiType == 1) { // BB Type
					mol.Epot += TypeBB_Potential(radius2, Lambda);
					force_abs = TypeBB_Force(radius2, Lambda);
					force = distance*force_abs;
					mol[i].Force -= force;
					mol[j].Force += force;
				}
				else { // AA Type
					mol.Epot += TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					if (force_abs > 0) count++;
					force = distance*force_abs;
					mol[i].Force -= force;
					mol[j].Force += force;
				}
			}
			for (auto& neighbor : mol[i].Neighbors) { //Fene bonds
				distance = relative_position(mol[i], *neighbor);
				radius2 = distance*distance;
				mol.Epot += 0.5*Fene_Potential(radius2);
				force_abs = Fene_Force(radius2);
				force = distance*force_abs;
				mol[i].Force -= force;

			}
		}
	}
	std::cout << count*2 << std::endl;
}

void Box::calculate_forces_verlet() {
	double radius2 { };
	double force_abs { };
	MatVec distance { };
	MatVec force { };
	for (auto& mol : Molecules) {
		mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
		for (auto& mono : mol.Monomers) {
			for (auto& other : mono.VerletList) {
				distance = relative_position(mono, *other);
				radius2 = distance*distance;
				if (mono.AmphiType == 1 && other -> AmphiType == 1) { // BB Type
					mol.Epot += TypeBB_Potential(radius2, Lambda);
					force_abs = TypeBB_Force(radius2, Lambda);
					force = distance*force_abs;
					mono.Force -= force;
				}
				else { // AA Type
					mol.Epot += TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					force = distance*force_abs;
					mono.Force -= force;
				}
			}
			for (auto& neighbor : mono.Neighbors) { //Fene bonds
				distance = relative_position(mono, *neighbor);
				radius2 = distance*distance;
				mol.Epot += 0.5*Fene_Potential(radius2);
				force_abs = Fene_Force(radius2);
				force = distance*force_abs;
				mono.Force -= force;
			}
		}
	}
}

void Box::update_VerletLists() {
	std::array<int,3> CellNumber { };
	double radius2 { };
	MatVec distance { };
	//clear all Lists;
	for (auto& sheet : CellList) {
		for (auto& row : sheet) {
			for (auto& list : row) {
				/*for (auto pointer : list) {
					delete(pointer);
				}*/
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
				CellNumber[i] = (int)(mono.Position[i]/CellSize[i]);
			}
			mono.VerletPosition = mono.Position;
			CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono);
		}
	}

	//make VerletLists
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			for (int i = 0; i < 3; i++) {
				CellNumber[i] = (int)(mono.Position[i]/CellSize[i]);
			}
			for (int j = CellNumber[0]-1; j < CellNumber[0]+2; j++) {
				for (int k = CellNumber[1]-1; k < CellNumber[1]+2; k++) {
					for (int l = CellNumber[2]-1; l < CellNumber[2]+2; l++) {
						for (auto& other : CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]]) {
							if (other == &mono) continue;
							distance = relative_position(mono, *other);
							radius2 = distance*distance;
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
	MatVec displacement { };
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			displacement = mono.Position - mono.VerletPosition;
			displacement -= round(displacement/Size) % Size;
			if (displacement.norm() > (VerletRadius - Cutoff)) {
				update_VerletLists();
				return;
			}
		}
	}
}

