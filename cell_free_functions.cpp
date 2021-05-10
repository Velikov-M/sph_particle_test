#define _USE_MATH_DEFINES

#include "free_functions.h"
#include "Particle.h"
#include "Global_consts.h"
#include "cell_free_functions.h"
#include <math.h>
#include <cmath>
#include <iostream> 
#include <stdexcept>

double cell_calculate_density(std::vector<Particle>& particles, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end) { //i_end is length - 1
	double density = 0;
	if (i_cell == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == 0) { // We need special treatment for left (i==0) and right (i==i_end) periodic BC
					int cell_size;
					try {
						cell_size = cells.at(i_end).at(j_cell + j - 1).size(); //Cheaky way of dealing with j==0 and j==j_end, in this cells we'll try to call non-existent cell
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[i_end][j_cell + j - 1][k];
						if (particles[j_part].get_type() == 'S') continue; // stone particles are not counted in density of fluid
						std::vector<double> j_pos = particles[j_part].get_position();
						std::vector<double> id_pos = particles[id_part].get_position();
						double distance = sqrt((id_pos[0] - j_pos[0] + length_of_area) * (id_pos[0] - j_pos[0] + length_of_area) + (id_pos[1] - j_pos[1]) * (id_pos[1] - j_pos[1])); //we can't use honest distance function!
						density += weight_function(distance, b);
					}
					continue;
				}
				int cell_size;
				try {
					cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
				}
				catch (const std::out_of_range& e) {
					continue;
				}
				if (cell_size == 0) continue;
				for (int k = 0; k < cell_size; k++) {
					int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
					if (particles[j_part].get_type() == 'S') continue;
					density += weight_function(calculate_distance(particles[id_part], particles[j_part]), b);
				}
			}
		}
		density *= mass;
		return density;
	}

	if (i_cell == i_end) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == 2) {
					int cell_size;
					try {
						cell_size = cells.at(0).at(j_cell + j - 1).size();
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[0][j_cell + j - 1][k];
						if (particles[j_part].get_type() == 'S') continue;
						std::vector<double> j_pos = particles[j_part].get_position();
						std::vector<double> id_pos = particles[id_part].get_position();
						double distance = sqrt((id_pos[0] - j_pos[0] - length_of_area) * (id_pos[0] - j_pos[0] - length_of_area) + (id_pos[1] - j_pos[1]) * (id_pos[1] - j_pos[1]));
						density += weight_function(distance, b);
					}
					continue;
				}
				int cell_size;
				try {
					cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
				}
				catch (const std::out_of_range& e) {
					continue;
				}
				if (cell_size == 0) continue;
				for (int k = 0; k < cell_size; k++) {
					int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
					if (particles[j_part].get_type() == 'S') continue;
					density += weight_function(calculate_distance(particles[id_part], particles[j_part]), b);
				}
			}
		}
		density *= mass;
		return density;
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			int cell_size;
			try {
				cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
			}
			catch (const std::out_of_range& e) {
				continue;
			}
			if (cell_size == 0) continue;
			for (int k = 0; k < cell_size; k++) {
				int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
				if (particles[j_part].get_type() == 'S') continue;
				density += weight_function(calculate_distance(particles[id_part], particles[j_part]), b);
			}
		}
	}
	density *= mass;
	return density;
}


std::vector<std::vector<double>> cell_calc_deform_tensor(std::vector<Particle>& particles, std::vector<double>& densities, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end) {
	std::vector < std::vector<double>> deform_tensor(2);
	deform_tensor[0] = { 0, 0 };
	deform_tensor[1] = { 0, 0 };
	std::vector<double> i_vel = particles[id_part].get_velosity();
	std::vector<double> i_pos = particles[id_part].get_position();
	if (i_cell == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == 0) { // We need special treatment for left (i==0) and right (i==i_end) periodic BC
					int cell_size;
					try {
						cell_size = cells.at(i_end).at(j_cell + j - 1).size(); //Cheaky way of dealing with j==0 and j==j_end, in this cells we'll try to call non-existent cell
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[i_end][j_cell + j - 1][k];
						if (particles[j_part].get_type() == 'S') continue; // stone particles are not counted in density of fluid
						std::vector<double> j_pos = particles[j_part].get_position();
						std::vector<double> j_vel = particles[j_part].get_velosity();
						std::vector<double> ji_vel = { j_vel[0] - i_vel[0], j_vel[1] - i_vel[1] };
						std::vector<double> ij_pos = { i_pos[0] - j_pos[0] + length_of_area, i_pos[1] - j_pos[1] };
						double ji_distance = sqrt(ij_pos[0] * ij_pos[0] + ij_pos[1] * ij_pos[1]); //we can't use honest distance function!

						deform_tensor[0][0] += (4.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[0][0] -= (2.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						// next, we'll calculate 0,1 component, as well as apply symmetry of linear deformation tensor
						deform_tensor[0][1] += mass * ji_vel[0] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[0][1] += mass * ji_vel[1] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[1][0] = deform_tensor[0][1];
						//and finish calculation for this ij-pair, by finding 1,1 component
						deform_tensor[1][1] += (4.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[1][1] -= (2.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					}
					continue;
				}
				int cell_size;
				try {
					cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
				}
				catch (const std::out_of_range& e) {
					continue;
				}
				if (cell_size == 0) continue;
				for (int k = 0; k < cell_size; k++) {
					int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
					if (id_part == j_part) continue;
					if (particles[j_part].get_type() == 'S') continue;
					std::vector<double> j_vel = particles[j_part].get_velosity();
					std::vector<double> j_pos = particles[j_part].get_position();
					std::vector<double> ji_vel = { j_vel[0] - i_vel[0], j_vel[1] - i_vel[1] };
					std::vector<double> ij_pos = { i_pos[0] - j_pos[0], i_pos[1] - j_pos[1] };
					double ji_distance = calculate_distance(particles[id_part], particles[j_part]);
					// let's calcualte add to 0,0 component of deformation tensor, u should ask question about formula (4.48), opening that delta part
					deform_tensor[0][0] += (4.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[0][0] -= (2.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					// next, we'll calculate 0,1 component, as well as apply symmetry of linear deformation tensor
					deform_tensor[0][1] += mass * ji_vel[0] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[0][1] += mass * ji_vel[1] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[1][0] = deform_tensor[0][1];
					//and finish calculation for this ij-pair, by finding 1,1 component
					deform_tensor[1][1] += (4.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[1][1] -= (2.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
				}
			}
		}
		return deform_tensor;
	}

	if (i_cell == i_end) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == 2) { // We need special treatment for left (i==0) and right (i==i_end) periodic BC
					int cell_size;
					try {
						cell_size = cells.at(0).at(j_cell + j - 1).size(); //Cheaky way of dealing with j==0 and j==j_end, in this cells we'll try to call non-existent cell
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[0][j_cell + j - 1][k];
						if (particles[j_part].get_type() == 'S') continue; // stone particles are not counted in density of fluid
						std::vector<double> j_pos = particles[j_part].get_position();
						std::vector<double> j_vel = particles[j_part].get_velosity();
						std::vector<double> ji_vel = { j_vel[0] - i_vel[0], j_vel[1] - i_vel[1] };
						std::vector<double> ij_pos = { i_pos[0] - j_pos[0] - length_of_area, i_pos[1] - j_pos[1] };
						double ji_distance = sqrt(ij_pos[0] * ij_pos[0] + ij_pos[1] * ij_pos[1]); //we can't use honest distance function!

						deform_tensor[0][0] += (4.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[0][0] -= (2.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						// next, we'll calculate 0,1 component, as well as apply symmetry of linear deformation tensor
						deform_tensor[0][1] += mass * ji_vel[0] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[0][1] += mass * ji_vel[1] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[1][0] = deform_tensor[0][1];
						//and finish calculation for this ij-pair, by finding 1,1 component
						deform_tensor[1][1] += (4.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
						deform_tensor[1][1] -= (2.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					}
					continue;
				}
				int cell_size;
				try {
					cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
				}
				catch (const std::out_of_range& e) {
					continue;
				}
				if (cell_size == 0) continue;
				for (int k = 0; k < cell_size; k++) {
					int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
					if (particles[j_part].get_type() == 'S') continue;
					if (id_part == j_part) continue;
					std::vector<double> j_vel = particles[j_part].get_velosity();
					std::vector<double> j_pos = particles[j_part].get_position();
					std::vector<double> ji_vel = { j_vel[0] - i_vel[0], j_vel[1] - i_vel[1] };
					std::vector<double> ij_pos = { i_pos[0] - j_pos[0], i_pos[1] - j_pos[1] };
					double ji_distance = calculate_distance(particles[id_part], particles[j_part]);
					// let's calcualte add to 0,0 component of deformation tensor, u should ask question about formula (4.48), opening that delta part
					deform_tensor[0][0] += (4.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[0][0] -= (2.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					// next, we'll calculate 0,1 component, as well as apply symmetry of linear deformation tensor
					deform_tensor[0][1] += mass * ji_vel[0] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[0][1] += mass * ji_vel[1] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[1][0] = deform_tensor[0][1];
					//and finish calculation for this ij-pair, by finding 1,1 component
					deform_tensor[1][1] += (4.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
					deform_tensor[1][1] -= (2.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
				}
			}
		}
		return deform_tensor;
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			int cell_size;
			try {
				cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
			}
			catch (const std::out_of_range& e) {
				continue;
			}
			if (cell_size == 0) continue;
			for (int k = 0; k < cell_size; k++) {
				int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
				if (id_part == j_part) continue;
				if (particles[j_part].get_type() == 'S') continue;
				std::vector<double> j_vel = particles[j_part].get_velosity();
				std::vector<double> j_pos = particles[j_part].get_position();
				std::vector<double> ji_vel = { j_vel[0] - i_vel[0], j_vel[1] - i_vel[1] }; // ji - velosity, i am not so sure
				std::vector<double> ij_pos = { i_pos[0] - j_pos[0], i_pos[1] - j_pos[1] };
				double ji_distance = calculate_distance(particles[id_part], particles[j_part]);
				// let's calcualte add to 0,0 component of deformation tensor, u should ask question about formula (4.48), opening that delta part
				deform_tensor[0][0] += (4.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
				deform_tensor[0][0] -= (2.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
				// next, we'll calculate 0,1 component, as well as apply symmetry of linear deformation tensor
				deform_tensor[0][1] += mass * ji_vel[0] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
				deform_tensor[0][1] += mass * ji_vel[1] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
				deform_tensor[1][0] = deform_tensor[0][1];
				//and finish calculation for this ij-pair, by finding 1,1 component
				deform_tensor[1][1] += (4.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
				deform_tensor[1][1] -= (2.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j_part] * ji_distance);
			}
		}
	}
	return deform_tensor;
}


std::vector<double> cell_calc_stone_accs(std::vector<Particle>& particles, std::vector<std::vector<std::vector<int>>>& cells, std::vector<double> stone_pos, int n_exterior, int n_stone, int i_end, bool gravity_flag) {
	std::vector <double> stone_acc = { 0, 0, 0 }; // x[0]'', x[1]'', \phi''
	for (int p = 0; p < n_stone; p++) {
		int id_part = n_exterior + p;
		std::vector<int> cell_of_prt = particles[id_part].get_cell_id();
		std::vector<double> r = { particles[id_part].get_position()[0] - stone_pos[0], particles[id_part].get_position()[1] - stone_pos[1] };
		int i_cell = cell_of_prt[0];
		int j_cell = cell_of_prt[1];
		if (i_cell == 0) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					if (i == 0) { // We need special treatment for left (i==0) and right (i==i_end) periodic BC
						int cell_size;
						try {
							cell_size = cells.at(i_end).at(j_cell + j - 1).size(); //Cheaky way of dealing with j==0 and j==j_end, in this cells we'll try to call non-existent cell
						}
						catch (const std::out_of_range& e) {
							continue;
						}
						if (cell_size == 0) continue;
						for (int k = 0; k < cell_size; k++) {
							int j_part = cells[i_end][j_cell + j - 1][k];
							if (particles[j_part].get_type() == 'S') continue;
							std::vector<double> ij_pos = { particles[id_part].get_position()[0] - particles[j_part].get_position()[0] + length_of_area, particles[id_part].get_position()[1] - particles[j_part].get_position()[1] };
							double distance = sqrt(ij_pos[0] * ij_pos[0] + ij_pos[1] * ij_pos[1]);
							double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
							double force_x0 = force_ij * ij_pos[0] / distance;
							double force_x1 = force_ij * ij_pos[1] / distance;
							double momentum = r[0] * force_x1 - r[1] * force_x0;
							stone_acc[0] += force_x0 / mass_stone;
							stone_acc[1] += force_x1 / mass_stone;
							stone_acc[2] += momentum / moment_inertia_stone;
						}
					}
					int cell_size;
					try {
						cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
						if (id_part == j_part) continue;
						if (particles[j_part].get_type() == 'S') continue;
						double distance = calculate_distance(particles[id_part], particles[j_part]);
						double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
						double force_x0 = force_ij * (particles[id_part].get_position()[0] - particles[j_part].get_position()[0]) / distance;
						double force_x1 = force_ij * (particles[id_part].get_position()[1] - particles[j_part].get_position()[1]) / distance;
						double momentum = r[0] * force_x1 - r[1] * force_x0;
						stone_acc[0] += force_x0 / mass_stone;
						stone_acc[1] += force_x1 / mass_stone;
						stone_acc[2] += momentum / moment_inertia_stone;
					}
				}
			}
			continue;
		}
		if (i_cell == i_end) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					if (i == 2) { // We need special treatment for left (i==0) and right (i==i_end) periodic BC
						int cell_size;
						try {
							cell_size = cells.at(0).at(j_cell + j - 1).size(); //Cheaky way of dealing with j==0 and j==j_end, in this cells we'll try to call non-existent cell
						}
						catch (const std::out_of_range& e) {
							continue;
						}
						if (cell_size == 0) continue;
						for (int k = 0; k < cell_size; k++) {
							int j_part = cells[0][j_cell + j - 1][k];
							if (particles[j_part].get_type() == 'S') continue;
							std::vector<double> ij_pos = { particles[id_part].get_position()[0] - particles[j_part].get_position()[0] - length_of_area, particles[id_part].get_position()[1] - particles[j_part].get_position()[1] };
							double distance = sqrt(ij_pos[0] * ij_pos[0] + ij_pos[1] * ij_pos[1]);
							double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
							double force_x0 = force_ij * ij_pos[0] / distance;
							double force_x1 = force_ij * ij_pos[1] / distance;
							double momentum = r[0] * force_x1 - r[1] * force_x0;
							stone_acc[0] += force_x0 / mass_stone;
							stone_acc[1] += force_x1 / mass_stone;
							stone_acc[2] += momentum / moment_inertia_stone;
						}
					}
					int cell_size;
					try {
						cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
						if (id_part == j_part) continue;
						if (particles[j_part].get_type() == 'S') continue;
						double distance = calculate_distance(particles[id_part], particles[j_part]);
						double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
						double force_x0 = force_ij * (particles[id_part].get_position()[0] - particles[j_part].get_position()[0]) / distance;
						double force_x1 = force_ij * (particles[id_part].get_position()[1] - particles[j_part].get_position()[1]) / distance;
						double momentum = r[0] * force_x1 - r[1] * force_x0;
						stone_acc[0] += force_x0 / mass_stone;
						stone_acc[1] += force_x1 / mass_stone;
						stone_acc[2] += momentum / moment_inertia_stone;
					}
				}
			}
			continue;
		}
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				int cell_size;
				try {
					cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
				}
				catch (const std::out_of_range& e) {
					continue;
				}
				if (cell_size == 0) continue;
				for (int k = 0; k < cell_size; k++) {
					int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
					if (id_part == j_part) continue;
					if (particles[j_part].get_type() == 'S') continue;
					double distance = calculate_distance(particles[id_part], particles[j_part]);
					double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
					double force_x0 = force_ij * (particles[id_part].get_position()[0] - particles[j_part].get_position()[0]) / distance;
					double force_x1 = force_ij * (particles[id_part].get_position()[1] - particles[j_part].get_position()[1]) / distance;
					double momentum = r[0] * force_x1 - r[1] * force_x0;
					stone_acc[0] += force_x0 / mass_stone;
					stone_acc[1] += force_x1 / mass_stone;
					stone_acc[2] += momentum / moment_inertia_stone;
				}
			}
		}
	}
	if (gravity_flag) {
		stone_acc[1] -= g;
	}
	return stone_acc;
}


std::vector<std::vector<double>> cell_calc_newton_pt(std::vector<Particle>& particles, std::vector<double>& densities, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end) {
	std::vector<std::vector<double>> pres_tensor(2);
	std::vector<std::vector<double>> eulier_part = calc_eulier_pres_tensor(densities, id_part);
	std::vector<std::vector<double>> deformations = cell_calc_deform_tensor(particles, densities, cells, id_part, i_cell, j_cell, i_end);

	pres_tensor[0] = { 0, 0 };
	pres_tensor[1] = { 0, 0 };

	pres_tensor[0][0] = eulier_part[0][0] + mu * deformations[0][0];
	pres_tensor[0][1] = eulier_part[0][1] + mu * deformations[0][1];
	pres_tensor[1][0] = eulier_part[1][0] + mu * deformations[1][0];
	pres_tensor[1][1] = eulier_part[1][1] + mu * deformations[1][1];

	return pres_tensor;
}


std::vector<double> cell_calc_acc(std::vector<Particle>& particles, std::vector<std::vector<std::vector<double>>>& pres_tensors, std::vector<double>& densities, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end, bool gravity_flag) {

	std::vector<double> acceleration = { 0, 0 }; // in this function, u should apply gravity
	std::vector<std::vector<double>> own_pres_tensor = pres_tensors[id_part];
	double own_density = densities[id_part];
	if (i_cell == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == 0) { // We need special treatment for left (i==0) and right (i==i_end) periodic BC
					int cell_size;
					try {
						cell_size = cells.at(i_end).at(j_cell + j - 1).size(); //Cheaky way of dealing with j==0 and j==j_end, in this cells we'll try to call non-existent cell
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[i_end][j_cell + j - 1][k];
						std::vector<double> ij_pos = { particles[id_part].get_position()[0] - particles[j_part].get_position()[0] + length_of_area, particles[id_part].get_position()[1] - particles[j_part].get_position()[1] };
						double distance = sqrt(ij_pos[0] * ij_pos[0] + ij_pos[1] * ij_pos[1]);
						if (particles[j_part].get_type() == 'S') {
							double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
							acceleration[0] += force_ij * ij_pos[0] / distance;
							acceleration[1] += force_ij * ij_pos[1] / distance;
							continue;
						}
						std::vector<std::vector<double>> j_pres_tensor = pres_tensors[j_part];
						std::vector <double> ij_grad_weigh = { 0, 0 };
						ij_grad_weigh[0] = ij_pos[0] * ddr_weight_fun(distance, b) / distance;
						ij_grad_weigh[1] = ij_pos[1] * ddr_weight_fun(distance, b) / distance;
						double comp00 = (own_pres_tensor[0][0] / (own_density * own_density) + j_pres_tensor[0][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
						double comp01 = (own_pres_tensor[0][1] / (own_density * own_density) + j_pres_tensor[0][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
						double comp10 = (own_pres_tensor[1][0] / (own_density * own_density) + j_pres_tensor[1][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
						double comp11 = (own_pres_tensor[1][1] / (own_density * own_density) + j_pres_tensor[1][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
						acceleration[0] += comp00 + comp01;
						acceleration[1] += comp10 + comp11;
					}
					continue;
				}
				int cell_size;
				try {
					cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
				}
				catch (const std::out_of_range& e) {
					continue;
				}
				if (cell_size == 0) continue;
				for (int k = 0; k < cell_size; k++) {
					int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
					if (id_part == j_part) continue;
					if (particles[j_part].get_type() == 'S') {
						double distance = calculate_distance(particles[id_part], particles[j_part]);
						double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
						acceleration[0] += force_ij * (particles[id_part].get_position()[0] - particles[j_part].get_position()[0]) / distance;
						acceleration[1] += force_ij * (particles[id_part].get_position()[1] - particles[j_part].get_position()[1]) / distance;
						continue;
					}
					std::vector<std::vector<double>> j_pres_tensor = pres_tensors[j_part];
					std::vector<double> ij_grad_weigh = calc_grad_weight_fun(particles, id_part, j_part);
					double comp00 = (own_pres_tensor[0][0] / (own_density * own_density) + j_pres_tensor[0][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
					double comp01 = (own_pres_tensor[0][1] / (own_density * own_density) + j_pres_tensor[0][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
					double comp10 = (own_pres_tensor[1][0] / (own_density * own_density) + j_pres_tensor[1][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
					double comp11 = (own_pres_tensor[1][1] / (own_density * own_density) + j_pres_tensor[1][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
					acceleration[0] += comp00 + comp01;
					acceleration[1] += comp10 + comp11;
				}
			}
		}
		if (gravity_flag) {
			acceleration[1] -= g;
		}
		return acceleration;
	}
	if (i_cell == i_end) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == 2) { // We need special treatment for left (i==0) and right (i==i_end) periodic BC
					int cell_size;
					try {
						cell_size = cells.at(0).at(j_cell + j - 1).size(); //Cheaky way of dealing with j==0 and j==j_end, in this cells we'll try to call non-existent cell
					}
					catch (const std::out_of_range& e) {
						continue;
					}
					if (cell_size == 0) continue;
					for (int k = 0; k < cell_size; k++) {
						int j_part = cells[0][j_cell + j - 1][k];
						std::vector<double> ij_pos = { particles[id_part].get_position()[0] - particles[j_part].get_position()[0] - length_of_area, particles[id_part].get_position()[1] - particles[j_part].get_position()[1] };
						double distance = sqrt(ij_pos[0] * ij_pos[0] + ij_pos[1] * ij_pos[1]);
						if (particles[j_part].get_type() == 'S') {
							double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
							acceleration[0] += force_ij * ij_pos[0] / distance;
							acceleration[1] += force_ij * ij_pos[1] / distance;
							continue;
						}
						std::vector<std::vector<double>> j_pres_tensor = pres_tensors[j_part];
						std::vector <double> ij_grad_weigh = { 0, 0 };
						ij_grad_weigh[0] = ij_pos[0] * ddr_weight_fun(distance, b) / distance;
						ij_grad_weigh[1] = ij_pos[1] * ddr_weight_fun(distance, b) / distance;
						double comp00 = (own_pres_tensor[0][0] / (own_density * own_density) + j_pres_tensor[0][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
						double comp01 = (own_pres_tensor[0][1] / (own_density * own_density) + j_pres_tensor[0][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
						double comp10 = (own_pres_tensor[1][0] / (own_density * own_density) + j_pres_tensor[1][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
						double comp11 = (own_pres_tensor[1][1] / (own_density * own_density) + j_pres_tensor[1][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
						acceleration[0] += comp00 + comp01;
						acceleration[1] += comp10 + comp11;
					}
					continue;
				}
				int cell_size;
				try {
					cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
				}
				catch (const std::out_of_range& e) {
					continue;
				}
				if (cell_size == 0) continue;
				for (int k = 0; k < cell_size; k++) {
					int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
					if (id_part == j_part) continue;
					if (particles[j_part].get_type() == 'S'){
						double distance = calculate_distance(particles[id_part], particles[j_part]);
						double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
						acceleration[0] += force_ij * (particles[id_part].get_position()[0] - particles[j_part].get_position()[0]) / distance;
						acceleration[1] += force_ij * (particles[id_part].get_position()[1] - particles[j_part].get_position()[1]) / distance;
						continue;
					}
					std::vector<std::vector<double>> j_pres_tensor = pres_tensors[j_part];
					std::vector<double> ij_grad_weigh = calc_grad_weight_fun(particles, id_part, j_part);
					double comp00 = (own_pres_tensor[0][0] / (own_density * own_density) + j_pres_tensor[0][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
					double comp01 = (own_pres_tensor[0][1] / (own_density * own_density) + j_pres_tensor[0][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
					double comp10 = (own_pres_tensor[1][0] / (own_density * own_density) + j_pres_tensor[1][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
					double comp11 = (own_pres_tensor[1][1] / (own_density * own_density) + j_pres_tensor[1][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
					acceleration[0] += comp00 + comp01;
					acceleration[1] += comp10 + comp11;
				}
			}
		}
		if (gravity_flag) {
			acceleration[1] -= g;
		}
		return acceleration;
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			int cell_size;
			try {
				cell_size = cells.at(i_cell + i - 1).at(j_cell + j - 1).size();
			}
			catch (const std::out_of_range& e) {
				continue;
			}
			if (cell_size == 0) continue;
			for (int k = 0; k < cell_size; k++) {
				int j_part = cells[i_cell + i - 1][j_cell + j - 1][k];
				if (id_part == j_part) continue;
				if (particles[j_part].get_type() == 'S') {
					double distance = calculate_distance(particles[id_part], particles[j_part]);
					double force_ij = force_by_Lennard_Jones(distance, Ene_bond, a_bond);
					acceleration[0] += force_ij * (particles[id_part].get_position()[0] - particles[j_part].get_position()[0]) / distance;
					acceleration[1] += force_ij * (particles[id_part].get_position()[1] - particles[j_part].get_position()[1]) / distance;
					continue;
				}
				std::vector<std::vector<double>> j_pres_tensor = pres_tensors[j_part];
				std::vector<double> ij_grad_weigh = calc_grad_weight_fun(particles, id_part, j_part);
				double comp00 = (own_pres_tensor[0][0] / (own_density * own_density) + j_pres_tensor[0][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
				double comp01 = (own_pres_tensor[0][1] / (own_density * own_density) + j_pres_tensor[0][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
				double comp10 = (own_pres_tensor[1][0] / (own_density * own_density) + j_pres_tensor[1][0] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[0];
				double comp11 = (own_pres_tensor[1][1] / (own_density * own_density) + j_pres_tensor[1][1] / (densities[j_part] * densities[j_part])) * ij_grad_weigh[1];
				acceleration[0] += comp00 + comp01;
				acceleration[1] += comp10 + comp11;
			}
		}
	}
	if (gravity_flag) {
		acceleration[1] -= g;
	}
	return acceleration;
}
