#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <cmath>
#include <ctime>
#include <algorithm>
#include "Particle.h"
#include "free_functions.h"
#include "Global_consts.h"
#include "Cell_free_functions.h"
#include "xy_work.h"

// All gobal variables contains in Global_consts.h
const double h = 0.0005; // timestep
const int N_timesteps = 1000;

void throw_a_stone() {

	const bool is_gravity = true; // we apply gravity by axe y
	const bool is_viscosity = true; //logic flag of viscosity
	const bool is_friction = false;

	int counter = 0;
	std::vector<Particle> vector_of_particles;

	for (int i = 0; i < 120; i++) { // floor of 60 prtcles
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_type('W');
		vector_of_particles[counter].set_position(0 + i / 4., 0.0);
		counter++;
	}

	int n_free_prt = read_initial_position(vector_of_particles, counter);

	int n_stone = 40;
	double a = 4.0;
	double b = 1.0;

	std::vector<double> stone_pos = { 5, 11, M_PI / 6 }; //initialising stone coordinates and it's derivatives (stone modeled as solid in 2D, therefore 3 kinematic coordinates (2 translational, 1 rotational))
	std::vector<double> stone_vel = { 500, 0, 0 };
	std::vector<std::vector<double>> initial_r_stone(n_stone);
	for (int i = 0; i < n_stone; i++) {
		vector_of_particles.push_back(Particle());
		double x_elipse = stone_pos[0] + a * cos(2.0 * M_PI * i / n_stone);
		double y_elipse = stone_pos[1] + b * sin(2.0 * M_PI * i / n_stone);
		initial_r_stone[i] = { x_elipse - stone_pos[0], y_elipse - stone_pos[1] };
		double x_new = stone_pos[0] + initial_r_stone[i][0] * cos(stone_pos[2]) + initial_r_stone[i][1] * sin(-stone_pos[2]);
		double y_new = stone_pos[1] - initial_r_stone[i][0] * sin(-stone_pos[2]) + initial_r_stone[i][1] * cos(stone_pos[2]);
		vector_of_particles[counter + n_free_prt + i].set_position(x_new, y_new);
		vector_of_particles[counter + n_free_prt + i].set_type('S');
		
	}

	const double x_0 = 0.0;
	const double x_end = 30;
	const double y_0 = 0.0;
	const double y_end = 30.0;

	std::vector<double> densities(counter + n_free_prt);
	//std::vector<double> test_densities(n + counter);
	create_XY(vector_of_particles, counter + n_free_prt + n_stone);
	std::vector<std::vector<std::vector<double>>> pres_tensors(counter + n_free_prt);
	std::vector<std::vector<std::vector<int>>> cells = create_cells(vector_of_particles, x_0, x_end, y_0, y_end, counter + n_free_prt + n_stone); // coordinates of mesh are custom, THINK ABOUT IT!
	const int n_horizontal = cells.size();

	for (int i = 0; i < counter + n_free_prt; i++) {
		std::vector<int> cell_of_prt = vector_of_particles[i].get_cell_id();
		densities[i] = cell_calculate_density(vector_of_particles, cells, i, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1);
		//std::cout << densities[i] - test_densities[i] << ", id = " << i << std::endl;
	}

	for (int i = 0; i < counter + n_free_prt; i++) {
		std::vector<int> cell_of_prt = vector_of_particles[i].get_cell_id();
		if (is_viscosity) {
			//std::cout << "id part = " << i << std::endl;
			//test_pt[i] = calc_newton_pres_tensor(vector_of_particles, densities, i, n + counter);
			pres_tensors[i] = cell_calc_newton_pt(vector_of_particles, densities, cells, i, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1);
			//std::cout << test_pt[i][0][0] - pres_tensors[i][0][0] << " " << test_pt[i][0][1] - pres_tensors[i][0][1] << std::endl;
			//std::cout << test_pt[i][1][0] - pres_tensors[i][1][0] << " " << test_pt[i][1][1] - pres_tensors[i][1][1] << std::endl;
		}
		else {
			pres_tensors[i] = calc_eulier_pres_tensor(densities, i);
		}
	}

	std::vector<std::vector<double>> vec_of_accs(n_free_prt); //accs of all free particles (fluid particles)
	std::vector<double> stone_accs; // we also need stone accs, first two translational, the third on is angular , {x0'', x1'', \phi''}
	for (int i = 0; i < n_free_prt; i++) {
		std::vector<int> cell_of_prt = vector_of_particles[i + counter].get_cell_id();
		vec_of_accs[i] = cell_calc_acc(vector_of_particles, pres_tensors, densities, cells, i + counter, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1, is_gravity);
	}
	stone_accs = cell_calc_stone_accs(vector_of_particles, cells, stone_pos, counter + n_free_prt, n_stone, n_horizontal - 1, is_gravity);

	std::vector < std::vector<double>> vec_of_hfvelosities(n_free_prt);
	for (int i = 0; i < n_free_prt; i++) {
		std::vector<double> prev_velosity = vector_of_particles[counter + i].get_velosity();
		vector_of_particles[counter + i].set_velosity(prev_velosity[0] + h * vec_of_accs[i][0] / 2.0, prev_velosity[1] + h * vec_of_accs[i][1] / 2.0); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
		if (is_friction) {
			apply_friction(vector_of_particles, counter + i);
		}
	}

	for (int i = 0; i < 3; i++) {
		stone_vel[i] = stone_vel[i] + h * stone_accs[i] / 2.0;
	}
	for (int i = 0; i < n_stone; i++) {
		std::vector <double> r_prt = { vector_of_particles[counter + n_free_prt + i].get_position()[0] - stone_pos[0], vector_of_particles[counter + n_free_prt + i].get_position()[1] - stone_pos[1] };
		double v0_from_omega = -stone_vel[2] * r_prt[1];
		double v1_from_omega = stone_vel[2] * r_prt[0];
		vector_of_particles[counter + n_free_prt + i].set_velosity(stone_vel[0] + v0_from_omega, stone_vel[1] + v1_from_omega); //aplying eulier formula for kinematics of solid
	}

	for (int i = 0; i < n_free_prt; i++) {
		std::vector<double> prev_position = vector_of_particles[counter + i].get_position();
		std::vector<double> cur_speed = vector_of_particles[counter + i].get_velosity();
		std::vector<int> cell_of_prt = vector_of_particles[counter + i].get_cell_id();
		cells[cell_of_prt[0]][cell_of_prt[1]].erase(std::remove(cells[cell_of_prt[0]][cell_of_prt[1]].begin(), cells[cell_of_prt[0]][cell_of_prt[1]].end(), counter + i), cells[cell_of_prt[0]][cell_of_prt[1]].end());
		//		std::vector<int>* my_cell = &cells[cell_of_prt[0]][cell_of_prt[1]];
		//		*my_cell->erase(std::remove(*my_cell->begin(), *my_cell->end(), i), *my_cell->end());
		double y_new = prev_position[1] + cur_speed[1] * h;
		if (y_new < y_0) {
			y_new = y_0 + b / 3.0;
		}
		if (y_new > y_end) {
			y_new = y_end - b / 3.0;
		}
		double x_new = prev_position[0] + cur_speed[0] * h;
		if (x_new < 0) {
			vector_of_particles[counter + i].set_position(x_new + x_end, y_new);
			continue;
		}
		if (x_new >= x_end) {
			vector_of_particles[counter + i].set_position(x_new - x_end, y_new);
			continue;
		}
		vector_of_particles[counter + i].set_position(x_new, y_new);
	}
	
	//now let's integrate positions of stone
	double x_new = stone_pos[0] + h * stone_vel[0];
	if (x_new < 0) x_new = x_new + x_end;
	if (x_new >= 30) x_new = x_new - x_end;
	stone_pos[0] = x_new;

	stone_pos[1] = stone_pos[1] + h * stone_vel[1];
	stone_pos[2] = stone_pos[2] + h * stone_vel[2];
	
	for (int i = 0; i < n_stone; i++) {
		std::vector<int> cell_of_prt = vector_of_particles[counter + n_free_prt + i].get_cell_id();
		cells[cell_of_prt[0]][cell_of_prt[1]].erase(std::remove(cells[cell_of_prt[0]][cell_of_prt[1]].begin(), cells[cell_of_prt[0]][cell_of_prt[1]].end(), counter + n_free_prt + i), cells[cell_of_prt[0]][cell_of_prt[1]].end());
		double x_new = stone_pos[0] + initial_r_stone[i][0] * cos(stone_pos[2]) + initial_r_stone[i][1] * sin(-stone_pos[2]);
		double y_new = stone_pos[1] -initial_r_stone[i][0] * sin(-stone_pos[2]) + initial_r_stone[i][1] * cos(stone_pos[2]);
		if (x_new < 0) x_new = x_new + x_end;
		if (x_new >= 30) x_new = x_new - x_end;
		vector_of_particles[counter + n_free_prt + i].set_position(x_new, y_new);
	}

	add_ts_XY(vector_of_particles, counter + n_free_prt + n_stone);
	refresh_cells(vector_of_particles, cells, x_0, y_0, counter, counter + n_free_prt + n_stone);

	for (int i = 0; i < N_timesteps; i++) {
		for (int j = 0; j < counter + n_free_prt; j++) {
			std::vector<int> cell_of_prt = vector_of_particles[j].get_cell_id();
			densities[j] = cell_calculate_density(vector_of_particles, cells, j, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1);
			//std::cout << densities[j] - test_densities[j] << ", id = " << j << std::endl;
		}

		for (int j = 0; j < counter + n_free_prt; j++) {
			std::vector<int> cell_of_prt = vector_of_particles[j].get_cell_id();
			if (is_viscosity) {
				//std::cout << "id part = " << j << std::endl;
				//test_pt[j] = calc_newton_pres_tensor(vector_of_particles, densities, j, n + counter);
				pres_tensors[j] = cell_calc_newton_pt(vector_of_particles, densities, cells, j, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1);
				//std::cout << test_pt[j][0][0] - pres_tensors[j][0][0] << " " << test_pt[j][0][1] - pres_tensors[j][0][1] << std::endl;
				//std::cout << test_pt[j][1][0] - pres_tensors[j][1][0] << " " << test_pt[j][1][1] - pres_tensors[j][1][1] << std::endl;
			}
			else {
				pres_tensors[j] = calc_eulier_pres_tensor(densities, j);
			}
		}

		std::vector<std::vector<double>> vec_of_accs(n_free_prt); //accs of all free particles (fluid particles)
		std::vector<double> stone_accs; // we also need stone accs, first two translational, the third on is angular , {x0'', x1'', \phi''}
		for (int j = 0; j < n_free_prt; j++) {
			std::vector<int> cell_of_prt = vector_of_particles[j + counter].get_cell_id();
			vec_of_accs[j] = cell_calc_acc(vector_of_particles, pres_tensors, densities, cells, j + counter, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1, is_gravity);
		}
		stone_accs = cell_calc_stone_accs(vector_of_particles, cells, stone_pos, counter + n_free_prt, n_stone, n_horizontal - 1, is_gravity);
		//std::vector<double> test_stone_accs;
		//test_stone_accs = calc_stone_accs(vector_of_particles, stone_pos, counter + n_free_prt, n_stone, is_gravity);
		//std::cout << stone_accs[0] - test_stone_accs[0] << " " << stone_accs[1] - test_stone_accs[1] << " " << stone_accs[2] - test_stone_accs[2] << std::endl;

		std::vector < std::vector<double>> vec_of_hfvelosities(n_free_prt);
		for (int j = 0; j < n_free_prt; j++) {
			std::vector<double> prev_velosity = vector_of_particles[counter + j].get_velosity();
			vector_of_particles[counter + j].set_velosity(prev_velosity[0] + h * vec_of_accs[j][0] / 2.0, prev_velosity[1] + h * vec_of_accs[j][1] / 2.0); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
			if (is_friction) {
				apply_friction(vector_of_particles, counter + j);
			}
		}

		for (int j = 0; j < 3; j++) {
			stone_vel[j] = stone_vel[j] + h * stone_accs[j] / 2.0;
		}

		for (int j = 0; j < n_stone; j++) {
			std::vector <double> r_prt = { vector_of_particles[counter + n_free_prt + j].get_position()[0] - stone_pos[0], vector_of_particles[counter + n_free_prt + j].get_position()[1] - stone_pos[1] };
			double v0_from_omega = -stone_vel[2] * r_prt[1];
			double v1_from_omega = stone_vel[2] * r_prt[0];
			vector_of_particles[counter + n_free_prt + j].set_velosity(stone_vel[0] + v0_from_omega, stone_vel[1] + v1_from_omega); //aplying eulier formula for kinematics of solid
		}

		for (int j = 0; j < n_free_prt; j++) {
			std::vector<double> prev_position = vector_of_particles[counter + j].get_position();
			std::vector<double> cur_speed = vector_of_particles[counter + j].get_velosity();
			std::vector<int> cell_of_prt = vector_of_particles[counter + j].get_cell_id();
			cells[cell_of_prt[0]][cell_of_prt[1]].erase(std::remove(cells[cell_of_prt[0]][cell_of_prt[1]].begin(), cells[cell_of_prt[0]][cell_of_prt[1]].end(), counter + j), cells[cell_of_prt[0]][cell_of_prt[1]].end());
			//		std::vector<int>* my_cell = &cells[cell_of_prt[0]][cell_of_prt[1]];
			//		*my_cell->erase(std::remove(*my_cell->begin(), *my_cell->end(), j), *my_cell->end());
			double y_new = prev_position[1] + cur_speed[1] * h;
			if (y_new < y_0) {
				y_new = y_0 + b / 3.0;
			}
			if (y_new > y_end) {
				y_new = y_end - b / 3.0;
			}
			double x_new = prev_position[0] + cur_speed[0] * h;
			if (x_new < 0) {
				vector_of_particles[counter + j].set_position(x_new + x_end, y_new);
				continue;
			}
			if (x_new >= x_end) {
				vector_of_particles[counter + j].set_position(x_new - x_end, y_new);
				continue;
			}
			vector_of_particles[counter + j].set_position(x_new, y_new);
		}

		//now let's integrate positions of stone
		double x_new = stone_pos[0] + h * stone_vel[0];
		if (x_new < 0) x_new = x_new + x_end;
		if (x_new >= 30) x_new = x_new - x_end;
		stone_pos[0] = x_new;

		stone_pos[1] = stone_pos[1] + h * stone_vel[1];
		stone_pos[2] = stone_pos[2] + h * stone_vel[2];

		for (int j = 0; j < n_stone; j++) {
			std::vector<int> cell_of_prt = vector_of_particles[counter + n_free_prt + j].get_cell_id();
			cells[cell_of_prt[0]][cell_of_prt[1]].erase(std::remove(cells[cell_of_prt[0]][cell_of_prt[1]].begin(), cells[cell_of_prt[0]][cell_of_prt[1]].end(), counter + n_free_prt + j), cells[cell_of_prt[0]][cell_of_prt[1]].end());
			double x_new = stone_pos[0] + initial_r_stone[j][0] * cos(stone_pos[2]) + initial_r_stone[j][1] * sin(-stone_pos[2]);
			double y_new = stone_pos[1] - initial_r_stone[j][0] * sin(-stone_pos[2]) + initial_r_stone[j][1] * cos(stone_pos[2]);
			if (x_new < 0) x_new = x_new + x_end;
			if (x_new >= 30) x_new = x_new - x_end;
			vector_of_particles[counter + n_free_prt + j].set_position(x_new, y_new);
		}

		add_ts_XY(vector_of_particles, counter + n_free_prt + n_stone);
		refresh_cells(vector_of_particles, cells, x_0, y_0, counter, counter + n_free_prt + n_stone);

	}
}