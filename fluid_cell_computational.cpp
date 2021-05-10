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
#include <omp.h>
// In this test, we will apply periodic BC on left and right corner on function level

// All gobal variables contains in Global_consts.h
const double h = 0.0005; // timestep
const int N_timesteps = 1000;



void set_down_with_cells() {
	omp_set_num_threads(16);
	const bool is_gravity = true; // we apply gravity by axe y
	const bool is_viscosity = true; //logic flag of viscosity
	const bool is_friction = false;

	// we'll create counter border particles
	int counter = 0;
	std::vector<Particle> vector_of_particles;
	for (int i = 0; i < 120; i++) { // floor of 60 prtcles
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_type('W');
		vector_of_particles[counter].set_position(0 + i / 2., 0.1);
		counter++;
	}
	int n = 1600;
	for (int i = 0; i < 40; i++) { //square of 400 particles 20x20, left lower corner in (0, 0.5)
		for (int j = 0; j < 40; j++) {
			vector_of_particles.push_back(Particle());
			vector_of_particles[counter + 40 * i + j].set_type('P');
			vector_of_particles[counter + 40 * i + j].set_position(j, 40.5 - i);
		}
	}

	// U CAN'T CHANGE X_END BEFORE REWRITING CELL_FREE_FUNCTIONS.CPP!!!! should be transfered to Global Constants?
	const double x_0 = 0.0;
	const double x_end = length_of_area;
	const double y_0 = 0.0;
	const double y_end = 60.0;


	std::vector<double> densities(n + counter);
	//std::vector<double> test_densities(n + counter);
	create_XY(vector_of_particles, n + counter);
	std::vector<std::vector<std::vector<double>>> pres_tensors(n + counter);

	//here we create cells, and fill them, step of cells equal to a domain of interst
	std::vector<std::vector<std::vector<int>>> cells = create_cells(vector_of_particles, x_0, x_end, y_0, y_end, counter + n); // coordinates of mesh are custom, THINK ABOUT IT!
	const int n_horizontal = cells.size();

	//unsigned int start_time_2 = clock();
	//for (int i = 0; i < counter + n; i++) {
	//	test_densities[i] = calculate_density(vector_of_particles, i, counter + n);
	//}
	//unsigned int end_time_2 = clock();
	//std::cout << "default computational time = " << end_time_2 - start_time_2 << std::endl;

	//unsigned int start_time = clock();
#pragma omp parallel for
	for (int i = 0; i < counter + n; i++) {
		std::vector<int> cell_of_prt = vector_of_particles[i].get_cell_id();
		densities[i] = cell_calculate_density(vector_of_particles, cells, i, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1);
		//std::cout << densities[i] - test_densities[i] << ", id = " << i << std::endl;
	}

	//unsigned int end_time = clock();
	//std::cout << "cell computational time = " << end_time - start_time << std::endl;
	//system("pause");
	// 
	// Now, we are ready for an initial step

	//std::vector<std::vector<std::vector<double>>> test_pt(n + counter);
#pragma omp parallel for
	for (int i = 0; i < counter + n; i++) {
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

	std::vector<std::vector<double>> vec_of_accs(n);
	std::vector<std::vector<double>> test_accs(n);
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		//std::cout << "id part = " << i + counter << std::endl;
		std::vector<int> cell_of_prt = vector_of_particles[i + counter].get_cell_id();
		//std::cout << "id cell = (" << cell_of_prt[0] << ", " << cell_of_prt[1] << ")" << std::endl;
		vec_of_accs[i] = cell_calc_acc(vector_of_particles, pres_tensors, densities, cells, i + counter, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1, is_gravity);
		//test_accs[i] = calculate_acceleration(vector_of_particles, pres_tensors, densities, i + counter, n + counter, 0, is_gravity);
		//std::cout << vec_of_accs[i][0] - test_accs[i][0] << " " << vec_of_accs[i][0] - test_accs[i][0] << std::endl;
	}
	std::vector < std::vector<double>> vec_of_hfvelosities(n);
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_velosity = vector_of_particles[counter + i].get_velosity();
		vector_of_particles[counter + i].set_velosity(prev_velosity[0] + h * vec_of_accs[i][0] / 2.0, prev_velosity[1] + h * vec_of_accs[i][1] / 2.0); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
		if (is_friction) {
			apply_friction(vector_of_particles, counter + i);
		}
	}
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_position = vector_of_particles[counter + i].get_position();
		std::vector<double> cur_speed = vector_of_particles[counter + i].get_velosity();
		std::vector<int> cell_of_prt = vector_of_particles[counter + i].get_cell_id();
		cells[cell_of_prt[0]][cell_of_prt[1]].erase(std::remove(cells[cell_of_prt[0]][cell_of_prt[1]].begin(), cells[cell_of_prt[0]][cell_of_prt[1]].end(), counter + i), cells[cell_of_prt[0]][cell_of_prt[1]].end());
//		std::vector<int>* my_cell = &cells[cell_of_prt[0]][cell_of_prt[1]];
//		*my_cell->erase(std::remove(*my_cell->begin(), *my_cell->end(), i), *my_cell->end());
		double x_new = prev_position[0] + cur_speed[0] * h;
		if (x_new < 0) {
			vector_of_particles[counter + i].set_position(x_new + x_end, prev_position[1] + cur_speed[1] * h);
			continue;
		}
		if (x_new >= x_end) {
			vector_of_particles[counter + i].set_position(x_new - x_end, prev_position[1] + cur_speed[1] * h);
			continue;
		}
		vector_of_particles[counter + i].set_position(prev_position[0] + cur_speed[0] * h, prev_position[1] + cur_speed[1] * h);
	}
	add_ts_XY(vector_of_particles, n + counter);
	refresh_cells(vector_of_particles, cells, x_0, y_0, counter, n + counter);

	for (int i = 0; i < N_timesteps; i++) {
	//std::cout << "step " << i << std::endl;
#pragma omp parallel for
		for (int j = 0; j < n + counter; j++) {
			std::vector<double> test_densities(counter + n);
			std::vector<int> cell_of_prt = vector_of_particles[j].get_cell_id();
			//test_densities[j] = calculate_density(vector_of_particles, j, counter + n);
			densities[j] = cell_calculate_density(vector_of_particles, cells, j, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1);
			//std::cout << test_densities[j] - densities[j] << std::endl;
		}
		if (i % 500 == 0) {
			std::cout << densities[200] << std::endl;
		}
#pragma omp parallel for
		for (int j = 0; j < n + counter; j++) {
			std::vector<int> cell_of_prt = vector_of_particles[j].get_cell_id();
			if (is_viscosity) {
				//test_pt[j] = calc_newton_pres_tensor(vector_of_particles, densities, j, n + counter);
				pres_tensors[j] = cell_calc_newton_pt(vector_of_particles, densities, cells, j, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1);
				//std::cout << "id cell = (" << cell_of_prt[0] << ", " << cell_of_prt[1] << ")" << std::endl;
				//std::cout << test_pt[j][0][0] - pres_tensors[j][0][0] << " " << test_pt[j][0][1] - pres_tensors[j][0][1] << std::endl;
				//std::cout << test_pt[j][1][0] - pres_tensors[j][1][0] << " " << test_pt[j][1][1] - pres_tensors[j][1][1] << std::endl;
			}
			else {
				pres_tensors[j] = calc_eulier_pres_tensor(densities, j);
			}
		}
#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			//if (i == 37) {
			//	std::cout << j << std::endl;
			//}
			std::vector<int> cell_of_prt = vector_of_particles[j + counter].get_cell_id();
			//std::cout << "id cell = (" << cell_of_prt[0] << ", " << cell_of_prt[1] << ")" << std::endl;
			vec_of_accs[j] = cell_calc_acc(vector_of_particles, pres_tensors, densities, cells, j + counter, cell_of_prt[0], cell_of_prt[1], n_horizontal - 1, is_gravity);
			//test_accs[j] = calculate_acceleration(vector_of_particles, pres_tensors, densities, j + counter, n + counter, 0, is_gravity);
			//std::cout << vec_of_accs[j][0] - test_accs[j][0] << " " << vec_of_accs[j][0] - test_accs[j][0] << std::endl;

		}
#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_velosity = vector_of_particles[counter + j].get_velosity();
			vector_of_particles[counter + j].set_velosity(prev_velosity[0] + h * vec_of_accs[j][0], prev_velosity[1] + h * vec_of_accs[j][1]);
			if (is_friction) {
				apply_friction(vector_of_particles, counter + j);
			}
		}
		for (int j = 0; j < n; j++) {
			std::vector<int> cell_of_prt = vector_of_particles[counter + j].get_cell_id();
			cells[cell_of_prt[0]][cell_of_prt[1]].erase(std::remove(cells[cell_of_prt[0]][cell_of_prt[1]].begin(), cells[cell_of_prt[0]][cell_of_prt[1]].end(),counter + j), cells[cell_of_prt[0]][cell_of_prt[1]].end());
			std::vector<double> prev_position = vector_of_particles[counter + j].get_position();
			std::vector<double> cur_speed = vector_of_particles[counter + j].get_velosity();
			double y_new = prev_position[1] + cur_speed[1] * h;
			if (y_new < y_0) {
				y_new = y_0 + b / 3.0;
			}
			if (y_new > y_end) {
				y_new = y_end - b / 3.0;
			}
			double x_new = prev_position[0] + cur_speed[0] * h;
			if (x_new <= 0) {
				vector_of_particles[counter + j].set_position(x_new + x_end, y_new);
				continue;
			}
			if (x_new > x_end) {
				vector_of_particles[counter + j].set_position(x_new - x_end, y_new);
				continue;
			}
			vector_of_particles[counter + j].set_position(prev_position[0] + cur_speed[0] * h, y_new);
		}
		add_ts_XY(vector_of_particles, n + counter);
		refresh_cells(vector_of_particles, cells, x_0, y_0, counter, n + counter);
	}
	system("pause");
	
}