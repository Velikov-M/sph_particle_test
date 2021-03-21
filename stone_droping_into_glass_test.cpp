#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "Particle.h"
#include "free_functions.h"
#include "xy_work.h"

// for now i'll define most of global variables, in "free_function.h", mb should be rewritten
const double h = 0.0005; // timestep
const double N_timesteps = 1000;


void droping_test() {
	bool is_gravity = true; // we apply gravity by axe y
	bool is_viscosity = true; //logic flag of viscosity
	bool is_friction = false; // old global dissipation, can be used only for tests

	int counter = 0;
	std::vector<Particle> vector_of_particles;

	// In the begining we're creating glass 
	for (int i = 0; i < 18; i++) {
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_position(-15.0, 9.0 - i / 2.);
		counter++;
	}
	for (int i = 0; i < 61; i++) {
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_position(-15.0 + i / 2., 0.0);
		counter++;
	}
	for (int i = 0; i < 18; i++) {
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_position(15.0, .5 + i / 2.);
		counter++;
	}

	// now let's create 100 free particles
	int n_free_prt = read_initial_position(vector_of_particles, counter); // returns the number of free prt, had been read from xy-file, vector_of_particles has been modified!!!
	// now we have to create elipse-shaped stone, a and b - are parametrs of elipse, h - is initial height
	int n_stone = 40;
	double a = 4.0;
	double b = 1.0;
	std::vector<double> stone_pos = { 0, 11, 0}; //initialising stone coordinates and it's derivatives (stone modeled as solid in 2D, therefore 3 kinematic coordinates (2 translational, 1 rotational))
	std::vector<double> stone_vel = { 0, 0, 0};
	std::vector<std::vector<double>> initial_r_stone(n_stone);
	for (int i = 0; i < n_stone; i++) {
		vector_of_particles.push_back(Particle());
		double x_elipse = stone_pos[0] + a * cos(2.0 * M_PI * i / n_stone);
		double y_elipse = stone_pos[1] + b * sin(2.0 * M_PI * i / n_stone);
		vector_of_particles[counter + n_free_prt + i].set_position(x_elipse, y_elipse);
		initial_r_stone[i] = { x_elipse - stone_pos[0], y_elipse - stone_pos[1] };
	}

	std::vector<double> densities(n_free_prt + counter);
	create_XY(vector_of_particles, n_free_prt + n_stone + counter);
	std::vector<std::vector<std::vector<double>>> pres_tensors(n_free_prt + counter);

	//Initial step
	for (int i = 0; i < n_free_prt + counter; i++) {
		densities[i] = calculate_density(vector_of_particles, i, n_free_prt + counter);
	}
	for (int i = 0; i < n_free_prt + counter; i++) {
		if (is_viscosity) {
			pres_tensors[i] = calc_newton_pres_tensor(vector_of_particles, densities, i, n_free_prt + counter);
		}
		else {
			pres_tensors[i] = calc_eulier_pres_tensor(densities, i, n_free_prt + counter);
		}
	}

	std::vector<std::vector<double>> vec_of_accs(n_free_prt); //accs of all free particles (fluid particles)
	std::vector<double> stone_accs; // we also need stone accs, first two translational, the third on is angular , {x0'', x1'', \phi''}

	// now we'll calculate accs in two steps, one for all free particles and one for stone (summing all forces and momentum from fluid prts to stone prts)
	for (int i = 0; i < n_free_prt; i++) {
		vec_of_accs[i] = calculate_acceleration(vector_of_particles, pres_tensors, densities, i + counter, n_free_prt + counter, n_stone, is_gravity);
	}
	stone_accs = calc_stone_accs(vector_of_particles, stone_pos, n_free_prt + counter, n_stone, is_gravity);
	// in this part we integrate velocities on initial step (not a simplex!)
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

	// now we ready to integrate positions

	for (int i = 0; i < n_free_prt; i++) {
		std::vector<double> prev_position = vector_of_particles[counter + i].get_position();
		std::vector<double> cur_speed = vector_of_particles[counter + i].get_velosity();
		vector_of_particles[counter + i].set_position(prev_position[0] + cur_speed[0] * h, prev_position[1] + cur_speed[1] * h);
	}

	for (int i = 0; i < 3; i++) {
		stone_pos[i] = stone_pos[i] + h * stone_vel[i];
	}
	for (int i = 0; i < n_stone; i++) {
		double x0_from_omega = initial_r_stone[i][0] * cos(stone_pos[2]) + initial_r_stone[i][1] * sin(-stone_pos[2]);
		double x1_from_omega = -initial_r_stone[i][0] * sin(-stone_pos[2]) + initial_r_stone[i][1] * cos(stone_pos[2]);
		vector_of_particles[counter + n_free_prt + i].set_position(stone_pos[0] + x0_from_omega, stone_pos[1] + x1_from_omega);
	}

	add_ts_XY(vector_of_particles, counter + n_free_prt + n_stone);

	// now we are ready for main integration process
	for (int i = 0; i < N_timesteps; i++) {
		for (int j = 0; j < counter + n_free_prt; j++) {
			densities[j] = calculate_density(vector_of_particles, j, counter +  n_free_prt);
		}
		for (int j = 0; j < counter + n_free_prt; j++) {
			if (is_viscosity) {
				pres_tensors[j] = calc_newton_pres_tensor(vector_of_particles, densities, j, counter + n_free_prt);
			}
			else {
				pres_tensors[j] = calc_eulier_pres_tensor(densities, j, counter + n_free_prt);
			}
		}
		if (i % 500 == 0) {
			std::cout << densities[49] << std::endl;
		}

		for (int j = 0; j < n_free_prt; j++) {
			vec_of_accs[j] = calculate_acceleration(vector_of_particles, pres_tensors, densities, j + counter, n_free_prt + counter, n_stone, is_gravity);
		}
		stone_accs = calc_stone_accs(vector_of_particles, stone_pos, n_free_prt + counter, n_stone, is_gravity);

		for (int j = 0; j < n_free_prt; j++) {
			std::vector<double> prev_velosity = vector_of_particles[j + counter].get_velosity();
			vector_of_particles[j + counter].set_velosity(prev_velosity[0] + h * vec_of_accs[j][0], prev_velosity[1] + h * vec_of_accs[j][1]);
			if (is_friction) {
				apply_friction(vector_of_particles, j + counter);
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
			std::vector<double> prev_position = vector_of_particles[j + counter].get_position();
			std::vector<double> cur_speed = vector_of_particles[j + counter].get_velosity();
			vector_of_particles[j + counter].set_position(prev_position[0] + h * cur_speed[0], prev_position[1] + h * cur_speed[1]);
		}
		for (int j = 0; j < 3; j++) {
			stone_pos[j] = stone_pos[j] + h * stone_vel[j];
		}
		for (int j = 0; j < n_stone; j++) {
			double x0_from_omega = initial_r_stone[j][0] * cos(stone_pos[2]) + initial_r_stone[j][1] * sin(-stone_pos[2]);
			double x1_from_omega = -initial_r_stone[j][0] * sin(-stone_pos[2]) + initial_r_stone[j][1] * cos(stone_pos[2]);
			vector_of_particles[counter + n_free_prt + j].set_position(stone_pos[0] + x0_from_omega, stone_pos[1] + x1_from_omega);
		}

		//std::cout << calc_norm_of_impulse(vector_of_particles, 0, n) << std::endl;
		//std::cout << calc_momentum(vector_of_particles, 0, n) << std::endl;
		add_ts_XY(vector_of_particles, counter + n_free_prt + n_stone);
	}




}