#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "Particle.h"
#include "free_functions.h"
#include "xy_work.h"

// for now i'll define some global variables, in "free_function.h", mb should be rewritten
const double h = 0.001; // timestep
const double N_timesteps = 1000;


void settling_down_test() {
	bool is_gravity = true; // we apply gravity by axe y
	bool is_viscosity = true; //logic flag of viscosity
	bool is_friction = false;

	// we'll create counter border particles
	int counter = 0;
	std::vector<Particle> vector_of_particles;
	for (int i = 0; i < 18; i++) {
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_position( -15.0, 9.0 - i / 2. );
		vector_of_particles[counter].set_type('W');
		counter++;
	}
	for (int i = 0; i < 61; i++) {
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_position(-15.0 + i / 2., 0.0 );
		vector_of_particles[counter].set_type('W');
		counter++;
	}
	for (int i = 0; i < 18; i++) {
		vector_of_particles.push_back(Particle());
		vector_of_particles[counter].set_position( 15.0, .5 + i / 2. );
		vector_of_particles[counter].set_type('W');
		counter++;
	}

	// now let's create 100 free particles and initialize densities
	int n = 100;
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			vector_of_particles.push_back(Particle());
			vector_of_particles[counter + 10 * i + j].set_position(-4.5 + j, 10.5 - i);
			vector_of_particles[counter + 10 * i + j].set_type('P');
		}
	}
	std::vector<double> densities(n + counter);
	create_XY(vector_of_particles, n + counter);
	std::vector<std::vector<std::vector<double>>> pres_tensors(n + counter);
	// Now, we are ready for an initial step
	for (int i = 0; i < n + counter; i++) {
		densities[i] = calculate_density(vector_of_particles, i, n + counter); // rewrite to calculate one density, not a vector
	}
	for (int i = 0; i < n + counter; i++) {
		if (is_viscosity) {
			pres_tensors[i] = calc_newton_pres_tensor(vector_of_particles, densities, i, n + counter);
		}
		else {
			pres_tensors[i] = calc_eulier_pres_tensor(densities, i);
		}
	}
	std::vector<std::vector<double>> vec_of_accs(n);
	for (int i = 0; i < n; i++) {
		vec_of_accs[i] = calculate_acceleration(vector_of_particles, pres_tensors, densities, i + counter, n + counter, 0, is_gravity);
	}
	std::vector < std::vector<double>> vec_of_hfvelosities(n);
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_velosity = vector_of_particles[counter + i].get_velosity();
		vector_of_particles[counter + i].set_velosity(prev_velosity[0] + h * vec_of_accs[i][0] / 2.0, prev_velosity[1] + h * vec_of_accs[i][1] / 2.0 ); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
		if (is_friction) {
			apply_friction(vector_of_particles, counter + i);
		}
	}
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_position = vector_of_particles[counter + i].get_position();
		std::vector<double> cur_speed = vector_of_particles[counter + i].get_velosity();
		vector_of_particles[counter + i].set_position(prev_position[0] + cur_speed[0] * h, prev_position[1] + cur_speed[1] * h);
	}

	add_ts_XY(vector_of_particles, n + counter);

	// now we are ready for main integration process
	for (int i = 0; i < N_timesteps; i++) {
		for (int j = 0; j < n + counter; j++) {
			densities[j] = calculate_density(vector_of_particles, j, n + counter);
		}
		if (i % 500 == 0) {
			std::cout << densities[49] << std::endl;
		}
		for (int j = 0; j < n + counter; j++) {
			if (is_viscosity) {
				pres_tensors[j] = calc_newton_pres_tensor(vector_of_particles, densities, j, n + counter);
			}
			else {
				pres_tensors[j] = calc_eulier_pres_tensor(densities, j);
			}
		}
		for (int j = 0; j < n; j++) {
			vec_of_accs[j] = calculate_acceleration(vector_of_particles, pres_tensors, densities, j + counter, n + counter, 0, is_gravity);
		}
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_velosity = vector_of_particles[counter + j].get_velosity();
			vector_of_particles[counter + j].set_velosity(prev_velosity[0] + h * vec_of_accs[j][0], prev_velosity[1] + h * vec_of_accs[j][1]);
			if (is_friction) {
				apply_friction(vector_of_particles, counter + j);
			}
		}
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_position = vector_of_particles[counter + j].get_position();
			std::vector<double> cur_speed = vector_of_particles[counter + j].get_velosity();
			vector_of_particles[counter + j].set_position(prev_position[0] + h * cur_speed[0], prev_position[1] + h * cur_speed[1]);
		}
		add_ts_XY(vector_of_particles, n + counter);
	}

	//std::cout << vector_of_particles[counter].get_position()[0] << std::endl;
}