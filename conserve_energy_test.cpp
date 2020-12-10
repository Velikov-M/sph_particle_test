#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "Particle.h"
#include "free_functions.h"
#include "xy_work.h"

const double h = 0.001; // timestep
const double N_timesteps = 100;
const double mass = 1.0; //wtf, i declare this const two times, IF U WANT TO CHANGE BE SURE TO MAKE IT IN 2 PLACES, "free_function.h" AS WELL

void conserve_energy_test() {
	bool is_gravity = false;
	bool is_viscosity = false;
	
	std::vector<double> energies_on_steps;
	std::vector<Particle> vector_of_particles;
	int n = 100;
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			vector_of_particles.push_back(Particle());
			vector_of_particles[10 * i + j].set_position(i, j);
		}
	}
	create_XY(vector_of_particles, n);
	std::vector<double> densities(n);

	//std::cout << calc_norm_of_impulse(vector_of_particles, 0, n) << std::endl;
	//std::cout << calc_momentum(vector_of_particles, 0, n) << std::endl;
	for (int i = 0; i < n; i++) {
		densities[i] = calculate_density(vector_of_particles, i, n);
	}
	std::vector<std::vector<double>> vec_of_accs(n);
	for (int i = 0; i < n; i++) {
		vec_of_accs[i] = calculate_acceleration(vector_of_particles, densities, i, n, is_gravity, is_viscosity);
	}
	std::vector < std::vector<double>> vec_of_hfvelosities(n);
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_velosity = vector_of_particles[i].get_velosity();
		vector_of_particles[i].set_velosity(prev_velosity[0] + h * vec_of_accs[i][0] / 2.0, prev_velosity[1] + h * vec_of_accs[i][1] / 2.0); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
	}
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_position = vector_of_particles[i].get_position();
		std::vector<double> cur_speed = vector_of_particles[i].get_velosity();
		vector_of_particles[i].set_position(prev_position[0] + cur_speed[0] * h, prev_position[1] + cur_speed[1] * h);
	}
	//std::cout << calc_norm_of_impulse(vector_of_particles, 0, n) << std::endl;
	//std::cout << calc_momentum(vector_of_particles, 0, n) << std::endl;
	add_ts_XY(vector_of_particles, n);

	// now we are ready for main integration process
	for (int i = 0; i < N_timesteps; i++) {
		for (int j = 0; j < n; j++) {
			densities[j] = calculate_density(vector_of_particles, j, n);
		}
		for (int j = 0; j < n; j++) {
			vec_of_accs[j] = calculate_acceleration(vector_of_particles, densities, j, n, is_gravity, is_viscosity);
		}
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_velosity = vector_of_particles[j].get_velosity();
			vector_of_particles[j].set_velosity(prev_velosity[0] + h * vec_of_accs[j][0], prev_velosity[1] + h * vec_of_accs[j][1]);  
		}
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_position = vector_of_particles[j].get_position();
			std::vector<double> cur_speed = vector_of_particles[j].get_velosity();
			vector_of_particles[j].set_position(prev_position[0] + h * cur_speed[0], prev_position[1] + h * cur_speed[1]);
		}
		//std::cout << calc_norm_of_impulse(vector_of_particles, 0, n) << std::endl;
		//std::cout << calc_momentum(vector_of_particles, 0, n) << std::endl;
		add_ts_XY(vector_of_particles, n);
	}

}