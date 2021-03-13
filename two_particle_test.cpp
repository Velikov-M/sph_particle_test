#include "Particle.h"
#include "free_functions.h"
#include "xy_work.h"
#include <iostream>

const double h = 0.001; // timestep
const double N_timesteps = 500;

void two_particle_test() {

	std::vector<Particle> vector_of_particles;
	int n = 3;
	bool is_gravity = false;
	bool is_viscosity = false;
	bool is_friction = true;

	// Let's  create two particles with coord (0, 0) and (1, 0). Initial velosities = 0
	vector_of_particles.push_back(Particle());
	vector_of_particles[0].set_position(-2, 0);
	vector_of_particles.push_back(Particle());
	vector_of_particles[1].set_position(3, 1);
	vector_of_particles.push_back(Particle());
	vector_of_particles[2].set_position(0, -1);

	create_XY(vector_of_particles, n);

	std::vector<double> densities(n);
	densities[0] = calculate_density(vector_of_particles, 0, n);
	densities[1] = calculate_density(vector_of_particles, 1, n);
	densities[2] = calculate_density(vector_of_particles, 2, n);

	std::cout << densities[0] << std::endl;
	std::cout << densities[1] << std::endl;

	std::vector<std::vector<double>> vec_of_accs(n);
	std::vector<std::vector<std::vector<double>>> pres_tensors(n);

	vec_of_accs[0] = calculate_acceleration(vector_of_particles, pres_tensors, densities, 0, n, is_gravity);
	vec_of_accs[1] = calculate_acceleration(vector_of_particles, pres_tensors, densities, 1, n, is_gravity);
	vec_of_accs[2] = calculate_acceleration(vector_of_particles, pres_tensors, densities, 2, n, is_gravity);

	for (int i = 0; i < n; i++) {
		std::vector<double> prev_velosity = vector_of_particles[i].get_velosity();
		vector_of_particles[i].set_velosity(prev_velosity[0] + h * vec_of_accs[i][0] / 2.0, prev_velosity[1] + h * vec_of_accs[i][1] / 2.0); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
		if (is_friction) {
			apply_friction(vector_of_particles, i);
		}
	}
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_position = vector_of_particles[i].get_position();
		std::vector<double> cur_speed = vector_of_particles[i].get_velosity();
		vector_of_particles[i].set_position(prev_position[0] + cur_speed[0] * h, prev_position[1] + cur_speed[1] * h);
	}

	add_ts_XY(vector_of_particles, n);

	for (int i = 0; i < N_timesteps; i++) {
		for (int j = 0; j < n; j++) {
			densities[j] = calculate_density(vector_of_particles, j, n);
		}
		for (int j = 0; j < n; j++) {
			vec_of_accs[j] = calculate_acceleration(vector_of_particles, pres_tensors, densities, j, n, is_gravity);
		}
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_velosity = vector_of_particles[j].get_velosity();
			vector_of_particles[j].set_velosity(prev_velosity[0] + h * vec_of_accs[j][0], prev_velosity[1] + h * vec_of_accs[j][1]);
			if (is_friction) {
				apply_friction(vector_of_particles, j);
			}
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