#include "Particle.h"
#include "free_functions.h"
#include "xy_work.h"
#include <iostream>

const double h = 0.01; // timestep
const double N_timesteps = 100;

void two_particle_test() {

	std::vector<Particle> vector_of_particles;
	int n = 2;
	bool is_gravity = false;

	// Let's  create two particles with coord (0, 0) and (1, 0). Initial velosities = 0
	vector_of_particles.push_back(Particle());
	vector_of_particles[0].set_position(0, 0);
	vector_of_particles.push_back(Particle());
	vector_of_particles[1].set_position(1, 0);

	std::vector<double> densities(n);
	densities[0] = calculate_density(vector_of_particles, 0, n);
	densities[1] = calculate_density(vector_of_particles, 1, n);

	std::cout << densities[0] << std::endl;
	std::cout << densities[1] << std::endl;

	std::vector<std::vector<double>> vec_of_accs(n);

	vec_of_accs[0] = calculate_acceleration(vector_of_particles, densities, 0, n, is_gravity);


}