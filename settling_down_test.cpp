#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "Particle.h"
#include "free_functions.h"

// for now i'll define some global variables, in "free_function.h", mb should be rewritten
const double h = 0.1; // timestep
const double N_timesteps = 10;
const double mass = 1.0; //wtf, i declare this const two times, IF U WANT TO CHANGE BE SURE TO MAKE IT IN 2 PLACES, "free_function.h" AS WELL


int settling_down_test() {
	bool is_gravity = true; // we apply gravity by axe y

	// we'll create 48 border particles
	std::vector<std::vector<double>> wall_prt_coord(48);
	for (int i = 0; i < 9; i++) {
		wall_prt_coord[i] = { -15.0, 9.0 - i };
	}
	for (int i = 0; i < 30; i++) {
		wall_prt_coord[i + 9] = { -15.0 + i, 0.0 };
	}
	for (int i = 0; i < 9; i++) {
		wall_prt_coord[i + 39] = { 15.0, 1.0 + i };
	}
	std::vector<Particle> vector_of_particles;
	for (int i = 0; i < wall_prt_coord.size(); i++) { // warning c4018, smth with signed and unsigned types, mb fix that later
		vector_of_particles.push_back(Particle());
		vector_of_particles[i].set_position(wall_prt_coord[i][0], wall_prt_coord[i][1]);
	}

	// now let's create 100 free particles and initialize densities
	int n = 100;
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			vector_of_particles.push_back(Particle());
			vector_of_particles[48 + 10 * i + j].set_position(-4.5 + j, 10.5 - i);
		}
	}
	std::vector<double> densities(n + 48);

	// Now, we are ready for an initial step
	for (int i = 0; i < n + 48; i++) {
		densities[i] = calculate_density(vector_of_particles, i, n + 48); // rewrite to calculate one density, not a vector
	}
	std::vector<std::vector<double>> vec_of_accs(n);
	for (int i = 0; i < n; i++) {
		vec_of_accs[i] = calculate_acceleration(vector_of_particles, densities, i + 48, n + 48, is_gravity);
	}
	std::vector < std::vector<double>> vec_of_hfvelosities(n);
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_velosity = vector_of_particles[48 + i].get_velosity();
		vector_of_particles[48 + i].set_velosity(prev_velosity[0] + h * vec_of_accs[i][0] / 2.0, prev_velosity[1] + h * vec_of_accs[i][1] / 2.0 ); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
	}
	for (int i = 0; i < n; i++) {
		std::vector<double> prev_position = vector_of_particles[48 + i].get_position();
		std::vector<double> cur_speed = vector_of_particles[48 + i].get_velosity();
		vector_of_particles[48 + i].set_position(prev_position[0] + cur_speed[0] * h, prev_position[1] + cur_speed[1] * h);
	}

	// now we are ready for main integration process
	for (int i = 0; i < N_timesteps; i++) {
		for (int j = 0; j < n + 48; j++) {
			densities[j] = calculate_density(vector_of_particles, j, n + 48);
		}
		for (int j = 0; j < n; j++) {
			vec_of_accs[j] = calculate_acceleration(vector_of_particles, densities, j + 48, n + 48, is_gravity);
		}
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_velosity = vector_of_particles[48 + j].get_velosity();
			vector_of_particles[48 + j].set_velosity(prev_velosity[0] + h * vec_of_accs[j][0], prev_velosity[1] + h * vec_of_accs[j][1]); //It's euler hf-step integration, not a simplex, hoewer we preserve accumalated h^2 
		}
		for (int j = 0; j < n; j++) {
			std::vector<double> prev_position = vector_of_particles[48 + j].get_position();
			std::vector<double> cur_speed = vector_of_particles[48 + j].get_velosity();
			vector_of_particles[48 + j].set_position(prev_position[0] + h * cur_speed[0], prev_position[1] + h * cur_speed[1]);
		}
	}

	std::cout << vector_of_particles[48].get_position()[0] << std::endl;
	system("pause");
	return 0;
}