#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "Particle.h"

//const double mass = 1.0;


//double weight_function(double r, double b) {
//	if (r >= b) {
//		return 0;
//	}
//	else {
//		return 5 / (M_PI * b * b) * (1 + 3 * r / b) * pow(1 - r / b, 3);
//	}
//}
//
//double calculate_distance(Particle& p1, Particle& p2) { //Probably should be rewritten without calling get_position() for every between-two-particle distance.
//	std::vector<double> r1 = p1.get_position();
//	std::vector<double> r2 = p2.get_position();
//	return sqrt((r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]));
//}
//
//std::vector<double> calculate_density(std::vector<Particle>& particles, int num) {
//	std::vector<double> vec_of_densties(num);
//	for (int i = 0; i < num; i++) {
//		double tmp = 0;
//		for (int j = 0; j < num; j++) {
//			tmp += weight_function(calculate_distance(particles[i], particles[j]), 2);
//		}
//		vec_of_densties[i] = tmp * mass;
//	}
//
//	return vec_of_densties;
//}
//
//int simple_tests() {
//	// Basic tests
//	Particle my_particle;
//	my_particle.set_velosity(3.2, 4.7);
//	my_particle.set_position(-1, -3 / 2.0);
//	std::cout << my_particle.get_velosity()[1] << std::endl;
//	std::cout << my_particle.get_position()[0] << std::endl;
//	my_particle.print_position();
//	std::cout << weight_function(1, 2) << std::endl;
//	
//
//	//small test on 4 particles, placed on square with
//	int n = 4;
//	std::vector<Particle> vector_of_particles(n);
//	vector_of_particles[0].set_position(1, 1);
//	vector_of_particles[0].set_velosity(-.1, 0);
//	vector_of_particles[1].set_position(-1, 1);
//	vector_of_particles[1].set_velosity(0, -.1);
//	vector_of_particles[2].set_position(-1, -1);
//	vector_of_particles[2].set_velosity(.1, 0);
//	vector_of_particles[3].set_position(1, -1);
//	vector_of_particles[3].set_velosity(0, .1); // in these 8 lines we define 4 particles positions and velosities
//
//	std::vector<double> densities = calculate_density(vector_of_particles, n);
//	for (int i = 0; i < n; i++) {
//		std::cout << densities[i] << ", ";
//	}
//	std::cout << std::endl;
//
//	system("pause");
//	return 0;
//};