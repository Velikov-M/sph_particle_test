#include "Particle.h"
#include <iostream>


Particle::Particle() {
	position = { 0, 0 };
	velosity = { 0, 0 };
}
void Particle::set_position(double x1, double x2) {
	position[0] = x1;
	position[1] = x2;
}

void Particle::set_velosity(double v1, double v2) {
	velosity[0] = v1;
	velosity[1] = v2;
}

void Particle::print_position() {
	std::cout << "(" << position[0] << ", " << position[1] << ")" << std::endl;
};

void Particle::print_velosity() {
	std::cout << "(" << velosity[0] << ", " << velosity[1] << ")" << std::endl;
};