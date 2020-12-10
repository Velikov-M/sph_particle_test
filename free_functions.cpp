#define _USE_MATH_DEFINES

#include "free_functions.h"
#include "Particle.h"
#include <math.h>
#include <cmath>

const double mass = 1; // mass of one particle, for now i am assuming be equal
const double B = 6428.6;  // next 3 global vars are consts in equation of state for liquid (Cole form, see Monaghan, SPH, 2005)
const double gamma = 7.0; //mb should think more about B and rho_0;
const double rho_0 = 0.28;
const double mu = 1; // dynamic viscosity of water by 25C temperature
const double b = 20.0; // domain of interest
const double g = 9.81; // gravity const

double weight_function(double r, double b) {
	if (r >= b) {
		return 0;
	}
	else {
		return 5 / (M_PI * b * b) * (1 + 3 * r / b) * pow(1 - r / b, 3);
	}
}

double ddr_weight_fun(double r, double b) {
	if (r >= b) {
		return 0;
	}
	else {
		return -60 * (b - r) * (b - r) * r / (M_PI * pow(b, 6));
	}
}

double calculate_distance(Particle& p1, Particle& p2) { //Probably should be rewritten without calling get_position() for every between-two-particle distance.
	std::vector<double> r1 = p1.get_position();
	std::vector<double> r2 = p2.get_position();
	return sqrt((r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]));
}

double calculate_density(std::vector<Particle>& particles, int id_part, int num) {
	double density = 0;
	for (int j = 0; j < num; j++) {
		density += weight_function(calculate_distance(particles[id_part], particles[j]), b);
	}
	density *= mass;
	return density;
}

std::vector<std::vector<double>> calc_eulier_pres_tensor(std::vector<double>& densities, int id_part, int num) {
	std::vector < std::vector<double>> pres_tensor(2);
	double pressure = B * (pow(densities[id_part] / rho_0, gamma) - 1);
	pres_tensor[0] = { -pressure, 0 };
	pres_tensor[1] = { 0, -pressure };
	return pres_tensor;
}

std::vector <double> calc_grad_weight_fun(std::vector<Particle>& particles, int id_part, int j_part) {
	std::vector <double> grad_weigh = { 0, 0 };
	double distance = calculate_distance(particles[id_part], particles[j_part]);
	grad_weigh[0] = (particles[id_part].get_position()[0] - particles[j_part].get_position()[0]) * ddr_weight_fun(distance, b) / distance;
	grad_weigh[1] = (particles[id_part].get_position()[1] - particles[j_part].get_position()[1]) * ddr_weight_fun(distance, b) / distance;
	return grad_weigh;
}

std::vector<double> calculate_acceleration(std::vector<Particle>& particles, std::vector<double>& densities, int id_part, int num, bool gravity_flag, bool viscosity_flag) {
	std::vector<double> acceleration = { 0, 0 }; // in this function, u should apply gravity
	std::vector<std::vector<double>> own_pres_tensor;
	if (viscosity_flag) {
		own_pres_tensor = calc_newton_pres_tensor(particles, densities, id_part, num);
	}
	else {
		own_pres_tensor = calc_eulier_pres_tensor(densities, id_part, num); // for now let's work with Eulier tensor
	}
	double own_density = densities[id_part];

	for (int j = 0; j < num; j++) {
		if (id_part == j) {
			continue;
		}

		std::vector<std::vector<double>> j_pres_tensor;
		if (viscosity_flag) {
			j_pres_tensor = calc_newton_pres_tensor(particles, densities, j, num);
		}
		else {
			j_pres_tensor = calc_eulier_pres_tensor(densities, j, num); // for now let's work with Eulier tensor
		}
		std::vector<double> ij_grad_weigh = calc_grad_weight_fun(particles, id_part, j);
		double comp00 = (own_pres_tensor[0][0] / (own_density * own_density) + j_pres_tensor[0][0] / (densities[j] * densities[j])) * ij_grad_weigh[0];
		double comp01 = (own_pres_tensor[0][1] / (own_density * own_density) + j_pres_tensor[0][1] / (densities[j] * densities[j])) * ij_grad_weigh[1];
		double comp10 = (own_pres_tensor[1][0] / (own_density * own_density) + j_pres_tensor[1][0] / (densities[j] * densities[j])) * ij_grad_weigh[0];
		double comp11 = (own_pres_tensor[1][1] / (own_density * own_density) + j_pres_tensor[1][1] / (densities[j] * densities[j])) * ij_grad_weigh[1];
		acceleration[0] += comp00 + comp01;
		acceleration[1] += comp10 + comp11;
	}


	if (gravity_flag) {
		acceleration[1] -= g;
	}
	return acceleration;
}

double calc_kinetic_nrg(std::vector<Particle>& particles, int start, int end) { // take as end an id of the last + 1 particle, so if u want to call energy of 0,1,2,3 prtcles, u should call (prt, 0,4)
	double sum_nrg = 0;
	for (int i = start; i < end; i++) {//check if name "i" raise a conflict
		std::vector<double> prt_velosity = particles[i].get_velosity();
		double abs_vel = sqrt(prt_velosity[0] * prt_velosity[0] + prt_velosity[1] * prt_velosity[1]);
		sum_nrg += mass * abs_vel * abs_vel / 2.0;
	}
	return sum_nrg;
}

std::vector<std::vector<double>> calc_deform_tensor(std::vector<Particle>& particles, std::vector<double>& densities, int id_part, int num) {
	std::vector < std::vector<double>> deform_tensor(2);
	deform_tensor[0] = { 0, 0 };
	deform_tensor[1] = { 0, 0 };
	std::vector<double> i_vel = particles[id_part].get_velosity();
	std::vector<double> i_pos = particles[id_part].get_position();
	for (int j = 0; j < num; j++) {
		if (j == id_part) continue;
		std::vector<double> j_vel = particles[j].get_velosity();
		std::vector<double> j_pos = particles[j].get_position();
		std::vector<double> ji_vel = { j_vel[0] - i_vel[0], j_vel[1] - i_vel[1]}; // ji - velosity, i am not so sure
		std::vector<double> ij_pos = { i_pos[0] - j_pos[0], i_pos[1] - j_pos[1]};
		double ji_distance = calculate_distance(particles[id_part], particles[j]);
		// let's calcualte add to 0,0 component of deformation tensor, u should ask question about formula (4.48), opening that delta part
		deform_tensor[0][0] += (4.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j] * ji_distance);
		deform_tensor[0][0] -= (2.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j] * ji_distance);
		// next, we'll calculate 0,1 component, as well as apply symmetry of linear deformation tensor
		deform_tensor[0][1] += mass * ji_vel[0] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j] * ji_distance);
		deform_tensor[0][1] += mass * ji_vel[1] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j] * ji_distance);
		deform_tensor[1][0] = deform_tensor[0][1];
		//and finish calculation for this ij-pair, by finding 1,1 component
		deform_tensor[1][1] += (4.0 / 3.0) * mass * ji_vel[1] * ij_pos[1] * ddr_weight_fun(ji_distance, b) / (densities[j] * ji_distance);
		deform_tensor[1][1] -= (2.0 / 3.0) * mass * ji_vel[0] * ij_pos[0] * ddr_weight_fun(ji_distance, b) / (densities[j] * ji_distance);
	}
	return deform_tensor;
}

std::vector<std::vector<double>> calc_newton_pres_tensor(std::vector<Particle>& particles, std::vector<double>& densities, int id_part, int num) {
	std::vector<std::vector<double>> pres_tensor(2);
	std::vector<std::vector<double>> eulier_part = calc_eulier_pres_tensor(densities, id_part, num);
	std::vector<std::vector<double>> deformations = calc_deform_tensor(particles, densities, id_part, num);

	pres_tensor[0] = {0, 0};
	pres_tensor[1] = {0, 0};
	
	pres_tensor[0][0] = eulier_part[0][0] + mu * deformations[0][0];
	pres_tensor[0][1] = eulier_part[0][1] + mu * deformations[0][1];
	pres_tensor[1][0] = eulier_part[1][0] + mu * deformations[1][0];
	pres_tensor[1][1] = eulier_part[1][1] + mu * deformations[1][1];

	return pres_tensor;
}

double calc_momentum(std::vector<Particle>& particles, int start, int end) {
	double sum_momentum = 0;
	for (int i = start; i < end; i++) {
		std::vector<double> prt_velosity = particles[i].get_velosity();
		std::vector<double> prt_position = particles[i].get_position();
		sum_momentum += (prt_position[0] * prt_velosity[1] - prt_position[1] * prt_velosity[0]) * mass;
	}
	return sum_momentum;
}

double calc_norm_of_impulse(std::vector<Particle>& particles, int start, int end) {
	std::vector<double> sum_impulse = { 0, 0 };
	for (int i = start; i < end; i++) {
		std::vector<double> prt_velosity = particles[i].get_velosity();
		sum_impulse[0] += prt_velosity[0];
		sum_impulse[1] += prt_velosity[1];
	}
	return sqrt(sum_impulse[0] * sum_impulse[0] + sum_impulse[1] * sum_impulse[1]);
}