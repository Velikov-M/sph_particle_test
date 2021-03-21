#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Particle.h"

void create_XY(std::vector<Particle>& particles, int num) {
	std::ofstream xyz_file;
	xyz_file.open("prt_positions.xyz");
	xyz_file << num << std::endl;
	xyz_file << std::endl;
	for (int i = 0; i < num; i++) {
		std::vector<double> cur_pos = particles[i].get_position();
		xyz_file << "P  " << cur_pos[0] << "  " << cur_pos[1] << std::endl;
	}
	xyz_file.close();
}

void add_ts_XY(std::vector<Particle>& particles, int num) {
	std::ofstream foutput;
	std::ifstream finput;

	finput.open("prt_positions.xyz");
	foutput.open("prt_positions.xyz",std::ios::app);
	
	if (finput.is_open()) {
		foutput << num << std::endl;
		foutput << std::endl;
		for (int i = 0; i < num; i++) {
			std::vector<double> cur_pos = particles[i].get_position();
			foutput << "P  " << cur_pos[0] << "  " << cur_pos[1] << std::endl;
		}
	}
	finput.close();
	foutput.close();
}

int read_initial_position(std::vector<Particle>& particles,int n) {
	std::ifstream finput;
	finput.open("initial_prt_pos.xyz");
	std::string line;
	int i = 0; //we have to know how many prts we read already
	while (std::getline(finput, line)) {
		particles.push_back(Particle());
		std::stringstream lineStream(line);
		std::string cell;
		std::getline(lineStream, cell, ' ');
		std::getline(lineStream, cell, ' ');
		double x_prt = stod(cell);
		std::getline(lineStream, cell, ' ');
		double y_prt = stod(cell);
		particles[n + i].set_position(x_prt, y_prt);
		i++;
	}
	finput.close();
	return i;
}