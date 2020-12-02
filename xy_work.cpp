#include <iostream>
#include <fstream>
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