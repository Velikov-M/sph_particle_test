#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

class Particle {
private:

	char type;
	std::vector<int> cell_id;
	std::vector<double> position;
	std::vector<double> velosity;

public:
	Particle();

	void set_position(double x1, double x2);

	inline std::vector<double> get_position() { return position;};

	void print_position();

	void set_velosity(double v1, double v2);

	inline std::vector<double> get_velosity() { return velosity;};

	void print_velosity();

	void set_type(char c);

	inline char get_type() { return type;};

	void set_cell_id(std::vector<int> id);

	inline std::vector<int> get_cell_id(){ return cell_id;};

};

#endif
