using namespace std;

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <omp.h>

#include "constants.h"
#include "vector_utils.h"
#include "biot_savart_routines.h"
#include "ParticleTrace.h"


void generate_uniform_grid(vector<PRECISION_TYPE>& grid_xyz,
							PRECISION_TYPE x_min, PRECISION_TYPE x_max,
							PRECISION_TYPE y_min, PRECISION_TYPE y_max,
							PRECISION_TYPE z_min, PRECISION_TYPE z_max,
							uint32_t num_points) {

	grid_xyz.clear();
	grid_xyz.reserve(3 * num_points * num_points * num_points);

	PRECISION_TYPE x_step = (x_max - x_min) / (num_points - 1);
	PRECISION_TYPE y_step = (y_max - y_min) / (num_points - 1);
	PRECISION_TYPE z_step = (z_max - z_min) / (num_points - 1);

	for (uint32_t i = 0; i < num_points; i++) {
		PRECISION_TYPE x = x_min + i * x_step;
		for (uint32_t j = 0; j < num_points; j++) {
			PRECISION_TYPE y = y_min + j * y_step;
			for (uint32_t k = 0; k < num_points; k++) {
				PRECISION_TYPE z = z_min + k * z_step;
				grid_xyz.push_back(x);
				grid_xyz.push_back(y);
				grid_xyz.push_back(z);
			}
		}
	}
}



//Function to read min/max values from CSV file
bool read_grid_bounds(const string& filename,
	PRECISION_TYPE& x_min, PRECISION_TYPE& x_max, PRECISION_TYPE& x_range, PRECISION_TYPE& x_midpoint,
	PRECISION_TYPE& y_min, PRECISION_TYPE& y_max, PRECISION_TYPE& y_range, PRECISION_TYPE& y_midpoint,
	PRECISION_TYPE& z_min, PRECISION_TYPE& z_max, PRECISION_TYPE& z_range, PRECISION_TYPE& z_midpoint) {


	ifstream file(filename);
	if (!file.is_open()) {
		cerr << "Error: Could not open file " << filename << endl;

		// Check if the file exists or is accessible
		if (errno == ENOENT) {
			cerr << "Error: The file does not exist." << endl;
		}
		else if (errno == EACCES) {
			cerr << "Error: Permission denied to open the file." << endl;
		}
		else {
			cerr << "Error: Unknown error occurred. Error code: " << errno << "" << endl;
		}

		return false;
	}

	string line;
	getline(file, line); // Read the first line (header)
	getline(file, line); // Read the data line

	stringstream ss(line);
	char comma; // For handling commas in CSV

	// Read values in order
	ss >> x_max >> comma >> y_max >> comma >> z_max >> comma >> x_min >> comma >> y_min >> comma >> z_min >> comma >> x_range >> comma >> y_range >> comma >> z_range >> comma >> x_midpoint >> comma >> y_midpoint >> comma >> z_midpoint;


	file.close();
	return true;
}




void save_magnetic_field(const string& filename,
						const PRECISION_TYPE* grid_xyz,
						const PRECISION_TYPE* bfield_xyz,
						uint32_t num_points) {

	// Open the file
	ofstream outfile(filename, ios::out | ios::trunc);
	if (!outfile.is_open()) {
		cerr << "Error: Could not open file " << filename << " for writing!" << endl;
		return;
	}

	// Set precision
	outfile << fixed << setprecision(OUTPUT_PRECISION);

	// Write the header
	outfile << "x,y,z,Bx,By,Bz\n";

	// Write data
	for (uint32_t i = 0; i < num_points; i++) {
		outfile << grid_xyz[3 * i] << "," << grid_xyz[3 * i + 1] << "," << grid_xyz[3 * i + 2] << ","
			<< bfield_xyz[3 * i] << "," << bfield_xyz[3 * i + 1] << "," << bfield_xyz[3 * i + 2] << "\n";
	}

	outfile.close();
	cout << "Magnetic field data saved to " << filename << endl;
}





void generate_straight_wire(uint32_t num_segments,
	PRECISION_TYPE zmin,
	PRECISION_TYPE dz,
	std::vector<PRECISION_TYPE>& segment_start,
	std::vector<PRECISION_TYPE>& segment_end) {

	segment_start.clear();
	segment_end.clear();
	segment_start.reserve(3 * num_segments);
	segment_end.reserve(3 * num_segments);

	for (uint32_t seg = 0; seg < num_segments; seg++) {
		segment_start.push_back(0.0);
		segment_start.push_back(0.0);
		segment_start.push_back(zmin + seg * dz);

		segment_end.push_back(0.0);
		segment_end.push_back(0.0);
		segment_end.push_back(zmin + (seg + 1) * dz);
	}
}



void generate_circular_loop(std::vector<PRECISION_TYPE>& segment_start,
	std::vector<PRECISION_TYPE>& segment_end,
	PRECISION_TYPE radius,
	uint32_t num_segments) {

	segment_start.clear();
	segment_end.clear();
	segment_start.reserve(3 * num_segments);
	segment_end.reserve(3 * num_segments);

	PRECISION_TYPE dtheta = 2.0 * PI / num_segments;

	for (uint32_t i = 0; i < num_segments; i++) {
		PRECISION_TYPE theta1 = i * dtheta;
		PRECISION_TYPE theta2 = (i + 1) * dtheta;

		segment_start.push_back(radius * cos(theta1));
		segment_start.push_back(radius * sin(theta1));
		segment_start.push_back(0.0);

		segment_end.push_back(radius * cos(theta2));
		segment_end.push_back(radius * sin(theta2));
		segment_end.push_back(0.0);
	}
}


void test_circular_loop(const uint32_t num_segments, PRECISION_TYPE loop_radius, const uint32_t grid_size) {
	// Create Loop
	vector<PRECISION_TYPE> segment_start(num_segments * 3);
	vector<PRECISION_TYPE> segment_end(num_segments * 3);
	generate_circular_loop(segment_start, segment_end, loop_radius, num_segments);

	// Grid
	vector<PRECISION_TYPE> grid(grid_size * grid_size * grid_size * 3);
	generate_uniform_grid(grid, -loop_radius * 5, loop_radius * 5,
		-loop_radius * 5, loop_radius * 5,
		-loop_radius * 5, loop_radius * 5, grid_size);

	// Initialize magnetic field vectors
	vector<PRECISION_TYPE> bfield_chukman(grid.size(), 0.0);
	vector<PRECISION_TYPE> bfield_se(grid.size(), 0.0);

	// Pass raw pointers using .data()
	calc_bfield_parallel_chukman(num_segments, segment_start.data(), segment_end.data(),
		grid.size() / 3, grid.data(), bfield_chukman.data());

	calc_bfield_parallel_se(num_segments, segment_start.data(), segment_end.data(),
		grid.size() / 3, grid.data(), bfield_se.data());


	// Save magnetic field data
	save_magnetic_field("B_field_loop_chukman.txt", grid.data(), bfield_chukman.data(), grid.size() / 3);
	save_magnetic_field("B_field_loop_se.txt", grid.data(), bfield_se.data(), grid.size() / 3);
}



void test_wire_z(const uint32_t num_segments, PRECISION_TYPE zmin, PRECISION_TYPE zmax, const uint32_t grid_size) {
	// Create Wire
	PRECISION_TYPE dz = (zmax - zmin) / num_segments;
	vector<PRECISION_TYPE> segment_start(num_segments * 3);
	vector<PRECISION_TYPE> segment_end(num_segments * 3);
	generate_straight_wire(num_segments, zmin, dz, segment_start, segment_end);

	// Grid
	vector<PRECISION_TYPE> grid(grid_size * grid_size * grid_size * 3);
	generate_uniform_grid(grid, zmin, zmax, zmin, zmax, zmin, zmax, grid_size);

	// Initialize magnetic field vectors
	vector<PRECISION_TYPE> bfield_se(grid.size(), 0.0);
	vector<PRECISION_TYPE> bfield_chukman(grid.size(), 0.0);

	// Pass raw pointers using .data()
	calc_bfield_parallel_chukman(num_segments, segment_start.data(), segment_end.data(),
		grid.size() / 3, grid.data(), bfield_chukman.data());
	calc_bfield_parallel_se(num_segments, segment_start.data(), segment_end.data(),
		grid.size() / 3, grid.data(), bfield_se.data());

	// Save magnetic field data
	save_magnetic_field("B_field_wire_se.txt", grid.data(), bfield_se.data(), grid.size() / 3);
	save_magnetic_field("B_field_wire_chukman.txt", grid.data(), bfield_chukman.data(), grid.size() / 3);
}









//void test_particle_trace(const uint32_t grid_size) {
//
//	// Read grid bounds
//	PRECISION_TYPE x_min, x_max, y_min, y_max, z_min, z_max, x_range, x_midpoint, y_range, y_midpoint, z_range, z_midpoint;
//
//	if (!read_grid_bounds("constants_test21.csv",	x_min, x_max, x_range, x_midpoint,
//													y_min, y_max, y_range, y_midpoint,
//													z_min, z_max, z_range, z_midpoint)) {
//
//		throw runtime_error("Failed to read grid bounds from constants_test21.csv");
//	}
//
//	cout << "Grid bounds: " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << endl;
//	
//	// Create grid
//	vector<vector<PRECISION_TYPE>> grid;
//	generate_uniform_grid(grid, x_min - x_range / 2, x_max + x_range / 2,
//								y_min - y_range / 2, y_max + y_range / 2,
//								z_min - z_range / 2, z_max + z_range / 2, grid_size);
//
//
//	string directory = TRACES_DIR;
//	string filepath = directory + "\\particle_1.csv";
//	ParticleTrace particle_trace(filepath);
//
//	// Retrieve segments
//	const auto& segment_start = particle_trace.getSegmentStart();
//	const auto& segment_end = particle_trace.getSegmentEnd();
//
//	// Initialize magnetic field vectors
//	vector<vector<PRECISION_TYPE>> bfield_chukman(grid.size(), vector<PRECISION_TYPE>(3, 0.0));
//	vector<vector<PRECISION_TYPE>> bfield_se(grid.size(), vector<PRECISION_TYPE>(3, 0.0));
//
//	//calc_bfield_parallel_chukman(particle_trace.getNumSegments(), segment_start, segment_end, grid.size(), grid, bfield_chukman);
//	calc_bfield_parallel_se(particle_trace.getNumSegments(), segment_start, segment_end, grid.size(), grid, bfield_se);
//
//	// Save magnetic field data
//	//save_magnetic_field("B_field_particle_trace_chukman.txt", grid, bfield_chukman, grid.size());
//	save_magnetic_field("B_field_particle_trace_se.txt", grid, bfield_se, grid.size());
//
//}
//
//void test_all_traces() {
//
//
//}



int main() {
	cout << fixed << setprecision(OUTPUT_PRECISION);
	omp_set_num_threads(NUM_THREADS); // Force 4 threads

	test_wire_z(100, -10.0, 10.0, 50);
	test_circular_loop(100, 2.0, 50);
	//test_particle_trace(15);


}

