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
#include "ParticleTrace.h"
#include "vector_utils.h"
#include "biot_savart_routines.h"

//constexpr auto OUTPUT_PRECISION = 20;
//constexpr auto TRACES_DIR = "D:\\OneDrive - University of Calgary\\PARA\\Projects\\Inventor\\Magnetic Quadruple\\comsol\\test_2\\particle_data";

void generate_uniform_grid(vector<vector<PRECISION_TYPE>>& grid_xyz,
	PRECISION_TYPE x_min, PRECISION_TYPE x_max,
	PRECISION_TYPE y_min, PRECISION_TYPE y_max,
	PRECISION_TYPE z_min, PRECISION_TYPE z_max,
	uint32_t num_points) {

	grid_xyz.clear();
	grid_xyz.reserve(num_points * num_points * num_points);

	PRECISION_TYPE x_step = (x_max - x_min) / (num_points - 1);
	PRECISION_TYPE y_step = (y_max - y_min) / (num_points - 1);
	PRECISION_TYPE z_step = (z_max - z_min) / (num_points - 1);

	for (uint32_t i = 0; i < num_points; i++) {
		PRECISION_TYPE x = x_min + i * x_step;
		for (uint32_t j = 0; j < num_points; j++) {
			PRECISION_TYPE y = y_min + j * y_step;
			for (uint32_t k = 0; k < num_points; k++) {
				PRECISION_TYPE z = z_min + k * z_step;
				grid_xyz.push_back({ x, y, z });
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


void save_magnetic_field(	const string& filename,
							const vector<vector<PRECISION_TYPE>>& grid_xyz,
							const vector<vector<PRECISION_TYPE>>& bfield_xyz,
							uint32_t num_grid) {

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

	// Use ostringstream for proper precision formatting
	ostringstream line;
	for (uint32_t i = 0; i < num_grid; i++) {
		line.str(""); // Clear previous contents
		line.clear(); // Reset error flags

		line << fixed << setprecision(OUTPUT_PRECISION)
			<< grid_xyz[i][0] << "," << grid_xyz[i][1] << "," << grid_xyz[i][2] << ","
			<< bfield_xyz[i][0] << "," << bfield_xyz[i][1] << "," << bfield_xyz[i][2] << "\n";

		outfile << line.str();  // Write formatted string
	}

	outfile.close();
	cout << "Magnetic field data saved to " << filename << endl;
}

// Function to generate a straight wire
void generate_straight_wire(uint32_t num_segments,
							PRECISION_TYPE zmin,
							PRECISION_TYPE dz,
							vector<vector<PRECISION_TYPE>>& segment_start,
							vector<vector<PRECISION_TYPE>>& segment_end) {

	segment_start.clear();
	segment_end.clear();
	for (uint32_t seg = 0; seg < num_segments; seg++) {
		segment_start.push_back({ 0.0, 0.0, zmin + seg * dz });
		segment_end.push_back({ 0.0, 0.0, zmin + (seg + 1) * dz });
	}
}

// Function to generate a circular loop
void generate_circular_loop(vector<vector<PRECISION_TYPE>>& segment_start,
	vector<vector<PRECISION_TYPE>>& segment_end,
	PRECISION_TYPE radius, int num_segments) {

	segment_start.clear();
	segment_end.clear();

	PRECISION_TYPE dtheta = 2.0 * PI / num_segments;

	for (int i = 0; i < num_segments; i++) {
		PRECISION_TYPE theta1 = i * dtheta;
		PRECISION_TYPE theta2 = (i + 1) * dtheta;

		vector<PRECISION_TYPE> s0 = { radius * cos(theta1), radius * sin(theta1), 0.0 }; // Start of segment
		vector<PRECISION_TYPE> s1 = { radius * cos(theta2), radius * sin(theta2), 0.0 }; // End of segment

		segment_start.push_back(s0);
		segment_end.push_back(s1);
	}
}


void test_circular_loop(const uint32_t num_segments, PRECISION_TYPE loop_radius, const uint32_t grid_size) {

	// Create Loop
	vector<vector<PRECISION_TYPE>> segment_start, segment_end;
	generate_circular_loop(segment_start, segment_end, loop_radius, num_segments);

	// Grid
	vector<vector<PRECISION_TYPE>> grid;
	generate_uniform_grid(grid, -loop_radius*5, loop_radius*5, -loop_radius*5, loop_radius*5, -loop_radius*5, loop_radius*5, grid_size);

	// Initialize magnetic field vectors
	vector<vector<PRECISION_TYPE>> bfield_chukman(grid.size(), vector<PRECISION_TYPE>(3, 0.0));
	vector<vector<PRECISION_TYPE>> bfield_se(grid.size(), vector<PRECISION_TYPE>(3, 0.0));

	calc_bfield_parallel_chukman(num_segments, segment_start, segment_end, grid.size(), grid, bfield_chukman);
	calc_bfield_parallel_se(num_segments, segment_start, segment_end, grid.size(), grid, bfield_se);

	// Save magnetic field data
	save_magnetic_field("B_field_loop_chukman.txt", grid, bfield_chukman, grid.size());
	save_magnetic_field("B_field_loop_se.txt", grid, bfield_se, grid.size());

}


void test_wire_z(const uint32_t num_segments, PRECISION_TYPE zmin, PRECISION_TYPE zmax, const uint32_t grid_size) {
	
	// Create Wire 
	PRECISION_TYPE dz = (zmax - zmin) / num_segments;
	vector<vector<PRECISION_TYPE>> segment_start, segment_end;
	generate_straight_wire(num_segments, zmin, dz, segment_start, segment_end);

	// Grid
	vector<vector<PRECISION_TYPE>> grid;
	generate_uniform_grid(grid, zmin, zmax, zmin, zmax, zmin, zmax, grid_size);

	// Initialize magnetic field vectors
	vector<vector<PRECISION_TYPE>> bfield_chukman(grid.size(), vector<PRECISION_TYPE>(3, 0.0));
	vector<vector<PRECISION_TYPE>> bfield_se(grid.size(), vector<PRECISION_TYPE>(3, 0.0));

	//calc_bfield_parallel_chukman(num_segments, segment_start, segment_end, grid.size(), grid, bfield_chukman);
	calc_bfield_parallel_se(num_segments, segment_start, segment_end, grid.size(), grid, bfield_se);

	// Save magnetic field data
	//save_magnetic_field("B_field_wire_chukman.txt", grid, bfield_chukman, grid.size());
	save_magnetic_field("B_field_wire_se.txt", grid, bfield_se, grid.size());


}


void test_particle_trace(const uint32_t grid_size) {

	// Read grid bounds
	PRECISION_TYPE x_min, x_max, y_min, y_max, z_min, z_max, x_range, x_midpoint, y_range, y_midpoint, z_range, z_midpoint;

	if (!read_grid_bounds("constants_test21.csv",	x_min, x_max, x_range, x_midpoint,
													y_min, y_max, y_range, y_midpoint,
													z_min, z_max, z_range, z_midpoint)) {

		throw runtime_error("Failed to read grid bounds from constants_test21.csv");
	}

	cout << "Grid bounds: " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << endl;
	
	// Create grid
	vector<vector<PRECISION_TYPE>> grid;
	generate_uniform_grid(grid, x_min - x_range / 2, x_max + x_range / 2,
								y_min - y_range / 2, y_max + y_range / 2,
								z_min - z_range / 2, z_max + z_range / 2, grid_size);


	string directory = TRACES_DIR;
	string filepath = directory + "\\particle_1.csv";
	ParticleTrace particle_trace(filepath);

	// Retrieve segments
	const auto& segment_start = particle_trace.getSegmentStart();
	const auto& segment_end = particle_trace.getSegmentEnd();

	// Initialize magnetic field vectors
	vector<vector<PRECISION_TYPE>> bfield_chukman(grid.size(), vector<PRECISION_TYPE>(3, 0.0));
	vector<vector<PRECISION_TYPE>> bfield_se(grid.size(), vector<PRECISION_TYPE>(3, 0.0));

	//calc_bfield_parallel_chukman(particle_trace.getNumSegments(), segment_start, segment_end, grid.size(), grid, bfield_chukman);
	calc_bfield_parallel_se(particle_trace.getNumSegments(), segment_start, segment_end, grid.size(), grid, bfield_se);

	// Save magnetic field data
	//save_magnetic_field("B_field_particle_trace_chukman.txt", grid, bfield_chukman, grid.size());
	save_magnetic_field("B_field_particle_trace_se.txt", grid, bfield_se, grid.size());

}

void test_all_traces() {


}

#ifdef _OPENMP
	#define OPENMP_ENABLED true
#else
	#define OPENMP_ENABLED false
#endif

int main() {
	cout << fixed << setprecision(OUTPUT_PRECISION);
	omp_set_num_threads(NUM_THREADS); // Force 4 threads

	test_wire_z(100, -10.0, 10.0, 30);
	//test_circular_loop(100, 2.0, 30);
	//test_particle_trace(15);

	/*#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		std::cout << "Hello from thread " << thread_id << std::endl;
	}
	return 0;*/

	//cout << "OpenMP enabled: " << (OPENMP_ENABLED ? "Yes" : "No") << endl;

	//#pragma omp parallel
	//{
	//	int thread_id = omp_get_thread_num();
	//	std::cout << "Hello from thread " << thread_id << std::endl;
	//}
	//return 0;



	//std::cout << "Max available threads: " << omp_get_max_threads() << std::endl;
	//return 0;


	//PRECISION_TYPE x_min, x_max, y_min, y_max, z_min, z_max, x_range, x_midpoint, y_range, y_midpoint, z_range, z_midpoint;

	//if (!read_grid_bounds("constants_test21.csv", x_min, x_max, x_range, x_midpoint,
	//	y_min, y_max, y_range, y_midpoint,
	//	z_min, z_max, z_range, z_midpoint)) {
	//	return 1; // Exit if reading fails
	//}

	//cout << "Grid bounds: " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << endl;

	//// Uniform Grid
	//int num_pts = 15;
	//vector<vector<PRECISION_TYPE>> grid;
	//generate_uniform_grid(grid, x_min - x_range / 2, x_max + x_range / 2,
	//							y_min - y_range / 2, y_max + y_range / 2,
	//							z_min - z_range / 2, z_max + z_range / 2, num_pts);










}

