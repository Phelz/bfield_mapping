using namespace std;

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "ParticleTrace.h"
#include "vector_utils.h"
#include "biot_savart_routines.h"


void generate_uniform_grid(vector<vector<long double>>& grid_xyz,
	long double x_min, long double x_max,
	long double y_min, long double y_max,
	long double z_min, long double z_max,
	uint32_t num_points) {

	grid_xyz.clear();
	grid_xyz.reserve(num_points * num_points * num_points);

	long double x_step = (x_max - x_min) / (num_points - 1);
	long double y_step = (y_max - y_min) / (num_points - 1);
	long double z_step = (z_max - z_min) / (num_points - 1);

	for (uint32_t i = 0; i < num_points; i++) {
		long double x = x_min + i * x_step;
		for (uint32_t j = 0; j < num_points; j++) {
			long double y = y_min + j * y_step;
			for (uint32_t k = 0; k < num_points; k++) {
				long double z = z_min + k * z_step;
				grid_xyz.push_back({ x, y, z });
			}
		}
	}
}



//Function to read min/max values from CSV file
bool read_grid_bounds(const string& filename,
	long double& x_min, long double& x_max, long double& x_range, long double& x_midpoint,
	long double& y_min, long double& y_max, long double& y_range, long double& y_midpoint,
	long double& z_min, long double& z_max, long double& z_range, long double& z_midpoint) {


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



// Function to generate a circular loop
void generate_circular_loop(vector<vector<long double>>& segment_start, vector<vector<long double>>& segment_end, long double radius, int num_segments) {
	segment_start.clear();
	segment_end.clear();

	long double dtheta = 2.0 * PI / num_segments;

	for (int i = 0; i < num_segments; i++) {
		long double theta1 = i * dtheta;
		long double theta2 = (i + 1) * dtheta;

		vector<long double> s0 = { radius * cos(theta1), radius * sin(theta1), 0.0 }; // Start of segment
		vector<long double> s1 = { radius * cos(theta2), radius * sin(theta2), 0.0 }; // End of segment

		segment_start.push_back(s0);
		segment_end.push_back(s1);
	}
}

int main() {
	cout << fixed << setprecision(20);

	long double x_min, x_max, y_min, y_max, z_min, z_max, x_range, x_midpoint, y_range, y_midpoint, z_range, z_midpoint;

	if (!read_grid_bounds("constants_test21.csv", x_min, x_max, x_range, x_midpoint,
		y_min, y_max, y_range, y_midpoint,
		z_min, z_max, z_range, z_midpoint)) {
		return 1; // Exit if reading fails
	}

	cout << "Grid bounds: " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << endl;

	// Uniform Grid
	int num_pts = 15;
	vector<vector<long double>> grid;
	generate_uniform_grid(grid, x_min - x_range / 2, x_max + x_range / 2,
		y_min - y_range / 2, y_max + y_range / 2,
		z_min - z_range / 2, z_max + z_range / 2, num_pts);

	// A simple wire in the z-axis
	//int num_segments = 10;
	//vector<vector<long double>> segment_start, segment_end;
	//long double zmin = -10.0, zmax = 10.0;
	//long double dz = (zmax - zmin) / num_segments;

	//for (int seg = 0; seg < num_segments; seg++) {
	//	segment_start.push_back({ 0.0, 0.0, zmin + seg * dz });
	//	segment_end.push_back({ 0.0, 0.0, zmin + (seg + 1) * dz });
	//}

	// Circular loop
	//long double loop_radius = 2.0;
	//int num_segments = 100;
	//vector<vector<long double>> segment_start, segment_end;
	//generate_circular_loop(segment_start, segment_end, loop_radius, num_segments);

	string directory = "D:\\OneDrive - University of Calgary\\PARA\\Projects\\Inventor\\Magnetic Quadruple\\comsol\\test_2\\particle_data";
	string filepath = directory + "\\particle_0.csv";
	ParticleTrace particle_trace(filepath);

	// Retrieve segments
	const auto& segment_start = particle_trace.getSegmentStart();
	const auto& segment_end = particle_trace.getSegmentEnd();
	const size_t num_segments = particle_trace.getNumSegments();



	//ofstream output_file_chukman("B_field_loop_chukman.txt");
	//ofstream output_file_se("B_field_loop_se.txt");
	//ofstream output_file_chukman("B_field_particle_trace_chukman.txt");
	ofstream output_file_se("B_field_particle_trace_se.txt");
	//ofstream output_file_chukman("B_field_wire_chukman.txt");
	//ofstream output_file_se("B_field_wire_se.txt");

	//output_file_chukman.precision(20);
	output_file_se.precision(20);

	//if (!output_file_chukman.is_open()) {
	//	cerr << "Error: Could not open file for writing!" << endl;
	//	return 1;
	//}

	if (!output_file_se.is_open()) {
		cerr << "Error: Could not open file for writing!" << endl;
		return 1;
	}

	//output_file_chukman << "x y z Bx By Bz\n";
	output_file_se << "x y z Bx By Bz\n";


	size_t total_points = grid.size();
	size_t count = 0;
	size_t update_interval = max(total_points / 10, size_t(1)); // Update every 10% or at least once

	// Call calculate_b0 for the entire grid
	vector<vector<long double>> bfield_chukman(grid.size(), vector<long double>(3, 0.0));
	vector<vector<long double>> bfield_se(grid.size(), vector<long double>(3, 0.0));

	//calc_bfield_parallel(num_segments, segment_start, segment_end, grid.size(), grid, bfield_chukman, biot_savart_chukman);
	calc_bfield_parallel(num_segments, segment_start, segment_end, grid.size(), grid, bfield_chukman, biot_savart_se);













	//for (const auto& r : grid) {
		//vector<long double> B_chukman = { 0.0, 0.0, 0.0 }; // Initialize B-field for Chukman routine

		//for (size_t i = 0; i < segment_start.size(); i++) {
		//	vector<long double> B_temp(3, 0.0);
		//	biot_savart_chukman(segment_start[i], segment_end[i], r, B_temp);
		//	B_chukman = add(B_chukman, B_temp);
		//}

		// Extract the corresponding magnetic field from the b1_xyz grid for calculate_b0
		//size_t grid_index = &r - &grid[0];  // Calculate the index of the current grid point
		//vector<long double> B_se = b1_xyz[grid_index];  // Extract the magnetic field at the current grid point

		//// Write to Chukman output file: x y z Bx By Bz
		//output_file_chukman << r[0] << " " << r[1] << " " << r[2] << " "
		//	<< B_chukman[0] << " " << B_chukman[1] << " " << B_chukman[2] << "\n";

		// Write to calculate_b0 output file: x y z Bx By Bz
		//output_file_se << r[0] << " " << r[1] << " " << r[2] << " "
			//<< B_se[0] << " " << B_se[1] << " " << B_se[2] << "\n";

		// Progress update
		//count++;
		//if (count % update_interval == 0) {
			//cout << "Progress: " << (100 * count / total_points) << "% ("
				//<< count << "/" << total_points << " points processed)" << endl;
		//}
	//}


	//output_file_chukman.close();
	output_file_se.close();

	cout << "Done!" << endl;

}

